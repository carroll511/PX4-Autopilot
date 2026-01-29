/****************************************************************************
 *
 *   Copyright (c) 2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file ControlAllocationStatic.cpp
 *
 * Static Control Allocation Algorithm Implementation
 * https://arxiv.org/abs/1801.04581
 */

#include "ControlAllocationStatic.hpp"
#include <matrix/matrix/math.hpp>
#include <mathlib/math/Limits.hpp>
#include <drivers/drv_hrt.h>
#include <px4_platform_common/log.h>
#include <cmath>

void
ControlAllocationStatic::setEffectivenessMatrix(
	const matrix::Matrix<float, ControlAllocation::NUM_AXES, ControlAllocation::NUM_ACTUATORS> &effectiveness,
	const ActuatorVector &actuator_trim, const ActuatorVector &linearization_point, int num_actuators,
	bool update_normalization_scale)
{
	ControlAllocation::setEffectivenessMatrix(effectiveness, actuator_trim, linearization_point, num_actuators,
			update_normalization_scale);
	_A_static_update_needed = true;
}

void
ControlAllocationStatic::setArmAngles(float theta1, float theta2, float theta3, float theta4)
{
	_arm_angles[0] = theta1;
	_arm_angles[1] = theta2;
	_arm_angles[2] = theta3;
	_arm_angles[3] = theta4;
	_A_static_update_needed = true;
}

void
ControlAllocationStatic::setParameters(float l, float k)
{
	_l = l;
	_k = k;
	_A_static_update_needed = true;
}

void
ControlAllocationStatic::buildStaticMatrix()
{
	// Clear A_static
	_A_static.setZero();
	
	// Build static allocation matrix for 8 motors (coaxial configuration) with fixed arm angles
	// Matrix structure:
	// [F; M] = A_static * F_dec
	// F_dec = [F_v1, F_l1, F_v2, F_l2, ..., F_v8, F_l8]  (16 columns total)
	//
	// Motor pairs (coaxial):
	// - Motors 1 & 6: theta1 (columns 0,1 and 10,11)
	// - Motors 2 & 5: theta2 (columns 2,3 and 8,9)
	// - Motors 3 & 7: theta3 (columns 4,5 and 12,13)
	// - Motors 4 & 8: theta4 (columns 6,7 and 14,15)
	//
	// For 8 motors, the matrix is 6x16:
	// Each row follows the 4-motor pattern extended to 8 motors with paired angles
	
	// Precompute sin and cos for each arm angle
	float sin_theta[4];
	float cos_theta[4];
	for (int i = 0; i < 4; i++) {
		sin_theta[i] = sinf(_arm_angles[i]);
		cos_theta[i] = cosf(_arm_angles[i]);
	}
	
	// Motor-to-theta mapping for 8 motors (coaxial pairs)
	// Motor pairs: (1,6)=theta1, (2,5)=theta2, (3,7)=theta3, (4,8)=theta4
	// Motors: 0=1, 1=2, 2=3, 3=4, 4=5, 5=6, 6=7, 7=8
	int theta_idx[8] = {0, 1, 2, 3, 1, 0, 2, 3};  // Which theta each motor uses
	
	// ROLL row coefficients: [motor0, motor1, motor2, motor3, motor4, motor5, motor6, motor7]
	float roll_coeff[8] = {-1.f, 1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f};  // Sign for l term
	float roll_k_sign[8] = {1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f};  // Sign for k*sin term
	
	// PITCH row: motor 0 uses theta[1] instead of theta[0] (per original spec)
	// PITCH theta indices: [motor0, motor1, motor2, motor3, motor4, motor5, motor6, motor7]
	int pitch_theta_idx[8] = {1, 1, 2, 3, 1, 1, 2, 3};  // Motor 0 and 5 use theta[1] for PITCH
	float pitch_coeff[8] = {1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, -1.f};  // Sign for l term
	float pitch_k_sign[8] = {-1.f, 1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f};  // Sign for k*cos term
	
	// YAW row coefficients
	float yaw_sign[8] = {-1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, 1.f};  // Overall sign
	float yaw_cos_sin_sign[8] = {1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f};  // Sign for (cos ± sin)
	float yaw_k_sign[8] = {1.f, -1.f, 1.f, -1.f, 1.f, -1.f, -1.f, -1.f};  // Sign for k term
	
	// Build matrix using loops
	for (int motor = 0; motor < 8; motor++) {
		int col_v = 2 * motor;      // Column for F_vi (vertical component)
		int col_l = 2 * motor + 1;   // Column for F_li (lateral component)
		int t_idx = theta_idx[motor];  // Which theta this motor uses
		
		// Row 0: THRUST_X
		_A_static(ControlAxis::THRUST_X, col_v) = 0.f;
		_A_static(ControlAxis::THRUST_X, col_l) = sin_theta[t_idx];
		
		// Row 1: THRUST_Y
		_A_static(ControlAxis::THRUST_Y, col_v) = 0.f;
		_A_static(ControlAxis::THRUST_Y, col_l) = -cos_theta[t_idx];
		
		// Row 2: THRUST_Z
		_A_static(ControlAxis::THRUST_Z, col_v) = -1.f;
		_A_static(ControlAxis::THRUST_Z, col_l) = 0.f;
		
		// Row 3: ROLL
		_A_static(ControlAxis::ROLL, col_v) = roll_coeff[motor] * _l + roll_k_sign[motor] * _k * sin_theta[t_idx];
		_A_static(ControlAxis::ROLL, col_l) = 0.f;
		
		// Row 4: PITCH
		int pitch_t_idx = pitch_theta_idx[motor];
		_A_static(ControlAxis::PITCH, col_v) = pitch_coeff[motor] * _l + pitch_k_sign[motor] * _k * cos_theta[pitch_t_idx];
		_A_static(ControlAxis::PITCH, col_l) = 0.f;
		
		// Row 5: YAW
		_A_static(ControlAxis::YAW, col_v) = 0.f;
		float cos_sin_term = yaw_cos_sin_sign[motor] > 0.f ? 
			(cos_theta[t_idx] + sin_theta[t_idx]) : 
			(cos_theta[t_idx] - sin_theta[t_idx]);
		_A_static(ControlAxis::YAW, col_l) = yaw_sign[motor] * _l * cos_sin_term + yaw_k_sign[motor] * _k;
	}
}

void
ControlAllocationStatic::updateStaticMatrix()
{
	if (_A_static_update_needed) {
		buildStaticMatrix();
		
		// Compute pseudo-inverse of A_static
		matrix::geninv(_A_static, _A_static_inv);
		
		// Debug: Print matrix once when it's built
		// static bool matrix_printed = false;
		// if (!matrix_printed) {
		// 	PX4_INFO("=== Static Allocation Matrix Built ===");
		// 	PX4_INFO("A_static (6x16):");
		// 	_A_static.print();
		// 	PX4_INFO("Arm angles: [%.3f, %.3f, %.3f, %.3f] rad",
		// 		(double)_arm_angles[0], (double)_arm_angles[1], 
		// 		(double)_arm_angles[2], (double)_arm_angles[3]);
		// 	PX4_INFO("Parameters: l=%.3f, k=%.3f", (double)_l, (double)_k);
		// 	PX4_INFO("Num actuators: %d", _num_actuators);
		// 	matrix_printed = true;
		// }
		
		_A_static_update_needed = false;
	}
}

void
ControlAllocationStatic::allocate()
{
	// Update A_static and its inverse if needed
	updateStaticMatrix();
	
	_prev_actuator_sp = _actuator_sp;
	
	// Compute control error: (F; M) - (F; M)_trim
	matrix::Vector<float, NUM_AXES> control_error = _control_sp - _control_trim;
	
	// Debug: Print control setpoints periodically (every ~100ms to avoid spam)
	// static hrt_abstime last_debug_print = 0;
	// hrt_abstime now = hrt_absolute_time();
	// bool should_debug = (now - last_debug_print > 100000); // 100ms
	
	// if (should_debug) {
	// 	PX4_INFO("Control SP: [T:%.3f,%.3f,%.3f] [F:%.3f,%.3f,%.3f]",
	// 		(double)_control_sp(0), (double)_control_sp(1), (double)_control_sp(2),
	// 		(double)_control_sp(3), (double)_control_sp(4), (double)_control_sp(5));
	// 	PX4_INFO("Control trim: [T:%.3f,%.3f,%.3f] [F:%.3f,%.3f,%.3f]",
	// 		(double)_control_trim(0), (double)_control_trim(1), (double)_control_trim(2),
	// 		(double)_control_trim(3), (double)_control_trim(4), (double)_control_trim(5));
	// 	PX4_INFO("Control error: [T:%.3f,%.3f,%.3f] [F:%.3f,%.3f,%.3f]",
	// 		(double)control_error(0), (double)control_error(1), (double)control_error(2),
	// 		(double)control_error(3), (double)control_error(4), (double)control_error(5));
	// 	last_debug_print = now;
	// }
	
	// Solve F_dec = A_static^{-1} * control_error
	matrix::Vector<float, 2 * NUM_ACTUATORS> F_dec = _A_static_inv * control_error;
	
	// if (should_debug) {
	// 	PX4_INFO("F_dec (first 4 motors): [%.3f,%.3f] [%.3f,%.3f] [%.3f,%.3f] [%.3f,%.3f]",
	// 		(double)F_dec(0), (double)F_dec(1), (double)F_dec(2), (double)F_dec(3),
	// 		(double)F_dec(4), (double)F_dec(5), (double)F_dec(6), (double)F_dec(7));
	// }
	
	// Convert F_dec back to actuator commands
	// Relationship: F_vi = n² * cos(alpha_i), F_li = n² * sin(alpha_i)
	// Therefore: F_vi² + F_li² = n⁴, so n² = sqrt(F_vi² + F_li²)
	// And: n = sqrt(n²) = (F_vi² + F_li²)^(1/4)
	
	// Enable debug output to diagnose erratic behavior
	static hrt_abstime last_debug_print = 0;
	hrt_abstime now = hrt_absolute_time();
	bool should_debug = (now - last_debug_print > 200000); // 200ms
	
	if (should_debug) {
		PX4_INFO("=== Static Allocation Debug ===");
		PX4_INFO("Control error: [T:%.3f,%.3f,%.3f] [F:%.3f,%.3f,%.3f]",
			(double)control_error(0), (double)control_error(1), (double)control_error(2),
			(double)control_error(3), (double)control_error(4), (double)control_error(5));
		last_debug_print = now;
	}
	
	float max_n_squared = 0.f;
	float min_n_squared = 1e10f;
	
	for (int i = 0; i < _num_actuators; i++) {
		float F_vi = F_dec(2 * i);
		float F_li = F_dec(2 * i + 1);
		
		// Compute n² from F_vi² + F_li² = n⁴
		// n² = sqrt(F_vi² + F_li²)
		float F_mag_sq = F_vi * F_vi + F_li * F_li;
		
		// Handle very small forces to avoid numerical issues
		if (F_mag_sq < 1e-6f) {
			F_mag_sq = 0.f;
		}
		
		float n_squared = sqrtf(F_mag_sq);
		
		if (n_squared > max_n_squared) {
			max_n_squared = n_squared;
		}
		if (n_squared < min_n_squared && n_squared > 1e-6f) {
			min_n_squared = n_squared;
		}
		
		// Compute n = sqrt(n²)
		float n = sqrtf(n_squared);
		
		// Compute alpha_i = atan2(F_vi, F_li)
		// Note: atan2(y, x) gives angle from x-axis
		// F_vi is the vertical component (y), F_li is the lateral component (x)
		float alpha_i = atan2f(F_vi, F_li);
		
		// Normalize alpha_i to [-π, π] (atan2 already returns in this range, but ensure it)
		const float pi_f = (float)M_PI;
		while (alpha_i > pi_f) alpha_i -= 2.0f * pi_f;
		while (alpha_i < -pi_f) alpha_i += 2.0f * pi_f;
		
		// Optional: Smooth servo angle changes to reduce erratic movement
		// This helps if servo angles are changing too rapidly
		static const float servo_smoothing = 0.8f;  // 0.0 = no smoothing, 1.0 = full smoothing
		if (servo_smoothing > 0.f && i < NUM_ACTUATORS) {
			// Compute shortest angular distance
			float angle_diff = alpha_i - _prev_motor_angles(i);
			// Normalize to [-pi, pi]
			while (angle_diff > pi_f) angle_diff -= 2.0f * pi_f;
			while (angle_diff < -pi_f) angle_diff += 2.0f * pi_f;
			// Apply smoothing
			alpha_i = _prev_motor_angles(i) + angle_diff * (1.0f - servo_smoothing);
			// Re-normalize after smoothing
			while (alpha_i > pi_f) alpha_i -= 2.0f * pi_f;
			while (alpha_i < -pi_f) alpha_i += 2.0f * pi_f;
		}
		
		// Store motor angle for servo command (normalized to [-π, π])
		_motor_angles(i) = alpha_i;
		_prev_motor_angles(i) = alpha_i;
		
		// Scale actuator command to [0, 1] range
		// The actuator setpoint should be in [0, 1] for normalized motor commands
		// Since F = n², we have n² = sqrt(F_vi² + F_li²)
		// The issue is that n² values can be very large when control errors are large
		// We need to either:
		// 1. Use n² directly and normalize it
		// 2. Use n and scale it appropriately
		// 
		// Given the large n values (3-4), it seems the scaling is off.
		// Let's try using n² normalized, which should be more stable.
		// Adjust n_squared_max based on your expected maximum force
		static const float n_squared_max = 10.0f;  // Maximum expected n² value (tune this!)
		float n_squared_normalized = fminf(n_squared / n_squared_max, 1.0f);
		
		// Set actuator setpoint using n² normalized
		// This maintains the F = n² relationship
		_actuator_sp(i) = _actuator_trim(i) + n_squared_normalized;
		
		if (should_debug && i < 4) {
			PX4_INFO("Motor %d: Fv=%.4f Fl=%.4f n2=%.4f n=%.4f alpha=%.3f rad sp=%.4f",
				i, (double)F_vi, (double)F_li, (double)n_squared, (double)n,
				(double)alpha_i, (double)_actuator_sp(i));
		}
	}
	
	if (should_debug) {
		float servo_angles[4];
		getServoAngles(servo_angles);
		PX4_INFO("Servo angles: [%.3f, %.3f, %.3f, %.3f] rad",
			(double)servo_angles[0], (double)servo_angles[1], 
			(double)servo_angles[2], (double)servo_angles[3]);
		PX4_INFO("n range: min=%.4f max=%.4f", (double)sqrtf(min_n_squared), (double)sqrtf(max_n_squared));
	}
	
	// Clip actuator setpoints to valid range
	clipActuatorSetpoint(_actuator_sp);
}

void
ControlAllocationStatic::getServoAngles(float servo_angles[4]) const
{	
	if (_num_actuators >= 8) {
		servo_angles[0] = _motor_angles(0);  // Motor 1 (pair with motor 6)
		servo_angles[1] = _motor_angles(1);  // Motor 2 (pair with motor 5)
		servo_angles[2] = _motor_angles(3);  // Motor 4 (pair with motor 8)
		servo_angles[3] = _motor_angles(2);  // Motor 3 (pair with motor 7)
	} else {
		// Fallback: if not 8 motors, fill with zeros
		for (int i = 0; i < 4; i++) {
			servo_angles[i] = 0.f;
		}
	}
}
