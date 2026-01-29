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
 * @file ControlAllocationStatic.hpp
 *
 * Static Control Allocation Algorithm
 *
 * This method uses (F; M) = A_static * F_dec, where F_dec is a vector of
 * decomposed forces (F_v1, F_l1, F_v2, F_l2, ...) where:
 * - F_vi = n^2 * cos(alpha_i) (vertical component)
 * - F_li = n^2 * sin(alpha_i) (lateral component)
 * 
 * Then F_dec = A_static^{-1} * (F; M) is solved, and converted back to
 * actuator commands: n^2 = sqrt(F_vi^2 + F_li^2) and alpha_i = atan2(F_vi, F_li)
 *
 * @author PX4 Development Team
 */

#pragma once

#include "ControlAllocation.hpp"

class ControlAllocationStatic: public ControlAllocation
{
public:
	ControlAllocationStatic() = default;
	virtual ~ControlAllocationStatic() = default;

	void allocate() override;
	void setEffectivenessMatrix(const matrix::Matrix<float, NUM_AXES, NUM_ACTUATORS> &effectiveness,
				    const ActuatorVector &actuator_trim, const ActuatorVector &linearization_point, int num_actuators,
				    bool update_normalization_scale) override;
	
	/**
	 * Set arm angles for motors (for 4 motors)
	 * @param theta1 Arm angle for motor 1 (radians)
	 * @param theta2 Arm angle for motor 2 (radians)
	 * @param theta3 Arm angle for motor 3 (radians)
	 * @param theta4 Arm angle for motor 4 (radians)
	 */
	void setArmAngles(float theta1, float theta2, float theta3, float theta4);
	
	/**
	 * Set parameters l (arm length) and k (moment coefficient)
	 * @param l Arm length parameter
	 * @param k Moment coefficient parameter
	 */
	void setParameters(float l, float k);
	
	/**
	 * Get the motor angles (alpha_i) for servo commands
	 * @return Vector of motor angles in radians
	 */
	const matrix::Vector<float, NUM_ACTUATORS> &getMotorAngles() const { return _motor_angles; }
	
	/**
	 * Get servo angles for motor pairs (4 servos total)
	 * Motor pairs: (1,6), (2,5), (4,8), (3,7) -> Servos 0,1,2,3
	 * @param servo_angles Output array for 4 servo angles in radians
	 */
	void getServoAngles(float servo_angles[4]) const;
	
	/**
	 * Get the static allocation matrix A_static for debugging
	 * @return Reference to the static allocation matrix (6x16 for 8 motors)
	 */
	const matrix::Matrix<float, NUM_AXES, 2 * NUM_ACTUATORS> &getStaticMatrix() const { return _A_static; }
	
	/**
	 * Get the pseudo-inverse of A_static for debugging
	 * @return Reference to the pseudo-inverse matrix (16x6 for 8 motors)
	 */
	const matrix::Matrix<float, 2 * NUM_ACTUATORS, NUM_AXES> &getStaticMatrixInv() const { return _A_static_inv; }
	
	/**
	 * Get arm angles for debugging
	 * @return Pointer to arm angles array [theta1, theta2, theta3, theta4]
	 */
	const float *getArmAngles() const { return _arm_angles; }
	
	/**
	 * Get parameters for debugging
	 * @param l Output arm length parameter
	 * @param k Output moment coefficient parameter
	 */
	void getParameters(float &l, float &k) const { l = _l; k = _k; }

protected:
	// Static allocation matrix: NUM_AXES x (2*NUM_ACTUATORS)
	matrix::Matrix<float, NUM_AXES, 2 * NUM_ACTUATORS> _A_static;
	
	// Pseudo-inverse of A_static: (2*NUM_ACTUATORS) x NUM_AXES
	matrix::Matrix<float, 2 * NUM_ACTUATORS, NUM_AXES> _A_static_inv;
	
	bool _A_static_update_needed{false};
	
	// Arm angles for each motor (fixed values)
	float _arm_angles[4]{0.785398f, -0.785398f, -2.35619f, 2.35619f};
	
	// Parameters: l (arm length) and k (moment coefficient)
	float _l{0.225f};  // Default arm length
	float _k{0.05f}; // Default moment coefficient
	
	// Motor angles (alpha_i) for each motor - used for servo commands
	matrix::Vector<float, NUM_ACTUATORS> _motor_angles;
	
	// Previous motor angles for smoothing (optional)
	matrix::Vector<float, NUM_ACTUATORS> _prev_motor_angles;
	
	/**
	 * Recalculate A_static and its inverse if required
	 */
	void updateStaticMatrix();
	
	/**
	 * Build A_static from fixed arm angles
	 * For 4 motors with fixed arm angles theta1-4, builds the static allocation matrix
	 */
	void buildStaticMatrix();
};
