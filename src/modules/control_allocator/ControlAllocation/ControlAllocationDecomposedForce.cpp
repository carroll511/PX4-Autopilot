/****************************************************************************
 *
 *   Copyright (c) 2024 PX4 Development Team. All rights reserved.
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
 * @file ControlAllocationDecomposedForce.cpp
 *
 * Control Allocation Algorithm that decomposes rotor forces into vertical
 * and lateral components using a static allocation matrix.
 * Rotor orientation is then calculated using arctan.
 *
 * @author Custom Implementation
 */

#include "ControlAllocationDecomposedForce.hpp"
#include <mathlib/math/Functions.hpp>
#include <px4_platform_common/log.h>

void
ControlAllocationDecomposedForce::setEffectivenessMatrix(
	const matrix::Matrix<float, ControlAllocation::NUM_AXES, ControlAllocation::NUM_ACTUATORS> &effectiveness,
	const ActuatorVector &actuator_trim, const ActuatorVector &linearization_point, int num_actuators,
	bool update_normalization_scale)
{
	ControlAllocation::setEffectivenessMatrix(effectiveness, actuator_trim, linearization_point, num_actuators,
			update_normalization_scale);
	_mix_update_needed = true;
	_num_rotors = num_actuators;
}

void
ControlAllocationDecomposedForce::updateStaticAllocationMatrix()
{
	if (!_mix_update_needed) {
		return;
	}

	// The static allocation matrix maps control setpoints to vertical and lateral
	// force components. We use the effectiveness matrix as the base and compute
	// a pseudo-inverse for the static allocation.
	// The key difference is that we treat each rotor as having two components:
	// vertical and lateral forces.

	// For now, we'll use a simplified approach where we compute the pseudo-inverse
	// of the effectiveness matrix. In a full implementation, you might want to
	// restructure this to explicitly separate vertical and lateral components.
	matrix::geninv(_effectiveness, _static_mix);

	// Set all the small elements to 0 to avoid issues
	for (int i = 0; i < _num_actuators; i++) {
		for (int j = 0; j < NUM_AXES; j++) {
			if (fabsf(_static_mix(i, j)) < 1e-3f) {
				_static_mix(i, j) = 0.f;
			}
		}
	}

	_mix_update_needed = false;
}

float
ControlAllocationDecomposedForce::computeRotorOrientation(float vertical, float lateral) const
{
	// Compute orientation using arctan2 for proper quadrant handling
	// This gives the angle of the rotor force vector in the vertical-lateral plane
	if (fabsf(vertical) < FLT_EPSILON && fabsf(lateral) < FLT_EPSILON) {
		return 0.0f; // No force, no orientation
	}

	return atan2f(lateral, vertical);
}

float
ControlAllocationDecomposedForce::computeActuatorMagnitude(float vertical, float lateral) const
{
	// Compute the magnitude of the force vector
	return sqrtf(vertical * vertical + lateral * lateral);
}

void
ControlAllocationDecomposedForce::allocate()
{
	// Update static allocation matrix if needed
	updateStaticAllocationMatrix();

	_prev_actuator_sp = _actuator_sp;

	// Compute vertical and lateral force components for each rotor
	// using the static allocation matrix
	matrix::Vector<float, NUM_AXES> control_demand = _control_sp - _control_trim;

	// Compute actuator forces using the static mix matrix
	ActuatorVector actuator_forces = _static_mix * control_demand;

	// Decompose forces into vertical and lateral components for each rotor
	// The effectiveness matrix encodes how each rotor contributes to control.
	// We decompose this into vertical (along rotor axis) and lateral (perpendicular) components.
	// For a coaxial quadcopter with tilting arms, this decomposition is done in the rotor's
	// local frame (vertical = along rotor axis, lateral = perpendicular to axis).

	for (int i = 0; i < _num_actuators; i++) {
		// Compute the total force vector from this rotor's contributions
		// Vertical component: force along the rotor's axis (primarily THRUST_Z)
		// For a tilting rotor, the vertical component is the projection along the rotor axis
		float vertical = _effectiveness(THRUST_Z, i) * control_demand(THRUST_Z);

		// Lateral components: forces perpendicular to the rotor axis
		// These come from THRUST_X, THRUST_Y (when rotor is tilted)
		float lateral_x = _effectiveness(THRUST_X, i) * control_demand(THRUST_X);
		float lateral_y = _effectiveness(THRUST_Y, i) * control_demand(THRUST_Y);

		// Torque contributions: roll and pitch torques are generated by lateral forces
		// when the rotor is tilted. Yaw torque is from propeller rotation.
		float torque_roll = _effectiveness(ROLL, i) * control_demand(ROLL);
		float torque_pitch = _effectiveness(PITCH, i) * control_demand(PITCH);
		float torque_yaw = _effectiveness(YAW, i) * control_demand(YAW);

		// For tilting rotors, torque contributions indicate lateral force generation
		// The lateral force magnitude is computed from the torque contributions
		// and the lateral thrust components
		float lateral = sqrtf(lateral_x * lateral_x + lateral_y * lateral_y);

		// Add torque contributions as additional lateral force indicators
		// (torque is generated by lateral force offset from center of mass)
		float lateral_from_torque = sqrtf(torque_roll * torque_roll + torque_pitch * torque_pitch);
		lateral = sqrtf(lateral * lateral + lateral_from_torque * lateral_from_torque);

		_vertical_forces(i) = vertical;
		_lateral_forces(i) = lateral;
	}

	// Compute actuator setpoints from vertical and lateral components
	// The actuator setpoint is the magnitude of the force vector
	// Rotor orientation is computed as arctan(lateral/vertical)
	for (int i = 0; i < _num_actuators; i++) {
		float vertical = _vertical_forces(i);
		float lateral = _lateral_forces(i);

		// Compute magnitude (actuator setpoint)
		_actuator_sp(i) = _actuator_trim(i) + computeActuatorMagnitude(vertical, lateral);

		// Rotor orientation can be computed as:
		// float orientation = computeRotorOrientation(vertical, lateral);
		// This gives the angle of the rotor force vector relative to vertical
		// orientation = 0 means purely vertical, orientation = π/2 means purely lateral
	}
}
