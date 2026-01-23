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
 * @file ControlAllocationDecomposedForce.hpp
 *
 * Control Allocation Algorithm that decomposes rotor forces into vertical
 * and lateral components using a static allocation matrix.
 * Rotor orientation is then calculated using arctan.
 *
 * @author Custom Implementation
 */

#pragma once

#include "ControlAllocation.hpp"
#include <matrix/matrix/math.hpp>

class ControlAllocationDecomposedForce: public ControlAllocation
{
public:
	ControlAllocationDecomposedForce() = default;
	virtual ~ControlAllocationDecomposedForce() = default;

	void allocate() override;
	void setEffectivenessMatrix(const matrix::Matrix<float, NUM_AXES, NUM_ACTUATORS> &effectiveness,
				    const ActuatorVector &actuator_trim, const ActuatorVector &linearization_point, int num_actuators,
				    bool update_normalization_scale) override;

protected:
	/**
	 * Static allocation matrix that maps control setpoints to vertical and lateral
	 * force components. This matrix doesn't depend on rotor orientation.
	 * For N rotors, the matrix has 2N columns (vertical, lateral for each rotor).
	 */
	matrix::Matrix<float, NUM_ACTUATORS, NUM_AXES> _static_mix;

	/**
	 * Number of rotors (actuators / 2, since each rotor has vertical and lateral components)
	 */
	int _num_rotors{0};

	/**
	 * Vertical force components for each rotor
	 */
	matrix::Vector<float, NUM_ACTUATORS> _vertical_forces;

	/**
	 * Lateral force components for each rotor
	 */
	matrix::Vector<float, NUM_ACTUATORS> _lateral_forces;

	/**
	 * Update the static allocation matrix from the effectiveness matrix
	 */
	void updateStaticAllocationMatrix();

	/**
	 * Compute rotor orientation from vertical and lateral force components
	 * @param vertical Vertical force component
	 * @param lateral Lateral force component
	 * @return Orientation angle in radians (arctan(lateral/vertical))
	 */
	float computeRotorOrientation(float vertical, float lateral) const;

	/**
	 * Compute actuator setpoint (magnitude) from vertical and lateral components
	 * @param vertical Vertical force component
	 * @param lateral Lateral force component
	 * @return Actuator setpoint magnitude
	 */
	float computeActuatorMagnitude(float vertical, float lateral) const;

	bool _mix_update_needed{false};
};
