/****************************************************************************
 *
 *   Copyright (C) 2019 PX4 Development Team. All rights reserved.
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
 * @file ControlAllocationTest.cpp
 *
 * Tests for Control Allocation Algorithms
 *
 * @author Julien Lecoeur <julien.lecoeur@gmail.com>
 */

#include <gtest/gtest.h>
#include <ControlAllocationPseudoInverse.hpp>
#include <ActuatorEffectiveness/ActuatorEffectivenessRotors.hpp>

using namespace matrix;

TEST(ControlAllocationTest, AllZeroCase)
{
	ControlAllocationPseudoInverse method;

	matrix::Vector<float, 6> control_sp;
	matrix::Vector<float, 6> control_allocated;
	matrix::Vector<float, 6> control_allocated_expected;
	matrix::Matrix<float, 6, 16> effectiveness;
	matrix::Vector<float, 16> actuator_sp;
	matrix::Vector<float, 16> actuator_trim;
	matrix::Vector<float, 16> linearization_point;
	matrix::Vector<float, 16> actuator_sp_expected;

	method.setEffectivenessMatrix(effectiveness, actuator_trim, linearization_point, 16, false);
	method.setControlSetpoint(control_sp);
	method.allocate();
	method.clipActuatorSetpoint();
	actuator_sp = method.getActuatorSetpoint();
	control_allocated_expected = method.getAllocatedControl();

	EXPECT_EQ(actuator_sp, actuator_sp_expected);
	EXPECT_EQ(control_allocated, control_allocated_expected);
}

TEST(ControlAllocationTest, CoaxialGeometryEffectiveness)
{
	ActuatorEffectivenessRotors::Geometry g{};
	const float L_x = 0.105f;                  // x arm length
	const float L_y = 0.105f;                  // y arm length
	const float ct = 1.0f;                     // thrust coef
	const float km = 0.05f;                    // torque ratio
	const float eta = 0.0f; // tilt angle in rad
	const Vector3f axis{-sinf(eta), 0.f, -cosf(eta)};
	const Vector3f axis_normalized = axis / axis.norm();

	const Vector3f positions[8] = {
		{ L_x,  L_y, 0.f},
		{ L_x, -L_y, 0.f},
		{-L_x, -L_y, 0.f},
		{-L_x,  L_y, 0.f},
		{ L_x, -L_y, 0.f},
		{ L_x,  L_y, 0.f},
		{-L_x,  L_y, 0.f},
		{-L_x, -L_y, 0.f},
	};

	// Alternate moment_ratio sign to match CW/CCW prop pairs
	for (int i = 0; i < 8; ++i) {
		const float moment_ratio = (i % 2 == 0) ? km : -km;
		g.rotors[i] = {positions[i], axis, ct, moment_ratio, -1};
	}

	g.num_rotors = 8;

	matrix::Matrix<float, 6, 16> eff{};
	int n = 0;

	// Build effectiveness matrix the same way ActuatorEffectivenessRotors does
	for (int i = 0; i < g.num_rotors; ++i) {
		const float moment_ratio = (i % 2 == 0) ? km : -km;
		const Vector3f thrust = ct * axis_normalized;
		const Vector3f moment = ct * positions[i].cross(axis_normalized) - ct * moment_ratio * axis_normalized;

		for (int j = 0; j < 3; ++j) {
			eff(j, i) = moment(j);
			eff(j + 3, i) = thrust(j);
		}

		++n;
	}

	ControlAllocationPseudoInverse alloc;
	Vector<float, 16> trim{}, lin{};
	alloc.setEffectivenessMatrix(eff, trim, lin, n, false);

	// Show the generated effectiveness matrix (matching CA_ROTOR* geometry above)
	printf("Effectiveness matrix (rows = actuators, cols = axes):\n");
	eff.print();

	ASSERT_EQ(n, 8);

	// Verify the effectiveness matrix matches the configured geometry
	for (int i = 0; i < n; ++i) {
		const float moment_ratio = (i % 2 == 0) ? km : -km;
		const Vector3f expected_thrust = ct * axis_normalized;
		const Vector3f expected_moment = ct * positions[i].cross(axis_normalized) - ct * moment_ratio * axis_normalized;

		for (int j = 0; j < 3; ++j) {
			EXPECT_EQ(eff(j, i), expected_moment(j));
			EXPECT_EQ(eff(j + 3, i), expected_thrust(j));
		}
	}
}
