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
#include <matrix/matrix/math.hpp>
#include <iostream>
#include <cmath>

using namespace matrix;

void updateControlAllocationMatrixScale(Matrix<float, 4, 6>& mix, Vector<float, 6>& control_allocation_scale, int num_actuators, bool normalize_rpy)
{
	if (normalize_rpy) {
		int num_non_zero_roll_torque = 0;
		int num_non_zero_pitch_torque = 0;

		for (int i = 0; i < num_actuators; i++) {
			if (fabsf(mix(i, 0)) > 1e-3f) {
				++num_non_zero_roll_torque;
			}
			if (fabsf(mix(i, 1)) > 1e-3f) {
				++num_non_zero_pitch_torque;
			}
		}

		float roll_norm_scale = 1.f;
		if (num_non_zero_roll_torque > 0) {
			roll_norm_scale = sqrtf(mix.col(0).norm_squared() / (num_non_zero_roll_torque / 2.f));
		}

		float pitch_norm_scale = 1.f;
		if (num_non_zero_pitch_torque > 0) {
			pitch_norm_scale = sqrtf(mix.col(1).norm_squared() / (num_non_zero_pitch_torque / 2.f));
		}

		control_allocation_scale(0) = fmaxf(roll_norm_scale, pitch_norm_scale);
		control_allocation_scale(1) = control_allocation_scale(0);

		// Scale yaw separately
		control_allocation_scale(2) = mix.col(2).max();
	} else {
		control_allocation_scale(0) = 1.f;
		control_allocation_scale(1) = 1.f;
		control_allocation_scale(2) = 1.f;
	}

	for (int axis_idx = 2; axis_idx >= 0; --axis_idx) {
		int num_non_zero_thrust = 0;
		float norm_sum = 0.f;

		for (int i = 0; i < num_actuators; i++) {
		float norm = fabsf(mix(i, 3 + axis_idx));
		norm_sum += norm;

			if (norm > FLT_EPSILON) {
				++num_non_zero_thrust;
			}
		}

		if (num_non_zero_thrust > 0) {
		control_allocation_scale(3 + axis_idx) = norm_sum / num_non_zero_thrust;
		} else {
		control_allocation_scale(3 + axis_idx) = control_allocation_scale(3);
		}
	}
}

void normalizeControlAllocationMatrix(Matrix<float, 4, 6>& mix, Vector<float, 6>& control_allocation_scale)
{
	if (control_allocation_scale(0) > FLT_EPSILON) {
		mix.col(0) /= control_allocation_scale(0);
		mix.col(1) /= control_allocation_scale(1);
	}
	if (control_allocation_scale(2) > FLT_EPSILON) {
		mix.col(2) /= control_allocation_scale(2);
	}
	if (control_allocation_scale(3) > FLT_EPSILON) {
		mix.col(3) /= control_allocation_scale(3);
		mix.col(4) /= control_allocation_scale(4);
		mix.col(5) /= control_allocation_scale(5);
	}

	// Set all the small elements to 0
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 6; j++) {
			if (fabsf(mix(i, j)) < 1e-3f) {
				mix(i, j) = 0.f;
			}
		}
	}
}

bool areMatricesEqual(const Matrix<float, 4, 6>& mat1, const Matrix<float, 4, 6>& mat2, float tolerance)
{
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            if (fabsf(mat1(i, j) - mat2(i, j)) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

TEST(ControlAllocationTest, PseudoInverseAndNormalization)
{
	constexpr size_t num_controls = 6;
	constexpr size_t num_actuators = 4;

	matrix::Matrix<float, num_controls, num_actuators> effectiveness;
	effectiveness(0, 0) = -1.0073246f; effectiveness(0, 1) = 1.0073246f; effectiveness(0, 2) = 1.0073246f;  effectiveness(0, 3) = -1.0073246f;
	effectiveness(1, 0) = 1.6425133f;  effectiveness(1, 1) = -2.29491f;  effectiveness(1, 2) = 1.6425133f;  effectiveness(1, 3) = -2.29491f;
	effectiveness(2, 0) = 0.53718674;  effectiveness(2, 1) = 0.53718674; effectiveness(2, 2) = -0.53718674; effectiveness(2, 3) = -0.53718674;
	effectiveness(3, 0) = -2.5382986;  effectiveness(3, 1) = 2.5382986;  effectiveness(3, 2) = -2.5382986;  effectiveness(3, 3) = 2.5382986;
	effectiveness(4, 0) = 0.0f;        effectiveness(4, 1) = 0.0f;       effectiveness(4, 2) = 0.0f;        effectiveness(4, 3) = 0.0f;
	effectiveness(5, 0) = -5.4366393;  effectiveness(5, 1) = -5.4366393; effectiveness(5, 2) = -5.4366393;  effectiveness(5, 3) = -5.4366393;

	matrix::Matrix<float, num_actuators, num_controls> ex_pseudo_inverse;
	ex_pseudo_inverse(0, 0) = -0.24818216f; ex_pseudo_inverse(0, 1) = 0.04587143f;  ex_pseudo_inverse(0, 2) = 0.46538751f;  ex_pseudo_inverse(0, 3) = -0.06291316f; ex_pseudo_inverse(0, 4) = 0.0f; ex_pseudo_inverse(0, 5) = -0.04873658f;
	ex_pseudo_inverse(1, 0) = 0.24818216f;  ex_pseudo_inverse(1, 1) = -0.04930917f; ex_pseudo_inverse(1, 2) = 0.46538751f;  ex_pseudo_inverse(1, 3) = 0.06024684f;  ex_pseudo_inverse(1, 4) = 0.0f; ex_pseudo_inverse(1, 5) = -0.04302574;
	ex_pseudo_inverse(2, 0) = 0.24818216f;  ex_pseudo_inverse(2, 1) = 0.04587143f;  ex_pseudo_inverse(2, 2) = -0.46538751f; ex_pseudo_inverse(2, 3) = -0.06291316f; ex_pseudo_inverse(2, 4) = 0.0f; ex_pseudo_inverse(2, 5) = -0.04873658f;
	ex_pseudo_inverse(3, 0) = -0.24818216f; ex_pseudo_inverse(3, 1) = -0.04930917f; ex_pseudo_inverse(3, 2) = -0.46538751f; ex_pseudo_inverse(3, 3) = 0.06024684f;  ex_pseudo_inverse(3, 4) = 0.0f; ex_pseudo_inverse(3, 5) = -0.04302574;

	// Step 1: Validation of pseudo-inverse computation
	matrix::Matrix<float, num_actuators, num_controls> pseudo_inverse;
	bool success = geninv(effectiveness, pseudo_inverse);

	ASSERT_TRUE(success) << "Psuedo-inverse computation failed";

	std::cout << "ex_Pesudo-inverse: " << ex_pseudo_inverse << std::endl;
	for (size_t i = 0; i < num_actuators; ++i) {
		for (size_t j = 0; j < num_controls; ++j) {
			std::cout << ex_pseudo_inverse(i, j) << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "Pesudo-inverse: " << pseudo_inverse << std::endl;
	for (size_t i = 0; i < num_actuators; ++i) {
		for (size_t j = 0; j < num_controls; ++j) {
			std::cout << pseudo_inverse(i, j) << " ";
		}
		std::cout << std::endl;
	}

	// Compare computed pseudo-inverse with expected pseudo-inverse
	EXPECT_TRUE(areMatricesEqual(pseudo_inverse, ex_pseudo_inverse, 1e-5))
        	<< "Computed pseudo-inverse does not match the expected pseudo-inverse.";

//	// Verify the Moore-Penrose condition G * G^+ * G ≈ G
//	matrix::Matrix<float, num_controls, num_controls> identity_approximation =
//        	effectiveness * pseudo_inverse * effectiveness;
//
//	EXPECT_TRUE(areMatricesEqual(identity_approximation, effectiveness, 1e-5))
//        	<< "Pseudo-inverse does not satisfy the Moore-Penrose condition.";


	// Output the difference between the two pseudo-inverse matrices for debugging
	std::cout << "Difference between computed and expected pseudo-inverse:\n"
		<< (pseudo_inverse - ex_pseudo_inverse) << std::endl;

	// Step 2: Validation of Mixer computation (Normaliztion)
	matrix::Vector<float, num_controls> control_allocation_scale;
	control_allocation_scale.setAll(1.0f);

	updateControlAllocationMatrixScale(pseudo_inverse, control_allocation_scale, num_actuators, true);
	normalizeControlAllocationMatrix(pseudo_inverse, control_allocation_scale);

	// Output normalized pseudo-inverse and scale
	std::cout << "Normalized Mixer Matrix:\n" << pseudo_inverse << std::endl;
	std::cout << "Control Allocation Scale:\n" << control_allocation_scale << std::endl;

	// Verify normalization
	for (size_t i = 0; i < num_controls; ++i) {
		float column_norm = pseudo_inverse.col(i).norm();
		EXPECT_NEAR(column_norm, 1.0f, 1e-5) << "Column " << i << " is not properly normalized";
	}
}
