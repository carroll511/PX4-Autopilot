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
 * @file ControlAllocationPseudoInverse.hpp
 *
 * Simple Control Allocation Algorithm
 *
 * @author Julien Lecoeur <julien.lecoeur@gmail.com>
 */

#include "ControlAllocationPseudoInverse.hpp"
#include <cstring>
#include <px4_log.h>
#include <uORB/uORB.h>
#include <cmath>  // std::fabs 사용
#include <drivers/drv_hrt.h> // For hrt_absolute_time()
// #include <iostream>

const float EPSILON = 1e-6; // 작은 차이를 허용하는 epsilon 값 정의
uORB::Publication <debug_value_s> _CAPI_debug_value_pub{ORB_ID(debug_value)};
debug_value_s CAPI_debug_msg = {};

void ControlAllocationPseudoInverse::initializeRCInput()
{
	// CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
	// CAPI_debug_msg.ind = 8;                         // Set the index (for example)
	// CAPI_debug_msg.value = 8.0f;                // Set the value to send as the debug output
	// _CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish the message

	_rc_sub = orb_subscribe(ORB_ID(input_rc));
}

// 변경할 파라미터 이름 배열
const char *param_names[12] = {
	"CA_ROTOR0_PX", "CA_ROTOR0_PY", "CA_ROTOR0_PZ",
	"CA_ROTOR1_PX", "CA_ROTOR1_PY", "CA_ROTOR1_PZ",
	"CA_ROTOR2_PX", "CA_ROTOR2_PY", "CA_ROTOR2_PZ",
	"CA_ROTOR3_PX", "CA_ROTOR3_PY", "CA_ROTOR3_PZ"
};

// 새로 설정할 값
const float values_30[12] = {
    0.06f , 0.19f , -0.47f, // ROTOR 1
    -0.16f, -0.19f, -0.47f, // ROTOR 2
    0.06f , -0.19f, -0.47f, // ROTOR 3
    -0.16f, 0.19f , -0.47f // ROTOR 4
};

const float values_60[12] = {
    -0.16f, 0.19f , -0.33f, // ROTOR 1
    -0.39f, -0.19f, -0.33f, // ROTOR 2
    -0.16f, -0.19f, -0.33f, // ROTOR 3
    -0.39f, 0.19f , -0.33f // ROTOR 4
};

const float values_90[12] = {
    -0.23f, 0.19f , -0.07f, // ROTOR 1
    -0.46f, -0.19f, -0.07f, // ROTOR 2
    -0.23f, -0.19f, -0.07f, // ROTOR 3
    -0.46f, 0.19f , -0.07f // ROTOR 4
};

float new_values[12] = {};

bool ControlAllocationPseudoInverse::updateRCInput()
{
	CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
	CAPI_debug_msg.ind = 7;                         // Set the index (for example)
	CAPI_debug_msg.value = 7.0f;                   // Set the value to send as the debug output
	_CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish the message

    if (_rc_sub < 0) {
        PX4_WARN("RC subscription not initialized!");
    }

    bool updated = false;
    orb_check(_rc_sub, &updated); // 업데이트 확인

    if (updated) {
        orb_copy(ORB_ID(input_rc), _rc_sub, &_rc_data); // 데이터 복사
    }


	if (std::abs(static_cast<int>(_rc_data.values[4])-_prev_rc) > EPSILON) {
		_prev_rc = _rc_data.values[4]; // 이전 값 저장
		return true;
	}
	else {
		_prev_rc = _rc_data.values[4]; // 이전 값 저장
		return false;
	}
}

void ControlAllocationPseudoInverse::enable_mix_update()
{
	// CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
	// CAPI_debug_msg.ind = 6;                         // Set the index (for example)
	// CAPI_debug_msg.value = 6.0f;                // Set the value to send as the debug output
	// _CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish the message

	_mix_update_needed = true;
}

void
ControlAllocationPseudoInverse::setEffectivenessMatrix(
	const matrix::Matrix<float, ControlAllocation::NUM_AXES, ControlAllocation::NUM_ACTUATORS> &effectiveness,
	const ActuatorVector &actuator_trim, const ActuatorVector &linearization_point, int num_actuators,
	bool update_normalization_scale)
{
	// CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
	// CAPI_debug_msg.ind = 5;                         // Set the index (for example)
	// CAPI_debug_msg.value = 5.0f;                // Set the value to send as the debug output
	// _CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish the message

	ControlAllocation::setEffectivenessMatrix(effectiveness, actuator_trim, linearization_point, num_actuators,
			update_normalization_scale);
	_mix_update_needed = true;
	_normalization_needs_update = update_normalization_scale;
}

void
ControlAllocationPseudoInverse::updatePseudoInverse()
{
	// CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
	// CAPI_debug_msg.ind = 4;                         // Set the index (for example)
	// CAPI_debug_msg.value = 4.0f;                // Set the value to send as the debug output
	// _CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish the message

	if (_mix_update_needed) {
		matrix::geninv(_effectiveness, _mix);

		// CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
		// CAPI_debug_msg.ind = 9;                         // Set the index (for example)
		// CAPI_debug_msg.value = 9.0f;                // Set the value to send as the debug output
		// _CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish the message

		if (_normalization_needs_update && !_had_actuator_failure) {
			updateControlAllocationMatrixScale();
			_normalization_needs_update = false;
		}

		normalizeControlAllocationMatrix();
		_mix_update_needed = false;
	}
}

void
ControlAllocationPseudoInverse::updateControlAllocationMatrixScale()
{
	// CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
	// CAPI_debug_msg.ind = 3;                         // Set the index (for example)
	// CAPI_debug_msg.value = 3.0f;                // Set the value to send as the debug output
	// _CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish the message

	// Same scale on roll and pitch
	if (_normalize_rpy) {

		int num_non_zero_roll_torque = 0;
		int num_non_zero_pitch_torque = 0;

		for (int i = 0; i < _num_actuators; i++) {

			if (fabsf(_mix(i, 0)) > 1e-3f) {
				++num_non_zero_roll_torque;
			}

			if (fabsf(_mix(i, 1)) > 1e-3f) {
				++num_non_zero_pitch_torque;
			}
		}

		float roll_norm_scale = 1.f;

		if (num_non_zero_roll_torque > 0) {
			roll_norm_scale = sqrtf(_mix.col(0).norm_squared() / (num_non_zero_roll_torque / 2.f));
		}

		float pitch_norm_scale = 1.f;

		if (num_non_zero_pitch_torque > 0) {
			pitch_norm_scale = sqrtf(_mix.col(1).norm_squared() / (num_non_zero_pitch_torque / 2.f));
		}

		_control_allocation_scale(0) = fmaxf(roll_norm_scale, pitch_norm_scale);
		_control_allocation_scale(1) = _control_allocation_scale(0);

		// Scale yaw separately
		_control_allocation_scale(2) = _mix.col(2).max();

	} else {
		_control_allocation_scale(0) = 1.f;
		_control_allocation_scale(1) = 1.f;
		_control_allocation_scale(2) = 1.f;
	}

	// Scale thrust by the sum of the individual thrust axes, and use the scaling for the Z axis if there's no actuators
	// (for tilted actuators)
	_control_allocation_scale(THRUST_Z) = 1.f;

	for (int axis_idx = 2; axis_idx >= 0; --axis_idx) {
		int num_non_zero_thrust = 0;
		float norm_sum = 0.f;

		for (int i = 0; i < _num_actuators; i++) {
			float norm = fabsf(_mix(i, 3 + axis_idx));
			norm_sum += norm;

			if (norm > FLT_EPSILON) {
				++num_non_zero_thrust;
			}
		}

		if (num_non_zero_thrust > 0) {
			_control_allocation_scale(3 + axis_idx) = norm_sum / num_non_zero_thrust;

		} else {
			_control_allocation_scale(3 + axis_idx) = _control_allocation_scale(THRUST_Z);
		}
	}
}

void
ControlAllocationPseudoInverse::normalizeControlAllocationMatrix()
{
	// CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
	// CAPI_debug_msg.ind = 2;                         // Set the index (for example)
	// CAPI_debug_msg.value = 2.0f;                // Set the value to send as the debug output
	// _CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish the message

	if (_control_allocation_scale(0) > FLT_EPSILON) {
		_mix.col(0) /= _control_allocation_scale(0);
		_mix.col(1) /= _control_allocation_scale(1);
	}

	if (_control_allocation_scale(2) > FLT_EPSILON) {
		_mix.col(2) /= _control_allocation_scale(2);
	}

	if (_control_allocation_scale(3) > FLT_EPSILON) {
		_mix.col(3) /= _control_allocation_scale(3);
		_mix.col(4) /= _control_allocation_scale(4);
		_mix.col(5) /= _control_allocation_scale(5);
	}

	// Set all the small elements to 0 to avoid issues
	// in the control allocation algorithms
	for (int i = 0; i < _num_actuators; i++) {
		for (int j = 0; j < NUM_AXES; j++) {
			if (fabsf(_mix(i, j)) < 1e-3f) {
				_mix(i, j) = 0.f;
			}
		}
	}
}

void
ControlAllocationPseudoInverse::allocate()
{
	CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
	CAPI_debug_msg.ind = 88;                         // Set the index (for example)
	CAPI_debug_msg.value = 88.0f;                // Set the value to send as the debug output
	_CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish the message

	// RC 입력 업데이트
	initializeRCInput();
	updateRCInput();
	// CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
	// CAPI_debug_msg.ind = 1;                         // Set the index (for example)
	// CAPI_debug_msg.value = 1.0f;                // Set the value to send as the debug output
	// _CAPI_debug_value_pub.publish(debug_msg);       // Publish the message

	// RC 채널 5 pwm 기준으로 조건문 추가
	// if (_rc_data.channel_count > 4 && std::abs(static_cast<int>(_rc_data.values[4])-_prev_rc) > EPSILON) {
		if(_rc_data.values[4] < 1200) {
		        std::memcpy(new_values, values_30, sizeof(values_30));
			CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
			CAPI_debug_msg.ind = static_cast<int8_t>(_rc_data.values[4]);                         // Set the index (for example)
			CAPI_debug_msg.value = 10.0f;                // Set the value to send as the debug output
			_CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish the message

		}

		else if(_rc_data.values[4] >= 1200 && _rc_data.values[4] < 1700) {
		        std::memcpy(new_values, values_60, sizeof(values_60));
			CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
			CAPI_debug_msg.ind = static_cast<int8_t>(_rc_data.values[4]);                         // Set the index (for example)
			CAPI_debug_msg.value = 11.0f;                // Set the value to send as the debug output
			_CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish
		}

		else if(_rc_data.values[4] >= 1700) {
		        std::memcpy(new_values, values_90, sizeof(values_90));
			CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
			CAPI_debug_msg.ind = static_cast<int8_t>(_rc_data.values[4]);                         // Set the index (for example)
			CAPI_debug_msg.value = 12.0f;                // Set the value to send as the debug output
			_CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish
		}

		PX4_INFO("RC channel 5 triggered geometry parameter update");

		enable_mix_update();

		for (int i = 0; i < 12; ++i) {
			// 파라미터 핸들 가져오기
			param_t param_handle = param_find(param_names[i]);

			if (param_handle != PARAM_INVALID) {
				// 파라미터 값 설정
				int result = param_set(param_handle, &new_values[i]);

				if (result == PX4_OK) {
					PX4_INFO("Parameter %s updated to %.2f", param_names[i], double(new_values[i]));
					// CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
					// CAPI_debug_msg.ind = 99;                         // Set the index (for example)
					// CAPI_debug_msg.value = 99.0f;                // Set the value to send as the debug output
					// _CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish
				}
				else {
					PX4_WARN("Failed to update parameter %s", param_names[i]);
					CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
					CAPI_debug_msg.ind = 22;                         // Set the index (for example)
					CAPI_debug_msg.value = 22.0f;                // Set the value to send as the debug output
					_CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish
				}
			}
			else {
				PX4_WARN("Parameter %s not found", param_names[i]);
				CAPI_debug_msg.timestamp = hrt_absolute_time();  // Get the current timestamp in microseconds
				CAPI_debug_msg.ind = 33;                         // Set the index (for example)
				CAPI_debug_msg.value = 33.0f;                // Set the value to send as the debug output
				_CAPI_debug_value_pub.publish(CAPI_debug_msg);       // Publish
			}
		}
	// }
	// _prev_rc = _rc_data.values[4];

	// 기존 로직 호출
	//Compute new gains if needed
	updatePseudoInverse();

	_prev_actuator_sp = _actuator_sp;

	// Allocate
	_actuator_sp = _actuator_trim + _mix * (_control_sp - _control_trim);
}
