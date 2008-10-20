// Copyright 2008 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@gmail.com>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 * accel.h: acceleration
 */

#ifndef __accel_h__
#define __accel_h__

extern "C" {
#include <ale.h>
}

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

class accel {
	static int use_gpu;
	static int _mask_gpu;
	static void *_accel_context;
	static void *_accel_backend;

public:

	static void mask_gpu() {
		_mask_gpu = 1;
	}

	static void unmask_gpu() {
		_mask_gpu = 0;
	}

	static void set_gpu() {
		use_gpu = 1;
	}

	static void set_none() {
		use_gpu = 0;
	}

	static void set_auto() {
		const char *accel_default = getenv("ALE_GPU_ACCEL_DEFAULT");
		if (!accel_default || strcmp(accel_default, "1")) {
			use_gpu = 0;
			return;
		}

		_accel_backend = ale_new_backend("glsl");

		if (!_accel_backend) {
			use_gpu = 0;
			return;
		}

		_accel_context = ale_new_context(_accel_backend);

		if (!_accel_context) {
			ale_delete_backend(_accel_backend);
			_accel_backend = NULL;
			use_gpu = 0;
			return;
		}

		use_gpu = 1;
	}

	static int is_gpu() {
		if (use_gpu == 2)
			set_auto();

		if (use_gpu == 0)
			return 0;

		if (_mask_gpu)
			return 0;

		return 1;
	}

	static void *context() {
		if (use_gpu == 2)
			set_auto();

		if (_accel_context)
			return _accel_context;

		if (use_gpu == 0) {
			_accel_backend = ale_new_backend("guile");
		} else {
			_accel_backend = ale_new_backend("glsl");
		}

		if (!_accel_backend) {
			fprintf(stderr, "Could not create backend.\n");
			exit(1);
		}

		_accel_context = ale_new_context(_accel_backend);

		if (!_accel_context) {
			fprintf(stderr, "Could not create context.\n");
			exit(1);
		}

		return _accel_context;

	}

};

#endif
