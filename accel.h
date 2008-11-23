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
	static int accel_type;
	static void *_accel_context;
	static void *_accel_backend;

public:

	static void set_gpu() {
		accel_type = 1;
	}

	static void set_cpu() {
		accel_type = 2;
	}

	static void set_none() {
		accel_type = 0;
	}

	static void set_auto() {
		const char *accel_default = getenv("ALE_ACCEL_DEFAULT");
		int attempted_accel_type = 0;

		if (!accel_default) {
			accel_type = 0;
			return;
		}

		if (!strcmp(accel_default, "none")) {
			accel_type = 0;
			return;
		} else if (!strcmp(accel_default, "gpu")) {
			_accel_backend = ale_new_backend("glsl");
			attempted_accel_type = 1;
		} else if (!strcmp(accel_default, "cpu")) {
			_accel_backend = ale_new_backend("dynlib");
			attempted_accel_type = 2;
		}

		if (!_accel_backend) {
			accel_type = 0;
			return;
		}

		_accel_context = ale_new_context(_accel_backend);

		if (!_accel_context) {
			ale_delete_backend(_accel_backend);
			_accel_backend = NULL;
			accel_type = 0;
			return;
		}

		accel_type = attempted_accel_type;
	}

	static void *context() {
		if (accel_type == -1)
			set_auto();

		if (_accel_context)
			return _accel_context;

		if (accel_type == 0) {
			_accel_backend = ale_new_backend("guile");
		} else if (accel_type == 1) {
			_accel_backend = ale_new_backend("glsl");
		} else {
			_accel_backend = ale_new_backend("dynlib");
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
