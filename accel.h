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
	static ale_context _accel_context;

public:

	static void set_accel() {
		accel_type = 0;
	}

	static void set_gpu() {
		accel_type = 1;
	}

	static void set_cpu() {
		accel_type = 2;
	}

	static void set_auto() {
		/*
		 * Set preference ACCELERATOR > GPU > DEFAULT > CPU
		 */

		cl_context cc;

		cc = clCreateContextFromType(0, CL_DEVICE_TYPE_ACCELERATOR, NULL, NULL, NULL);

	        if (cc == ((cl_context) 0))
	                cc = clCreateContextFromType(0, CL_DEVICE_TYPE_GPU, NULL, NULL, NULL);

	        if (cc == ((cl_context) 0))
	                cc = clCreateContextFromType(0, CL_DEVICE_TYPE_DEFAULT, NULL, NULL, NULL);

	        if (cc == ((cl_context) 0))
	                cc = clCreateContextFromType(0, CL_DEVICE_TYPE_CPU, NULL, NULL, NULL);

	        if (cc == ((cl_context) 0)) {
			fprintf(stderr, "Could not create an OpenCL context.\n");
			exit(1);
		}

		_accel_context = ale_new_context(cc);
		clReleaseContext(cc);

		if (!_accel_context) {
			fprintf(stderr, "Could not create an ALE context.\n");
			exit(1);
		}
	}

	static ale_context context() {
		if (_accel_context)
			return _accel_context;

		if (accel_type == -1) {
			set_auto();
			return _accel_context;
		}

		cl_context cc;

		if (accel_type == 0) {
			cc = clCreateContextFromType(0, CL_DEVICE_TYPE_ACCELERATOR, NULL, NULL, NULL);
		} else if (accel_type == 1) {
			cc = clCreateContextFromType(0, CL_DEVICE_TYPE_GPU, NULL, NULL, NULL);
		} else {
			cc = clCreateContextFromType(0, CL_DEVICE_TYPE_CPU, NULL, NULL, NULL);
		}

	        if (cc == ((cl_context) 0)) {
			fprintf(stderr, "Could not create an OpenCL context.\n");
			exit(1);
		}

		_accel_context = ale_new_context(cc);
		clReleaseContext(cc);

		if (!_accel_context) {
			fprintf(stderr, "Could not create an ALE context.\n");
			exit(1);
		}

		return _accel_context;
	}
};

#endif
