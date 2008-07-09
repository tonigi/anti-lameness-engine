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

#include "gpu.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

class accel {
	static int use_gpu;
	static int mask_gpu;

public:

	static void mask_gpu() {
		mask_gpu = 1;
	}

	static void unmask_gpu() {
		mask_gpu = 0;
	}

	static void set_gpu() {
		use_gpu = 1;
	}

	static void set_none() {
		use_gpu = 0;
	}

	static void set_auto() {
#ifdef USE_GLEW
		const char *accel_default = getenv("ALE_GPU_ACCEL_DEFAULT");
		if (accel_default && !strcmp(accel_default, "1") && gpu::is_ok())
			use_gpu = 1;
		else 
			use_gpu = 0;
#else
		use_gpu = 0;
#endif
	}

	static int is_gpu() {
		if (use_gpu == 2)
			set_auto();

		if (use_gpu == 0)
			return 0;

		if (mask_gpu)
			return 0;

		if (!gpu::is_ok()) {
			fprintf(stderr, "GPU acceleration error.\n");
			exit(1);
		}

#ifdef USE_GLEW
		const char *extensions[] = {
			"GL_ARB_fragment_shader",
			"GL_ARB_texture_float",
			"GL_ARB_texture_rectangle",
			"GL_EXT_framebuffer_object",
			NULL
		};

		int unsupported = 0;
		for (const char **c = extensions; *c; c++) {
			if (!glewIsSupported(*c)) {
				fprintf(stderr, "GL feature %s not supported.\n", *c);
				unsupported = 1;
			}
		}
		if (unsupported)
			exit(1);
#else
		fprintf(stderr, "GLEW linkage is required for acceleration.\n");
		exit(1);
#endif

		return 1;
	}
};

#endif
