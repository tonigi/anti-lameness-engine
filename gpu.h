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
 * gpu.h: Graphics processor management
 */

#ifndef __gpu_h__
#define __gpu_h__

#include <config.h>

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef USE_GLEW
#include <GL/glew.h>
#endif

#ifdef USE_FREEGLUT
#include <GL/freeglut.h>
#elif USE_GLUT
#include <GL/glut.h>
#endif

#include "thread.h"


class gpu {
	static int gpu_initialized;
	static thread::lock_t _gpu_lock;

	static void try_init_gpu() {
		assert(!gpu_initialized);

#ifdef USE_GLUT
		char program_name[20];
		int fake_argc = 1;
		char *fake_argv[] = { program_name };

		strcpy(program_name, "ale");

		glutInit(&fake_argc, fake_argv);
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE); 
		glutCreateWindow(PACKAGE_NAME " " VERSION);
		gpu_initialized = 1;
#else
		gpu_initialized = 0;
#endif

#ifdef USE_GLEW
		if (gpu_initialized) {
			if (glewInit() != GLEW_OK) {
				fprintf(stderr, "glewInit() returned an error.\n");
				exit(1);
			}
		}
#endif
	}

public:

	static int is_ok() {
		if (!gpu_initialized) {
			try_init_gpu();
		}

		return gpu_initialized;
	}

	static void lock() {
		_gpu_lock.lock();
	}

	static void unlock() {
		_gpu_lock.unlock();
	}
};

#endif
