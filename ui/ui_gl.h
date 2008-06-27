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

#ifndef __ui_gl_h__
#define __ui_gl_h__

#ifdef USE_FREEGLUT
#include <GL/freeglut.h>
#elif USE_GLUT
#include <GL/glut.h>
#endif

#include "ui.h"
#include "gpu.h"

/*
 * OpenGL user interface.
 */

class ui_gl : public ui {
private:
	void printf(const char *format, ...) {
		va_list ap;
		va_start(ap, format);
		vfprintf(ui_stream, format, ap);
		va_end(ap);
	}

	void update() {
#ifdef HAVE_GLUT_MAIN_LOOP_EVENT
		gpu::lock();
		glutMainLoopEvent();
		gpu::unlock();
#endif
	}

public:
	ui_gl() {
#ifndef HAVE_GLUT_MAIN_LOOP_EVENT
		const char *freeglut_error = "need extension glutMainLoopEvent()";
		throw freeglut_error;
#endif

		if (!gpu::is_ok()) {
			const char *gpu_init_error = "GPU initialization error";
			throw gpu_init_error;
		}
	}
};

#endif
