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

#ifdef USE_PTHREAD
	pthread_t update_thread;
#endif

	void printf(const char *format, ...) {
		va_list ap;
		va_start(ap, format);
		vfprintf(ui_stream, format, ap);
		va_end(ap);
	}

	void update() {
		/*
		 * Do nothing.  We only update via a thread loop.
		 */
	}

	static void *thread_update_loop(void *vu) {
#if defined HAVE_GLUT_MAIN_LOOP_EVENT && defined HAVE_NANOSLEEP
		for(;;) {
			struct timespec t;
			t.tv_sec = 0;
			t.tv_nsec = 10000000;
			nanosleep(&t, NULL);

			gpu::lock();
			glutMainLoopEvent();
			gpu::unlock();
		}
#endif

		return NULL;
	}
	
	/*
	 * GLUT callback logic.
	 */

	static ui_gl *instance;

	/*
	 * Reshape callback.
	 */
	static void glutReshape(int width, int height) {
		/*
		 * Adjust projection parameters.
		 */
		
		glViewport(0, 0, width, height);

		/*
		 * Correct projection
		 */
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		/*
		 * Normalize dimensions
		 */
		if (width > height) {
			GLdouble aspect = (GLdouble) width / height;
			gluOrtho2D(-aspect, aspect, -1., 1.);
		} else {
			GLdouble aspect = (GLdouble) height / width;
			gluOrtho2D(-1., 1., -aspect, aspect);
		}
			
		/*
		 * Back to ModelView
		 */
		glMatrixMode( GL_MODELVIEW );
	}

	/*
	 * Display callback
	 */
	static void glutDisplay() {

		/*
		 * Clear buffer
		 */
		glClearColor( 0.f, 0.f, 0.f, 0.f );
		glClear( GL_COLOR_BUFFER_BIT );

		/*
		 * Show result.
		 */
		glutSwapBuffers();
	}

public:
	ui_gl() {
#ifndef USE_GLUT
		const char *glut_error = "this build was not configured for GLUT";
		throw glut_error;
#endif
#ifndef HAVE_GLUT_MAIN_LOOP_EVENT
		const char *freeglut_error = "need extension glutMainLoopEvent()";
		throw freeglut_error;
#endif
#ifndef USE_PTHREAD
		const char *pthread_error = "need POSIX threads";
		throw pthread_error;
#endif
#ifndef HAVE_NANOSLEEP
		const char *nanosleep_error = "need nanosleep()";
		throw nanosleep_error;
#endif


		if (!gpu::is_ok()) {
			const char *gpu_init_error = "GPU initialization error";
			throw gpu_init_error;
		}

		instance = this;

		/*
		 * Initialization
		 */

		gpu::lock();
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		gluOrtho2D( -1., 1., -1., 1. );
		glMatrixMode( GL_MODELVIEW );
		gpu::unlock();
		
		/*
		 * XXX: We don't lock the GPU for these, which could cause a
		 * problem.
		 */

		glutReshapeFunc(glutReshape);
		glutDisplayFunc(glutDisplay);
		glutReshapeWindow(640, 480);

		/*
		 * Start an update thread.
		 */
#ifdef USE_PTHREAD
		pthread_attr_t pattr;
		pthread_attr_init(&pattr);
		pthread_attr_setdetachstate(&pattr, PTHREAD_CREATE_JOINABLE);
		pthread_create(&update_thread, NULL, thread_update_loop, (void *) this);
#endif
	}

	~ui_gl() {
#ifdef USE_PTHREAD
		pthread_cancel(update_thread);
		pthread_join(update_thread, NULL);
#endif
	}
};

#endif
