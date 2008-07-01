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

#include "../d2.h"
#include "ui.h"
#include "accel.h"
#include "gpu.h"

/*
 * OpenGL user interface.
 */

class ui_gl : public ui {
private:
	int osd_mode;
	int width, height;
	int cursor_line;
	int cursor_period;
	int dirty;

#ifdef USE_PTHREAD
	pthread_t update_thread;
#endif

	void printf(const char *format, ...) {
		update();

		/*
		 * Reject messages that aren't loud.
		 */

		if (!strstr(format, "\n***"))
			return;

		va_list ap;
		va_start(ap, format);
		vfprintf(ui_stream, format, ap);
		va_end(ap);
	}

	void update() {
		dirty = 1;
	}

	static void *glut_thread(void *vu) {
#if defined HAVE_GLUT_MAIN_LOOP_EVENT && defined HAVE_NANOSLEEP
		/*
		 * Initialization
		 */

		if (!gpu::is_ok()) {
			const char *gpu_init_error = "GPU initialization error";
			throw gpu_init_error;
		}

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
		glutKeyboardFunc(glutKeyboard);
		glutReshapeWindow(640, 480);

		for(;;) {
			struct timespec t;
			t.tv_sec = 0;
			t.tv_nsec = 10000000;
			nanosleep(&t, NULL);

			if (instance->dirty) {
				instance->dirty = 0;
				glutPostRedisplay();
			}

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
		 * Correct projection and window parameters.
		 */
		glViewport(0, 0, width, height);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(0, width, height, 0);
		glMatrixMode( GL_MODELVIEW );

		instance->width = width;
		instance->height = height;
	}

	static void bitmap_string(const char *s) {
		while (*s)
			glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *s++);
	}

	static void outline_string(int i, int j, const char *s) {
		for (int ii = -1; ii <= 1; ii++)
		for (int jj = -1; jj <= 1; jj++) {
			glColor3f(0, 0, 0);
			glRasterPos2f(i + ii, j + jj);
			bitmap_string(s);
		}

		glColor3f(1, 1, 1);
		glRasterPos2f(i, j);
		bitmap_string(s);
	}

	static void printf_glut(const char *format, ...) {
		char line_buffer[1000];
		line_buffer[999] = '\0';
		int line_buffer_len = 999;

		va_list ap;
		va_start(ap, format);
		vsnprintf(line_buffer, line_buffer_len, format, ap);
		va_end(ap);

		outline_string(10, instance->cursor_line, line_buffer);
		instance->cursor_line += instance->cursor_period;
	}

	/*
	 * Keyboard callback
	 */
	static void glutKeyboard(unsigned char c, int x, int y) {
		if (c == 'o') {
			instance->osd_mode = (instance->osd_mode + 1) % 3;
			glutPostRedisplay();
		} else if (c == 'O') {
			instance->osd_mode = (instance->osd_mode + 2) % 3;
			glutPostRedisplay();
		}
	}


	/*
	 * Display callback
	 */
	static void glutDisplay() {
		/*
		 * Clear buffer
		 */
		glClearColor(0, 0, 0, 0);
		glClear(GL_COLOR_BUFFER_BIT);

		/*
		 * Display text.
		 */

		instance->cursor_line = 23;
		instance->cursor_period = 18;

		int baseline = instance->height - 15;

		if (instance->osd_mode >= 1) {
			const char *action = "";
			int current_frame = instance->status.frame_num;
			int total_frames = d2::image_rw::count();
			
			switch(instance->status.code) {
			case status_type::LOAD_FILE:
				action = "Loading frame";
				break;
			case status_type::EXPOSURE_PASS_1:
			case status_type::PREMATCH:
			case status_type::ALIGN:
			case status_type::GLOBAL_ALIGN:
			case status_type::POSTMATCH:
			case status_type::EXPOSURE_PASS_2:
			case status_type::MULTI:
				action = "Aligning";
				break;
			case status_type::RENDERA:
			case status_type::RENDERD:
			case status_type::RENDERO:
				action = "Rendering";
				break;
			case status_type::WRITED:
			case status_type::WRITEO:
			case status_type::IP_WRITE:
				action = "Writing";
				break;
			case status_type::IP_RENDER:
			case status_type::IP_UPDATE:
				action = "Irani-Peleg";
				break;
			case status_type::D3_CONTROL_POINT_SOLVE:
			case status_type::D3_SUBDIVIDING_SPACE:
			case status_type::D3_UPDATING_OCCUPANCY:
			case status_type::D3_RENDER:
				action = "3D operations";
				break;
			default:
				break;
			}
			if (strcmp(action, "")) {
				printf_glut("%s (%u/%u)", action, current_frame + 1, total_frames);
			}
		}

		if (instance->osd_mode >= 2) {
			switch(instance->status.code) {
			case status_type::EXPOSURE_PASS_1:
			case status_type::EXPOSURE_PASS_2:
			case status_type::PREMATCH:
			case status_type::ALIGN:
			case status_type::GLOBAL_ALIGN:
			case status_type::POSTMATCH:
			case status_type::MULTI:
				printf_glut("  degree=%u/%u", instance->status.element_degree + 1, instance->status.total_degree + 1);
				printf_glut(" element=%u/%u", instance->status.element_index + 1, instance->status.num_elements);
				printf_glut("   match=%.6f%%", instance->status.match_value);
				printf_glut(" perturb= %.3f", instance->status.perturb_size);
				printf_glut("     lod= %.3f", pow(2, -instance->status.align_lod));
				printf_glut("exp_mult= %.3f %.3f %.3f",
						instance->status.exp_multiplier[0],
						instance->status.exp_multiplier[1],
						instance->status.exp_multiplier[2]);
				break;
			case status_type::WRITED:
			case status_type::IP_WRITE:
				printf_glut(d2::image_rw::output_name());
				break;

			case status_type::WRITEO:
				break;

			case status_type::IP_RENDER:
			case status_type::IP_UPDATE:
				printf_glut("Step %u/%u", instance->status.irani_peleg_step + 1, instance->status.irani_peleg_steps);
				if (instance->status.irani_peleg_stage == 1) {
					printf_glut("Simulate");
				} else if (instance->status.irani_peleg_stage == 2) {
					printf_glut("Backproject");
				}
				break;
			case status_type::D3_CONTROL_POINT_SOLVE:
			case status_type::D3_SUBDIVIDING_SPACE:
			case status_type::D3_UPDATING_OCCUPANCY:
			case status_type::D3_RENDER:
				break;
			default:
				break;
			}
		}

		if (instance->cursor_line < baseline - instance->cursor_period)
			instance->cursor_line = baseline - instance->cursor_period;

		if (accel::is_gpu() == 0)
			printf_glut("No preview available (No GPU acceleration)");
		else 
			printf_glut("No preview available");

		printf_glut("'o' cycles text");

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


		instance = this;

		osd_mode = 2;

		/*
		 * Start a GLUT thread.
		 */
#ifdef USE_PTHREAD
		pthread_attr_t pattr;
		pthread_attr_init(&pattr);
		pthread_attr_setdetachstate(&pattr, PTHREAD_CREATE_JOINABLE);
		pthread_create(&update_thread, NULL, glut_thread, (void *) this);
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
