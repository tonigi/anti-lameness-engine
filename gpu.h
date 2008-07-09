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

#define ALE_GLSL_ASSERT_INCLUDE \
"#define assert(x) { if ((x)); }\n"

class gpu {
	static int gpu_initialized;

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

	static int is_enabled();

	class image {
#ifdef USE_GLUT
		/*
		 * For now, assume that the image fits in a single texture.
		 */
		GLuint texture;
#endif
	};

	class program {

		static const char **program_buffer;
		static int program_buffer_size;

		static void size_program_buffer() {
			program_buffer = (const char **) realloc(program_buffer, sizeof(const char *) * program_buffer_size);
			assert(program_buffer);
		}

		static void check_log(const char *description, const char *log, int program_lines = 0, const char **program = NULL) {
			if (strcmp(log, "")) {
				fprintf(stderr, "GLSL %s error log:\n", description);
				fprintf(stderr, "%s", log);
				if (program) {
					fprintf(stderr, "Program text:\n");
					for (int i = 0; i < program_lines; i++)
						fprintf(stderr, "%s", program[i]);
				}
				exit(1);
			}
		}

#ifdef USE_GLEW
		static void check_id(GLuint id) {
			if (id == 0) {
				fprintf(stderr, "GLSL object creation failed.\n");
				exit(1);
			}
		}

		GLuint id;
#endif

	public:
		class shader {
#ifdef USE_GLEW
			GLuint id;
#endif

		public:
			void attach_to_program(const program *p) const {
#ifdef USE_GLEW
				glAttachShader(p->id, id);
#endif
			}
			shader(const char *source) {
#ifdef USE_GLEW

				assert(is_enabled());

				program_buffer_size++;
				size_program_buffer();
				program_buffer[program_buffer_size - 1] = source;

				char log[1000] = "";
				id = glCreateShader(GL_FRAGMENT_SHADER_ARB);
				check_id(id);
				glShaderSource(id, program_buffer_size, program_buffer, NULL);
				glCompileShader(id);
				glGetShaderInfoLog(id, 1000, NULL, log);
				check_log("shader compilation", log, program_buffer_size, program_buffer);

				program_buffer_size--;
				size_program_buffer();
#endif
			}
			shader (const shader &p) {
				assert(0);
			}

			shader &operator=(const shader &p) {
				assert(0);
			}

			~shader() {
#ifdef USE_GLEW
				glDeleteShader(id);
#endif
			}
		};

		class library {
		protected:
			shader *gpu_shader;
		public:
			library() {
				gpu_shader = NULL;
			}

			/*
			 * For libraries containing multiple shaders, or
			 * referencing shaders not stored locally, this
			 * method can be overridden.
			 */
			virtual void attach_shaders(program *p) const {
				assert(p);
				if (gpu_shader)
					p->attach(gpu_shader);
			}
		};

		static void set_constant(const char *name, int value) {
			program_buffer_size++;
			size_program_buffer();

			const int line_length = 1000;

			char *program_line = (char *) malloc(line_length * sizeof(char));
			assert(program_line);

#if 0
			/*
			 * XXX: This seems to generate link errors, for some reason.
			 */
			snprintf(program_line, line_length, "const int %s = %d;\n", name, value);
#else
			snprintf(program_line, line_length, "#define %s %d\n", name, value);
#endif

			program_buffer[program_buffer_size - 1] = program_line;
		}

		program() {
#ifdef USE_GLEW

			assert(is_enabled());

			id = glCreateProgram();
			check_id(id);
#endif
		}

		program (const program &p) {
			assert(0);
		}

		program &operator=(const program &p) {
			assert(0);
		}

		~program() {
#ifdef USE_GLEW
			glDeleteProgram(id);
#endif
		}

		void attach(const shader *s) {
			s->attach_to_program(this);
		}

		void attach(const library *l) {
			l->attach_shaders(this);
		}

		void link() {
#ifdef USE_GLEW
			char log[1000] = "";
			glLinkProgram(id);
			glGetProgramInfoLog(id, 1000, NULL, log);
			check_log("program linking", log);
#endif
		}

	};
};

#endif
