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

extern "C" {
#include <ale.h>
}

#include "accel.h"
#include "thread.h"

#define ALE_GPU_ASSERT_INCLUDE \
"(ale-define-macro (assert expr) `(if (not ,expr) (ale-error \"Assertion failed: \" ,(object->string expr))))"

class gpu {
public:
	class program {

		/*
		 * Libale kernel.
		 */

		void *lk;

	public:
		class shader {

			/*
			 * Libale kernel module.
			 */

			void *lkm;

		public:
			void attach_to_program(const program *p) const {
				// ale_add_kernel_module(p->lk, lkm);
			}
			shader(const char *source) {
				// lkm = ale_new_kernel_module(accel::context(), source);
			}
			shader (const shader &p) {
				assert(0);
			}
			shader &operator=(const shader &p) {
				assert(0);
			}

			~shader() {
				// ale_delete_kernel_module(lkm);
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

			/*
			 * XXX: if possible, it would probably be better not to
			 * use this method, as it would require addition of
			 * constant handling to Libale.
			 */

			assert(0);

#if 0
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
#endif
		}

		program() {
			// lk = ale_new_kernel(accel::context());
		}

		program (const program &p) {
			assert(0);
		}

		program &operator=(const program &p) {
			assert(0);
		}

		~program() {
			// ale_delete_kernel(lk);
		}

		void attach(const shader *s) {
			s->attach_to_program(this);
		}

		void attach(const library *l) {
			l->attach_shaders(this);
		}

		void link() {
			// ale_link_kernel(lk);
		}

	};
};

#endif
