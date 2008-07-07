// Copyright 2002, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

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
 * render.h: A superclass for all rendering classes.
 */

#ifndef __render_h__
#define __render_h__

#include "transformation.h"
#include "image.h"
#include "point.h"

#define ACTIVE_RENDERER_COUNT 30

#define ALE_GLSL_RENDER_INCLUDE \
"struct exclusion {\n"\
"	bool is_render;\n"\
"	float x[6];\n"\
"};\n"\
"struct render {\n"\
"	int rx_count;\n"\
"	exclusion rx_parameters[EXCLUSION_ARRAY_SIZE];\n"\
"};\n"\
"uniform render render_static;\n"\
"bool render_is_excluded_r(vec2 offset, vec4 position, int frame);\n"

/*
 * Class render accepts messages synchronizing rendering steps through the
 * methods sync(n) and sync(), and returns information about the currently
 * rendered image via methods get_image() and get_defined().  This class is
 * abstract, and must be subclassed to be instantiated.
 */

class render : public gpu::program::library {
private:
	static unsigned int rx_count;
	static exclusion *rx_parameters;
	static int rx_show;
	static render *directory[ACTIVE_RENDERER_COUNT];
	static int directory_length;
	static int extend;
	static ale_pos scale_factor;
	static ale_real wt;

	image **queue;
	unsigned int queue_size;
	int step_num;
	int entry_number;

	static int strpfix(const char *a, const char *b) {
		return strncmp(a, b, strlen(a));
	}

protected:
	/*
	 * Constructor
	 */
	render() {
		const char *shader_code =
			ALE_GLSL_RENDER_INCLUDE
			"bool render_is_excluded_r(vec2 offset, vec4 p, int f) {\n"
			"	for (int param = 0; param < render_static.rx_count; param++)\n"
			"		if (render_static.rx_parameters[param].is_render\n"
			"		 && p.x + offset.x >= render_static.rx_parameters[param].x[0]\n"
			"		 && p.x + offset.x <= render_static.rx_parameters[param].x[1]\n"
			"		 && p.y + offset.y >= render_static.rx_parameters[param].x[2]\n"
			"		 && p.y + offset.y <= render_static.rx_parameters[param].x[3]\n"
			"		 && float(f) >= render_static.rx_parameters[param].x[4]\n"
			"		 && float(f) <= render_static.rx_parameters[param].x[5])\n"
			"			return true;\n"
			"	return false;\n"
			"}\n";

		if (accel::is_gpu()) {
			gpu_shader = new gpu::program::shader(shader_code);
		}

		if (directory_length >= ACTIVE_RENDERER_COUNT) {
			fprintf(stderr, "\n\n*** Too many renderers in d2::render::render() ***\n\n");
			exit(1);
		}

		directory[directory_length] = this;
		entry_number = directory_length;

		directory_length++;

		step_num = -1;
		queue = NULL;
		queue_size = 0;
	}

	/*
	 * Get extension state
	 */
	int is_extend() {
		return extend;
	}

	/*
	 * Get the scale factor
	 */
	ale_pos get_scale_factor() {
		return scale_factor;
	}

	/*
	 * Get the current step number
	 */

	int get_step() {
		return step_num;
	}

	/*
	 * Perform the current rendering step.
	 */

	virtual void step() = 0;

public:

	/*
	 * Check for render-coordinate excluded regions.  (Applies an offset to
	 * spatial coordinates internally.)
	 */
	static int is_excluded_r(point offset, point p, int f) {

		for (unsigned int param = 0; param < rx_count; param++)
			if (rx_parameters[param].type == exclusion::RENDER
			 && p[0] + offset[0] >= rx_parameters[param].x[0]
			 && p[0] + offset[0] <= rx_parameters[param].x[1]
			 && p[1] + offset[1] >= rx_parameters[param].x[2]
			 && p[1] + offset[1] <= rx_parameters[param].x[3]
			 && f >= rx_parameters[param].x[4]
			 && f <= rx_parameters[param].x[5])
				return 1;

		return 0;
	}
	static int is_excluded_r(point offset, int i, int j, int f) {

		for (unsigned int param = 0; param < rx_count; param++)
			if (rx_parameters[param].type == exclusion::RENDER
			 && i + offset[0] >= rx_parameters[param].x[0]
			 && i + offset[0] <= rx_parameters[param].x[1]
			 && j + offset[1] >= rx_parameters[param].x[2]
			 && j + offset[1] <= rx_parameters[param].x[3]
			 && f >= rx_parameters[param].x[4]
			 && f <= rx_parameters[param].x[5])
				return 1;

		return 0;
	}
	int is_excluded_r(int i, int j, int f) {
		return is_excluded_r(get_image()->offset(), i, j, f);
	}

	/*
	 * Check for frame-coordinate excluded regions.
	 */
	static int is_excluded_f(point p, int f) {

		for (unsigned int param = 0; param < rx_count; param++)
			if (rx_parameters[param].type == exclusion::FRAME
			 && p[0] >= rx_parameters[param].x[0]
			 && p[0] <= rx_parameters[param].x[1]
			 && p[1] >= rx_parameters[param].x[2]
			 && p[1] <= rx_parameters[param].x[3]
			 && f >= rx_parameters[param].x[4]
			 && f <= rx_parameters[param].x[5])
				return 1;

		return 0;
	}
	static int is_excluded_f(int i, int j, int f) {

		for (unsigned int param = 0; param < rx_count; param++)
			if (rx_parameters[param].type == exclusion::FRAME
			 && i >= rx_parameters[param].x[0]
			 && i <= rx_parameters[param].x[1]
			 && j >= rx_parameters[param].x[2]
			 && j <= rx_parameters[param].x[3]
			 && f >= rx_parameters[param].x[4]
			 && f <= rx_parameters[param].x[5])
				return 1;

		return 0;
	}

	static int render_count() {
		return directory_length;
	}
	static render *render_num(int n) {
		assert (n < directory_length);
		return directory[n];
	}

	static void render_init(unsigned int _rx_count, exclusion *_rx_parameters, 
			int _rx_show, int _extend, ale_pos _scale_factor) {
		rx_count = _rx_count;
		rx_show = _rx_show;
		extend = _extend;
		scale_factor = _scale_factor;

		rx_parameters = (exclusion *) malloc(rx_count * sizeof(exclusion));

		for (unsigned int region = 0; region < rx_count; region++) {

			rx_parameters[region] = _rx_parameters[region];

			/*
			 * Scale spatial rendering coordinates
			 */

			if (rx_parameters[region].type == exclusion::RENDER)
				for (int p = 0; p < 4; p++)
					rx_parameters[region].x[p] *= scale_factor;
		}
	}

	static void set_wt(ale_real _wt) {
		wt = _wt;
	}

	static ale_real get_wt() {
		return wt;
	}

	static int is_rx_show() {
		return rx_show;
	}

	static unsigned int get_rx_count() {
		return rx_count;
	}

	static const exclusion *get_rx_parameters() {
		return rx_parameters;
	}

	/*
	 * Current rendering result.
	 */

	virtual const image *get_image() const = 0;

	/*
	 * Result of rendering at the given frame.
	 */

	const image *get_image(unsigned int n) {
		sync(n);

		if (n == (unsigned int) step_num)
			return get_image();

		n = step_num - n - 1;

		assert (n < queue_size);

		return queue[n];
	}

	/*
	 * Extend the rendering queue.
	 */

	void extend_queue(unsigned int n) {
		/*
		 * Increase the size of the queue, if necessary, to
		 * accommodate the given lag.
		 */
		if (n > queue_size) {
			unsigned int new_size = n;
			queue = (image **) realloc(queue, new_size * sizeof(image *));
			assert(queue);
			if (queue == NULL) {
				fprintf(stderr, "\n\n*** VISE: Unable to allocate memory ***\n\n\n");
				exit(1);
			}
			memset(queue + queue_size, 0, (new_size - queue_size) * sizeof(image *));
			queue_size = new_size;
		}
	}

	/*
	 * Definition map.  Unit-depth image whose pixels are nonzero where
	 * the image is defined.
	 */

	virtual const image *get_defined() const = 0;

	/*
	 * Sync.
	 */

	virtual void sync(int n) {
		assert (step_num >= -1);
		for (int i = step_num + 1; i <= n; i++) {
			if (queue_size > 0 && step_num >= 0) {
				/*
				 * Shift the current queue so that the new head remains at the
				 * zero index.  There are more time-efficient ways to handle
				 * queues, but the benefits are not clear in this case.
				 */
				delete queue[queue_size - 1];
				for (int i = queue_size - 1; i > 0; i--) {
					queue[i] = queue[i - 1];
				}
				queue[0] = get_image()->clone("Render queue clone");
			}

			step_num++;
			step();
		}
	}

	/*
	 * Perform any final rendering steps.  Return a non-zero value if
	 * anything changed.
	 */

	virtual int sync() {
		return 0;
	}

	/*
	 * Set point rendering bounds, if possible.
	 */

	virtual void init_point_renderer(unsigned int h, unsigned int w, unsigned int d) {
		assert(0);
		fprintf(stderr, "Error: init_point_renderer() not supported by this renderer\n");
		exit(1);
	}

	/*
	 * Point render.
	 */

	virtual void point_render(unsigned int i, unsigned int j, unsigned int f, transformation t) {
		assert(0);
		fprintf(stderr, "Error: point_render() not supported by this renderer\n");
		exit(1);
	}

	/*
	 * Finish point rendering.
	 */

	virtual void finish_point_rendering() {
		assert(0);
		fprintf(stderr, "Error: finish_point_rendering() not supported by this renderer\n");
		exit(1);
	}

	virtual ~render() {
		directory[entry_number] = NULL;
	}

	int entry() {
		return entry_number;
	}

	virtual void free_memory() = 0;

	static void free_entry(int n) {
		if (directory[n] != NULL) {
			directory[n]->free_memory();
			delete directory[n];
			directory[n] = NULL;
		}
	}

	static void free_all_memory() {
		for (int i = 0; i < ACTIVE_RENDERER_COUNT; i++)
			free_entry(i);

		directory_length = 0;
	}

	static void reset() {
		free_all_memory();
	}
};

#endif
