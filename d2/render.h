// Copyright 2002, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
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

/*
 * Class render accepts messages synchronizing rendering steps through the
 * methods sync(n) and sync(), and returns information about the currently
 * rendered image via methods get_image() and get_defined().  This class is
 * abstract, and must be subclassed to be instantiated.
 */

class render {
private:
	static unsigned int rx_count;
	static ale_pos *rx_parameters;
	static int rx_show;
	static render *directory[ACTIVE_RENDERER_COUNT];
	static int directory_length;
	static int extend;
	static double scale_factor;
	static double wt;

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
	double get_scale_factor() {
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
	 * Check for excluded regions.  (Applies an offset to spatial
	 * coordinates internally.)
	 */
	static int is_excluded(point offset, point p, int f) {

		for (unsigned int param = 0; param < rx_count; param++)
			if (p[0] + offset[0] >= rx_parameters[6 * param + 0]
			 && p[0] + offset[0] <= rx_parameters[6 * param + 1]
			 && p[1] + offset[1] >= rx_parameters[6 * param + 2]
			 && p[1] + offset[1] <= rx_parameters[6 * param + 3]
			 && f >= rx_parameters[6 * param + 4]
			 && f <= rx_parameters[6 * param + 5])
				return 1;

		return 0;
	}
	static int is_excluded(point offset, int i, int j, int f) {

		for (unsigned int param = 0; param < rx_count; param++)
			if (i + offset[0] >= rx_parameters[6 * param + 0]
			 && i + offset[0] <= rx_parameters[6 * param + 1]
			 && j + offset[1] >= rx_parameters[6 * param + 2]
			 && j + offset[1] <= rx_parameters[6 * param + 3]
			 && f >= rx_parameters[6 * param + 4]
			 && f <= rx_parameters[6 * param + 5])
				return 1;

		return 0;
	}

	int is_excluded(int i, int j, int f) {
		return is_excluded(get_image()->offset(), i, j, f);
	}

	static int render_count() {
		return directory_length;
	}
	static render *render_num(int n) {
		assert (n < directory_length);
		return directory[n];
	}

	static void render_init(unsigned int _rx_count, int *_rx_parameters, 
			int _rx_show, int _extend, double _scale_factor) {
		rx_count = _rx_count;
		rx_show = _rx_show;
		extend = _extend;
		scale_factor = _scale_factor;

		rx_parameters = (ale_pos *) malloc(rx_count * 6 * sizeof(ale_pos));

		for (unsigned int param = 0; param < rx_count; param++) {
			rx_parameters[param * 6 + 0] = _rx_parameters[param * 6 + 0] * scale_factor;
			rx_parameters[param * 6 + 1] = _rx_parameters[param * 6 + 1] * scale_factor;
			rx_parameters[param * 6 + 2] = _rx_parameters[param * 6 + 2] * scale_factor;
			rx_parameters[param * 6 + 3] = _rx_parameters[param * 6 + 3] * scale_factor;
			rx_parameters[param * 6 + 4] = _rx_parameters[param * 6 + 4];
			rx_parameters[param * 6 + 5] = _rx_parameters[param * 6 + 5];
		}
	}

	static void set_wt(double _wt) {
		wt = _wt;
	}

	static double get_wt() {
		return wt;
	}

	static int is_rx_show() {
		return rx_show;
	}

	static unsigned int get_rx_count() {
		return rx_count;
	}

	static const ale_pos *get_rx_parameters() {
		return rx_parameters;
	}

	/*
	 * Current rendering result.
	 */

	virtual const image *get_image() = 0;

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

	virtual const image *get_defined() = 0;

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

	virtual void free_memory() = 0;

	static void free_all_memory() {

		for (int i = 0; i < ACTIVE_RENDERER_COUNT; i++)
			if (directory[i] != NULL) {
				directory[i]->free_memory();
				delete directory[i];
				directory[i] = NULL;
			}

		directory_length = 0;

	}

	static void reset() {
		free_all_memory();
	}
};

#endif
