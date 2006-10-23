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

#ifndef __incremental_h__
#define __incremental_h__

#include "invariant.h"
#include "../render.h"
#include "../transformation.h"
#include "../image.h"
#include "../point.h"

/*
 * Class for incremental renderers.
 */

class incremental : public render {
protected:
	image_weighted_avg *accum_image;
	invariant *inv;

	/*
	 * Set extents of image and weight according to a new image to be
	 * merged.  This function should remove only superfluous undefined
	 * areas.
	 */

	void set_extents_by_map(unsigned int frame_num, transformation t) {

		assert (accum_image  != NULL);

		ale_pos extend_offset_i = accum_image->offset()[0];
		ale_pos extend_offset_j = accum_image->offset()[1];

		int extend_top = 0;
		int extend_bottom = 0;
		int extend_left = 0;
		int extend_right = 0;

		double zero = 0;
		double infinity = 1 / zero;
		
		assert (!finite(infinity));
		assert (!isnan(infinity));
		assert (infinity > 0);

		point min, max;

		min[0] = min[1] = infinity;
		max[0] = max[1] = -infinity;

		for (unsigned int i = 0; i < t.unscaled_height(); i++)
		for (unsigned int j = 0; j < t.unscaled_width(); j++) {

			if (is_excluded_f(i, j, frame_num))
				continue;
			
			point p = t.transform_unscaled(point(i, j));

			if (is_excluded_r(accum_image->offset(), p, frame_num))
				continue;

			if (p[0] < min[0]) {
				min[0] = p[0];
			}
			if (p[0] > max[0]) {
				max[0] = p[0];
			}
			if (p[1] < min[1]) {
				min[1] = p[1];
			}
			if (p[1] > max[1]) {
				max[1] = p[1];
			}
		}

		if (!finite(max[0])
		 || !finite(max[1])
		 || !finite(min[0])
		 || !finite(min[1]))
			return;

		extend_top = (int) ceil(extend_offset_i - floor(min[0]));
		extend_left = (int) ceil(extend_offset_j - floor(min[1]));
		extend_bottom = (int) ceil(ceil(max[0]) - (accum_image->height() - 1 + extend_offset_i));
		extend_right = (int) ceil(ceil(max[1]) - (accum_image->width() - 1 + extend_offset_j));

		accum_image->extend(extend_top, extend_bottom, 
			extend_left, extend_right);
	}	

	void increase_extents_by_map(unsigned int frame_num, transformation t) {

		assert (accum_image  != NULL);

		ale_pos extend_offset_i = accum_image->offset()[0];
		ale_pos extend_offset_j = accum_image->offset()[1];

		int extend_top = 0;
		int extend_bottom = 0;
		int extend_left = 0;
		int extend_right = 0;

		double zero = 0;
		double infinity = 1 / zero;
		
		assert (!finite(infinity));
		assert (!isnan(infinity));
		assert (infinity > 0);

		point min, max;

		min[0] = min[1] = infinity;
		max[0] = max[1] = -infinity;

		for (unsigned int i = 0; i < t.unscaled_height(); i++)
		for (unsigned int j = 0; j < t.unscaled_width(); j++) {

			if (is_excluded_f(i, j, frame_num))
				continue;
			
			point p = t.transform_unscaled(point(i, j));

			if (is_excluded_r(point(0, 0), p, frame_num))
				continue;

			if (p[0] < min[0]) {
				min[0] = p[0];
			}
			if (p[0] > max[0]) {
				max[0] = p[0];
			}
			if (p[1] < min[1]) {
				min[1] = p[1];
			}
			if (p[1] > max[1]) {
				max[1] = p[1];
			}
		}

		if (!finite(max[0])
		 || !finite(max[1])
		 || !finite(min[0])
		 || !finite(min[1]))
			return;

		if (ceil(min[0]) < extend_offset_i)
			extend_top = (int) ceil(extend_offset_i - floor(min[0]));
		if (ceil(min[1]) < extend_offset_j)
			extend_left = (int) ceil(extend_offset_j - floor(min[1]));
		if (floor(max[0]) > accum_image->height() - 1 + extend_offset_i)
			extend_bottom = (int) ceil(ceil(max[0]) - (accum_image->height() - 1 + extend_offset_i));
		if (floor(max[1]) > accum_image->width() - 1 + extend_offset_j)
			extend_right = (int) ceil(ceil(max[1]) - (accum_image->width() - 1 + extend_offset_j));

		accum_image->extend(extend_top, extend_bottom, 
			extend_left, extend_right);
	}	

	/*
	 * Merge operation for a single pixel in the accumulated image.
	 */
	void _merge_pixel(int frame, const image *delta, transformation t, int i, int j, const filter::ssfe *_ssfe) {

		if (_ssfe->ex_is_honored() && is_excluded_r(i, j, frame))
			return;

		if (accum_image->accumulate_norender(i, j))
			return;
		
		/*
		 * Pixel value to be merged, and the associated
		 * confidence
		 */

		pixel value, confidence;

		_ssfe->filtered(i, j, frame, &value, &confidence);

		accum_image->accumulate(i, j, frame, value, confidence);
	}

	/*
	 * Merge part of a delta frame with part of the accumulated image using
	 * the specified transformation.
	 */

	struct subdomain_args {
		incremental *instance;
		int frame;
		const image *delta;
		transformation t;
		unsigned int i_min, i_max, j_min, j_max;
	};

	void _merge_subdomain(void *args) {
		subdomain_args *sargs = (subdomain_args *) args;

		int frame = sargs->frame;
		const image *delta = sargs->delta;
		transformation t = sargs->t;
		unsigned int i_min = sargs->i_min;
		unsigned int i_max = sargs->i_max;
		unsigned int j_min = sargs->j_min;
		unsigned int j_max = sargs->j_max;

		point offset = accum_image->offset();

		assert (accum_image != NULL);
		assert (delta != NULL);

		const filter::ssfe *_ssfe = inv->ssfe();

		for (unsigned int i = i_min; i < i_max; i++)
		for (unsigned int j = j_min; j < j_max; j++) {

#if 0
			/*
			 * This is untested, but it should work, and is less
			 * verbose than what follows.
			 */

			_merge_pixel(frame, delta, t, i, j, _ssfe);
#else

			if (_ssfe->ex_is_honored() && is_excluded_r(i, j, frame))
				continue;

			if (accum_image->accumulate_norender(i, j))
				continue;
			
			/*
			 * Pixel value to be merged, and the associated
			 * confidence
			 */

			pixel value, confidence;

			_ssfe->filtered(i, j, frame, &value, &confidence);

			accum_image->accumulate(i, j, frame, value, confidence);
#endif
		}
	}

	static void *_merge_run_subdomain(void *args) {
		subdomain_args *sargs = (subdomain_args *) args;
		sargs->instance->_merge_subdomain(args);
#ifdef USE_PTHREAD
		pthread_exit(0);
#else
		return NULL;
#endif
	}

	void
	_merge(int frame, const image *delta, transformation t) {

		point offset = accum_image->offset();

		assert (accum_image != NULL);
		assert (delta != NULL);

		const filter::ssfe *_ssfe = inv->ssfe();

		_ssfe->set_parameters(t, delta, offset);

		int N;

#ifdef USE_PTHREAD
		N = thread::count();

		pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t) * N);
		pthread_attr_t *thread_attr = (pthread_attr_t *) malloc(sizeof(pthread_attr_t) * N);

#else
		N = 1;
#endif

		subdomain_args *args = (subdomain_args *) malloc(sizeof(subdomain_args) * N);

		for (int ti = 0; ti < N; ti++) {
			args[ti].instance = this;
			args[ti].frame = frame;
			args[ti].delta = delta;
			args[ti].t = t;
			args[ti].i_min = (accum_image->height() * ti) / N;
			args[ti].i_max = (accum_image->height() * (ti + 1)) / N;
			args[ti].j_min = 0;
			args[ti].j_max = accum_image->width();

#ifdef USE_PTHREAD
			pthread_attr_init(&thread_attr[ti]);
			pthread_attr_setdetachstate(&thread_attr[ti], PTHREAD_CREATE_JOINABLE);
			pthread_create(&threads[ti], &thread_attr[ti], _merge_run_subdomain, &args[ti]);
#else
			_merge_subdomain(&args[ti]);
#endif
		}

#ifdef USE_PTHREAD
		for (int ti = 0; ti < N; ti++) {
			pthread_join(threads[ti], NULL);
		}
#endif

		free(args);
	}

public:

	/*
	 * Constructor
	 */
	incremental(invariant *inv) {
		this->inv = inv;
		accum_image = NULL;
	}
	
	/*
	 * Invariant
	 */
	const invariant *get_invariant() const {
		return inv;
	}

	/*
	 * Result of rendering.
	 */

	virtual const image *get_image() {
		assert (accum_image != NULL);
		return accum_image->get_colors();
	}

	/*
	 * Definition map.  Unit-depth image whose pixels are nonzero where
	 * the image is defined.
	 */

	virtual const image *get_defined() {
		assert (accum_image != NULL);
		return accum_image->get_weights();
	}

	/*
	 * Perform the current rendering step.
	 */
	virtual void step() {
		assert (get_step() >= -1);
		if (get_step() == 0) {
			transformation t = align::of(0);

			const image *im = image_rw::open(0);

			ui::get()->rendering();

			if (inv->is_median())
				accum_image = new image_weighted_median(1, 1, 3);
			else
				accum_image = new image_weighted_simple(1, 1, 3, inv);

			set_extents_by_map(0, t);

			_merge(0, im, t);

			image_rw::close(0);
		} else if (align::match(get_step())) {
			transformation t = align::of(get_step());
			ui::get()->rendering();
			if (is_extend())
				increase_extents_by_map(get_step(), t);
			const image *im = image_rw::open(get_step());
			_merge(get_step(), im, t);
			image_rw::close(get_step());
		}
	}


	virtual void init_point_renderer(unsigned int h, unsigned int w, unsigned int d) {
		assert(accum_image == NULL);

		if (inv->is_median())
			accum_image = new image_weighted_median(h, w, d);
		else
			accum_image = new image_weighted_simple(h, w, d, inv);

		assert(accum_image);
	}
	
	virtual void point_render(unsigned int i, unsigned int j, unsigned int f, transformation t) {
		const image *im = d2::image_rw::get_open(f);

		const filter::ssfe *_ssfe = inv->ssfe();

		_ssfe->set_parameters(t, im, accum_image->offset());
		_merge_pixel(f, im, t, i, j, _ssfe);
	}

	virtual void finish_point_rendering() {
		return;
	}

	void free_memory() {
		delete accum_image;
		accum_image = NULL;
	}
};

#endif
