// Copyright 2002, 2004, 2007 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                          <dhilvert@ugcs.caltech.edu>

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
	gpu::program *gpu_program;
	gpu::program::shader *gpu_shader;
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

		ale_pos zero = 0;
		ale_pos infinity = 1 / zero;
		
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
		extend_bottom = (int) ceil(ceil(max[0]) - (ale_pos) (accum_image->height() - 1 + extend_offset_i));
		extend_right = (int) ceil(ceil(max[1]) - (ale_pos) (accum_image->width() - 1 + extend_offset_j));

		accum_image->_extend(extend_top, extend_bottom, 
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
			extend_bottom = (int) ceil(ceil(max[0]) - (ale_pos) (accum_image->height() - 1 + extend_offset_i));
		if (floor(max[1]) > accum_image->width() - 1 + extend_offset_j)
			extend_right = (int) ceil(ceil(max[1]) - (ale_pos) (accum_image->width() - 1 + extend_offset_j));

		accum_image->_extend(extend_top, extend_bottom, 
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

		if (exposure::get_confidence() != 0) {
			_ssfe->filtered(i, j, frame, &value, &confidence, ((pixel) accum_image->get_pixel(i, j)), accum_image->get_weights()->get_pixel(i, j));
		} else {
			_ssfe->filtered(i, j, frame, &value, &confidence);
		}

		accum_image->accumulate(i, j, frame, value, confidence);
	}

	/*
	 * Merge part of a delta frame with part of the accumulated image using
	 * the specified transformation.
	 */

	class merge : public thread::decompose_domain {
		incremental *instance;
		int frame;
		const image *delta;
		transformation t;
		invariant *inv;
		image_weighted_avg *accum_image;
	protected:

		/*
		 * Texture implementation.
		 */
		void accel_run() {
			point offset = accum_image->offset();

			assert (accum_image != NULL);
			assert (delta != NULL);

#if 0
			const filter::ssfe *_ssfe = inv->ssfe();
#endif

			accum_image->accumulate_accel(instance->gpu_program);
		}

		/*
		 * Iterative domain decomposition implementation.
		 */
		void prepare_subdomains(unsigned int N) {
			ale_pos_disable_casting();
			ale_real_disable_casting();
			ale_accum_disable_casting();
		}
		void subdomain_algorithm(unsigned int thread,
				int i_min, int i_max, int j_min, int j_max) {

			point offset = accum_image->offset();

			assert (accum_image != NULL);
			assert (delta != NULL);

			const filter::ssfe *_ssfe = inv->ssfe();

			for (int i = i_min; i < i_max; i++)
			for (int j = j_min; j < j_max; j++) {

#if 0
				/*
				 * This is untested, but it should work, and is less
				 * verbose than what follows.
				 */

				instance->_merge_pixel(frame, delta, t, i, j, _ssfe);
#else

				if (_ssfe->ex_is_honored() && instance->is_excluded_r(i, j, frame))
					continue;

				if (accum_image->accumulate_norender(i, j))
					continue;
				
				/*
				 * Pixel value to be merged, and the associated
				 * confidence
				 */

				pixel value, confidence;

				if (exposure::get_confidence() != 0) {
					_ssfe->filtered(i, j, frame, &value, &confidence, 
							((pixel) accum_image->get_pixel(i, j)), 
							accum_image->get_weights()->get_pixel(i, j));
				} else {
					_ssfe->filtered(i, j, frame, &value, &confidence);
				}

				accum_image->accumulate(i, j, frame, value, confidence);
#endif
			}
		}
		void finish_subdomains(unsigned int N) {
			ale_pos_enable_casting();
			ale_real_enable_casting();
			ale_accum_enable_casting();
		}

	public:
		merge(incremental *_instance,
		      int _frame,
		      const image *_delta,
		      transformation _t) : decompose_domain(0, _instance->accum_image->height(),
		                                            0, _instance->accum_image->width()),
					   t(_t) {

			instance = _instance;
			frame = _frame;
			delta = _delta;
			t = _t;

			inv = instance->inv;
			accum_image = instance->accum_image;
		}

#if 0
		void run() {
			accel_run();
		}
#endif
	};

	void
	_merge(int frame, const image *delta, transformation t) {

		ui::get()->d2_incremental_start();

		point offset = accum_image->offset();

		assert (accum_image != NULL);
		assert (delta != NULL);

		const filter::ssfe *_ssfe = inv->ssfe();

		_ssfe->set_parameters(t, delta, offset);

		merge m(this, frame, delta, t);
		m.run();

		ui::get()->d2_incremental_stop();
	}

public:

	/*
	 * Constructor
	 */
	incremental(invariant *inv) {
		const char *shader_main =
			ALE_GLSL_RENDER_INCLUDE
			ALE_GLSL_SSFE_INCLUDE
			ALE_GLSL_IMAGE_WEIGHTED_AVG_INCLUDE

			"uniform int frame;\n"
			"uniform ssfe _ssfe;\n"
			"uniform image_weighted_avg accum_image;\n"
			"uniform bool use_certainty; /* == exposure->get_confidence() */ \n"

			"void main() {\n"
			"	vec3 value = vec3(0, 0, 0);\n"
			"	vec3 confidence = vec3(0, 0, 0);\n"
			"	if (ssfe_ex_is_honored(_ssfe) && render_is_excluded_r(accum_image.offset, gl_TexCoord[0], frame));\n"
			"	else if (image_weighted_avg_accumulate_norender(accum_image, gl_TexCoord[0]));\n"
			"	else if (use_certainty) {\n"
			"		ssfe_filtered(_ssfe, gl_TexCoord[0], frame, value, confidence,\n"
			"				image_weighted_avg_get_pixel(accum_image, gl_TexCoord[0]),\n"
			"				image_get_pixel(image_weighted_avg_get_weights(accum_image), gl_TexCoord[0]));\n"
			"	} else {\n"
			"		ssfe_filtered(_ssfe, gl_TexCoord[0], frame, value, confidence);\n"
			"	}\n"
			"	image_weighted_avg_accumulate(accum_image, gl_TexCoord[0], frame, value, confidence);\n"
			"}\n";

		this->inv = inv;

		if (inv->is_median())
			accum_image = new image_weighted_median(1, 1, 3);
		else
			accum_image = new image_weighted_simple(1, 1, 3, inv);

		assert (accum_image);

		gpu_program = new gpu::program();
		gpu_shader = new gpu::program::shader(shader_main);
		gpu_program->attach(gpu_shader);
		gpu_program->attach(inv->ssfe());
		gpu_program->attach(accum_image);
		gpu_program->attach((render *) this);
		gpu_program->link();
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

	virtual const image *get_image() const {
		assert (accum_image != NULL);
		return accum_image;
	}

	/*
	 * Definition map.  Unit-depth image whose pixels are nonzero where
	 * the image is defined.
	 */

	virtual const image *get_defined() const {
		assert (accum_image != NULL);
		return accum_image->get_weights();
	}

	/*
	 * Perform the current rendering step.
	 */
	virtual void step() {

		/*
		 * Dynamic invariants are not incrementally updated.
		 */
		if (inv->ssfe()->get_scaled_filter()->is_dynamic())
			return;

		assert (get_step() >= -1);
		if (get_step() == 0) {
			transformation t = align::of(0);

			const image *im = image_rw::open(0);

			ui::get()->rendering();

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
		accum_image->_extend(0, h - 1, 0, w - 1);
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
