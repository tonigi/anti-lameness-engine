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
			point p = t.transform_unscaled(point(i, j));

			if (is_excluded(accum_image->offset(), p, frame_num))
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
			point p = t.transform_unscaled(point(i, j));

			if (is_excluded(point(0, 0), p, frame_num))
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
	 * Merge part of a delta frame with part of the accumulated image using
	 * the specified transformation.
	 */
	void
	_merge(int frame, const image *delta, transformation t) {

		point offset = accum_image->offset();

		assert (accum_image != NULL);
		assert (delta != NULL);

		const filter::ssfe *_ssfe = inv->ssfe();

		_ssfe->set_parameters(t, delta, offset);

		for (unsigned int i = 0; i < accum_image->height(); i++)
		for (unsigned int j = 0; j < accum_image->width(); j++) {

			if (_ssfe->ex_is_honored() && is_excluded(i, j, frame))
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
		}
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
		assert (get_step() >= 0);
		return accum_image->get_colors();
	}

	/*
	 * Definition map.  Unit-depth image whose pixels are nonzero where
	 * the image is defined.
	 */

	virtual const image *get_defined() {
		assert (get_step() >= 0);
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

	void free_memory() {
		delete accum_image;
	}
};

#endif
