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
	image *accum_image;
	image *accum_weight;
	invariant *inv;

	/*
	 * Set extents of image and weight according to a new image to be
	 * merged.  This function should remove only superfluous undefined
	 * areas.
	 */

	void set_extents_by_map(unsigned int frame_num, transformation t) {

		assert (accum_weight != NULL);
		assert (accum_image  != NULL);

		ale_pos extend_offset_i = accum_image->offset()[0];
		ale_pos extend_offset_j = accum_image->offset()[1];

		assert (accum_weight->offset()[0] == extend_offset_i);
		assert (accum_weight->offset()[1] == extend_offset_j);

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
		accum_weight->extend(extend_top, extend_bottom, 
			extend_left, extend_right);
	}	

	void increase_extents_by_map(unsigned int frame_num, transformation t) {

		assert (accum_weight != NULL);
		assert (accum_image  != NULL);

		ale_pos extend_offset_i = accum_image->offset()[0];
		ale_pos extend_offset_j = accum_image->offset()[1];

		assert (accum_weight->offset()[0] == extend_offset_i);
		assert (accum_weight->offset()[1] == extend_offset_j);

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
		accum_weight->extend(extend_top, extend_bottom, 
			extend_left, extend_right);
	}	

	/*
	 * Pre-transformation check for whether an area should be skipped.
	 * Takes image weights as an argument.
	 */
	int accumulate_norender(int i, int j, const image *weight) {

		/*
		 * Initial value
		 */
		if (inv->is_first() && weight->get_pixel(i, j)[0] != 0)
			return 1;

		return 0;
	}

	/*
	 * Accumulate pixels into the accumulated image.
	 */

	void accumulate(image *accum, image *weight, int i, int j, int f, pixel new_value, pixel new_weight) {

		/*
		 * Perform operations separately for each channel
		 */
		for (unsigned int k = 0; k < 3; k++) {

			/*
			 * Cases independent of the old pixel value and weight
			 * for which the update can be ignored.   XXX: the
			 * first bound should probably be less than or equal to
			 * the certainty floor.
			 */

			if (fabs(new_weight[k]) < 0.001
			 || (!inv->is_avg()
			  && new_weight[k] < render::get_wt())) {
				continue;
			}

			/*
			 * Cases independent of the old pixel value and weight for which
			 * previous pixel values can be ignored. 
			 */

			if (inv->is_last() && new_weight[k] >= render::get_wt()) {
				accum->chan(i, j, k) = new_value[k];
				weight->chan(i, j, k) = new_weight[k];
				continue;
			}

			/*
			 * Obtain the old pixel weight.
			 */

			ale_real old_weight = weight->chan(i, j, k);

			/*
			 * Cases independent of the old pixel value for which the
			 * update can be ignored.
			 */

			if (old_weight >= render::get_wt()
			 && inv->is_first())
				continue;

			/*
			 * Cases independent of the old pixel value for which previous
			 * pixel values can be ignored.
			 */

			if (old_weight == 0
			 || (old_weight < render::get_wt()
			  && !inv->is_avg())) {
				weight->chan(i, j, k) = new_weight[k];
				accum ->chan(i, j, k) = new_value[k];
				continue;
			}

			/*
			 * Obtain the old pixel value
			 */

			ale_real old_value = accum->chan(i, j, k);

			/*
			 * Cases in which the old pixel value can be ignored
			 */

			if ((inv->is_max()
			  && new_value[k] > old_value)
			 || (inv->is_min()
			  && new_value[k] < old_value)) {
				weight->chan(i, j, k) = new_weight[k];
				accum-> chan(i, j, k) = new_value[k];
				continue;
			}

			/*
			 * Cases in which the new pixel value can be ignored
			 */

			if ((inv->is_max()
			  && old_value > new_value[k])
			 || (inv->is_min()
			  && old_value < new_value[k])) {
				continue;
			}

			/*
			 * Update the weight
			 */

			ale_real updated_weight = old_weight + new_weight[k];

			if (fabs(updated_weight) < 0.001) {
				/*
				 * XXX: Give up.  Because of the way we
				 * represent weights and values, we can't
				 * handle small weights (although we could
				 * probably handle smaller weights than we
				 * currently handle).  This could be fixed
				 * by always storing unnormalized values
				 * separately from weights, and dividing only
				 * once, just prior to image output.
				 */
				continue;
			}

			weight->chan(i, j, k) = updated_weight;

			/*
			 * Update the channel
			 */

			accum->chan(i, j, k) = (old_weight * accum->chan(i, j, k)
					      + new_weight[k] * new_value[k])
				             / updated_weight;
		}
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
		assert (accum_weight != NULL);

		const filter::ssfe *_ssfe = inv->ssfe();

		_ssfe->set_parameters(t, delta, offset);

		for (unsigned int i = 0; i < accum_image->height(); i++)
		for (unsigned int j = 0; j < accum_image->width(); j++) {

			if (_ssfe->ex_is_honored() && is_excluded(i, j, frame))
				continue;

			if (accumulate_norender(i, j, accum_weight))
				continue;
			
			/*
			 * Pixel value to be merged, and the associated
			 * confidence
			 */

			pixel value, confidence;

			_ssfe->filtered(i, j, frame, &value, &confidence);

			accumulate(accum_image, accum_weight, i, j, frame, value, confidence);
		}
	}

public:

	/*
	 * Constructor
	 */
	incremental(invariant *inv) {
		this->inv = inv;
		accum_weight = NULL;
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
		return accum_image;
	}

	/*
	 * Definition map.  Unit-depth image whose pixels are nonzero where
	 * the image is defined.
	 */

	virtual const image *get_defined() {
		assert (get_step() >= 0);
		return accum_weight;
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

			accum_image = new image_ale_real(1, 1, 3);

			accum_weight = new image_ale_real(1, 1, 3);

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
		delete accum_weight;
	}
};

#endif
