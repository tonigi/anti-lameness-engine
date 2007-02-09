// Copyright 2002, 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                          <dhilvert@ugcs.caltech.edu>

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
 * image_weighted_simple.h: Image representing a weighted average of inputs.
 * Covers simple cases that require space constant with frame count.
 */

#ifndef __image_weighted_simple_h__
#define __image_weighted_simple_h__

#include "image_weighted_avg.h"
#include "render/invariant.h"

class image_weighted_simple : public image_weighted_avg {
private:
	invariant *inv;
	image_ale_real *colors;
	image_ale_real *weights;
	
public:
	image_weighted_simple (unsigned int dimy, unsigned int dimx, unsigned int
			depth, invariant *inv, char *name = "anonymous") 
			: image_weighted_avg(dimy, dimx, depth, name) {
		colors = new image_ale_real(dimy, dimx, depth);
		weights = new image_ale_real(dimy, dimx, depth);
		this->inv = inv;
	}

	virtual ~image_weighted_simple() {
		delete colors;
		delete weights;
	}

	/*
	 * Extend the image area to the top, bottom, left, and right,
	 * initializing the new image areas with black pixels.  Negative values
	 * shrink the image.
	 */
	void extend(int top, int bottom, int left, int right) {
		colors->extend(top, bottom, left, right);
		weights->extend(top, bottom, left, right);

		_dimx = colors->width();
		_dimy = colors->height();
		_offset = colors->offset();
	}

	/*
	 * Pre-transformation check for whether an area should be skipped.
	 * Takes image weights as an argument.
	 */
	int accumulate_norender(int i, int j) {
		/*
		 * Initial value
		 */
		if (inv->is_first() && weights->get_pixel(i, j)[0] != 0)
			return 1;

		return 0;
	}

	/*
	 * Accumulate pixels 
	 */
	void accumulate(int i, int j, int f, pixel new_value, pixel new_weight) {

		/*
		 * Perform operations separately for each channel
		 */
		for (unsigned int k = 0; k < 3; k++) {

			const double min_weight = 0.001;

			/*
			 * Cases independent of the old pixel value and weight
			 * for which the update can be ignored.   XXX: the
			 * minimum weight used here should probably be less
			 * than or equal to the certainty floor.
			 */

			if (fabs(new_weight[k]) < min_weight
			 || (!inv->is_avg()
			  && new_weight[k] < render::get_wt())) {
				continue;
			}

			/*
			 * Cases independent of the old pixel value and weight for which
			 * previous pixel values can be ignored. 
			 */

			if (inv->is_last() && new_weight[k] >= render::get_wt()) {
				colors->chan(i, j, k) = new_value[k];
				weights->chan(i, j, k) = new_weight[k];
				continue;
			}

			/*
			 * Obtain the old pixel weight.
			 */

			ale_real old_weight = weights->chan(i, j, k);

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
				weights->chan(i, j, k) = new_weight[k];
				colors ->chan(i, j, k) = new_value[k];
				continue;
			}

			/*
			 * Obtain the old pixel value
			 */

			ale_real old_value = colors->chan(i, j, k);

			/*
			 * Cases in which the old pixel value can be ignored
			 */

			if ((inv->is_max()
			  && new_value[k] > old_value)
			 || (inv->is_min()
			  && new_value[k] < old_value)) {
				weights->chan(i, j, k) = new_weight[k];
				colors-> chan(i, j, k) = new_value[k];
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

			if (fabs(updated_weight) < min_weight) {
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

			weights->chan(i, j, k) = updated_weight;

			/*
			 * Update the channel
			 */

			colors->chan(i, j, k) = (old_weight * colors->chan(i, j, k)
					      + new_weight[k] * new_value[k])
				             / updated_weight;
		}
	}


	pixel get_pixel(unsigned int y, unsigned int x) const {
		return colors->get_pixel(y, x);
	}

	image *get_weights() {
		return weights;
	}

	image *get_colors() {
		assert(0);
		return colors;
	}
};

#endif
