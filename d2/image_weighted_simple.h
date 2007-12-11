// Copyright 2002, 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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
	image *colors;
	image *weights;
	
public:
	image_weighted_simple (unsigned int dimy, unsigned int dimx, unsigned int
			depth, invariant *inv, const char *name = "anonymous") 
			: image_weighted_avg(dimy, dimx, depth, name) {
		colors = new_image_ale_real(dimy, dimx, depth);
		weights = new_image_ale_real(dimy, dimx, depth);
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
		if (inv->is_first() && weights->pix(i, j)[0] != 0)
			return 1;
		/*
		 * Weight limit satisfied
		 */
		if (inv->is_avgf() 
		 && weights->pix(i, j)[0] > inv->get_param()
		 && weights->pix(i, j)[1] > inv->get_param()
		 && weights->pix(i, j)[2] > inv->get_param())
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

			/*
			 * Cases independent of the old pixel value and weight
			 * for which the update can be ignored.
			 */

			if (!inv->is_avgx() 
			 && new_weight[k] < render::get_wt())
				continue;

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

			if (inv->is_avgf()
			 && old_weight >= inv->get_param())
				continue;

			/*
			 * Cases independent of the old pixel value for which previous
			 * pixel values can be ignored.
			 */

			if (old_weight == 0
			 || (old_weight < render::get_wt()
			  && !inv->is_avgx())) {
				weights->chan(i, j, k) = new_weight[k];
				colors ->chan(i, j, k) = new_value[k];
				continue;
			}

			if (inv->is_max()) {

				/*
				 * Cases in which the old pixel value can be ignored
				 */

				if (new_value[k] * (ale_real) weights->chan(i, j, k)
				  > (ale_real) colors->chan(i, j, k) * new_weight[k]) {
					weights->chan(i, j, k) = new_weight[k];
					colors-> chan(i, j, k) = new_value[k];
				}
				
				continue;

			} else if (inv->is_min()) {
				/*
				 * Cases in which the old pixel value can be ignored
				 */

				if (new_value[k] * (ale_real) weights->chan(i, j, k)
				  < (ale_real) colors->chan(i, j, k) * new_weight[k]) {
					weights->chan(i, j, k) = new_weight[k];
					colors-> chan(i, j, k) = new_value[k];
				}
				
				continue;
			}


			/*
			 * Update weight and color values.
			 */

			weights->chan(i, j, k) += new_weight[k];
			colors->chan(i, j, k) += new_value[k];
		}
	}


	pixel get_pixel(unsigned int y, unsigned int x) const {
		return colors->get_pixel(y, x) / weights->get_pixel(y, x);
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
