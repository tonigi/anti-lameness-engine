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
 * image_weighted_median.h: Image representing a weighted median of inputs.
 */

#ifndef __image_weighted_median_h__
#define __image_weighted_median_h__

#include "exposure/exposure.h"
#include "point.h"
#include "image.h"

class image_weighted_median : public image_weighted_avg {
private:

	/*
	 * XXX: This storage approach may have poor cache characteristics.
	 * It might be better to localize elements having identical spatial
	 * coordinates.
	 */
	image_ale_real **colors;
	image_ale_real **weights;

public:
	image_weighted_median (unsigned int dimy, unsigned int dimx, unsigned int
			depth, char *name = "anonymous") 
			: image_weighted_avg(dimy, dimx, depth, name) {

		colors  = (image_ale_real **) malloc(image_rw::count() * sizeof(image_ale_real *));
		weights = (image_ale_real **) malloc(image_rw::count() * sizeof(image_ale_real *));

		assert(colors);
		assert(weights);

		if (!colors || !weights) {
			fprintf(stderr, "Could not allocate memory for image data.\n");
			exit(1);
		}

		for (unsigned int f = 0; f < image_rw::count(); f++) {
			colors[f] = new image_ale_real(dimy, dimx, depth);
			weights[f] = new image_ale_real(dimy, dimx, depth);

			assert(colors[f]);
			assert(weights[f]);

			if (!colors[f] || !weights[f]) {
				fprintf(stderr, "Could not allocate memory for image data.\n");
				exit(1);
			}
		}
	}

	virtual ~image_weighted_median() {
		for (unsigned int f = 0; f < image_rw::count(); f++) {
			delete colors[f];
			delete weights[f];
		}

		free(colors);
		free(weights);
	}

	/*
	 * Extend the image area to the top, bottom, left, and right,
	 * initializing the new image areas with black pixels.  Negative values
	 * shrink the image.
	 */
	void extend(int top, int bottom, int left, int right) {

		for (unsigned int f = 0; f < image_rw::count(); f++) {
			colors[f]->extend(top, bottom, left, right);
			weights[f]->extend(top, bottom, left, right);
		}

		_dimx = colors[0]->width();
		_dimy = colors[0]->height();
		_offset = colors[0]->offset();
	}

	int accumulate_norender(int i, int j) {
		return 0;
	}

	/*
	 * Perform insertion sort on the arrays, where sort is by color.
	 */
	void accumulate(int i, int j, int f, pixel new_value, pixel new_weight) {
		for (unsigned int k = 0; k < 3; k++) {

			if (new_weight[k] <= 0)
				continue;

			for (unsigned int ff = 0; ff < image_rw::count(); ff++) {
				if (ff == image_rw::count() - 1) {
					colors[ff]->pix(i, j)[k] = new_value[k];
					weights[ff]->pix(i, j)[k] += new_weight[k];
					return;
				}
				if (colors[ff]->pix(i, j)[k] > new_value[k]) {
					for (unsigned int fff = image_rw::count() - 1; fff > ff; fff--) {
						weights[fff]->pix(i, j)[k] = weights[fff - 1]->pix(i, j)[k] + new_weight[k];
						colors[fff]->pix(i, j)[k]  = colors[fff - 1]->pix(i, j)[k];
					}
					colors[ff]->pix(i, j)[k] = new_value[k];
					weights[ff]->pix(i, j)[k] = new_weight[k];
					if (ff > 0)
						weights[ff]->pix(i, j)[k] += weights[ff - 1]->pix(i, j)[k];

					return;
				}
				if (colors[ff]->pix(i, j)[k] == new_value[k]) {
					for (unsigned int fff = ff; fff < image_rw::count(); fff++)
						weights[fff]->pix(i, j)[k] += new_weight[k];
					return;
				}
				if ((ff == 0 && weights[ff] == 0)
				 || (ff >  0 && weights[ff] == weights[ff - 1])) {
					colors[ff]->pix(i, j)[k] = new_value[k];
					for (unsigned int fff = ff; fff < image_rw::count(); fff++)
						weights[fff]->pix(i, j)[k] += new_weight[k];
				}
			}
		}
	}

	/*
	 * XXX: This is inefficient in cases where only one channel is desired.
	 */
	pixel get_pixel(unsigned int y, unsigned int x) const {
		pixel result;

		for (int k = 0; k < 3; k++) {
			ale_real total = weights[image_rw::count() - 1]->chan(y, x, k);
			ale_real accum = 0;

			if (total == 0)
				return pixel::zero();

			for (unsigned int f = 0; f < image_rw::count(); f++) {
				accum += weights[f]->chan(y, x, k);
				if (accum > total / 2) {
					result[k] = colors[f]->chan(y, x, k);
					break;
				} else if (accum == total / 2) {
					result[k] = (colors[f]->chan(y, x, k)
					           + colors[f + 1]->chan(y, x, k)) / 2;
					break;
				}
			}
		}

		return result;
	}

	image *get_weights() {
		return weights[image_rw::count() - 1];
	}

	image *get_colors() {
		return this;
	}

};

#endif
