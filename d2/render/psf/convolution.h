// Copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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

#ifndef __psf_convolution_h__
#define __psf_convolution_h__

#include "../../point.h"
#include "psf.h"

/*
 * XXX: This doesn't work yet.
 */

/*
 * Point-spread function module.
 *
 * This module implements the convolution (f1 * f2) of point-spread functions f1 and
 * f2.
 */

class convolution : public psf {
	ale_pos _radius;
	psf *f1, *f2;
	ale_real _min_i, _max_i, _min_j, _max_j;

public:
	/*
	 * The following four functions indicate filter boundaries.  Filter
	 * support may include everything up to and including the boundaries
	 * specified here.
	 */
	ale_real min_i() const { return _min_i; }
	ale_real max_i() const { return _max_i; }
	ale_real min_j() const { return _min_j; }
	ale_real max_j() const { return _max_j; }

	/*
	 * Get the number of varieties supported by this PSF.  These usually
	 * correspond to different points in the sensor array.
	 */
	virtual unsigned int varieties() {
		return f1->varieties() * f2->varieties();
	}

	/*
	 * Select the variety appropriate for a given position in the sensor
	 * array.
	 */
	virtual unsigned int select(unsigned int i, unsigned int j) {
		return (f1->select(i, j) * f2->varieties() + f2->select(i, j));
	}

	/*
	 * Response function
	 *
	 * Get the response to the rectangle bounded by (top, bot, lef, rig).
	 * This function must correctly handle points which fall outside of the
	 * filter support.  The variety of the responding pixel is provided, in
	 * case response is not uniform for all pixels (e.g. some sensor arrays
	 * stagger red, green, and blue sensors).
	 */
	psf_result operator()(ale_real top, ale_real bot, ale_real lef, ale_real rig, 
			unsigned int variety) const {
		psf_result result;
		psf_result r1, r2;

		unsigned int v1 = variety / f2->varieties();
		unsigned int v2 = variety % f2->varieties();

		/*
		 * This code uses a rasterized approximation of the filters involved.
		 */

		ale_real vertical_center = (top + bot) / 2;
		ale_real horizontal_center = (lef + rig) / 2;
		ale_real vertical_resolution = bot - top;
		ale_real horizontal_resolution = rig - lef;

		if (!(vertical_resolution > 0
		   && horizontal_resolution > 0))
			return result;  /* zero */

		for (ale_real i = f1->min_i() + (vertical_resolution / 2); 
			   i < f1->max_i() - (vertical_resolution / 2);
			   i += vertical_resolution)
		for (ale_real j = f1->min_j() + (horizontal_resolution / 2);
		           j < f1->max_j() - (horizontal_resolution / 2);
			   j += horizontal_resolution) {

			ale_real t = i - (vertical_resolution / 2);
			ale_real b = i + (vertical_resolution / 2);
			ale_real l = j - (horizontal_resolution / 2);
			ale_real r = j + (horizontal_resolution / 2);
			ale_real vc = vertical_center;
			ale_real hc = horizontal_center;

			r1 = (*f1)(t, b, l, r, v1);
			r2 = (*f2)(vc - b, vc - t, hc - r, hc - l, v2);

			for (int k1 = 0; k1 < 3; k1++)
			for (int k2 = 0; k2 < 3; k2++)
				result.set_matrix(k1, k2, result.get_matrix(k1, k2)
						        + r1.get_matrix(k1, k2) 
							* r2.get_matrix(k1, k2));
		}

		return result;
	}

	convolution(psf *f1, psf *f2) {

		this->f1 = f1;
		this->f2 = f2;

		/*
		 * XXX: I'm fairly sure that this is correct for filters with
		 * zero-centered bounding boxes, and I _think_ it's correct for
		 * other filters also, but I haven't formally proven this.
		 */

		_min_i = f1->min_i() + f2->min_i();
		_min_j = f1->min_j() + f2->min_j();
		_max_i = f1->max_i() + f2->max_i();
		_max_j = f1->max_j() + f2->max_j();
	}
};

#endif
