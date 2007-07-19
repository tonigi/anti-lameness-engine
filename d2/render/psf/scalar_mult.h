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

#ifndef __psf_scalar_mult_h__
#define __psf_scalar_mult_h__

#include "../../point.h"
#include "psf.h"

/*
 * Point-spread function module.
 *
 * This module implements the scalar_mult (f1 * f2) of point-spread functions f1 and
 * f2.
 */

class scalar_mult : public psf {
	ale_pos _radius;
	psf *f;
	ale_real scalar;
	float _min_i, _max_i, _min_j, _max_j;

public:
	/*
	 * The following four functions indicate filter boundaries.  Filter
	 * support may include everything up to and including the boundaries
	 * specified here.
	 */
	float min_i() const { return _min_i; }
	float max_i() const { return _max_i; }
	float min_j() const { return _min_j; }
	float max_j() const { return _max_j; }

	/*
	 * Get the number of varieties supported by this PSF.  These usually
	 * correspond to different points in the sensor array.
	 */
	virtual unsigned int varieties() {
		return f->varieties();
	}

	/*
	 * Select the variety appropriate for a given position in the sensor
	 * array.
	 */
	virtual unsigned int select(unsigned int i, unsigned int j) {
		return f->select(i, j);
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
	psf_result operator()(float top, float bot, float lef, float rig, 
			unsigned int variety) const {
		psf_result result;
		psf_result r;

		r = (*f)(top, bot, lef, rig, variety);

		for (int k1 = 0; k1 < 3; k1++)
		for (int k2 = 0; k2 < 3; k2++)
			result.set_matrix(k1, k2, scalar * r.get_matrix(k1, k2));

		return result;
	}

	scalar_mult(ale_real s, psf *f) {
		this->scalar = s;
		this->f = f;

		_min_i = f->min_i();
		_min_j = f->min_j();
		_max_i = f->max_i();
		_max_j = f->max_j();
	}
};

#endif
