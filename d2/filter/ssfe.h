// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>

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

#ifndef __ssfe_h__
#define __ssfe_h__

#include "scaled_filter.h"
#include "filter.h"

/*
 * Scaled filter class with exclusion.
 */

class ssfe {
private:
	/*
	 * Honor exclusion?
	 */
	int honor_exclusion;
	scaled_filter *f;

public:

	ssfe(scaled_filter *f, int honor_exclusion) {
		this->honor_exclusion = honor_exclusion;
		this->f = f;
	}

	const scaled_filter *get_scaled_filter() const {
		return f;
	}

	int equals(const ssfe *s) {
		return (honor_exclusion == s->honor_exclusion
		     && f->equals(s->f));
	}
	
	int ex_is_honored() const {
		return honor_exclusion;
	}

	/*
	 * Set the parameters for filtering.
	 */
	void set_parameters(transformation t, const image *im, point offset) const {
		f->set_parameters(t, im, offset);
	}
	void set_parameters(transformation t, transformation s, const image *im) const {
		f->set_parameters(t, s, im);
	}

	/*
	 * Return filtered RESULT and WEIGHT at point P in a coordinate system
	 * specified by the inverse of transformation T based on data taken
	 * from image IM.
	 */
	void filtered(int i, int j, int frame, pixel *result, pixel *weight) const {

		*result = pixel(0, 0, 0);
		*weight = pixel(0, 0, 0);
		
		if (honor_exclusion && render::is_excluded(i, j, frame))
			return;

		f->filtered(i, j, result, weight);
	}
};
#endif
