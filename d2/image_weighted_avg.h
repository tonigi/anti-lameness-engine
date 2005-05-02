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
 * image_weighted_avg.h: Image representing a weighted average of inputs.
 */

#ifndef __image_weighted_avg_h__
#define __image_weighted_avg_h__

#include "image_ale_real.h"
#include "exposure/exposure.h"
#include "point.h"
#include "image.h"

class image_weighted_avg : public image {
private:
	void trigger(pixel multiplier) {
		assert(0);
	}

public:
	image_weighted_avg (unsigned int dimy, unsigned int dimx, unsigned int
			depth, char *name = "anonymous") 
			: image(dimy, dimx, depth, name, NULL) {
	}

	virtual ~image_weighted_avg() {
	}

	spixel &pix(unsigned int y, unsigned int x) {
		assert(0);
		static spixel foo;
		return foo;
	}

	const spixel &pix(unsigned int y, unsigned int x) const {
		assert(0);
		static spixel foo;
		return foo;
	}

	ale_real &chan(unsigned int y, unsigned int x, unsigned int k) {
		assert(0);
		static ale_real foo;
		return foo;
	}

	const ale_real &chan(unsigned int y, unsigned int x, unsigned int k) const {
		assert(0);
		static ale_real foo;
		return foo;
	}

	/*
	 * Make a new image suitable for receiving scaled values.
	 */
	virtual image *scale_generator(int height, int width, int depth, char *name) const {
		return new image_ale_real(height, width, depth, name, _exp);
	}

	/*
	 * Pre-transformation check for whether an area should be skipped.
	 * Takes image weights as an argument.
	 */
	virtual int accumulate_norender(int i, int j) = 0;

	/*
	 * Accumulate pixels 
	 */
	virtual void accumulate(int i, int j, int f, pixel new_value, pixel new_weight) = 0;

	/*
	 * Get color map
	 */
	virtual image *get_colors() = 0;

	/*
	 * Get weight map
	 */
	virtual image *get_weights() = 0;
};

#endif
