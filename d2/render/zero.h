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

/*
 * zero.h: A renderer for the zero filter.
 */

#ifndef __render_zero_h__
#define __render_zero_h__

#include "incremental.h"
#include "../image_zero.h"

class zero : public incremental {
public:
	zero(invariant *inv) : incremental(inv) {
		assert (typeid(*inv->ssfe()->get_scaled_filter()->get_filter()) == typeid(filter::zero));
	}

	/*
	 * Perform the current rendering step.
	 */
	virtual void step() {
		assert (get_step() >= -1);
		if (get_step() == 0) {
			transformation t = align::of(0);

			const image *im = image_rw::open(0);

			/*
			 * XXX: This approach to trimming pixels is probably a
			 * bit too aggressive.  If it is changed, be sure to
			 * also change the corresponding lines in
			 * incremental.h.
			 */

			unsigned int trim_size = (int) ceil(get_scale_factor()) - 1;

			accum_image = new image_zero((int) floor(im->height() * get_scale_factor()) - trim_size,
						     (int) floor(im->width()  * get_scale_factor()) - trim_size, 3);

			accum_weight = new image_zero(accum_image->height(),
					accum_image->width(), accum_image->depth());

			if (is_extend())
				increase_extents(t);

			image_rw::close(0);
		} else if (align::match(get_step())) {
			transformation t = align::of(get_step());
			if (is_extend())
				increase_extents(t);
		}
	}

};

#endif
