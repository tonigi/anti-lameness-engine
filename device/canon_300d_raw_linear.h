// Copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the License,
    or (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef __canon_300d_raw_linear_h__
#define __canon_300d_raw_linear_h__

#include "../d2.h"

/*
 * Device module for the Canon 300D (Digital Rebel).
 *
 * CRW files are not directly supported; use Dave Coffin's dcraw program
 * with arguments -d -4 to produce a raw linear PPM file.
 *
 * http://www.cybercom.net/~dcoffin/dcraw/
 */

class canon_300d_raw_linear {
public:

	/*
	 * Linear colorspace PSF
	 */

	class lpsf : public d2::box {
	public:
		/*
		 * A box filter of diameter 1 results in an optical fill factor
		 * of 100%.  This probably isn't exactly right; there is lots
		 * of room for experimentation.
		 */

		lpsf() : d2::box (0.5) {
		}
	};

	/*
	 * Exposure
	 */

	class exposure : public d2::exposure_default {
		d2::pixel linearize(d2::pixel input) const {
			return input * get_multiplier();
		}
		d2::pixel unlinearize(d2::pixel input) const {
			return input / get_multiplier();
		}
	};

	/*
	 * View Angle
	 */

	static ale_pos view_angle() {
		fprintf(stderr, "\n\n*** Error: tried to obtain view angle for a sensor device. ***\n\n");
		exit(1);
		// return 30;
	}
};

#undef LPSF_ROWS
#undef LPSF_COLS
#undef NLPSF_ROWS
#undef NLPSF_COLS

#endif
