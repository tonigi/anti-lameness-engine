// Copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 3 of the License,
    or (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef __ov7620_raw_linear_h__
#define __ov7620_raw_linear_h__

#include "../d2.h"

/*
 * Device module for the OmniView OV7620 chip.
 *
 * This module is designed for use with the Linux 2.6.x driver patched with
 * code available at:
 *
 * http://auricle.dyndns.org/ov7620/
 *
 */

class ov7620_raw_linear {
public:

	/*
	 * Linear colorspace PSF
	 */

	class lpsf : public d2::box {
	public:
		/*
		 * A box filter of diameter ~0.6325 results in an optical fill
		 * factor of 40%.  Experimentation shows that filters of about 
		 * this size perform well with this sensor, and the fill factor
		 * value matches data available on the web about this sensor:
		 *
		 * http://www.stanford.edu/class/ee109/reference/camera/DS-OV7620-1.3.pdf
		 */

		lpsf() : d2::box (0.3162) {
		}
	};

	/*
	 * Exposure
	 */

	class exposure : public d2::exposure_default {
		d2::pixel linearize(d2::pixel input) const {
			ale_real lo = 16 / (ale_real) 255;
			ale_real scale = 255 / (ale_real) 224;
			return (input - d2::pixel(lo, lo, lo)) * scale * get_multiplier();
		}
		d2::pixel unlinearize(d2::pixel input) const {
			ale_real lo = 16 / (ale_real) 255;
			ale_real scale = 224 / (ale_real) 255;
			return (input / get_multiplier()) * scale + d2::pixel(lo, lo, lo);
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
