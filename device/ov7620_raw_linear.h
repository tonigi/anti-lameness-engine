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

#ifndef __ov7620_raw_linear_h__
#define __ov7620_raw_linear_h__

#include "../d2.h"

/*
 * Device module for the IBM PC Camera Pro (IBM VGA Camera Model XVP610; FCC ID
 * KSX-X9902).
 *
 * This module is designed for use with the Linux 2.6.x driver patched with
 * code available at:
 *
 * http://auricle.dyndns.org/xvp610/
 *
 */

/*
 * PSF Data
 */

#define LPSF_ROWS 3
#define LPSF_COLS 3

static const ale_real ov7620_raw_linear_lpsf_calibrated_response[LPSF_ROWS][LPSF_COLS][3] = {
	{
		{ 1.6, 4.8, 3.2 },
		{ 6.4, 6.4, 8 },
		{ 0, 0, 0 }
	}, {
		{ 1.6, -6.4, -1.6 },
		{ 58.4, 61.6, 80.8 },
		{ 0, 0, 0 }
	}, {
		{ -1.6, -1.6, -1.6 },
		{ 0, 0, 0 },
		{ 0, 0, 0 }
	}
};

#define NLPSF_ROWS 1
#define NLPSF_COLS 5

static const ale_real ov7620_raw_linear_nlpsf_calibrated_response[NLPSF_ROWS][NLPSF_COLS][3] = {
	{
#if 0
		/* Filter 12 */
		{ -2.9, -4.9, -2.7 },
		{ -11.4, -5.4, -9.1 },
		{ 40.9, 39.4, 34.2 },
		{ -3.8, -3, -3.6 },
		{ -4.1, -4.5, -3.2 }
#endif
		/* Filter 11 */
		{ -6.1, -4.9, -5.9 },
		{ -8.2, -5.4, -5.9 },
		{ 36.1, 34.6, 31 },
		{ -3.8, -3,  -3.6 },
		{ -4.1, -4.5, -3.2 }
	}
};

class ov7620_raw_linear {
public:

	/*
	 * Linear colorspace PSF
	 */

	class lpsf : public d2::box {
	public:
		lpsf() : d2::box (0.5) {
			fprintf(stderr, "\n\n*** Error: Undefined device-specific response. ***\n\n");
			exit(1);
		}
	};

	/*
	 * Non-linear colorspace PSF
	 */

	class nlpsf : public d2::psf_template<NLPSF_ROWS, NLPSF_COLS> {
	public:
		nlpsf() : d2::psf_template<NLPSF_ROWS, NLPSF_COLS> (1, 10, ov7620_raw_linear_nlpsf_calibrated_response) {
			fprintf(stderr, "\n\n*** Error: Undefined device-specific response. ***\n\n");
			exit(1);
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
