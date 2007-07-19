// code by HJ Hornbeck, based on code copyright David Hilvert


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

#ifndef __nikon_d50_h__
#define __nikon_d50_h__

#include "../d2.h"
#include "../d2/render/psf/psf.h"

/*
 * Device module for the Nikon D50.
 *
 * Much like the Canon, alignment works best with RAW input created by
 *  dcraw, with the -D/-d and -4 flags enabled. If you get bayer patterns
 *  in the output image, and have Irani-Peleg iteration on, try swapping
 *  the inputs from RAW to processed (but lossless) and ditch the device
 *  setting.
 *
 */

class nikon_d50 {
public:

	/*
	 * Linear colorspace PSF
	 *
	 * This PSF isn't an exact match, instead balancing quality, speed, and
	 *  robustness. While some blur will probably remain, at scale factors
	 *  less than 3 it's easily corrected for in the GIMP (or similar).
	 *  Larger diameters perform more sharpening, but run time increases
	 *  by the square of the diameter. Finally, the PSF should vary by
	 *  lens and focal length, but this mild one works relatively well
	 *  across a wide range of lenses.
	 */

	static d2::psf *lpsf() {

		return d2::psf_parse::get( 1, "circle=1.3^circle=1.3" );
		}

	/*
	 * Non-linear colorspace PSF
	 * 
	 * The linear PSF works perfectly with RAW input, and even performs
	 *  well if I substitute images with non-linear colourspace. I
	 *  haven't found a need for this, yet.
	 */

	static d2::psf *nlpsf() {

		return NULL;
		}

	/*
	 * Exposure
	 *
	 * The defaults for the Canon work well enough.
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
	 *
	 * Copy/Paste from the Canon again.
	 */

	static ale_pos view_angle() {
		fprintf(stderr, "\n\n*** Error: tried to obtain view angle for a sensor device. ***\n\n");
		exit(1);
		// return 30;
	}

	/*
	 * Bayer pattern
	 *
	 * Adding this function because this info should be here, instead
	 *  of hard-wired into the argument parsing code.
	 */

	static unsigned int bayer() {
		return IMAGE_BAYER_BGRG;
	}
};

#undef LPSF_ROWS
#undef LPSF_COLS
#undef NLPSF_ROWS
#undef NLPSF_COLS

#endif
