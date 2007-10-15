// code by HJ Hornbeck, based on code copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//				    <dhilvert@ugcs.caltech.edu>

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

#ifndef __psf_gauss_h__
#define __psf_gauss_h__

#include "../../point.h"
#include "psf.h"

/*
 * A Gaussian point-spread function. It's aimed at duplicating the most common type
 *  of blurring in many optical systems. It is uniform across the entire image, so
 *  it can't correct for poor focus at the edges.
 */

#define D2_GAUSS_CUTOFF ((ale_real) 2.0)

class gauss : public psf {
	ale_real sigma;		  // radius, in pixels per standard deviation

	/*
	 * Disabled the following definition because some compilers may not be
	 * able to handle static const definitions within a class (and because
	 * the C++ specification may disallow such for non-integral types,
	 * anyway).
	 *
	 * -- dhilvert@auricle.dyndns.org  18-May-2007
	 */

//	static const ale_pos cutoff = 2;	// standard deviations before we cut off

	// helper variables
	ale_real radius;
	ale_pos sigma_premult;
public:

	/*
	 * The following four functions indicate filter boundaries.  Filter
	 * support may include everything up to and including the boundaries
	 * specified here.
	 */
	inline ale_real max_i() const { return radius; }
	inline ale_real min_i() const { return -max_i(); } // we're symmetrical, so it works!
	inline ale_real min_j() const { return -max_i(); }
	inline ale_real max_j() const { return max_i(); }

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

		// calculate some needed values
		ale_pos area_premult = (bot - top) * (rig - lef) / 25;
		float vert_step = (bot - top) / 4;
		float horiz_step = (rig - lef) / 4;
		float total = 0;


		// determine the final value by simple sampling:
		for (float i = top; i < bot + vert_step / 2; i += vert_step)
		for (float j = lef; j < rig + horiz_step / 2; j += horiz_step) {

			// calculate radius for given sample
			ale_pos r = sqrt( i*i + j*j );

			if ( r < radius ) // calculate gaussian falloff
				total += exp( -r * r * sigma_premult ) ;
			// outside our radius? must be 0...
			}

		// adjust for point sampling and area
		total *= area_premult;

		// pre-fill the colour result matrix
		for (int k = 0; k < 3; k++)
			result.matrix(k, k) = 0;

		// fill in the results
		for (int k = 0; k < 3; k++)
			result.matrix(k, k) = total;
	       
		return result;
	}

	/*
	 * Our glorious constructor
	 */
	gauss(ale_real sig) {

		sigma = sig;

		// fill in our helper variables
		radius = sigma * D2_GAUSS_CUTOFF;
		sigma_premult = 1 / (sigma * sigma);
		}
};

#undef D2_GAUSS_CUTOFF

#endif

