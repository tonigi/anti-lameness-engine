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

#ifndef __psf_template_h__
#define __psf_template_h__

#include "psf.h"

/*
 * Point-spread function template.
 */

template <unsigned int rows, unsigned int cols>
class psf_template : public psf {
	const ale_real (&response)[rows][cols][3];
	ale_pos height, width;
public:
	psf_template(ale_pos h, ale_pos w, const ale_real (&_response)[rows][cols][3]) : response(_response) {
		height = h / 2;
		width = w / 2;
	}

	/*
	 * The following four functions indicate filter boundaries.  Filter
	 * support may include everything up to and including the boundaries
	 * specified here.
	 */
	float min_i() const { return -height; }
	float max_i() const { return  height; }
	float min_j() const { return -width; }
	float max_j() const { return  width; }

	/*
	 * Response function
	 *
	 * Get the response to the rectangle bounded by (top, bot, lef, rig).
	 * This function must correctly handle points which fall outside of the
	 * filter support.  The variety of the responding pixel is
	 * provided, in case response is not uniform for all pixels
	 * (e.g. some sensor arrays stagger red, green, and blue
	 * sensors).
	 */
	psf_result operator()(float top, float bot, float lef, float rig, 
			unsigned int variety) const {
		psf_result result;

		if (top < min_i())
			top = min_i();
		if (bot > max_i())
			bot = max_i();
		if (lef < min_j())
			lef = min_j();
		if (rig > max_j())
			rig = max_j();

		int il = (int) floor((top - min_i()) / (max_i() - min_i()) * rows);
		int ih = (int) floor((bot - min_i()) / (max_i() - min_i()) * (rows - 0.001));
		int jl = (int) floor((lef - min_j()) / (max_j() - min_j()) * cols);
		int jh = (int) floor((rig - min_j()) / (max_j() - min_j()) * (cols - 0.001));

		for (int ii = il; ii <= ih; ii++)
		for (int jj = jl; jj <= jh; jj++) {

			float ltop = ((float) ii) / rows * (max_i() - min_i()) + min_i();
			float lbot = ((float) ii + 1) / rows * (max_i() - min_i()) + min_i();
			float llef = ((float) jj) / cols * (max_j() - min_j()) + min_j();
			float lrig = ((float) jj + 1) / cols * (max_j() - min_j()) + min_j();

			if (ltop < top)
				ltop = top;
			if (lbot > bot)
				lbot = bot;
			if (llef < lef)
				llef = lef;
			if (lrig > rig)
				lrig = rig;

			assert (ii >= 0);
			assert (ii < (int) rows);
			assert (jj >= 0);
			assert (jj < (int) cols);

			for (int k = 0; k < 3; k++)
				result.matrix(k, k) += ((lbot - ltop) * (lrig - llef)
					      * response[ii][jj][k]);
		}

		return result;
	}
};

#endif
