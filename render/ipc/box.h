// Copyright 2003 David Hilvert <dhilvert@ugcs.caltech.edu>

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

#ifndef __ipc_box_h__
#define __ipc_box_h__

#include <assert.h>
#include "../../point.h"

/*
 * Module for Irani-Peleg image reconstruction.  
 *
 * This module implements the box filter.
 */

class box {
	double _radius;
public:
	/*
	 * Result type.
	 */
	class ipc_result {
		friend class box;
		double response;
	public:

		/*
		 * Response intensity on channel K2 resulting from
		 * stimulus on channel K1.  Channel R=0, G=1, B=2.
		 */
		double operator()(int k1, int k2) {
			if (k1 == k2)
				return response;
			else
				return 0;
		}
	};

	/*
	 * The following four functions indicate filter boundaries.  Filter
	 * support may include everything up to and including the boundaries
	 * specified here.
	 */
	float min_i() { return -_radius; }
	float max_i() { return  _radius; }
	float min_j() { return -_radius; }
	float max_j() { return  _radius; }

	/*
	 * Get the response to the rectangle bounded by (top, bot, lef, rig).
	 * This function must correctly handle points which fall outside of the
	 * filter support.  The position (i, j) of the responding pixel is
	 * provided, in case response is not uniform for all pixels (e.g. some
	 * sensor arrays stagger red, green, and blue sensors).
	 */
	struct ipc_result operator()(float top, float bot, float lef, float
			rig, int i, int j) {
		struct ipc_result result;

		if (top < min_i())
			top = min_i();
		if (bot > max_i())
			bot = max_i();
		if (lef < min_j())
			lef = min_j();
		if (rig > max_j())
			rig = max_j();

		if (bot > top && rig > lef)
			result.response = (bot - top) * (rig - lef);
		else
			result.response = 0;

		return result;
	}

	void radius(double r) {
		this->_radius = r;
	}
};

#endif
