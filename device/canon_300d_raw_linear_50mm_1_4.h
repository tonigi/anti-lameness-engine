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

#ifndef __canon_300d_raw_linear_50mm_1_4_h__
#define __canon_300d_raw_linear_50mm_1_4_h__

#include "../d2.h"

/*
 * Extend the Canon 300D (Digital Rebel) device with an EF 50mm f/1.4 lens.
 */

class canon_300d_raw_linear_50mm_1_4 : public canon_300d_raw_linear {
public:
	/*
	 * View Angle
	 *
	 * According to:
	 *
	 * http://www.acapixus.dk/photography/angle_of_view.htm
	 */

	static ale_pos view_angle() {
		return 31 * M_PI / 180;
	}
};

#endif
