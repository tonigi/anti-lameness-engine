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
 * view.h: A render subclass implementing views
 */

#ifndef __view_h__
#define __view_h__

#include "render.h"

class view : public render {

	/*
	 * Z-buffer rendering data
	 *
	 * XXX: there is no immediately obvious reason for depth to use RGB
	 * triples.
	 */

	d2::image *result;
	d2::image *depth;

public:

	/*
	 * Describe a sphere at position p with radius r.  Non-negative 'frame'
	 * less than the total number of frames indicates that the position is
	 * specified in the local space of that frame number.  Otherwise, the
	 * position is specified in world coordinates.
	 */

	virtual void describe(int frame, point p, ale_real r) {
	}

	d2::image *output() {
	}
};

#endif
