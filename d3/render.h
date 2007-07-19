// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>

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

/*
 * render.h: A superclass for all rendering classes.
 */

#ifndef __3drender_h__
#define __3drender_h__

#include "point.h"

/*
 * Class render accepts descriptions of scenes.  Subclasses may produce various
 * types of output.
 */

class render {
public:

	/*
	 * Describe a sphere at position p with radius r.  Non-negative 'frame'
	 * less than the total number of frames indicates that the position is
	 * specified in the local space of that frame number.  Otherwise, the
	 * position is specified in world coordinates.
	 */

	virtual void describe(int frame, point p, ale_real r) = 0;

	virtual ~render() {
	}
};

#endif
