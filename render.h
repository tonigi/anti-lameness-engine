// Copyright 2002 David Hilvert <dhilvert@ugcs.caltech.edu>

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

#ifndef __render_h__
#define __render_h__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "gpt.h"
#include "image.h"
#include "point.h"

class render {
public:

	/*
	 * Result of rendering.
	 */

	virtual const image *get_image() = 0;

	/*
	 * Definition map.  Unit-depth image whose pixels are nonzero where
	 * the image is defined.
	 */

	virtual const image_weights *get_defined() = 0;

	/*
	 * Perform rendering steps requiring no frames beyond frame N.
	 */

	virtual void operator()(int n) = 0;

	/*
	 * Perform any final rendering steps.  Return a non-zero value if
	 * anything changed.
	 */

	virtual int operator()() {
		return 0;
	}

};

#endif
