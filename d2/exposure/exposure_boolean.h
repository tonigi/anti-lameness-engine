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
 * exposure_boolean.h: Boolean exposure properties.
 */

#ifndef __exposure_boolean_h__
#define __exposure_boolean_h__

/*
 * Boolean exposure.
 */

class exposure_boolean : public exposure {
public:
	pixel linearize(pixel input) const {
		for (int k = 0; k < 3; k++)
			input[k] = (input[k] == 0) ? 0 : 1;
		return input;
	}
	pixel unlinearize(pixel input) const {
		for (int k = 0; k < 3; k++)
			input[k] = (input[k] == 0) ? 0 : 1;
		return input;
	}
};

#endif
