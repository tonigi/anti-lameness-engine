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

#ifndef __filterclass_h__
#define __filterclass_h__

/*
 * A superclass for all filtering classes.
 */

class filter {
public:

	/*
	 * Size of filter support, in number of half-cycles to each side of the
	 * filter center.
	 */
	virtual ale_real support() const = 0;

	/*
	 * Response of filter at point p
	 */
	virtual ale_real response(point p) const = 0;

	virtual int equals(const filter *f) const = 0;

	virtual ~filter(){
	}
};
#endif
