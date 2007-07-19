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

#ifndef __zero_h__
#define __zero_h__

/*
 * A zero filter.
 */

class zero : public filter {
public:

	/*
	 * Size of filter support, in number of half-cycles to each side of the
	 * filter center.
	 */
	ale_real support() const {
		return 0;
	}

	virtual int equals(const filter *f) const {
		if (typeid(*f) == typeid(*this))
			return 1;
		return 0;
	}

	/*
	 * Response of filter at point p
	 */
	virtual ale_real response(point p) const {
		return 0;
	}
};
#endif
