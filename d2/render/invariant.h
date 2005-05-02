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

#ifndef __invariant_h__
#define __invariant_h__

#include "../filter.h"

/*
 * Class for incremental renderer invariants.
 */

#define min 0
#define max 1
#define avg 2
#define first 3
#define last 4
#define median 5

class invariant {
public:
	int type;
	filter::ssfe *s;

	invariant(filter::ssfe *s) {
		this->s = s;
		type = 2;
	}
	int equals(const invariant *i) const {
		return (i->type == type
		     && s->equals(i->ssfe()));
	}
	const filter::ssfe *ssfe() const {
		return s;
	}
	int is_max() const {
		return type == max;
	}
	int is_min() const {
		return type == min;
	}
	int is_avg() const {
		return type == avg;
	}
	int is_first() const {
		return type == first;
	}
	int is_last() const {
		return type == last;
	}
	int is_median() const {
		return type == median;
	}
	void set_max() {
		type = max;
	}
	void set_min() {
		type = min;
	}
	void set_avg() {
		type = avg;
	}
	void set_first() {
		type = first;
	}
	void set_last() {
		type = last;
	}
	void set_median() {
		type = median;
	}
	~invariant() {
		delete s;
	}
};

#undef min
#undef max
#undef avg
#undef first
#undef last

#endif
