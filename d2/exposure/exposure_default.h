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
 * exposure_default.h: Default exposure properties.
 */

#ifndef __exposure_default_h__
#define __exposure_default_h__

/*
 * The default exposure is modeled after the simple power transfer function
 * described in
 *
 * http://netpbm.sourceforge.net/doc/pnmgamma.html
 *
 * Note: optimizations in d2/image_rw.h depend on the details of this function.
 */

class exposure_default : public exposure {
public:
	pixel linearize(pixel input) const {
		return ppow(input, 1/0.45) * get_multiplier();
	}
	pixel unlinearize(pixel input) const {
		return ppow(input / get_multiplier(), 0.45);
	}
};

#endif
