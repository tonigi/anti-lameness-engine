// Copyright 2002, 2004 David Hilvert <dhilvert@auricle.dyndns.org>, 
//                                    <dhilvert@ugcs.caltech.edu>

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
 * transformation.h: Represent transformations of the kind q = c(b^-1(p)),
 * where p is a point in the source coordinate system, q is a point in the
 * target coordinate system, b^-1 is a transformation correcting barrel
 * distortion, and c is a transformation of projective or Euclidean type.
 * (Note that ^-1 in this context indicates the function inverse rather than
 * the exponential.)
 */

#ifndef __transformation_h__
#define __transformation_h__

#include "trans_single.h"
#include "trans_multi.h"

typedef trans_multi transformation;

#endif
