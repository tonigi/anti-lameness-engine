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

#ifndef __ale_pos_h__
#define __ale_pos_h__

#define SINGLE 1
#define DOUBLE 2

/*
 * Real-valued type used to represent coordinates in an image domain.
 */

#if ALE_COORDINATES == SINGLE

typedef float ale_pos;

#define ALE_POS_PRECISION_STRING "SINGLE"

#elif ALE_COORDINATES == DOUBLE

typedef double ale_pos;

#define ALE_POS_PRECISION_STRING "DOUBLE"

#else

#warning Unknown positional precision in ale_pos.h: Choosing PPRECISION=SINGLE.

typedef float ale_pos;

#define ALE_POS_PRECISION_STRING "SINGLE"

#endif

#undef SINGLE
#undef DOUBLE
#undef HALF

#endif
