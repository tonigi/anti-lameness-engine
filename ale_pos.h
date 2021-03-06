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

#include "ale_fixed.h"

#define SINGLE 1
#define DOUBLE 2
#define FIXED16 3
#define FIXED32 4

/*
 * Real-valued type used to represent coordinates in an image domain.
 */

#if ALE_COORDINATES == SINGLE

typedef float ale_pos;

#define ALE_POS_PRECISION_STRING "SINGLE"

#define ale_pos_disable_casting()
#define ale_pos_enable_casting() 
#define ale_pos_casting_status() 1

#elif ALE_COORDINATES == DOUBLE

typedef double ale_pos;

#define ALE_POS_PRECISION_STRING "DOUBLE"

#define ale_pos_disable_casting()
#define ale_pos_enable_casting() 
#define ale_pos_casting_status() 1

#elif ALE_COORDINATES == FIXED32

typedef ale_fixed<ale_fixed_32,16> ale_pos;

#define ALE_POS_PRECISION_STRING "FIXED32"

#define ale_pos_disable_casting() ale_pos::disable_casting()
#define ale_pos_enable_casting() ale_pos::enable_casting()
#define ale_pos_casting_status() ale_pos::casting_status()

#elif ALE_COORDINATES == FIXED16

typedef ale_fixed<ale_fixed_16_calc,10> ale_pos;

#define ALE_POS_PRECISION_STRING "FIXED16"

#define ale_pos_disable_casting() ale_pos::disable_casting()
#define ale_pos_enable_casting() ale_pos::enable_casting()
#define ale_pos_casting_status() ale_pos::casting_status()

#else

#warning Unknown positional precision in ale_pos.h: Choosing PPRECISION=SINGLE.

typedef float ale_pos;

#define ALE_POS_PRECISION_STRING "SINGLE"

#define ale_pos_disable_casting()
#define ale_pos_enable_casting() 
#define ale_pos_casting_status() 1

#endif

const ale_pos ale_pos_0 = (ale_pos) 0;

#if ALE_COLORS == FIXED && ALE_COORDINATES == FIXED

#define ale_pos_to_real(x) (x)
#define ale_real_to_pos(x) (x)

#else

#define ale_pos_to_real(x) (x)
#define ale_real_to_pos(x) (x)

#endif

#undef SINGLE
#undef DOUBLE
#undef HALF
#undef FIXED16
#undef FIXED32

#endif
