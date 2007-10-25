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

#ifndef __ale_accum_h__
#define __ale_accum_h__

#include "ale_fixed.h"

#define FIXED16 4
#define FIXED32 5

/*
 * Real-valued type used when accumulating over the domain of an image.
 */

#if ALE_COLORS == FIXED16

typedef ale_fixed<ale_fixed_16_accum,12> ale_accum;

#define ale_accum_enable_casting() ale_accum::enable_casting()
#define ale_accum_disable_casting() ale_accum::disable_casting()

#elif ALE_COLORS == FIXED32

typedef ale_fixed<ale_fixed_32_accum,15> ale_accum;

#define ale_accum_enable_casting() ale_accum::enable_casting()
#define ale_accum_disable_casting() ale_accum::disable_casting()

#else

typedef double ale_accum;

#define ale_accum_disable_casting()
#define ale_accum_enable_casting()

#endif

#undef FIXED32
#undef FIXED16

#endif
