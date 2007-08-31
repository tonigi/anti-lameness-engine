// Copyright 2002, 2003, 2004, 2005, 2006 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                                      <dhilvert@ugcs.caltech.edu>

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

#ifndef __ale_math_h__
#define __ale_math_h__

/*
 * Certain versions of Mac OSX may disable math.h definitions of is*() macros
 * when iostream is included, so we include the latter here.
 */

#include <math.h>
#include <iostream>

/*
 * isnan() and isinf() code logic are based on those noted in the GNU Autoconf
 * manual, by the Free Software Foundation:
 *
 * http://www.gnu.org/software/autoconf/manual/html_node/Function-Portability.html
 * 
 * As C++ is available here, we use C++ overloading instead of sizeof()
 * switches to handle different types.
 */

#ifndef isnan
# define isnan(x) ale_isnan(x)
static inline int ale_isnan(float       x) { return x != x; }
static inline int ale_isnan(double      x) { return x != x; }
static inline int ale_isnan(long double x) { return x != x; }
#endif

#ifndef isinf
# define isinf(x) ale_isinf(x)
static inline int ale_isinf(float       x) { return isnan (x - x); }
static inline int ale_isinf(double      x) { return isnan (x - x); }
static inline int ale_isinf(long double x) { return isnan (x - x); }
#endif

#endif
