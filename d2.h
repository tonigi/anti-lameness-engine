// Copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>, 
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
 * Top-level header file for classes treating scenes as two-dimensional data.
 * autoconf 'config.h' should be included after this file, as we undefine
 * various autoconf defines herein.
 */

#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include <vector>
#include <stack>
#include <algorithm>

#include "ale_real.h"
#include "ale_accum.h"
#include "ale_pos.h"
#include "ui/ui.h"


#include <iostream>
#include <ostream>
#include <typeinfo>

#ifdef USE_MAGICK
#include <magick/api.h>

/*
 * ImageMagick defines these, for reasons unclear.
 * Since they clash with autotools names, undefine
 * them here.
 */

#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef PACKAGE_STRING
#undef PACKAGE_NAME
#undef PACKAGE_BUGREPORT
#endif

#ifdef USE_FFTW
#include <fftw3.h>
#endif

#ifdef USE_UNIX
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef USE_PTHREAD
#include <pthread.h>
#include "thread.h"
#endif

#include "rand.h"



/*
 * All header files in the d2 namespace.
 */

namespace d2 {
#include "d2/exclusion.h"
#include "d2/pixel.h"
#include "d2/spixel.h"
#include "d2/pixel_accum.h"
#include "d2/exposure/exposure.h"
#include "d2/exposure/exposure_default.h"
#include "d2/exposure/exposure_linear.h"
#include "d2/exposure/exposure_boolean.h"
#include "d2/align.h"
#include "d2/transformation.h"
#include "d2/image.h"
}

/*
 * XXX: The placement of this file is somewhat of a hack.  What should
 * be done about this?
 */
#include "optimizations.h"

namespace d2 {
#include "d2/image_ale_real.h"
#include "d2/image_weighted_avg.h"
#include "d2/image_weighted_simple.h"
#include "d2/image_weighted_median.h"
#include "d2/image_zero.h"
#include "d2/image_rw.h"
#include "d2/point.h"
#include "d2/ppm.h"
#include "d2/render.h"
#include "d2/render_parse.h"
#include "d2/tfile.h"
#include "d2/filter.h"
#include "d2/render/combine.h"
// #include "d2/render/drizzle.h"
// #include "d2/render/usm.h"
#include "d2/render/ipc.h"
// #include "d2/render/merge.h"
#include "d2/render/psf/psf.h"
#include "d2/render/psf/psf_template.h"
#include "d2/render/psf/box.h"	
#include "d2/render/psf/circle.h"	
#include "d2/render/psf/sum.h"
#include "d2/render/psf/stdin.h"	
#include "d2/render/psf/stdin_vg.h"	
#include "d2/render/psf/convolution.h"
#include "d2/render/psf/scalar_mult.h"
#include "d2/render/psf/psf_parse.h"
#include "d2/render/psf/psf_calibrate.h"
#include "d2/vise_core.h"

}
