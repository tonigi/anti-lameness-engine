// Copyright 2002 David Hilvert <dhilvert@auricle.dyndns.org>, 
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

#include "align.h"

/*
 * See align.h for details on these variables.
 */

int align::_exp_register = 1;

ale_pos align::scale_factor;

int align::_keep = 0;
transformation *align::kept_t = NULL;
int *align::kept_ok = NULL;

tload_t *align::tload = NULL;
tsave_t *align::tsave = NULL;
render *align::reference = NULL;
const image *align::reference_image = NULL;
const image *align::reference_defined = NULL;
transformation align::latest_t;
int align::latest_ok;
int align::latest = 0;
int align::extend = 0;

int align::alignment_class = 1;
int align::default_initial_alignment_type = 0;
transformation align::default_initial_alignment;
int align::is_default = 1;
int align::old_is_default;
transformation align::old_initial_alignment;
transformation align::old_final_alignment;
int align::channel_alignment_type = 2;
float align::metric_exponent = 2;
float align::match_threshold = 0;
ale_pos align::perturb_lower = 0.125;
ale_pos align::perturb_upper = 32;

/*
 * Empirically, it's okay to use a level-of-detail equal to twice the
 * resolution of the perturbation, so we set the default lod_max to 1, as
 * 2^1==2.  lod_max of zero seems okay also, but lower values seem to cause
 * problems.
 */

int align::lod_max = 1;

ale_pos align::rot_max = 32.0;
ale_accum align::match_sum = 0;
int align::match_count = 0;

ale_pos align::_mc = 0;
