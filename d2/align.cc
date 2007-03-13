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

transformation align::orig_t;
int align::_keep = 0;
transformation *align::kept_t = NULL;
int *align::kept_ok = NULL;

tload_t *align::tload = NULL;
tsave_t *align::tsave = NULL;
render *align::reference = NULL;
filter::scaled_filter *align::interpolant = NULL;
const image *align::reference_image = NULL;
const image *align::reference_defined = NULL;
const image *align::weight_map = NULL;
image *align::alignment_weights = NULL;
const image *align::alignment_weights_const = NULL;
const char *align::wmx_file = NULL;
const char *align::wmx_exec = NULL;
const char *align::wmx_defs = NULL;
const char *align::fw_output = NULL;
double align::horiz_freq_cut = 0;
double align::vert_freq_cut = 0;
double align::avg_freq_cut = 0;
transformation align::latest_t;
int align::latest_ok;
int align::latest = -1;

int align::alignment_class = 1;
int align::default_initial_alignment_type = 1;
int align::perturb_type = 0;
transformation align::default_initial_alignment;
int align::is_fail_default = 0;
int align::channel_alignment_type = 2;
float align::metric_exponent = 2;
float align::match_threshold = 0;

/*
 * Upper/lower bounds
 */

ale_pos align::perturb_lower = 0.125;
int align::perturb_lower_percent = 0;
ale_pos align::perturb_upper = 14;
int align::perturb_upper_percent = 1;

/*
 * XXX: These empirical observations may have been based on buggy code.
 * 
 * Empirically, it's okay to use a level-of-detail equal to twice the
 * resolution of the perturbation, so we set the default lod_max to 1, as
 * 2^1==2.  lod_max of zero seems okay also, but lower values seem to cause
 * problems.
 */

int align::lod_max = 1;

ale_pos align::rot_max = 32.0;
ale_pos align::bda_mult = 0.0001;
ale_pos align::bda_rate = 0.0004;
ale_accum align::match_sum = 0;
int align::match_count = 0;

ale_pos align::_mc = 0;
int align::_mcd_limit = 100;
int align::certainty_weights = 0;
int align::_gs = 0;
int align::_ma_card = 0;
double align::_ma_cont = 100;
unsigned int align::_gs_mo = 16;
exclusion *align::ax_parameters = NULL;
int align::ax_count = 0;

const point **align::cp_array = NULL;
unsigned int align::cp_count = 0;
