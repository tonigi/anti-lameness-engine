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
ale_pos align::bda_mult = 2;
ale_pos align::bda_rate = 8;
ale_accum align::match_sum = 0;
int align::match_count = 0;

ale_pos align::_mc = 0;
int align::_mcd_limit = 1;
int align::certainty_weights = 0;
int align::_gs = 6;
unsigned int align::_ma_card = 1;
double align::_ma_cont = 100;

ale_pos align::_gs_mo = 67;
int align::gs_mo_percent = 1;

exclusion *align::ax_parameters = NULL;
int align::ax_count = 0;

const point **align::cp_array = NULL;
unsigned int align::cp_count = 0;

std::vector<align::diff_stat_t::memo_entry> align::diff_stat_t::memos;

void *d2::align::diff_stat_t::diff_subdomain(void *args) {

	subdomain_args *sargs = (subdomain_args *) args;

	struct scale_cluster c = sargs->c;
	std::vector<run> runs = sargs->runs;
	ale_pos _mc_arg = sargs->_mc_arg;
	int ax_count = sargs->ax_count;
	int f = sargs->f;
	int i_min = sargs->i_min;
	int i_max = sargs->i_max;
	int j_min = sargs->j_min;
	int j_max = sargs->j_max;
	int subdomain = sargs->subdomain;

	assert (reference_image);

	point offset = c.accum->offset();

	int i, j;

	/*
	 * We always use the same code for exhaustive and Monte Carlo
	 * pixel sampling, setting _mc_arg = 1 when all pixels are to
	 * be sampled.
	 */

	if (_mc_arg <= 0 || _mc_arg >= 1)
		_mc_arg = 1;

	int index;

	int index_max = (i_max - i_min) * (j_max - j_min);
	
	/*
	 * We use a random process for which the expected
	 * number of sampled pixels is +/- .000003 from mc_arg
	 * in the range [.005,.995] for an image with 100,000
	 * pixels.  (The actual number may still deviate from
	 * the expected number by more than this amount,
	 * however.)  The method is as follows:
	 *
	 * We have mc_arg == USE/ALL, or (expected # pixels to
	 * use)/(# total pixels).  We derive from this
	 * SKIP/USE.
	 *
	 * SKIP/USE == (SKIP/ALL)/(USE/ALL) == (1 - (USE/ALL))/(USE/ALL)
	 *
	 * Once we have SKIP/USE, we know the expected number
	 * of pixels to skip in each iteration.  We use a random
	 * selection process that provides SKIP/USE close to
	 * this calculated value.
	 *
	 * If we can draw uniformly to select the number of
	 * pixels to skip, we do.  In this case, the maximum
	 * number of pixels to skip is twice the expected
	 * number.
	 *
	 * If we cannot draw uniformly, we still assign equal
	 * probability to each of the integer values in the
	 * interval [0, 2 * (SKIP/USE)], but assign an unequal
	 * amount to the integer value ceil(2 * SKIP/USE) + 1.
	 */

	ale_pos su = (1 - _mc_arg) / _mc_arg;

	ale_pos mc_max = (floor(2*su) * (1 + floor(2*su)) + 2*su)
		      / (2 + 2 * floor(2*su) - 2*su);

	/*
	 * Reseed the random number generator;  we want the
	 * same set of pixels to be used when comparing two
	 * alignment options.  If we wanted to avoid bias from
	 * repeatedly utilizing the same seed, we could seed
	 * with the number of the frame most recently aligned:
	 *
	 * 	srand(latest);
	 *
	 * However, in cursory tests, it seems okay to just use
	 * the default seed of 1, and so we do this, since it
	 * is simpler; both of these approaches to reseeding
	 * achieve better results than not reseeding.  (1 is
	 * the default seed according to the GNU Manual Page
	 * for rand(3).)
	 * 
	 * For subdomain calculations, we vary the seed by subdomain.
	 */

	rng_t rng;

	rng.seed(1 + subdomain);

	for(index = -1 + (int) ceil((mc_max+1) 
				  * ( (1 + ((ale_pos) (rng.get())) ) 
				    / (1 + ((ale_pos) RAND_MAX)) ));
	    index < index_max;
	    index += (int) ceil((mc_max+1) 
			      * ( (1 + ((ale_pos) (rng.get())) ) 
				/ (1 + ((ale_pos) RAND_MAX)) ))){

		i = index / (j_max - j_min) + i_min;
		j = index % (j_max - j_min) + j_min;

		/*
		 * Check for exclusion in render coordinates.
		 */

		if (ref_excluded(i, j, offset, c.ax_parameters, ax_count))
			continue;

		/*
		 * Check transformation support.
		 */

		if (!runs.back().offset.supported((int) (i + offset[0]), (int) (j + offset[1])))
			continue;

		/*
		 * Transform
		 */

		struct point q, r = point::undefined();

		q = (c.input_scale < 1.0 && interpolant == NULL)
		  ? runs.back().offset.scaled_inverse_transform(
			point(i + offset[0], j + offset[1]))
		  : runs.back().offset.unscaled_inverse_transform(
			point(i + offset[0], j + offset[1]));

		if (runs.size() == 2) {
			r = (c.input_scale < 1.0)
			  ? runs.front().offset.scaled_inverse_transform(
				point(i + offset[0], j + offset[1]))
			  : runs.front().offset.unscaled_inverse_transform(
				point(i + offset[0], j + offset[1]));
		}

		ale_pos ti = q[0];
		ale_pos tj = q[1];

		/*
		 * Check that the transformed coordinates are within
		 * the boundaries of array c.input and that they
		 * are not subject to exclusion.
		 *
		 * Also, check that the weight value in the accumulated array
		 * is nonzero, unless we know it is nonzero by virtue of the
		 * fact that it falls within the region of the original frame
		 * (e.g. when we're not increasing image extents).
		 */

		if (input_excluded(ti, tj, c.ax_parameters, ax_count))
			continue;

		if (input_excluded(r[0], r[1], c.ax_parameters, ax_count))
			r = point::undefined();

		/*
		 * Check the boundaries of the input frame
		 */

		if (!(ti >= 0
		   && ti <= c.input->height() - 1
		   && tj >= 0
		   && tj <= c.input->width() - 1
		   && c.defined->get_pixel(i, j)[0] != 0))
			continue;

		if (!(r[0] >= 0
		   && r[0] <= c.input->height() - 1
		   && r[1] >= 0
		   && r[1] <= c.input->width() - 1))
			r = point::undefined();

		sargs->runs.back().sample(f, c, i, j, q, r);
	}

	return NULL;
}

void d2::align::diff_stat_t::try_decrease() {

	assert (_mc <= 0);
	
	diff_stat_t half(*this);

	half.mcd_decrease();

	if (!half.reliable())
		return;

	*this = half;

	ui::get()->alignment_monte_carlo_parameter(get_mc());
}

