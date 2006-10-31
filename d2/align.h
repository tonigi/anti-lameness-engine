// Copyright 2002, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

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

/*
 * align.h: Handle alignment of frames.
 */

#ifndef __d2align_h__
#define __d2align_h__

#include "filter.h"
#include "transformation.h"
#include "image.h"
#include "point.h"
#include "render.h"
#include "tfile.h"
#include "image_rw.h"

class align {
private:

	/*
	 * Private data members
	 */

	static ale_pos scale_factor;

	/*
	 * Original frame transformation
	 */
	static transformation orig_t;

	/*
	 * Keep data older than latest
	 */
	static int _keep;
	static transformation *kept_t;
	static int *kept_ok;

	/*
	 * Transformation file handlers
	 */

	static tload_t *tload;
	static tsave_t *tsave;

	/*
	 * Control point variables
	 */

	static const point **cp_array;
	static unsigned int cp_count;

	/*
	 * Reference rendering to align against
	 */

	static render *reference;
	static filter::scaled_filter *interpolant;
	static const image *reference_image;
	static const image *reference_defined;

 	/*
	 * Per-pixel alignment weight map
	 */

	static const image *weight_map;

	/*
	 * Frequency-dependent alignment weights
	 */

	static double horiz_freq_cut;
	static double vert_freq_cut;
	static double avg_freq_cut;
	static const char *fw_output;

	/*
	 * Algorithmic alignment weighting
	 */

	static const char *wmx_exec;
	static const char *wmx_file;
	static const char *wmx_defs;

	/*
	 * Consolidated alignment weights
	 *
	 * The final value of alignment_weights_const is the canonical weight
	 * set.  If alignment_weights_const is set but alignment_weights is
	 * not, then the memory is not ours, and the object should not be
	 * modified or deleted.
	 */

	static image *alignment_weights;
	static const image *alignment_weights_const;

	/*
	 * Latest transformation.
	 */

	static transformation latest_t;

	/*
	 * Flag indicating whether the latest transformation 
	 * resulted in a match.
	 */

	static int latest_ok;

	/*
	 * Frame number most recently aligned.
	 */

	static int latest;

	/*
	 * Exposure registration
	 *
	 * 0. Preserve the original exposure of images.
	 *
	 * 1. Match exposure between images.
	 *
	 * 2. Use only image metadata for registering exposure.
	 */

	static int _exp_register;

	/*
	 * Alignment class.
	 *
	 * 0. Translation only.  Only adjust the x and y position of images.
	 * Do not rotate input images or perform projective transformations.
	 * 
	 * 1. Euclidean transformations only.  Adjust the x and y position
	 * of images and the orientation of the image about the image center.
	 *
	 * 2. Perform general projective transformations.  See the file gpt.h
	 * for more information about general projective transformations.
	 */

	static int alignment_class;

	/*
	 * Default initial alignment.
	 *
	 * 0. Identity transformation.
	 *
	 * 1. Most recently accepted frame's final transformation.
	 */

	static int default_initial_alignment_type;
	static transformation default_initial_alignment;

	/*
	 * Projective group behavior
	 *
	 * 0. Perturb in output coordinates.
	 *
	 * 1. Perturb in source coordinates
	 */

	static int perturb_type;

	/*
	 * Helper variables for new --follow semantics (as of 0.5.0).
	 *
	 * These variables implement delta following.  The change between the
	 * non-default old initial alignment and old final alignment is used to
	 * adjust the non-default current initial alignment.  If either the old
	 * or new initial alignment is a default alignment, the old --follow
	 * semantics are preserved.
	 */

	static int is_default, old_is_default;
	static int old_lod;
	static transformation old_initial_alignment;
	static transformation old_final_alignment;

	/*
	 * Alignment for failed frames -- default or optimal?
	 *
	 * A frame that does not meet the match threshold can be assigned the
	 * best alignment found, or can be assigned its alignment default.
	 */

	static int is_fail_default;

	/*
	 * Alignment code.
	 * 
	 * 0. Align images with an error contribution from each color channel.
	 * 
	 * 1. Align images with an error contribution only from the green channel.
	 * Other color channels do not affect alignment.
	 *
	 * 2. Align images using a summation of channels.  May be useful when dealing
	 * with images that have high frequency color ripples due to color aliasing.
	 */

	static int channel_alignment_type;

	/*
	 * Error metric exponent
	 */

	static float metric_exponent;

	/*
	 * Match threshold
	 */

	static float match_threshold;

	/*
	 * Perturbation lower and upper bounds.
	 */

	static ale_pos perturb_lower;
	static int perturb_lower_percent;
	static ale_pos perturb_upper;
	static int perturb_upper_percent;

	/*
	 * Maximum level-of-detail scale factor is 2^lod_max/perturb.
	 */

	static int lod_max;

	/*
	 * Maximum rotational perturbation
	 */

	static ale_pos rot_max;

	/*
	 * Barrel distortion alignment multiplier
	 */

	static ale_pos bda_mult;

	/*
	 * Barrel distortion maximum adjustment rate
	 */

	static ale_pos bda_rate;

	/*
	 * Alignment match sum
	 */

	static ale_accum match_sum;

	/*
	 * Alignment match count.
	 */

	static int match_count;

	/*
	 * Monte Carlo parameter
	 *
	 * 0.  Don't use monte carlo alignment sampling.
	 *
	 * (0,1].  Select, on average, a number of pixels which is the
	 * specified fraction of the number of pixels in the accumulated image.
	 */

	static ale_pos _mc;
	static int _mcd_limit;

	/*
	 * Certainty weight flag
	 *
	 * 0. Don't use certainty weights for alignment.
	 *
	 * 1. Use certainty weights for alignment.
	 */

	static int certainty_weights;

	/*
	 * Global search parameter
	 *
	 * 0.  Local:   Local search only.
	 * 1.  Inner:   Alignment reference image inner region
	 * 2.  Outer:   Alignment reference image outer region
	 * 3.  All:     Alignment reference image inner and outer regions.
	 * 4.  Central: Inner if possible; else, best of inner and outer.
	 * 5.  Points:  Align by control points.
	 */

	static int _gs;

	/*
	 * Minimum overlap for global searches
	 */

	static unsigned int _gs_mo;

	/*
	 * Exclusion regions
	 */

	static exclusion *ax_parameters;
	static int ax_count;

	/*
	 * Type for scale cluster.
	 */

	struct scale_cluster {
		const image *accum;
		const image *defined;
		const image *aweight;
		exclusion *ax_parameters;

		ale_pos input_scale;
		const image *input;
	};

	/*
	 * Check for exclusion region coverage in the reference 
	 * array.
	 */
	static int ref_excluded(int i, int j, point offset, exclusion *params, int param_count) {
		for (int idx = 0; idx < param_count; idx++)
			if (params[idx].type == exclusion::RENDER
			 && i + offset[0] >= params[idx].x[0]
			 && i + offset[0] <= params[idx].x[1]
			 && j + offset[1] >= params[idx].x[2]
			 && j + offset[1] <= params[idx].x[3])
				return 1;

		return 0;
	}

	/*
	 * Check for exclusion region coverage in the input 
	 * array.
	 */
	static int input_excluded(ale_pos ti, ale_pos tj, exclusion *params, int param_count) {
		for (int idx = 0; idx < param_count; idx++)
			if (params[idx].type == exclusion::FRAME
			 && ti >= params[idx].x[0]
			 && ti <= params[idx].x[1]
			 && tj >= params[idx].x[2]
			 && tj <= params[idx].x[3])
				return 1;

		return 0;
	}

	/*
	 * Overlap function.  Determines the number of pixels in areas where
	 * the arrays overlap.  Uses the reference array's notion of pixel
	 * positions.
	 */
	static unsigned int overlap(struct scale_cluster c, transformation t, int ax_count) {
		assert (reference_image);

		unsigned int result = 0;

		point offset = c.accum->offset();

		for (unsigned int i = 0; i < c.accum->height(); i++)
		for (unsigned int j = 0; j < c.accum->width();  j++) {

			if (ref_excluded(i, j, offset, c.ax_parameters, ax_count))
				continue;

			/*
			 * Transform
			 */

			struct point q;

			q = (c.input_scale < 1.0 && interpolant == NULL)
			  ? t.scaled_inverse_transform(
				point(i + offset[0], j + offset[1]))
			  : t.unscaled_inverse_transform(
				point(i + offset[0], j + offset[1]));

			ale_pos ti = q[0];
			ale_pos tj = q[1];

			/*
			 * Check that the transformed coordinates are within
			 * the boundaries of array c.input, and check that the
			 * weight value in the accumulated array is nonzero,
			 * unless we know it is nonzero by virtue of the fact
			 * that it falls within the region of the original
			 * frame (e.g. when we're not increasing image
			 * extents).  Also check for frame exclusion.
			 */

			if (input_excluded(ti, tj, c.ax_parameters, ax_count))
				continue;

			if (ti >= 0
			 && ti <= c.input->height() - 1
			 && tj >= 0
			 && tj <= c.input->width() - 1
			 && c.defined->get_pixel(i, j)[0] != 0)
				result++;
		}

		return result;
	}

	/*
	 * Not-quite-symmetric difference function.  Determines the difference in areas
	 * where the arrays overlap.  Uses the reference array's notion of pixel positions.
	 *
	 * For the purposes of determining the difference, this function divides each
	 * pixel value by the corresponding image's average pixel magnitude, unless we
	 * are:
	 *
	 * a) Extending the boundaries of the image, or
	 * 
	 * b) following the previous frame's transform
	 *
	 * If we are doing monte-carlo pixel sampling for alignment, we
	 * typically sample a subset of available pixels; otherwise, we sample
	 * all pixels.
	 *
	 */

	class diff_stat_t {
		ale_accum result;
		ale_accum divisor;

		typedef unsigned int hist_bin;

		int hist_min_r;
		int hist_min_d;

		hist_bin *histogram;
		hist_bin *histogram_integral;
		hist_bin hist_total;
		int hist_size;

		void init_histogram_integral() {
			histogram_integral = (hist_bin *)malloc(sizeof(hist_bin) * hist_size);
			hist_bin total = 0;
			for (int i = 0; i < hist_size; i++) {
				total += histogram[i];
				histogram_integral[i] = total;
			}
		}

		int histogram_integral_inverse(unsigned int i) {
			assert (hist_size > 0);

			int min = -1;
			int max = hist_size;

			while (min < max - 1) {
				int mid = (min + max) / 2;
				if (histogram_integral[mid] <= i) {
					min = mid;
				} else {
					max = mid;
				}
			}

			return min + 1;
		}

		ale_accum simulated_error(rng_t rng, ale_pos mc) {
			ale_accum result = 0;
			ale_accum divisor = 0;

			int samples = (int) floor(mc * hist_total);

			ui::get()->d2_align_sim_start();

			for (int s = 0; s < samples; s++) {
				int index = rng.get() % hist_total;
				int histogram_index = histogram_integral_inverse((unsigned int) index);

				result += pow(2, hist_min_r + histogram_index / hist_size);
				divisor += pow(2, hist_min_d + histogram_index % hist_size);
			}

			ui::get()->d2_align_sim_stop();

			return pow(result / divisor, 1 / metric_exponent);
		}

		struct consensus_subdomain {
			diff_stat_t *better;
			diff_stat_t *worse;
			ale_pos mc;
			int run_count;
			int success;
		};

		static void *consensus_thread (void *args) {
			consensus_subdomain *cs = (consensus_subdomain *) args;

			cs->better->init_histogram_integral();
			cs->worse->init_histogram_integral();

			for (int run = 0; run < cs->run_count; run++) {
				rng_t rng;

				rng.seed(run);

				ale_accum b = cs->better->simulated_error(rng, cs->mc);
				ale_accum w = cs->worse->simulated_error(rng, cs->mc);

				if (b < w)
					cs->success++;
			}

			free(cs->better->histogram_integral);
			free(cs->worse->histogram_integral);

			return NULL;
		}

		void add_hist(int r, int d, int count) {

			hist_total += count;

			int r_shift = 0, d_shift = 0;

			if (r - hist_min_r >= hist_size) {
				r_shift = (r - hist_min_r) - hist_size + 1;
				hist_min_r += r_shift;
			}

			if (d - hist_min_d >= hist_size) {
				d_shift = (d - hist_min_d) - hist_size + 1;
				hist_min_d += d_shift;
			}

			assert (r_shift >= 0);
			assert (d_shift >= 0);

			if (r_shift || d_shift) {
				for (int rr = 0; rr < hist_size; rr++)
				for (int dd = 0; dd < hist_size; dd++) {

					hist_bin value = histogram[rr * hist_size + dd];

					histogram[rr * hist_size + dd] = 0;

					int rrr = rr - r_shift;
					int ddd = dd - d_shift;

					if (rrr < 0)
						rrr = 0;
					if (ddd < 0)
						ddd = 0;

					histogram[rrr * hist_size + ddd] += value;
				}
			}

			r -= hist_min_r;
			d -= hist_min_d;

			if (r < 0)
				r = 0;
			if (d < 0)
				d = 0;

			histogram[r * hist_size + d] += count;
		}

		void add_hist(ale_accum result, ale_accum divisor) {
			ale_accum rbin = log(result) / log(2);
			ale_accum dbin = log(divisor) / log(2);

			if (!(rbin > INT_MIN))
				rbin = INT_MIN;
			if (!(dbin > INT_MIN))
				dbin = INT_MIN;

			add_hist((int) floor(rbin), (int) floor(dbin), 1);
		}

	public:
		diff_stat_t() {
			result = 0;
			divisor = 0;
			hist_min_r = INT_MIN;
			hist_min_d = INT_MIN;
			hist_size = 20;
			hist_total = 0;
			histogram = (hist_bin *) calloc(hist_size * hist_size, sizeof(hist_bin));
		}

		void clear() {
			free(histogram);

			result = 0;
			divisor = 0;
			hist_min_r = INT_MIN;
			hist_min_d = INT_MIN;
			hist_total = 0;
			histogram = (hist_bin *) calloc(hist_size * hist_size, sizeof(hist_bin));
		}

		~diff_stat_t() {
			free(histogram);
		}

		int check_removal(diff_stat_t *with) {
			ale_accum bresult, bdivisor, wresult, wdivisor;
			hist_bin *bhist, *whist;

			bhist = (hist_bin *)malloc(sizeof(hist_bin) * hist_size);
			whist = (hist_bin *)malloc(sizeof(hist_bin) * with->hist_size);

			bresult = result;
			bdivisor = divisor;
			wresult = with->result;
			wdivisor = with->divisor;

			for (int i = 0; i < hist_size; i++)
				bhist[i] = histogram[i];

			for (int i = 0; i < with->hist_size; i++)
				whist[i] = with->histogram[i];

			for (int r = 0; r < _mcd_limit; r++) {
				int max_gradient_bin = -1;
				int max_gradient_hist = -1;
				ale_accum max_gradient = 0;

				for (int i = 0; i < hist_size; i++) {
					if (bhist[i] <= 0)
						continue;

					ale_accum br = pow(2, hist_min_r + i / hist_size);
					ale_accum bd = pow(2, hist_min_d + i % hist_size);
					ale_accum b_test_gradient = 
						bresult / bdivisor - (bresult - br) / (bdivisor - bd);
					
					if (b_test_gradient > max_gradient) {
						max_gradient_bin = i;
						max_gradient_hist = 0;
						max_gradient = b_test_gradient;
					}
				}

				for (int i = 0; i < with->hist_size; i++) {
					if (whist[i] <= 0)
						continue;

					ale_accum wr = pow(2, with->hist_min_r + i / with->hist_size);
					ale_accum wd = pow(2, with->hist_min_d + i % with->hist_size);
					ale_accum w_test_gradient = 
						(wresult - wr) / (wdivisor - wd) - wresult / wdivisor;
					
					if (w_test_gradient > max_gradient) {
						max_gradient_bin = i;
						max_gradient_hist = 1;
						max_gradient = w_test_gradient;
					}
				}

				if (max_gradient_hist == 0) {
					bhist[max_gradient_bin]--;
					bresult -= pow(2, hist_min_r + max_gradient_bin / hist_size);
					bdivisor -= pow(2, hist_min_d + max_gradient_bin / hist_size);
				} else if (max_gradient_hist == 1) {
					whist[max_gradient_bin]--;
					wresult -= pow(2, with->hist_min_r + max_gradient_bin / with->hist_size);
					wdivisor -= pow(2, with->hist_min_d + max_gradient_bin / with->hist_size);
				}
			}

			free(bhist);
			free(whist);

			if (bresult / bdivisor < wresult / wdivisor)
				return 1;
			
			return 0;
		}

		ale_pos consensus(diff_stat_t *with, ale_pos mc) {

			int _mcd_runs = 10;

			consensus_subdomain cs = {
				this, with, mc, _mcd_runs, 0
			};

			consensus_thread(&cs);

			return ((double) cs.success / (double) _mcd_runs);
		}

		int reliable(diff_stat_t *with, ale_pos mc) {
#if 0
			return consensus(with, mc) >= _mcd_lower;
#else
			return check_removal(with);
#endif
		}

		void add(const diff_stat_t *ds) {
			result += ds->result;
			divisor += ds->divisor;

			for (int r = 0; r < ds->hist_size; r++)
			for (int d = 0; d < ds->hist_size; d++)
				add_hist(r + ds->hist_min_r, d + ds->hist_min_d, 
						ds->histogram[r * hist_size + d]);
		}

		ale_accum get_result() {
			return result;
		}

		ale_accum get_divisor() {
			return divisor;
		}

		ale_accum get_error() {
			return pow(result / divisor, 1/metric_exponent);
		}

		void sample(int f, scale_cluster c, int i, int j, ale_pos ti, ale_pos tj) {
			pixel pa = c.accum->get_pixel(i, j);
			pixel pb;
			pixel weight;
			ale_accum this_result = 0;
			ale_accum this_divisor = 0;

			if (interpolant != NULL)
				interpolant->filtered(i, j, &pb, &weight, 1, f);
			else {
				pixel result[2];
				c.input->get_bl(point(ti, tj), result);
				pb = result[0];
				weight = result[1];
			}

			/*
			 * Handle certainty.
			 */

			if (certainty_weights == 0)
				weight = pixel(1, 1, 1);

			if (c.aweight != NULL)
				weight *= c.aweight->get_pixel(i, j);

			/*
			 * Determine alignment type.
			 */

			if (channel_alignment_type == 0) {
				/*
				 * Align based on all channels.
				 */


				for (int k = 0; k < 3; k++) {
					ale_real achan = pa[k];
					ale_real bchan = pb[k];

					this_result += weight[k] * pow(fabs(achan - bchan), metric_exponent);
					this_divisor += weight[k] * pow(achan > bchan ? achan : bchan, metric_exponent);
				}
			} else if (channel_alignment_type == 1) {
				/*
				 * Align based on the green channel.
				 */

				ale_real achan = pa[1];
				ale_real bchan = pb[1];

				this_result = weight[1] * pow(fabs(achan - bchan), metric_exponent);
				this_divisor = weight[1] * pow(achan > bchan ? achan : bchan, metric_exponent);
			} else if (channel_alignment_type == 2) {
				/*
				 * Align based on the sum of all channels.
				 */

				ale_real asum = 0;
				ale_real bsum = 0;
				ale_real wsum = 0;

				for (int k = 0; k < 3; k++) {
					asum += pa[k];
					bsum += pb[k];
					wsum += weight[k] / 3;
				}

				this_result = wsum * pow(fabs(asum - bsum), metric_exponent);
				this_divisor = wsum * pow(asum > bsum ? asum : bsum, metric_exponent);
			}

			result += this_result;
			divisor += this_divisor;

			add_hist(this_result, this_divisor);
		}

		void print_hist() {
			fprintf(stderr, "\n");
			fprintf(stderr, "hist_min_r = %d\n", hist_min_r);
			fprintf(stderr, "hist_min_d = %d\n", hist_min_d);
			fprintf(stderr, "hist_size = %d\n", hist_size);
			fprintf(stderr, "hist_total = %d\n", hist_total);
			fprintf(stderr, "\n");

			hist_bin recalc_total = 0;

			for (int r = 0; r < hist_size; r++) {
				for (int d = 0; d < hist_size; d++) {
					recalc_total += histogram[r * hist_size + d];
					fprintf(stderr, "\t%d", histogram[r * hist_size + d]);
				}
				fprintf(stderr, "\n");
			}

			fprintf(stderr, "\n");
			fprintf(stderr, "recalc_total = %d\n", recalc_total);
		}
	};

	struct subdomain_args {
		struct scale_cluster c;
		transformation t;
		ale_pos _mc_arg;
		int ax_count;
		int f;
		diff_stat_t* diff_stat;
		int i_min, i_max, j_min, j_max;
		int subdomain;
	};

	static void *diff_subdomain(void *args) {

		subdomain_args *sargs = (subdomain_args *) args;

		struct scale_cluster c = sargs->c;
		transformation t = sargs->t;
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

		ale_pos u = (1 - _mc_arg) / _mc_arg;

		ale_pos mc_max = (floor(2*u) * (1 + floor(2*u)) + 2*u)
			      / (2 + 2 * floor(2*u) - 2*u);

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
			 * Transform
			 */

			struct point q;

			q = (c.input_scale < 1.0 && interpolant == NULL)
			  ? t.scaled_inverse_transform(
				point(i + offset[0], j + offset[1]))
			  : t.unscaled_inverse_transform(
				point(i + offset[0], j + offset[1]));

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

			if (ti >= 0
			 && ti <= c.input->height() - 1
			 && tj >= 0
			 && tj <= c.input->width() - 1
			 && c.defined->get_pixel(i, j)[0] != 0)

				sargs->diff_stat->sample(f, c, i, j, ti, tj);

		}

		return NULL;
	}

	static ale_accum diff(struct scale_cluster c, transformation t,
			ale_pos _mc_arg, int ax_count, int f, diff_stat_t *diff_stat_in = NULL) {

		diff_stat_t *diff_stat = diff_stat_in;

		if (diff_stat == NULL)
			diff_stat = new diff_stat_t();

		diff_stat->clear();

		ui::get()->d2_align_sample_start();

		if (interpolant != NULL) 
			interpolant->set_parameters(t, c.input, c.accum->offset());

		int N;
#ifdef USE_PTHREAD
		N = thread::count();

		pthread_t *threads = (pthread_t *) malloc(sizeof(pthread_t) * N);
		pthread_attr_t *thread_attr = (pthread_attr_t *) malloc(sizeof(pthread_attr_t) * N);

#else
		N = 1;
#endif

		subdomain_args *args = new subdomain_args[N];

		for (int ti = 0; ti < N; ti++) {
			args[ti].c = c;
			args[ti].t = t;
			args[ti]._mc_arg = _mc_arg;
			args[ti].ax_count = ax_count;
			args[ti].f = f;
			args[ti].diff_stat = new diff_stat_t();
			args[ti].i_min = (c.accum->height() * ti) / N;
			args[ti].i_max = (c.accum->height() * (ti + 1)) / N;
			args[ti].j_min = 0;
			args[ti].j_max = c.accum->width();
			args[ti].subdomain = ti;

#ifdef USE_PTHREAD
			pthread_attr_init(&thread_attr[ti]);
			pthread_attr_setdetachstate(&thread_attr[ti], PTHREAD_CREATE_JOINABLE);
			pthread_create(&threads[ti], &thread_attr[ti], diff_subdomain, &args[ti]);
#else
			diff_subdomain(&args[ti]);
#endif
		}

		for (int ti = 0; ti < N; ti++) {
#ifdef USE_PTHREAD
			pthread_join(threads[ti], NULL);
#endif
			diff_stat->add(args[ti].diff_stat);

			delete args[ti].diff_stat;
		}

		delete[] args;

		ui::get()->d2_align_sample_stop();

		ale_accum result = diff_stat->get_error();

		if (diff_stat_in == NULL)
			delete diff_stat;

		return result;
	}


	/*
	 * Adjust exposure for an aligned frame B against reference A.
	 *
	 * Expects full-LOD images.
	 *
	 * This function is a bit of a mess, as it reflects rather ad-hoc rules
	 * regarding what seems to work w.r.t. certainty.  Using certainty in the
	 * first pass seems to result in worse alignment, while not using certainty
	 * in the second pass results in incorrect determination of exposure.
	 *
	 * [Note that this may have been due to a bug in certainty determination
	 * within this function.]
	 */
	static void set_exposure_ratio(unsigned int m, struct scale_cluster c,
			transformation t, int ax_count, int pass_number) {

		if (_exp_register == 2) {
			/*
			 * Use metadata only.
			 */
			ale_real gain_multiplier = image_rw::exp(m).get_gain_multiplier();
			pixel multiplier = pixel(gain_multiplier, gain_multiplier, gain_multiplier);

			image_rw::exp(m).set_multiplier(multiplier);
			ui::get()->exp_multiplier(multiplier[0],
					          multiplier[1],
						  multiplier[2]);

			return;
		}

		pixel_accum asum, bsum;

		point offset = c.accum->offset();

		for (unsigned int i = 0; i < c.accum->height(); i++)
		for (unsigned int j = 0; j < c.accum->width(); j++) {

			if (ref_excluded(i, j, offset, c.ax_parameters, ax_count))
				continue;

			/*
			 * Transform
			 */

			struct point q;

			q = (c.input_scale < 1.0 && interpolant == NULL)
			  ? t.scaled_inverse_transform(
				point(i + offset[0], j + offset[1]))
			  : t.unscaled_inverse_transform(
				point(i + offset[0], j + offset[1]));

			/*
			 * Check that the transformed coordinates are within
			 * the boundaries of array c.input, that they are not
			 * subject to exclusion, and that the weight value in
			 * the accumulated array is nonzero.
			 */

			if (input_excluded(q[0], q[1], c.ax_parameters, ax_count))
				continue;

			if (q[0] >= 0
			 && q[0] <= c.input->height() - 1
			 && q[1] >= 0
			 && q[1] <= c.input->width() - 1
			 && c.defined->get_pixel(i, j).minabs_norm() != 0) { 
				pixel a = c.accum->get_pixel(i, j);
				pixel b[2];

				c.input->get_bl(q, b);

#if 1
				pixel weight = (c.aweight
					      ? c.aweight->get_pixel(i, j)
					      : pixel(1, 1, 1))
					     * ((!certainty_weights && pass_number)
					      ? c.defined->get_pixel(i, j)
					      : pixel(1, 1, 1))
					     * (pass_number
					      ? b[1]
					      : pixel(1, 1, 1));
#else
				pixel weight = pixel(1, 1, 1);
#endif

				asum += a    * weight;
				bsum += b[0] * weight;
			}
		}

		// std::cerr << (asum / bsum) << " ";
		
		pixel_accum new_multiplier;

		new_multiplier = asum / bsum * image_rw::exp(m).get_multiplier();

		if (finite(new_multiplier[0])
		 && finite(new_multiplier[1])
		 && finite(new_multiplier[2])) {
			image_rw::exp(m).set_multiplier(new_multiplier);
			ui::get()->exp_multiplier(new_multiplier[0],
					          new_multiplier[1],
						  new_multiplier[2]);
		}
	}

	/*
	 * Copy all ax parameters.
	 */
	static exclusion *copy_ax_parameters(int local_ax_count, exclusion *source) {

		exclusion *dest = (exclusion *) malloc(local_ax_count * sizeof(exclusion));

		assert (dest);

		if (!dest)
			ui::get()->memory_error("exclusion regions");

		for (int idx = 0; idx < local_ax_count; idx++)
			dest[idx] = source[idx];

		return dest;
	}

	/*
	 * Copy ax parameters according to frame.
	 */
	static exclusion *filter_ax_parameters(int frame, int *local_ax_count) {

		exclusion *dest = (exclusion *) malloc(ax_count * sizeof(exclusion));

		assert (dest);

		if (!dest)
			ui::get()->memory_error("exclusion regions");

		*local_ax_count = 0;

		for (int idx = 0; idx < ax_count; idx++) {
			if (ax_parameters[idx].x[4] > frame 
			 || ax_parameters[idx].x[5] < frame)
				continue;

			dest[*local_ax_count] = ax_parameters[idx];

			(*local_ax_count)++;
		}

		return dest;
	}

	static void scale_ax_parameters(int local_ax_count, exclusion *ax_parameters, 
			ale_pos ref_scale, ale_pos input_scale) {
		for (int i = 0; i < local_ax_count; i++) {
			ale_pos scale = (ax_parameters[i].type == exclusion::RENDER)
				      ? ref_scale
				      : input_scale;

			for (int n = 0; n < 6; n++) {
				ax_parameters[i].x[n] = (int) floor(ax_parameters[i].x[n]
						                  * scale);
			}
		}
	}

	/*
	 * Prepare the next level of detail for ordinary images.
	 */
	static const image *prepare_lod(const image *current) {
		if (current == NULL)
			return NULL;

		return current->scale_by_half("prepare_lod");
	}

	/*
	 * Prepare the next level of detail for weighted images.
	 */
	static const image *prepare_lod(const image *current, const image *weights) {
		if (current == NULL)
			return NULL;

		return current->scale_by_half(weights, "prepare_lod");
	}

	/*
	 * Prepare the next level of detail for definition maps.
	 */
	static const image *prepare_lod_def(const image *current) {
		assert(current);

		return current->defined_scale_by_half("prepare_lod_def"); 
	}

	/*
	 * Initialize the scale cluster data structure.
	 */
	static struct scale_cluster *init_clusters(int frame, ale_real scale_factor,
			const image *input_frame, unsigned int steps,
			int *local_ax_count) {

		/*
		 * Allocate memory for the array.
		 */

		struct scale_cluster *scale_clusters = 
			(struct scale_cluster *) malloc(steps * sizeof(struct scale_cluster));

		assert (scale_clusters);

		if (!scale_clusters)
			ui::get()->memory_error("alignment");

		/*
		 * Prepare images and exclusion regions for the highest level
		 * of detail.  
		 */

		scale_clusters[0].accum = reference_image;

		ui::get()->constructing_lod_clusters(0.0);
		scale_clusters[0].input_scale = scale_factor;
		if (scale_factor < 1.0 && interpolant == NULL)
			scale_clusters[0].input = input_frame->scale(scale_factor, "alignment");
		else
			scale_clusters[0].input = input_frame;

		scale_clusters[0].defined = reference_defined;
		scale_clusters[0].aweight = alignment_weights_const;
		scale_clusters[0].ax_parameters = filter_ax_parameters(frame, local_ax_count);

		scale_ax_parameters(*local_ax_count, scale_clusters[0].ax_parameters, scale_factor, 
				(scale_factor < 1.0 && interpolant == NULL) ? scale_factor : 1);

		/*
		 * Prepare reduced-detail images and exclusion
		 * regions.
		 */

		for (unsigned int step = 1; step < steps; step++) {
			ui::get()->constructing_lod_clusters(step);
			scale_clusters[step].accum = prepare_lod(scale_clusters[step - 1].accum, scale_clusters[step - 1].aweight);
			scale_clusters[step].defined = prepare_lod_def(scale_clusters[step - 1].defined);
			scale_clusters[step].aweight = prepare_lod(scale_clusters[step - 1].aweight);
			scale_clusters[step].ax_parameters 
				= copy_ax_parameters(*local_ax_count, scale_clusters[step - 1].ax_parameters);

			double sf = scale_clusters[step - 1].input_scale / 2;
			scale_clusters[step].input_scale = sf;

			if (sf >= 1.0 || interpolant != NULL) {
				scale_clusters[step].input = scale_clusters[step - 1].input;
				scale_ax_parameters(*local_ax_count, scale_clusters[step].ax_parameters, 0.5, 1);
			} else if (sf > 0.5) {
				scale_clusters[step].input = scale_clusters[step - 1].input->scale(sf, "alignment");
				scale_ax_parameters(*local_ax_count, scale_clusters[step].ax_parameters, 0.5, sf);
			} else {
				scale_clusters[step].input = scale_clusters[step - 1].input->scale(0.5, "alignment");
				scale_ax_parameters(*local_ax_count, scale_clusters[step].ax_parameters, 0.5, 0.5);
			}
		}

		return scale_clusters;
	}

	/*
	 * Destroy the first element in the scale cluster data structure.
	 */
	static void final_clusters(struct scale_cluster *scale_clusters, ale_real scale_factor,
			unsigned int steps) {

		if (scale_clusters[0].input_scale < 1.0)
			delete scale_clusters[0].input;

		free((void *)scale_clusters[0].ax_parameters);

		for (unsigned int step = 1; step < steps; step++) {
			delete scale_clusters[step].accum;
			delete scale_clusters[step].defined;
			delete scale_clusters[step].aweight;

			if (scale_clusters[step].input_scale < 1.0)
				delete scale_clusters[step].input;

			free((void *)scale_clusters[step].ax_parameters);
		}

		free(scale_clusters);
	}

	/*
	 * Calculate the centroid of a control point for the set of frames
	 * having index lower than m.  Divide by any scaling of the output.
	 */
	static point unscaled_centroid(unsigned int m, unsigned int p) {
		assert(_keep);

		point point_sum(0, 0);
		ale_accum divisor = 0;

		for(unsigned int j = 0; j < m; j++) {
			point pp = cp_array[p][j];

			if (pp.defined()) {
				point_sum += kept_t[j].transform_unscaled(pp) 
					   / kept_t[j].scale();
				divisor += 1;
			}
		}

		if (divisor == 0)
			return point::undefined();

		return point_sum / divisor;
	}

	/*
	 * Calculate centroid of this frame, and of all previous frames,
	 * from points common to both sets.
	 */
	static void centroids(unsigned int m, point *current, point *previous) {
		/*
		 * Calculate the translation
		 */
		point other_centroid(0, 0);
		point this_centroid(0, 0);
		ale_pos divisor = 0;

		for (unsigned int i = 0; i < cp_count; i++) {
			point other_c = unscaled_centroid(m, i);
			point this_c  = cp_array[i][m];

			if (!other_c.defined() || !this_c.defined())
				continue;

			other_centroid += other_c;
			this_centroid  += this_c;
			divisor        += 1;

		}

		if (divisor == 0) {
			*current = point::undefined();
			*previous = point::undefined();
			return;
		}

		*current = this_centroid / divisor;
		*previous = other_centroid / divisor;
	}

	/*
	 * Calculate the RMS error of control points for frame m, with
	 * transformation t, against control points for earlier frames.
	 */
	static ale_accum cp_rms_error(unsigned int m, transformation t) {
		assert (_keep);

		ale_accum err = 0;
		ale_accum divisor = 0;

		for (unsigned int i = 0; i < cp_count; i++)
		for (unsigned int j = 0; j < m;        j++) {
			const point *p = cp_array[i];
			point p_ref = kept_t[j].transform_unscaled(p[j]);
			point p_cur = t.transform_unscaled(p[m]);

			if (!p_ref.defined() || !p_cur.defined())
				continue;

			err += p_ref.lengthtosq(p_cur);
			divisor += 1;
		}

		return sqrt(err / divisor);
	}

	/*
	 * Align frame m against the reference.
	 *
	 * XXX: the transformation class currently combines ordinary
	 * transformations with scaling.  This is somewhat convenient for
	 * some things, but can also be confusing.  This method, _align(), is
	 * one case where special care must be taken to ensure that the scale
	 * is always set correctly (by using the 'rescale' method).
	 */
	static ale_accum _align(int m, int local_gs) {

		/*
		 * Open the input frame.
		 */


		ui::get()->loading_file();
		const image *input_frame = image_rw::open(m);

		/*
		 * Local upper/lower data, possibly dependent on image
		 * dimensions.
		 */

		ale_pos local_lower, local_upper;

		/*
		 * Select the minimum dimension as the reference.
		 */

		ale_pos reference_size = input_frame->height();
		if (input_frame->width() < reference_size)
			reference_size = input_frame->width();

		if (perturb_lower_percent)
			local_lower = perturb_lower
				    * reference_size
				    * 0.01
				    * scale_factor;
		else
			local_lower = perturb_lower;

		if (perturb_upper_percent)
			local_upper = perturb_upper
				    * reference_size
				    * 0.01
				    * scale_factor;
		else
			local_upper = perturb_upper;

		local_upper = pow(2, floor(log(local_upper) / log(2)));

		/*
		 * Logarithms aren't exact, so we divide repeatedly to discover
		 * how many steps will occur, and pass this information to the
		 * user interface.
		 */

		int step_count = 0;
		double step_variable = local_upper;
		while (step_variable >= local_lower) {
			step_variable /= 2;
			step_count++;
		}

		ui::get()->set_steps(step_count);

		ale_pos perturb = local_upper;
		transformation offset;
		ale_accum here;
		diff_stat_t *here_diff_stat = new diff_stat_t();
		int lod;

		if (_keep) {
			kept_t[latest] = latest_t;
			kept_ok[latest] = latest_ok;
		}

		/*
		 * Maximum level-of-detail.  Use a level of detail at most
		 * 2^lod_diff finer than the adjustment resolution.  lod_diff
		 * is a synonym for lod_max.
		 */

		const int lod_diff = lod_max;

		/*
		 * Determine how many levels of detail should be prepared.
		 */

		unsigned int steps = (perturb > pow(2, lod_max)) 
			           ? (unsigned int) (log(perturb) / log(2)) - lod_max + 1 : 1;


		/*
		 * Prepare multiple levels of detail.
		 */

		int local_ax_count;
		struct scale_cluster *scale_clusters = init_clusters(m,
				scale_factor, input_frame, steps,
				&local_ax_count);

		/*
		 * Initialize variables used in the main loop.
		 */

		lod = (steps - 1);

		/*
		 * Initialize the default initial transform
		 */

		if (default_initial_alignment_type == 0) {
			
			/*
			 * Follow the transformation of the original frame,
			 * setting new image dimensions.
			 */

			default_initial_alignment = orig_t;
			default_initial_alignment.set_dimensions(input_frame);

		} else if (default_initial_alignment_type == 1)

			/*
			 * Follow previous transformation, setting new image
			 * dimensions.
			 */

			default_initial_alignment.set_dimensions(input_frame);

		else
			assert(0);

		old_is_default = is_default;

		/*
		 * Scale default initial transform for lod
		 */

		default_initial_alignment.rescale(1 / pow(2, lod));

		/*
		 * Load any file-specified transformation
		 */

		offset = tload_next(tload, alignment_class == 2, default_initial_alignment, &is_default);

		transformation new_offset = offset;
		
		if (!old_is_default && !is_default && default_initial_alignment_type == 1) {

			/*
			 * Implement new delta --follow semantics.
			 *
			 * If we have a transformation T such that
			 *
			 * 	prev_final == T(prev_init)
			 *
			 * Then we also have
			 *
			 * 	current_init_follow == T(current_init)
			 *
			 * We can calculate T as follows:
			 *
			 * 	T == prev_final(prev_init^-1)
			 *
			 * Where ^-1 is the inverse operator.
			 */

			/*
			 * Ensure that the lod for the old initial and final
			 * alignments are equal to the lod for the new initial
			 * alignment.
			 */

			old_final_alignment.rescale (1 / pow(2, lod));
			old_initial_alignment.rescale(1 / pow(2, lod - old_lod));

			if (alignment_class == 0) {
				/*
				 * Translational transformations
				 */
	
				ale_pos t0 = -old_initial_alignment.eu_get(0) + old_final_alignment.eu_get(0);
				ale_pos t1 = -old_initial_alignment.eu_get(1) + old_final_alignment.eu_get(1);

				new_offset.eu_modify(0, t0);
				new_offset.eu_modify(1, t1);

			} else if (alignment_class == 1) {
				/*
				 * Euclidean transformations
				 */

				ale_pos t2 = -old_initial_alignment.eu_get(2) + old_final_alignment.eu_get(2);

				new_offset.eu_modify(2, t2);

				point p( offset.scaled_height()/2 + offset.eu_get(0) - old_initial_alignment.eu_get(0),
					 offset.scaled_width()/2 + offset.eu_get(1) - old_initial_alignment.eu_get(1) );

				p = old_final_alignment.transform_scaled(p);

				new_offset.eu_modify(0, p[0] - offset.scaled_height()/2 - offset.eu_get(0));
				new_offset.eu_modify(1, p[1] - offset.scaled_width()/2 - offset.eu_get(1));
				
			} else if (alignment_class == 2) {
				/*
				 * Projective transformations
				 */

				point p[4];

				p[0] = old_final_alignment.transform_scaled(old_initial_alignment
				     . scaled_inverse_transform(offset.transform_scaled(point(      0        ,       0       ))));
				p[1] = old_final_alignment.transform_scaled(old_initial_alignment
				     . scaled_inverse_transform(offset.transform_scaled(point(offset.scaled_height(),       0       ))));
				p[2] = old_final_alignment.transform_scaled(old_initial_alignment
				     . scaled_inverse_transform(offset.transform_scaled(point(offset.scaled_height(), offset.scaled_width()))));
				p[3] = old_final_alignment.transform_scaled(old_initial_alignment
				     . scaled_inverse_transform(offset.transform_scaled(point(      0        , offset.scaled_width()))));

				new_offset.gpt_set(p);
			}
		}

		old_initial_alignment = offset;
		old_lod = lod;
		offset = new_offset;

		struct scale_cluster si = scale_clusters[lod];
		ale_pos _mc_arg = (_mc > 0) ? (_mc * pow(2, 2 * lod))
			                    /* : ((double)_mcd_min / (si.accum->height() * si.accum->width())); */
			                    : 1;

		/*
		 * Projective adjustment value
		 */

		ale_pos adj_p = (perturb >= pow(2, lod_diff))
			     ? pow(2, lod_diff) : (double) perturb;

		/*
		 * Pre-alignment exposure adjustment
		 */

		if (_exp_register) {
			ui::get()->exposure_1();
			transformation o = offset;
			for (int k = lod; k > 0; k--)
				o.rescale(2);
			set_exposure_ratio(m, scale_clusters[0], o, local_ax_count, 0);
		}

		/*
		 * Current difference (error) value
		 */

		ui::get()->prematching();
		here = diff(si, offset, _mc_arg, local_ax_count, m, here_diff_stat);
		ui::get()->set_match(here);

		/*
		 * Current and modified barrel distortion parameters
		 */

		ale_pos current_bd[BARREL_DEGREE];
		ale_pos modified_bd[BARREL_DEGREE];
		offset.bd_get(current_bd);
		offset.bd_get(modified_bd);

		/*
		 * Translational global search step
		 */

		if (perturb >= local_lower && local_gs != 0 && local_gs != 5) {
			
			ui::get()->aligning(perturb, lod);
			
			transformation lowest_t = offset;
			ale_accum lowest_v = here;

			point min, max;

			transformation offset_p = offset;

			if (!offset_p.is_projective())
				offset_p.eu_to_gpt();

			min = max = offset_p.gpt_get(0);
			for (int p_index = 1; p_index < 4; p_index++) {
				point p = offset_p.gpt_get(p_index);
				if (p[0] < min[0])
					min[0] = p[0];
				if (p[1] < min[1])
					min[1] = p[1];
				if (p[0] > max[0])
					max[0] = p[0];
				if (p[1] > max[1])
					max[1] = p[1];
			}

			point inner_min_t = -min;
			point inner_max_t = -max + point(si.accum->height(), si.accum->width());
			point outer_min_t = -max + point(adj_p - 1, adj_p - 1);
			point outer_max_t = point(si.accum->height(), si.accum->width()) - point(adj_p, adj_p);

			if (local_gs == 1 || local_gs == 3 || local_gs == 4) {

				/*
				 * Inner
				 */

				for (ale_pos i = inner_min_t[0]; i <= inner_max_t[0]; i += adj_p)
				for (ale_pos j = inner_min_t[1]; j <= inner_max_t[1]; j += adj_p) {
					transformation t = offset;
					t.translate(point(i, j));
					ale_accum v = diff(si, t, _mc_arg, local_ax_count, m);
					unsigned int ovl = overlap(si, t, local_ax_count);

					if ((v < lowest_v && ovl >= _gs_mo) 
					 || (!finite(lowest_v) && finite(v))) {
						lowest_v = v;
						lowest_t = t;
					}
				}
			} 
			
			if (local_gs == 2 || local_gs == 3 || local_gs == -1) {

				/*
				 * Outer
				 */

				for (ale_pos i = outer_min_t[0]; i <= outer_max_t[0]; i += adj_p)
				for (ale_pos j = outer_min_t[1]; j <  inner_min_t[1]; j += adj_p) {
					transformation t = offset;
					t.translate(point(i, j));
					ale_accum v = diff(si, t, _mc_arg, local_ax_count, m);
					unsigned int ovl = overlap(si, t, local_ax_count);

					if ((v < lowest_v && ovl >= _gs_mo)
					 || (!finite(lowest_v) && finite(v))) {
						lowest_v = v;
						lowest_t = t;
					}
				}
				for (ale_pos i = outer_min_t[0]; i <= outer_max_t[0]; i += adj_p)
				for (ale_pos j = outer_max_t[1]; j >  inner_max_t[1]; j -= adj_p) {
					transformation t = offset;
					t.translate(point(i, j));
					ale_accum v = diff(si, t, _mc_arg, local_ax_count, m);
					unsigned int ovl = overlap(si, t, local_ax_count);

					if ((v < lowest_v && ovl >= _gs_mo)
					 || (!finite(lowest_v) && finite(v))) {
						lowest_v = v;
						lowest_t = t;
					}
				}
				for (ale_pos i = outer_min_t[0]; i <  inner_min_t[0]; i += adj_p)
				for (ale_pos j = outer_min_t[1]; j <= outer_max_t[1]; j += adj_p) {
					transformation t = offset;
					t.translate(point(i, j));
					ale_accum v = diff(si, t, _mc_arg, local_ax_count, m);
					unsigned int ovl = overlap(si, t, local_ax_count);

					if ((v < lowest_v && ovl >= _gs_mo)
					 || (!finite(lowest_v) && finite(v))) {
						lowest_v = v;
						lowest_t = t;
					}
				}
				for (ale_pos i = outer_max_t[0]; i >  inner_max_t[0]; i -= adj_p)
				for (ale_pos j = outer_min_t[1]; j <= outer_max_t[1]; j += adj_p) {
					transformation t = offset;
					t.translate(point(i, j));
					ale_accum v = diff(si, t, _mc_arg, local_ax_count, m);
					unsigned int ovl = overlap(si, t, local_ax_count);

					if ((v < lowest_v && ovl >= _gs_mo)
					 || (!finite(lowest_v) && finite(v))) {
						lowest_v = v;
						lowest_t = t;
					}
				}
			}

			offset = lowest_t;
			here = lowest_v;

			ui::get()->set_match(here);
		}

		/*
		 * Control point alignment
		 */

		if (local_gs == 5) {

			transformation o = offset;
			for (int k = lod; k > 0; k--)
				o.rescale(2);

			/*
			 * Determine centroid data
			 */

			point current, previous;
			centroids(m, &current, &previous);

			if (current.defined() && previous.defined()) {
				o = orig_t;
				o.set_dimensions(input_frame);
				o.translate((previous - current) * o.scale());
				current = previous;
			}

			/*
			 * Determine rotation for alignment classes other than translation.
			 */

			ale_accum lowest_error = cp_rms_error(m, o);

			ale_pos rot_lower = 2 * local_lower
					  / sqrt(pow(scale_clusters[0].input->height(), 2)
					       + pow(scale_clusters[0].input->width(),  2))
					  * 180
					  / M_PI;

			if  (alignment_class > 0)
			for (ale_pos rot = 30; rot > rot_lower; rot /= 2) 
			for (ale_pos srot = -rot; srot <= rot; srot += rot * 2) {
				int is_improved = 1;
				while (is_improved) {
					is_improved = 0;
					transformation test_t = o;
					test_t.rotate(current * o.scale(), srot);
					ale_pos test_v = cp_rms_error(m, test_t);

					if (test_v < lowest_error) {
						lowest_error = test_v;
						o = test_t;
						srot += 3 * rot;
						is_improved = 1;
					}
				}
			}

			/*
			 * Determine projective parameters through a local
			 * minimum search.
			 */

			if (alignment_class == 2) {
				ale_accum adj_p = lowest_error;

				if (adj_p < local_lower)
					adj_p = local_lower;

				while (adj_p >= local_lower) {
					transformation test_t = o;
					int is_improved = 1;
					ale_accum test_v;
					ale_accum adj_s;

					while (is_improved) {
						is_improved = 0;

						for (int i = 0; i < 4; i++)
						for (int j = 0; j < 2; j++)
						for (adj_s = -adj_p; adj_s <= adj_p; adj_s += 2 * adj_p) {

							test_t = o;

							if (perturb_type == 0)
								test_t.gpt_modify(j, i, adj_s);
							else if (perturb_type == 1)
								test_t.gr_modify(j, i, adj_s);
							else
								assert(0);

							test_v = cp_rms_error(m, test_t);

							if (test_v < lowest_error) {
								lowest_error = test_v;
								o = test_t;
								adj_s += 3 * adj_p;
								is_improved = 1;
							}
						}
					}
					adj_p /= 2;
				}
			}

			if (_exp_register)
				set_exposure_ratio(m, scale_clusters[0], o, local_ax_count, 0);

			for (int k = lod; k > 0; k--)
				o.rescale(0.5);

			offset = o;
		}

		/*
		 * Perturbation adjustment loop.  
		 */

		while (perturb >= local_lower) {

			ui::get()->aligning(perturb, lod);

			ale_pos adj_s;

			/*
			 * Orientational adjustment value in degrees.
			 *
			 * Since rotational perturbation is now specified as an
			 * arclength, we have to convert.
			 */

			ale_pos adj_o = 2 * perturb 
				          / sqrt(pow(scale_clusters[0].input->height(), 2)
					       + pow(scale_clusters[0].input->width(),  2))
					  * 180
					  / M_PI;

			/*
			 * Barrel distortion adjustment value
			 */

			ale_pos adj_b = perturb * bda_mult;

			transformation test_t;
			transformation old_offset = offset;
			ale_accum test_d;
			ale_accum old_here = here;
			diff_stat_t *old_here_diff_stat = here_diff_stat;
			here_diff_stat = new diff_stat_t();
			int found_better = 0;
			int found_reliable_better = 0;
			int found_reliable_worse = 0;
			int found_unreliable_worse = 0;

			if (alignment_class < 2 && alignment_class >= 0) {

				/* 
				 * Translational or euclidean transformation
				 */

				for (int i = 0; i < 2; i++)
				for (adj_s = -adj_p; adj_s <= adj_p; adj_s += 2 * adj_p) {

					test_t = offset;

					test_t.eu_modify(i, adj_s);

					test_d = diff(si, test_t, _mc_arg, local_ax_count, m, here_diff_stat);

					if (test_d < here || (!finite(here) && finite(test_d))) {
						found_better = 1;
						if (_mc > 0
						 || _mc_arg >= 1
						 || here_diff_stat->reliable(old_here_diff_stat, _mc_arg)) {
							found_reliable_better = 1;
							here = test_d;
							offset = test_t;
							goto done;
						}
					}

					if (_mc <= 0 && test_d > here) {
						if (_mc_arg >= 1
						 || old_here_diff_stat->reliable(here_diff_stat, _mc_arg))
							found_reliable_worse = 1;
						else
							found_unreliable_worse = 1;
					}
				}

				if (alignment_class == 1 && adj_o < rot_max)
				for (adj_s = -adj_o; adj_s <= adj_o; adj_s += 2 * adj_o) {

					test_t = offset;

					test_t.eu_modify(2, adj_s);

					test_d = diff(si, test_t, _mc_arg, local_ax_count, m, here_diff_stat);

					if (test_d < here || (!finite(here) && finite(test_d))) {
						found_better = 1;
						if (_mc > 0
						 || _mc_arg >= 1
						 || here_diff_stat->reliable(old_here_diff_stat, _mc_arg)) {
							found_reliable_better = 1;
							here = test_d;
							offset = test_t;
							goto done;
						}
					}

					if (_mc <= 0 && test_d > here) {
						if (_mc_arg >= 1
						 || old_here_diff_stat->reliable(here_diff_stat, _mc_arg))
							found_reliable_worse = 1;
						else
							found_unreliable_worse = 1;
					}
				}

			} else if (alignment_class == 2) {

				/*
				 * Projective transformation
				 */

				for (int i = 0; i < 4; i++)
				for (int j = 0; j < 2; j++)
				for (adj_s = -adj_p; adj_s <= adj_p; adj_s += 2 * adj_p) {

					test_t = offset;

					if (perturb_type == 0)
						test_t.gpt_modify(j, i, adj_s);
					else if (perturb_type == 1)
						test_t.gr_modify(j, i, adj_s);
					else
						assert(0);

					test_d = diff(si, test_t, _mc_arg, local_ax_count, m, here_diff_stat);

#if 0
					fprintf(stderr, "old\n");
					old_here_diff_stat->print_hist();
					fprintf(stderr, "new\n");
					here_diff_stat->print_hist();
#endif

					if (test_d < here || (!finite(here) && finite(test_d))) {
						// fprintf(stderr, "found better\n");
						found_better = 1;
						if (_mc > 0
						 || _mc_arg >= 1
						 || here_diff_stat->reliable(old_here_diff_stat, _mc_arg)) {
							// fprintf(stderr, "found reliable better\n");
							found_reliable_better = 1;
							here = test_d;
							offset = test_t;
							adj_s += 3 * adj_p;
						}
					}

					if (_mc <= 0 && test_d > here) {
						if (_mc_arg >= 1
						 || old_here_diff_stat->reliable(here_diff_stat, _mc_arg))
							found_reliable_worse = 1;
						else
							found_unreliable_worse = 1;
					}
				}

			} else assert(0);

			/*
			 * Barrel distortion
			 */

			if (bda_mult != 0) {

				static unsigned int d_rotation = 0;

				for (unsigned int d = 0; d < offset.bd_count(); d++)
				for (adj_s = -adj_b; adj_s <= adj_b; adj_s += 2 * adj_b) {

					unsigned int rd = (d + d_rotation) % offset.bd_count();
					d_rotation = (d_rotation + 1) % offset.bd_count();

					if (bda_rate > 0 && fabs(modified_bd[rd] + adj_s - current_bd[rd]) > bda_rate)
						continue;
				
					test_t = offset;

					test_t.bd_modify(rd, adj_s);

					test_d = diff(si, test_t, _mc_arg, local_ax_count, m, here_diff_stat);

					if (test_d < here || (!finite(here) && finite(test_d))) {
						found_better = 1;
						if (_mc > 0
						 || _mc_arg >= 1
						 || here_diff_stat->reliable(old_here_diff_stat, _mc_arg)) {
							found_reliable_better = 1;
							here = test_d;
							offset = test_t;
							modified_bd[rd] += adj_s;
							goto done;
						}
					}

					if (_mc <= 0 && test_d > here) {
						if (_mc_arg >= 1
						 || old_here_diff_stat->reliable(here_diff_stat, _mc_arg))
							found_reliable_worse = 1;
						else
							found_unreliable_worse = 1;
					}
				}
			}
			
	done:

			if (_mc_arg < 1 && _mc <= 0 
			 && (!found_reliable_worse || found_better) 
			 /* && found_unreliable_worse */
			 && !found_reliable_better) {
				_mc_arg *= 2;
				fprintf(stderr, "increasing mc to %f\n", _mc_arg);
				continue;
			}

			if (_mc <= 0
			 && found_reliable_better) {
				diff_stat_t *here_diff_stat_half = new diff_stat_t();
				diff_stat_t *old_here_diff_stat_half = new diff_stat_t();

				diff(si, offset, _mc_arg / 2, local_ax_count, m, here_diff_stat_half);
				diff(si, old_offset, _mc_arg / 2, local_ax_count, m, old_here_diff_stat_half);

				if (here_diff_stat_half->reliable(old_here_diff_stat_half, _mc_arg / 2)) {
					_mc_arg /= 2;
					fprintf(stderr, "decreasing mc to %f\n", _mc_arg);
				}
				delete here_diff_stat;
				delete old_here_diff_stat;
				here_diff_stat = here_diff_stat_half;
				old_here_diff_stat = old_here_diff_stat_half;
			}

			if (!(here < old_here) && !(!finite(old_here) && finite(here))) {
				fprintf(stderr, "increasing perturbation\n");
				perturb *= 0.5;

				if (lod > 0) {

					/* 
					 * We're about to work with images
					 * twice as large, so rescale the 
					 * transforms.
					 */

					offset.rescale(2);
					default_initial_alignment.rescale(2);

					/*
					 * Work with images twice as large
					 */

					lod--;
					si = scale_clusters[lod];

					if (_mc > 0)
						_mc_arg /= 4;
					else
						_mc_arg /= 2;

					here = diff(si, offset, _mc_arg, local_ax_count, m, here_diff_stat);
					delete old_here_diff_stat;

				} else {
					delete here_diff_stat;
					here_diff_stat = old_here_diff_stat;
					adj_p = perturb;
				}

				/*
				 * Announce that we've dropped a perturbation level.
				 */

				ui::get()->alignment_perturbation_level(perturb, lod);

			} else {
				delete old_here_diff_stat;
			}

			ui::get()->set_match(here);
		}

		if (lod > 0) {
			offset.rescale(pow(2, lod));
			default_initial_alignment.rescale(pow(2, lod));
		}

		/*
		 * Post-alignment exposure adjustment
		 */

		if (_exp_register == 1) {
			ui::get()->exposure_2();
			set_exposure_ratio(m, scale_clusters[0], offset, local_ax_count, 1);
		}

		/*
		 * Recalculate error
		 */

		ui::get()->postmatching();
		diff_stat_t *new_diff_stat = new diff_stat_t();
		here = diff(scale_clusters[0], offset, _mc_arg, local_ax_count, m, new_diff_stat);
		new_diff_stat->print_hist();
		delete new_diff_stat;
		ui::get()->set_match(here);

		/*
		 * Free the level-of-detail structures
		 */

		final_clusters(scale_clusters, scale_factor, steps);

		/*
		 * Ensure that the match meets the threshold.
		 */

		if ((1 - here) * 100 > match_threshold) {
			/*
			 * Update alignment variables
			 */
			latest_ok = 1;
			default_initial_alignment = offset;
			old_final_alignment = offset;
			ui::get()->alignment_match_ok();
		} else if (local_gs == 4) {

			/*
			 * Align with outer starting points.
			 */

			/*
			 * XXX: This probably isn't exactly the right thing to do,
			 * since variables like old_initial_value have been overwritten.
			 */

			ale_accum nested_result = _align(m, -1);

			if ((1 - nested_result) * 100 > match_threshold) {
				return nested_result;
			} else if (nested_result < here) {
				here = nested_result;
				offset = latest_t;
			}

			if (is_fail_default)
				offset = default_initial_alignment;

			ui::get()->set_match(here);
			ui::get()->alignment_no_match();
		
		} else if (local_gs == -1) {
			
			latest_ok = 0;
			latest_t = offset;
			return here;

		} else {
			if (is_fail_default)
				offset = default_initial_alignment;
			latest_ok = 0;
			ui::get()->alignment_no_match();
		}
 
                /*
                 * Write the tonal registration multiplier as a comment.
                 */
 
                pixel trm = image_rw::exp(m).get_multiplier();
                tsave_trm(tsave, trm[0], trm[1], trm[2]);

		/*
		 * Save the transformation information
		 */

		tsave_next(tsave, offset, alignment_class == 2);

		latest_t = offset;

		/*
		 * Update match statistics.
		 */

		match_sum += (1 - here) * 100;
		match_count++;
		latest = m;

		/*
		 * Close the image file.
		 */

		image_rw::close(m);

		return here;
	}

#ifdef USE_FFTW
	/*
	 * High-pass filter for frequency weights
	 */
	static void hipass(int rows, int cols, fftw_complex *inout) {
		for (int i = 0; i < rows * vert_freq_cut; i++)
		for (int j = 0; j < cols; j++)
		for (int k = 0; k < 2; k++)
			inout[i * cols + j][k] = 0;
		for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols * horiz_freq_cut; j++)
		for (int k = 0; k < 2; k++)
			inout[i * cols + j][k] = 0;
		for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
		for (int k = 0; k < 2;    k++)
			if (i / (double) rows + j / (double) cols < 2 * avg_freq_cut)
				inout[i * cols + j][k] = 0;
	}
#endif


	/*
	 * Reset alignment weights
	 */
	static void reset_weights() {

		alignment_weights_const = NULL;

		if (alignment_weights != NULL)
			delete alignment_weights;

		alignment_weights = NULL;
	}

	/*
	 * Initialize alignment weights
	 */
	static void init_weights() {
		if (alignment_weights != NULL) 
			return;

		int rows = reference_image->height();
		int cols = reference_image->width();
		int colors = reference_image->depth();

		alignment_weights = new image_ale_real(rows, cols,
				colors, "alignment_weights");

		assert (alignment_weights);

		for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			alignment_weights->set_pixel(i, j, pixel(1, 1, 1));
	}
	
	/*
	 * Update alignment weights with weight map
	 */
	static void map_update() {

		if (weight_map == NULL)
			return;

		init_weights();

		point map_offset = reference_image->offset() - weight_map->offset();

		int rows = reference_image->height();
		int cols = reference_image->width();

		for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++) {
			point map_weight_position = map_offset + point(i, j);
			if (map_weight_position[0] >= 0
			 && map_weight_position[1] >= 0
			 && map_weight_position[0] <= weight_map->height() - 1
			 && map_weight_position[1] <= weight_map->width() - 1)
				alignment_weights->pix(i, j) *= weight_map->get_bl(map_weight_position);
		}
	}

	/*
	 * Update alignment weights with an internal weight map, reflecting a
	 * summation of certainty values.  Use existing memory structures if
	 * possible.  This function updates alignment_weights_const; hence, it
	 * should not be called prior to any functions that modify the
	 * alignment_weights structure.
	 */
	static void imap_update() {
		if (alignment_weights == NULL) {
			alignment_weights_const = reference_defined;
		} else {
			int rows = reference_image->height();
			int cols = reference_image->width();

			for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				alignment_weights->pix(i, j) *= reference_defined->get_pixel(i, j);

			alignment_weights_const = alignment_weights;
		}
	}

	/*
	 * Update alignment weights with algorithmic weights
	 */
	static void wmx_update() {
#ifdef USE_UNIX

		static exposure *exp_def = new exposure_default();
		static exposure *exp_bool = new exposure_boolean();

		if (wmx_file == NULL || wmx_exec == NULL || wmx_defs == NULL)
			return;

		unsigned int rows = reference_image->height();
		unsigned int cols = reference_image->width();

		image_rw::write_image(wmx_file, reference_image);
		image_rw::write_image(wmx_defs, reference_defined, exp_bool);

		/* execute ... */
		int exit_status = 1;
		if (!fork()) {
			execlp(wmx_exec, wmx_exec, wmx_file, wmx_defs, NULL);
			ui::get()->exec_failure(wmx_exec, wmx_file, wmx_defs);
		}

		wait(&exit_status);

		if (exit_status)
			ui::get()->fork_failure("d2::align");

		image *wmx_weights = image_rw::read_image(wmx_file, exp_def);

		if (wmx_weights->height() != rows || wmx_weights->width() != cols)
			ui::get()->error("algorithmic weighting must not change image size");

		if (alignment_weights == NULL)
			alignment_weights = wmx_weights;
		else 
			for (unsigned int i = 0; i < rows; i++)
			for (unsigned int j = 0; j < cols; j++)
				alignment_weights->pix(i, j) *= wmx_weights->pix(i, j);
#endif
	}

	/*
	 * Update alignment weights with frequency weights
	 */
	static void fw_update() {
#ifdef USE_FFTW
		if (horiz_freq_cut == 0
		 && vert_freq_cut  == 0
		 && avg_freq_cut   == 0)
			return;

		/*
		 * Required for correct operation of --fwshow
		 */

		assert (alignment_weights == NULL);

		int rows = reference_image->height();
		int cols = reference_image->width();
		int colors = reference_image->depth();

		alignment_weights = new image_ale_real(rows, cols,
				colors, "alignment_weights");

		fftw_complex *inout = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * rows * cols);

		assert (inout);

		fftw_plan pf = fftw_plan_dft_2d(rows, cols,
						inout, inout, 
						FFTW_FORWARD, FFTW_ESTIMATE);

		fftw_plan pb = fftw_plan_dft_2d(rows, cols,
						inout, inout, 
						FFTW_BACKWARD, FFTW_ESTIMATE);

		for (int k = 0; k < colors; k++) {
			for (int i = 0; i < rows * cols; i++) {
				inout[i][0] = reference_image->get_pixel(i / cols, i % cols)[k];
				inout[i][1] = 0;
			}

			fftw_execute(pf);
			hipass(rows, cols, inout);
			fftw_execute(pb);

			for (int i = 0; i < rows * cols; i++) {
#if 0
				alignment_weights->pix(i / cols, i % cols)[k] = fabs(inout[i][0] / (rows * cols));
#else
				alignment_weights->pix(i / cols, i % cols)[k] = 
					sqrt(pow(inout[i][0] / (rows * cols), 2)
					   + pow(inout[i][1] / (rows * cols), 2));
#endif
			}
		}

		fftw_destroy_plan(pf);
		fftw_destroy_plan(pb);
		fftw_free(inout);

		if (fw_output != NULL)
			image_rw::write_image(fw_output, alignment_weights);
#endif
	}

	/*
	 * Update alignment to frame N.
	 */
	static void update_to(int n) {
		assert (n <= latest + 1);

		if (latest < 0) {
			const image *i = image_rw::open(n);
			int is_default;
			transformation result = alignment_class == 2
				? transformation::gpt_identity(i, scale_factor)
				: transformation::eu_identity(i, scale_factor);
			result = tload_first(tload, alignment_class == 2, result, &is_default);
			tsave_first(tsave, result, alignment_class == 2);
			image_rw::close(n);

			if (_keep > 0) {
				kept_t = (transformation *) malloc(image_rw::count()
						* sizeof(transformation));
				kept_ok = (int *) malloc(image_rw::count()
						* sizeof(int));
				assert (kept_t);
				assert (kept_ok);

				if (!kept_t || !kept_ok)
					ui::get()->memory_error("alignment");

				kept_ok[0] = 1;
				kept_t[0] = result;
			}

			latest = 0;
			latest_ok = 1;
			latest_t = result;

			default_initial_alignment = result;
			orig_t = result;
		}

		for (int i = latest + 1; i <= n; i++) {
			assert (reference != NULL);

			ui::get()->set_arender_current();
			reference->sync(i - 1);
			ui::get()->clear_arender_current();
			reference_image = reference->get_image();
			reference_defined = reference->get_defined();

			reset_weights();
			fw_update();
			wmx_update();
			map_update();

			if (certainty_weights)
				imap_update();  /* Must be called after all other _updates */

			assert (reference_image != NULL);
			assert (reference_defined != NULL);

			_align(i, _gs);
		}
	}

public:

	/*
	 * Set the control point count
	 */
	static void set_cp_count(unsigned int n) {
		assert (cp_array == NULL);

		cp_count = n;
		cp_array = (const point **) malloc(n * sizeof(const point *));
	}

	/*
	 * Set control points.
	 */
	static void set_cp(unsigned int i, const point *p) {
		cp_array[i] = p;
	}

	/*
	 * Register exposure
	 */
	static void exp_register() {
		_exp_register = 1;
	}

	/*
	 * Register exposure only based on metadata
	 */
	static void exp_meta_only() {
		_exp_register = 2;
	}

	/*
	 * Don't register exposure
	 */
	static void exp_noregister() {
		_exp_register = 0;
	}

	/*
	 * Set alignment class to translation only.  Only adjust the x and y
	 * position of images.  Do not rotate input images or perform
	 * projective transformations.
	 */
	static void class_translation() {
		alignment_class = 0;
	}

	/* 
	 * Set alignment class to Euclidean transformations only.  Adjust the x
	 * and y position of images and the orientation of the image about the
	 * image center.
	 */
	static void class_euclidean() {
		alignment_class = 1;
	}

	/*
	 * Set alignment class to perform general projective transformations.
	 * See the file gpt.h for more information about general projective
	 * transformations.
	 */
	static void class_projective() {
		alignment_class = 2;
	}

	/*
	 * Set the default initial alignment to the identity transformation.
	 */
	static void initial_default_identity() {
		default_initial_alignment_type = 0;
	}

	/*
	 * Set the default initial alignment to the most recently matched
	 * frame's final transformation.
	 */
	static void initial_default_follow() {
		default_initial_alignment_type = 1;
	}

	/*
	 * Perturb output coordinates.
	 */
	static void perturb_output() {
		perturb_type = 0;
	}

	/*
	 * Perturb source coordinates.
	 */
	static void perturb_source() {
		perturb_type = 1;
	}

	/*
	 * Frames under threshold align optimally
	 */
	static void fail_optimal() {
		is_fail_default = 0;
	}

	/*
	 * Frames under threshold keep their default alignments.
	 */
	static void fail_default() {
		is_fail_default = 1;
	}

	/*
	 * Align images with an error contribution from each color channel.
	 */
	static void all() {
		channel_alignment_type = 0;
	}

	/*
	 * Align images with an error contribution only from the green channel.
	 * Other color channels do not affect alignment.
	 */
	static void green() {
		channel_alignment_type = 1;
	}

	/*
	 * Align images using a summation of channels.  May be useful when
	 * dealing with images that have high frequency color ripples due to
	 * color aliasing.
	 */
	static void sum() {
		channel_alignment_type = 2;
	}

	/*
	 * Error metric exponent
	 */

	static void set_metric_exponent(float me) {
		metric_exponent = me;
	}

	/*
	 * Match threshold
	 */

	static void set_match_threshold(float mt) {
		match_threshold = mt;
	}

	/*
	 * Perturbation lower and upper bounds.
	 */

	static void set_perturb_lower(ale_pos pl, int plp) {
		perturb_lower = pl;
		perturb_lower_percent = plp;
	}

	static void set_perturb_upper(ale_pos pu, int pup) {
		perturb_upper = pu;
		perturb_upper_percent = pup;
	}

	/*
	 * Maximum rotational perturbation.
	 */

	static void set_rot_max(int rm) {

		/*
		 * Obtain the largest power of two not larger than rm.
		 */

		rot_max = pow(2, floor(log(rm) / log(2)));
	}

	/*
	 * Barrel distortion adjustment multiplier
	 */

	static void set_bda_mult(ale_pos m) {
		bda_mult = m;
	}

	/*
	 * Barrel distortion maximum rate of change
	 */

	static void set_bda_rate(ale_pos m) {
		bda_rate = m;
	}

	/*
	 * Level-of-detail
	 */

	static void set_lod_max(int lm) {
		lod_max = lm;
	}

	/*
	 * Set the scale factor
	 */
	static void set_scale(ale_pos s) {
		scale_factor = s;
	}

	/*
	 * Set reference rendering to align against
	 */
	static void set_reference(render *r) {
		reference = r;
	}

	/*
	 * Set the interpolant
	 */
	static void set_interpolant(filter::scaled_filter *f) {
		interpolant = f;
	}

 	/*
	 * Set alignment weights image
	 */
	static void set_weight_map(const image *i) {
		weight_map = i;
	}

	/*
	 * Set frequency cuts
	 */
	static void set_frequency_cut(double h, double v, double a) {
		horiz_freq_cut = h;
		vert_freq_cut  = v;
		avg_freq_cut   = a;
	}

	/*
	 * Set algorithmic alignment weighting
	 */
	static void set_wmx(const char *e, const char *f, const char *d) {
		wmx_exec = e;
		wmx_file = f;
		wmx_defs = d;
	}

	/*
	 * Show frequency weights
	 */
	static void set_fl_show(const char *filename) {
		fw_output = filename;
	}

	/*
	 * Set transformation file loader.
	 */
	static void set_tload(tload_t *tl) {
		tload = tl;
	}

	/*
	 * Set transformation file saver.
	 */
	static void set_tsave(tsave_t *ts) {
		tsave = ts;
	}
	
	/*
	 * Get match statistics for frame N.
	 */
	static int match(int n) {
		update_to(n);

		if (n == latest)
			return latest_ok;
		else if (_keep)
			return kept_ok[n];
		else {
			assert(0);
			exit(1);
		}
	}

	/*
	 * Message that old alignment data should be kept.
	 */
	static void keep() {
		assert (latest == -1);
		_keep = 1;
	}

	/*
	 * Get alignment for frame N.
	 */
	static transformation of(int n) {
		update_to(n);
		if (n == latest)
			return latest_t;
		else if (_keep)
			return kept_t[n];
		else {
			assert(0);
			exit(1);
		}
	}

	/*
	 * Use Monte Carlo alignment sampling with argument N.
	 */
	static void mc(ale_pos n) {
		_mc = n;
	}

	static void mcd_limit(int n) {
		_mcd_limit = n;
	}

	/*
	 * Set the certainty-weighted flag.
	 */
	static void certainty_weighted(int flag) {
		certainty_weights = flag;
	}

	/*
	 * Set the global search type.
	 */
	static void gs(const char *type) {
		if (!strcmp(type, "local")) {
			_gs = 0;
		} else if (!strcmp(type, "inner")) {
			_gs = 1;
		} else if (!strcmp(type, "outer")) {
			_gs = 2;
		} else if (!strcmp(type, "all")) {
			_gs = 3;
		} else if (!strcmp(type, "central")) {
			_gs = 4;
		} else if (!strcmp(type, "points")) {
			_gs = 5;
			keep();
		} else {
			ui::get()->error("bad global search type");
		}
	}

	/*
	 * Set the minimum overlap for global searching
	 */
	static void gs_mo(unsigned int value) {
		_gs_mo = value;
	}

	/*
	 * Don't use Monte Carlo alignment sampling.
	 */
	static void no_mc() {
		_mc = 0;
	}

	/*
	 * Set alignment exclusion regions
	 */
	static void set_exclusion(exclusion *_ax_parameters, int _ax_count) {
		ax_count = _ax_count;
		ax_parameters = _ax_parameters;
	}

	/*
	 * Get match summary statistics.
	 */
	static ale_accum match_summary() {
		return match_sum / match_count;
	}
};

#endif
