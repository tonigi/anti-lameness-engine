// Copyright 2002, 2004, 2007 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                          <dhilvert@ugcs.caltech.edu>

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
	 * Non-certainty alignment weights
	 */

	static image *alignment_weights;

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
	 * Default initial alignment type.
	 *
	 * 0. Identity transformation.
	 *
	 * 1. Most recently accepted frame's final transformation.
	 */

	static int default_initial_alignment_type;

	/*
	 * Projective group behavior
	 *
	 * 0. Perturb in output coordinates.
	 *
	 * 1. Perturb in source coordinates
	 */

	static int perturb_type;

	/*
	 * Alignment element
	 *
	 * This structure contains variables necessary for handling a
	 * multi-alignment element.  The change between the non-default old
	 * initial alignment and old final alignment is used to adjust the
	 * non-default current initial alignment.  If either the old or new
	 * initial alignment is a default alignment, the old --follow semantics
	 * are preserved.
	 */

	struct element_t {
		int is_default, old_is_default;
		int is_primary;
		transformation old_initial_alignment;
		transformation old_final_alignment;
		transformation default_initial_alignment;
		const image *input_frame;

	public:
		element_t() {
			is_default = 1;
			input_frame = NULL;
		}
	};

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

	static ale_real metric_exponent;

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
	 */

	static ale_pos _mc;

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
	 * Multi-alignment cardinality.
	 */

	static unsigned int _ma_card;

	/*
	 * Multi-alignment contiguity.
	 */

	static double _ma_cont;

	/*
	 * Minimum overlap for global searches
	 */

	static ale_accum _gs_mo;
	static int gs_mo_percent;

	/*
	 * Exclusion regions
	 */

	static exclusion *ax_parameters;
	static int ax_count;

	/*
	 * Types for scale clusters.
	 */

	struct nl_scale_cluster {
		const image *accum_max;
		const image *accum_min;
		const image *certainty_max;
		const image *certainty_min;
		const image *aweight_max;
		const image *aweight_min;
		exclusion *ax_parameters;

		ale_pos input_scale;
		const image *input_certainty_max;
		const image *input_certainty_min;
		const image *input_max;
		const image *input_min;
	};

	struct scale_cluster {
		const image *accum;
		const image *certainty;
		const image *aweight;
		exclusion *ax_parameters;

		ale_pos input_scale;
		const image *input_certainty;
		const image *input;

		nl_scale_cluster *nl_scale_clusters;
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
			 && c.certainty->get_pixel(i, j)[0] != 0)
				result++;
		}

		return result;
	}

	/*
	 * Calculate the region associated with the current multi-alignment
	 * element.
	 */
	static void calculate_element_region(transformation *t, scale_cluster si, 
			int local_ax_count) {

		unsigned int i_max = si.accum->height();
		unsigned int j_max = si.accum->width();
		point offset = si.accum->offset();

		if (si.input_scale < 1.0 && interpolant == NULL)
			t->begin_calculate_scaled_region(i_max, j_max, offset);
		else
			t->begin_calculate_unscaled_region(i_max, j_max, offset);

		for (unsigned int i = 0; i < i_max; i++) 
		for (unsigned int j = 0; j < j_max; j++) {

			if (ref_excluded(i, j, offset, si.ax_parameters, local_ax_count))
				continue;

			point q;

			while ((q = t->get_query_point((int) (i + offset[0]), 
							(int) (j + offset[1]))).defined()) {

				ale_pos ti = q[0];
				ale_pos tj = q[1];

				if (input_excluded(ti, tj, si.ax_parameters, ax_count))
					continue;

				if (ti >= 0
				 && ti <= si.input->height() - 1
				 && tj >= 0
				 && tj <= si.input->width() - 1
				 && si.certainty->get_pixel(i, j)[0] != 0) {

					assert(0);
				}
			}
		}

		t->end_calculate_region();
	}

	/*
	 * Monte carlo iteration class.
	 *
	 * Monte Carlo alignment has been used for statistical comparisons in
	 * spatial registration, and is now used for tonal registration
	 * and final match calculation.
	 */

	/*
	 * We use a random process for which the expected number of sampled
	 * pixels is +/- .000003 from the coverage in the range [.005,.995] for
	 * an image with 100,000 pixels.  (The actual number may still deviate
	 * from the expected number by more than this amount, however.)  The
	 * method is as follows:
	 *
	 * We have coverage == USE/ALL, or (expected # pixels to use)/(# total
	 * pixels).  We derive from this SKIP/USE.
	 *
	 * SKIP/USE == (SKIP/ALL)/(USE/ALL) == (1 - (USE/ALL))/(USE/ALL)
	 *
	 * Once we have SKIP/USE, we know the expected number of pixels to skip
	 * in each iteration.  We use a random selection process that provides
	 * SKIP/USE close to this calculated value.
	 *
	 * If we can draw uniformly to select the number of pixels to skip, we
	 * do.  In this case, the maximum number of pixels to skip is twice the
	 * expected number.
	 *
	 * If we cannot draw uniformly, we still assign equal probability to
	 * each of the integer values in the interval [0, 2 * (SKIP/USE)], but
	 * assign an unequal amount to the integer value ceil(2 * SKIP/USE) +
	 * 1.
	 */

	/*
	 * When reseeding the random number generator, we want the same set of
	 * pixels to be used in cases where two alignment options are compared.
	 * If we wanted to avoid bias from repeatedly utilizing the same seed,
	 * we could seed with the number of the frame most recently aligned:
	 *
	 * 	srand(latest);
	 *
	 * However, in cursory tests, it seems okay to just use the default
	 * seed of 1, and so we do this, since it is simpler; both of these
	 * approaches to reseeding achieve better results than not reseeding.
	 * (1 is the default seed according to the GNU Manual Page for
	 * rand(3).)
	 * 
	 * For subdomain calculations, we vary the seed by adding the subdomain
	 * index.
	 */

	class mc_iterate {
		ale_pos mc_max;
		unsigned int index;
		unsigned int index_max;
		int i_min;
		int i_max;
		int j_min;
		int j_max;

		rng_t rng;

	public:
		mc_iterate(int _i_min, int _i_max, int _j_min, int _j_max, unsigned int subdomain) 
				: rng() {

			ale_pos coverage;

			i_min = _i_min;
			i_max = _i_max;
			j_min = _j_min;
			j_max = _j_max;

			index_max = (i_max - i_min) * (j_max - j_min);

			if (index_max < 500 || _mc > 100 || _mc <= 0)
				coverage = 1;
			else
				coverage = _mc / 100;

			double su = (1 - coverage) / coverage;

			mc_max = (floor(2*su) * (1 + floor(2*su)) + 2*su)
			       / (2 + 2 * floor(2*su) - 2*su);

			rng.seed(1 + subdomain);

#define FIXED16 3
#if ALE_COORDINATES == FIXED16
			/*
			 * XXX: This calculation might not yield the correct
			 * expected value.
			 */
			index = -1 + (int) ceil(((ale_pos) mc_max+1) 
				   / (ale_pos) ( (1 + 0xffffff)
				                 / (1 + (rng.get() & 0xffffff))));
#else
                        index = -1 + (int) ceil((ale_accum) (mc_max+1) 
                                   * ( (1 + ((ale_accum) (rng.get())) ) 
                                     / (1 + ((ale_accum) RAND_MAX)) ));
#endif
#undef FIXED16
		}

		int get_i() {
			return index / (j_max - j_min) + i_min;
		}

		int get_j() {
			return index % (j_max - j_min) + j_min;
		}

		void operator++(int whats_this_for) {
#define FIXED16 3
#if ALE_COORDINATES == FIXED16
			index += (int) ceil(((ale_pos) mc_max+1) 
				   / (ale_pos) ( (1 + 0xffffff)
				                 / (1 + (rng.get() & 0xffffff))));
#else
                        index += (int) ceil((ale_accum) (mc_max+1) 
                               * ( (1 + ((ale_accum) (rng.get())) ) 
                                 / (1 + ((ale_accum) RAND_MAX)) ));
#endif
#undef FIXED16
		}

		int done() {
			return (index >= index_max);
		}
	};

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
	
		struct run {

			transformation offset;
			ale_pos perturb;

			ale_accum result;
			ale_accum divisor;

			point max, min;
			ale_accum centroid[2], centroid_divisor;
			ale_accum de_centroid[2], de_centroid_v, de_sum;

			void init() {

				result = 0;
				divisor = 0;

				min = point::posinf();
				max = point::neginf();

				centroid[0] = 0;
				centroid[1] = 0;
				centroid_divisor = 0;

				de_centroid[0] = 0;
				de_centroid[1] = 0;

				de_centroid_v = 0;

				de_sum = 0;
			}

			void init(transformation _offset, ale_pos _perturb) {
				offset = _offset;
				perturb = _perturb;
				init();
			}

			/*
			 * Required for STL sanity.
			 */
			run() : offset() {
				init();
			}

			run(transformation _offset, ale_pos _perturb) : offset() {
				init(_offset, _perturb);
			}

			void add(const run &_run) {
				result += _run.result;
				divisor += _run.divisor;

				for (int d = 0; d < 2; d++) {
					if (min[d] > _run.min[d])
						min[d] = _run.min[d];
					if (max[d] < _run.max[d])
						max[d] = _run.max[d];
					centroid[d] += _run.centroid[d];
					de_centroid[d] += _run.de_centroid[d];
				}

				centroid_divisor += _run.centroid_divisor;
				de_centroid_v += _run.de_centroid_v;
				de_sum += _run.de_sum;
			}

			run(const run &_run) : offset() {

				/*
				 * Initialize
				 */
				init(_run.offset, _run.perturb);

				/*
				 * Add
				 */
				add(_run);
			}

			run &operator=(const run &_run) {

				/*
				 * Initialize
				 */
				init(_run.offset, _run.perturb);

				/*
				 * Add
				 */
				add(_run);

				return *this;
			}

			~run() {
			}

			ale_accum get_error() const {
				return pow(result / divisor, 1/(ale_accum) metric_exponent);
			}

			void sample(int f, scale_cluster c, int i, int j, point t, point u,
					const run &comparison) {

				pixel pa = c.accum->get_pixel(i, j);

				ale_real this_result[2] = { 0, 0 };
				ale_real this_divisor[2] = { 0, 0 };

				pixel p[2];
				pixel weight[2];
				weight[0] = pixel(1, 1, 1);
				weight[1] = pixel(1, 1, 1);

				if (interpolant != NULL) {
					interpolant->filtered(i, j, &p[0], &weight[0], 1, f);

					if (weight[0].min_norm() > ale_real_weight_floor) {
						p[0] /= weight[0];
					} else {
						return;
					}

				} else {
                                        p[0] = c.input->get_bl(t);
				}

				if (u.defined()) {
					p[1] = c.input->get_bl(u);
				}


				/*
				 * Handle certainty.
				 */

				if (certainty_weights == 1) {

					/*
					 * For speed, use arithmetic interpolation (get_bl(.))
					 * instead of geometric (get_bl(., 1))
					 */

					weight[0] *= c.input_certainty->get_bl(t);
					if (u.defined())
						weight[1] *= c.input_certainty->get_bl(u);
					weight[0] *= c.certainty->get_pixel(i, j);
					weight[1] *= c.certainty->get_pixel(i, j);
				}

				if (c.aweight != NULL) {
					weight[0] *= c.aweight->get_pixel(i, j);
					weight[1] *= c.aweight->get_pixel(i, j);
				}

				/*
				 * Update sampling area statistics
				 */

				if (min[0] > i)
					min[0] = i;
				if (min[1] > j)
					min[1] = j;
				if (max[0] < i)
					max[0] = i;
				if (max[1] < j)
					max[1] = j;

				centroid[0] += (weight[0][0] + weight[0][1] + weight[0][2]) * i;
				centroid[1] += (weight[0][0] + weight[0][1] + weight[0][2]) * j;
				centroid_divisor += (weight[0][0] + weight[0][1] + weight[0][2]);

				/*
				 * Determine alignment type.
				 */

				for (int m = 0; m < (u.defined() ? 2 : 1); m++)
				if (channel_alignment_type == 0) {
					/*
					 * Align based on all channels.
					 */


					for (int k = 0; k < 3; k++) {
						ale_real achan = pa[k];
						ale_real bchan = p[m][k];

						this_result[m] += weight[m][k] * pow(fabs(achan - bchan), metric_exponent);
						this_divisor[m] += weight[m][k] * pow(achan > bchan ? achan : bchan, metric_exponent);
					}
				} else if (channel_alignment_type == 1) {
					/*
					 * Align based on the green channel.
					 */

					ale_real achan = pa[1];
					ale_real bchan = p[m][1];

					this_result[m] = weight[m][1] * pow(fabs(achan - bchan), metric_exponent);
					this_divisor[m] = weight[m][1] * pow(achan > bchan ? achan : bchan, metric_exponent);
				} else if (channel_alignment_type == 2) {
					/*
					 * Align based on the sum of all channels.
					 */

					ale_real asum = 0;
					ale_real bsum = 0;
					ale_real wsum = 0;

					for (int k = 0; k < 3; k++) {
						asum += pa[k];
						bsum += p[m][k];
						wsum += weight[m][k] / 3;
					}

					this_result[m] = wsum * pow(fabs(asum - bsum), metric_exponent);
					this_divisor[m] = wsum * pow(asum > bsum ? asum : bsum, metric_exponent);
				}

				if (u.defined()) {
//					ale_real de = fabs(this_result[0] / this_divisor[0]
//						         - this_result[1] / this_divisor[1]);
					ale_real de = fabs(this_result[0] - this_result[1]);

					de_centroid[0] += de * (ale_real) i;
					de_centroid[1] += de * (ale_real) j;

					de_centroid_v += de * (ale_real) t.lengthto(u);

					de_sum += de;
				}

				result += (this_result[0]);
				divisor += (this_divisor[0]);
			}

			void rescale(ale_pos scale) {
				offset.rescale(scale);

				de_centroid[0] *= scale;
				de_centroid[1] *= scale;
				de_centroid_v *= scale;
			}

			point get_centroid() {
				point result = point(centroid[0] / centroid_divisor, centroid[1] / centroid_divisor);

				assert (finite(centroid[0]) 
				     && finite(centroid[1]) 
				     && (result.defined() || centroid_divisor == 0));

				return result;
			}

			point get_error_centroid() {
				point result = point(de_centroid[0] / de_sum, de_centroid[1] / de_sum);
				return result;
			}


			ale_pos get_error_perturb() {
				ale_pos result = de_centroid_v / de_sum;

				return result;
			}

		};

		/*
		 * When non-empty, runs.front() is best, runs.back() is
		 * testing.
		 */

		std::vector<run> runs;

		/*
		 * old_runs stores the latest available perturbation set for
		 * each multi-alignment element.
		 */

		typedef std::pair<unsigned int, unsigned int> run_index;
		std::map<run_index, run> old_runs;

		static void *diff_subdomain(void *args);

		struct subdomain_args {
			struct scale_cluster c;
			std::vector<run> runs;
			int ax_count;
			int f;
			int i_min, i_max, j_min, j_max;
			int subdomain;
		};

		int get_current_index() const {
			assert(runs.size());
			return runs[0].offset.get_current_index();
		}

		struct scale_cluster si;
		int ax_count;
		int frame;

		std::vector<ale_pos> perturb_multipliers;

public:
		void diff(struct scale_cluster c, ale_pos perturb, 
				transformation t, 
				int _ax_count, int f) {

			if (runs.size() == 2)
				runs.pop_back();

			runs.push_back(run(t, perturb));

			si = c;
			ax_count = _ax_count;
			frame = f;

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
				args[ti].runs = runs;
				args[ti].ax_count = ax_count;
				args[ti].f = f;
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
				runs.back().add(args[ti].runs.back());
			}

			delete[] args;

			ui::get()->d2_align_sample_stop();

		}

	private:
		void rediff() {
			std::vector<transformation> t_array;
			std::vector<ale_pos> p_array;

			for (unsigned int r = 0; r < runs.size(); r++) {
				t_array.push_back(runs[r].offset);
				p_array.push_back(runs[r].perturb);
			}

			runs.clear();

			for (unsigned int r = 0; r < t_array.size(); r++)
				diff(si, p_array[r], t_array[r], ax_count, frame);
		}


	public:
		int better() {
			assert(runs.size() >= 2);
			assert(runs[0].offset.scale() == runs[1].offset.scale());

			return (runs[1].get_error() < runs[0].get_error() 
			     || (!finite(runs[0].get_error()) && finite(runs[1].get_error())));
		}

		diff_stat_t() : runs(), old_runs(), perturb_multipliers() {
		}

		run_index get_run_index(unsigned int perturb_index) {
			return run_index(get_current_index(), perturb_index);
		}

		run &get_run(unsigned int perturb_index) {
			run_index index = get_run_index(perturb_index);

			assert(old_runs.count(index));
			return old_runs[index];
		}

		void rescale(ale_pos scale, scale_cluster _si) {
			assert(runs.size() == 1);

			si = _si;

			runs[0].rescale(scale);
			
			rediff();
		}

		void push_element() {
			assert(runs.size() == 1);

			runs[0].offset.push_element();

			rediff();
		}

		unsigned int get_current_index() {
			assert (runs.size() > 0);

			return runs[0].offset.get_current_index();
		}

		void set_current_index(unsigned int i) {
			assert(runs.size() == 1);
			runs[0].offset.set_current_index(i);
			rediff();
		}

		void calculate_element_region() {
			assert(runs.size() == 1);

			if (get_offset().get_current_index() > 0
			 && get_offset().is_nontrivial()) 
				align::calculate_element_region(&runs[0].offset, si, ax_count);
		}

		~diff_stat_t() {
		}

		diff_stat_t &operator=(const diff_stat_t &dst) {
			/*
			 * Copy run information.
			 */
			runs = dst.runs;
			old_runs = dst.old_runs;

			/*
			 * Copy diff variables
			 */
			si = dst.si;
			ax_count = dst.ax_count;
			frame = dst.frame;
			perturb_multipliers = dst.perturb_multipliers;

			return *this;
		}

		diff_stat_t(const diff_stat_t &dst) : runs(), old_runs(),
		                                      perturb_multipliers() {
			operator=(dst);
		}

		ale_accum get_result() {
			assert(runs.size() == 1);
			return runs[0].result;
		}

		ale_accum get_divisor() {
			assert(runs.size() == 1);
			return runs[0].divisor;
		}

		transformation get_offset() {
			assert(runs.size() == 1);
			return runs[0].offset;
		}

		int operator!=(diff_stat_t &param) {
			return (get_error() != param.get_error());
		}

		int operator==(diff_stat_t &param) {
			return !(operator!=(param));
		}

		ale_pos get_error_perturb() {
			assert(runs.size() == 1);
			return runs[0].get_error_perturb();
		}

		ale_accum get_error() const {
			assert(runs.size() == 1);
			return runs[0].get_error();
		}

	public:
		/*
		 * Get the set of transformations produced by a given perturbation
		 */
		void get_perturb_set(std::vector<transformation> *set, 
				ale_pos adj_p, ale_pos adj_o, ale_pos adj_b, 
				ale_pos *current_bd, ale_pos *modified_bd,
				std::vector<ale_pos> multipliers = std::vector<ale_pos>()) {

			assert(runs.size() == 1);

			transformation test_t;

			/* 
			 * Translational or euclidean transformation
			 */

			for (unsigned int i = 0; i < 2; i++)
			for (unsigned int s = 0; s < 2; s++) {

				if (!multipliers.size())
					multipliers.push_back(1);

				assert(finite(multipliers[0]));

				test_t = get_offset();

				// test_t.eu_modify(i, (s ? -adj_p : adj_p) * multipliers[0]);
				test_t.translate((i ? point(1, 0) : point(0, 1))
				               * (s ? -adj_p : adj_p)
					       * multipliers[0]);
				
				test_t.snap(adj_p / 2);

				set->push_back(test_t);
				multipliers.erase(multipliers.begin());

			}

			if (alignment_class > 0)
			for (unsigned int s = 0; s < 2; s++) {

				if (!multipliers.size())
					multipliers.push_back(1);

				assert(multipliers.size());
				assert(finite(multipliers[0]));

				if (!(adj_o * multipliers[0] < rot_max))
					return;

				ale_pos adj_s = (s ? 1 : -1) * adj_o * multipliers[0];

				test_t = get_offset();

				run_index ori = get_run_index(set->size());
				point centroid = point::undefined();

				if (!old_runs.count(ori))
					ori = get_run_index(0);

				if (!centroid.finite() && old_runs.count(ori)) {
					centroid = old_runs[ori].get_error_centroid();

					if (!centroid.finite())
						centroid = old_runs[ori].get_centroid();

					centroid *= test_t.scale() 
					          / old_runs[ori].offset.scale();
				}

				if (!centroid.finite() && !test_t.is_projective()) {
					test_t.eu_modify(2, adj_s);
				} else if (!centroid.finite()) {
					centroid = point(si.input->height() / 2, 
					                 si.input->width() / 2);

					test_t.rotate(centroid + si.accum->offset(),
							adj_s);
				} else {
					test_t.rotate(centroid + si.accum->offset(), 
							adj_s);
				}

				test_t.snap(adj_p / 2);

				set->push_back(test_t);
				multipliers.erase(multipliers.begin());
			}

			if (alignment_class == 2) {

				/*
				 * Projective transformation
				 */

				for (unsigned int i = 0; i < 4; i++)
				for (unsigned int j = 0; j < 2; j++)
				for (unsigned int s = 0; s < 2; s++) {

					if (!multipliers.size())
						multipliers.push_back(1);

					assert(multipliers.size());
					assert(finite(multipliers[0]));

					ale_pos adj_s = (s ? -1 : 1) * adj_p * multipliers [0];

					test_t = get_offset();

					if (perturb_type == 0)
						test_t.gpt_modify(j, i, adj_s);
					else if (perturb_type == 1)
						test_t.gr_modify(j, i, adj_s);
					else
						assert(0);

					test_t.snap(adj_p / 2);

					set->push_back(test_t);
					multipliers.erase(multipliers.begin());
				}

			}

			/*
			 * Barrel distortion
			 */

			if (bda_mult != 0 && adj_b != 0) {

				for (unsigned int d = 0; d < get_offset().bd_count(); d++)
				for (unsigned int s = 0; s < 2; s++) {

					if (!multipliers.size())
						multipliers.push_back(1);

					assert (multipliers.size());
					assert (finite(multipliers[0]));

					ale_pos adj_s = (s ? -1 : 1) * adj_b * multipliers[0];

					if (bda_rate > 0 && fabs(modified_bd[d] + adj_s - current_bd[d]) > bda_rate)
						continue;
				
					transformation test_t = get_offset();

					test_t.bd_modify(d, adj_s);

					set->push_back(test_t);
				}
			}
		}

		void confirm() {
			assert(runs.size() == 2);
			runs[0] = runs[1];
			runs.pop_back();
		}

		void discard() {
			assert(runs.size() == 2);
			runs.pop_back();
		}

		void perturb_test(ale_pos perturb, ale_pos adj_p, ale_pos adj_o, ale_pos adj_b, 
				ale_pos *current_bd, ale_pos *modified_bd, int stable) {

			assert(runs.size() == 1);

			std::vector<transformation> t_set;

			if (perturb_multipliers.size() == 0) {
				get_perturb_set(&t_set, adj_p, adj_o, adj_b, 
						current_bd, modified_bd);

				for (unsigned int i = 0; i < t_set.size(); i++) {
					diff_stat_t test = *this;

					test.diff(si, perturb, t_set[i], ax_count, frame);

					test.confirm();

					if (finite(test.get_error_perturb()))
						perturb_multipliers.push_back(adj_p / test.get_error_perturb());
					else
						perturb_multipliers.push_back(1);

				}

				t_set.clear();
			}

			get_perturb_set(&t_set, adj_p, adj_o, adj_b, current_bd, modified_bd,
					perturb_multipliers);

			int found_unreliable = 1;
			std::vector<int> tested(t_set.size(), 0);

			for (unsigned int i = 0; i < t_set.size(); i++) {
				run_index ori = get_run_index(i);

				/*
				 * Check for stability
				 */
				if (stable
				 && old_runs.count(ori)
				 && old_runs[ori].offset == t_set[i])
					tested[i] = 1;
			}

			std::vector<ale_pos> perturb_multipliers_original = perturb_multipliers;

			while (found_unreliable) {

				found_unreliable = 0;

				for (unsigned int i = 0; i < t_set.size(); i++) {

					if (tested[i])
						continue;

					diff(si, perturb, t_set[i], ax_count, frame);

					if (!(i < perturb_multipliers.size())
					 || !finite(perturb_multipliers[i])) {

						perturb_multipliers.resize(i + 1);

						if (finite(perturb_multipliers[i])
						 && finite(adj_p)
						 && finite(adj_p / runs[1].get_error_perturb())) {
							perturb_multipliers[i] = 
								adj_p / runs[1].get_error_perturb();

							found_unreliable = 1;
						} else
							perturb_multipliers[i] = 1;

						continue;
					}

					run_index ori = get_run_index(i);

					if (old_runs.count(ori) == 0)
						old_runs.insert(std::pair<run_index, run>(ori, runs[1]));
					else 
						old_runs[ori] = runs[1];

					if (finite(perturb_multipliers_original[i])
					 && finite(runs[1].get_error_perturb())
					 && finite(adj_p)
					 && finite(perturb_multipliers_original[i] * adj_p / runs[1].get_error_perturb()))
						perturb_multipliers[i] = perturb_multipliers_original[i]
							* adj_p / runs[1].get_error_perturb();
					else
						perturb_multipliers[i] = 1;

					tested[i] = 1;

					if (better()
					 && runs[1].get_error() < runs[0].get_error()
					 && perturb_multipliers[i] 
					  / perturb_multipliers_original[i] < 2) {
						runs[0] = runs[1];
						runs.pop_back();
						return;
					}

				}

				if (runs.size() > 1)
					runs.pop_back();

				if (!found_unreliable)
					return;
			}
		}

		/*
		 * Attempt to make the current element non-trivial, by finding a nearby
		 * alignment admitting a non-empty element region.
		 */
		void make_element_nontrivial(ale_pos adj_p, ale_pos adj_o) {
			assert(runs.size() == 1);

			transformation *t = &runs[0].offset;

			if (t->is_nontrivial())
				return;

			calculate_element_region();

			if (t->is_nontrivial())
				return;

			std::vector<transformation> t_set;
			get_perturb_set(&t_set, adj_p, adj_o, 0, NULL, NULL);

			for (unsigned int i = 0; i < t_set.size(); i++) {

				align::calculate_element_region(&t_set[i], si, ax_count);

				if (t_set[i].is_nontrivial()) {
					*t = t_set[i];
					return;
				}
			}
		}

	};


	/*
	 * Adjust exposure for an aligned frame B against reference A.
	 *
	 * Expects full-LOD images.
	 *
	 * Note: This method does not use any weighting, by certainty or
	 * otherwise, in the first exposure registration pass, as any bias of
	 * weighting according to color may also bias the exposure registration
	 * result; it does use weighting, including weighting by certainty
	 * (even if certainty weighting is not specified), in the second pass,
	 * under the assumption that weighting by certainty improves handling
	 * of out-of-range highlights, and that bias of exposure measurements
	 * according to color may generally be less harmful after spatial
	 * registration has been performed.
	 */
	class exposure_ratio_iterate : public thread::decompose_domain {
		pixel_accum *asums;
		pixel_accum *bsums;
		pixel_accum *asum;
		pixel_accum *bsum;
		struct scale_cluster c;
		transformation t;
		int ax_count;
		int pass_number;
	protected:
		void prepare_subdomains(unsigned int N) {
			asums = new pixel_accum[N];
			bsums = new pixel_accum[N];
		}
		void subdomain_algorithm(unsigned int thread,
				int i_min, int i_max, int j_min, int j_max) {

			point offset = c.accum->offset();

			for (mc_iterate m(i_min, i_max, j_min, j_max, thread); !m.done(); m++) {
				
				unsigned int i = (unsigned int) m.get_i();
				unsigned int j = (unsigned int) m.get_j();

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
				 && ((pixel) c.certainty->get_pixel(i, j)).minabs_norm() != 0) { 
					pixel a = c.accum->get_pixel(i, j);
					pixel b;

					b = c.input->get_bl(q);

					pixel weight = ((c.aweight && pass_number)
						      ? (pixel) c.aweight->get_pixel(i, j)
						      : pixel(1, 1, 1))
						     * (pass_number
						      ? psqrt(c.certainty->get_pixel(i, j)
						            * c.input_certainty->get_bl(q, 1))
						      : pixel(1, 1, 1));

					asums[thread] += a * weight;
					bsums[thread] += b * weight;
				}
			}
		}

		void finish_subdomains(unsigned int N) {
			for (unsigned int n = 0; n < N; n++) {
				*asum += asums[n];
				*bsum += bsums[n];
			}
			delete[] asums;
			delete[] bsums;
		}
	public:
		exposure_ratio_iterate(pixel_accum *_asum,
				       pixel_accum *_bsum,
				       struct scale_cluster _c,
				       transformation _t,
				       int _ax_count,
				       int _pass_number) : decompose_domain(0, _c.accum->height(), 
				                                            0, _c.accum->width()){

			asum = _asum;
			bsum = _bsum;
			c = _c;
			t = _t;
			ax_count = _ax_count;
			pass_number = _pass_number;
		}
	};

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

		pixel_accum asum(0, 0, 0), bsum(0, 0, 0);

		exposure_ratio_iterate eri(&asum, &bsum, c, t, ax_count, pass_number);
		eri.run();

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
				ax_parameters[i].x[n] = ax_parameters[i].x[n] * scale;
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
	 * Prepare the next level of detail for definition maps.
	 */
	static const image *prepare_lod_def(const image *current) {
		if (current == NULL)
			return NULL;

		return current->defined_scale_by_half("prepare_lod_def"); 
	}

	/*
	 * Initialize scale cluster data structures.
	 */

	static void init_nl_cluster(struct scale_cluster *sc) {
	}

	static struct scale_cluster *init_clusters(int frame, ale_pos scale_factor,
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

		scale_clusters[0].certainty = reference_defined;
		scale_clusters[0].aweight = alignment_weights;
		scale_clusters[0].ax_parameters = filter_ax_parameters(frame, local_ax_count);

		/*
		 * Allocate and determine input frame certainty.
		 */

		if (scale_clusters[0].input->get_bayer() != IMAGE_BAYER_NONE) {
			scale_clusters[0].input_certainty = new_image_bayer_ale_real(
					scale_clusters[0].input->height(),
					scale_clusters[0].input->width(),
					scale_clusters[0].input->depth(),
					scale_clusters[0].input->get_bayer());
		} else {
			scale_clusters[0].input_certainty = scale_clusters[0].input->clone("certainty");
		}

		for (unsigned int i = 0; i < scale_clusters[0].input_certainty->height(); i++)
		for (unsigned int j = 0; j < scale_clusters[0].input_certainty->width(); j++)
		for (unsigned int k = 0; k < 3; k++)
		if (scale_clusters[0].input->get_channels(i, j) & (1 << k))
			((image *) scale_clusters[0].input_certainty)->set_chan(i, j, k,
				scale_clusters[0].input->
					exp().confidence(scale_clusters[0].input->get_pixel(i, j))[k]);

		scale_ax_parameters(*local_ax_count, scale_clusters[0].ax_parameters, scale_factor, 
				(scale_factor < 1.0 && interpolant == NULL) ? scale_factor : (ale_pos) 1);

		init_nl_cluster(&(scale_clusters[0]));

		/*
		 * Prepare reduced-detail images and exclusion
		 * regions.
		 */

		for (unsigned int step = 1; step < steps; step++) {
			ui::get()->constructing_lod_clusters(step);
			scale_clusters[step].accum = prepare_lod(scale_clusters[step - 1].accum);
			scale_clusters[step].certainty = prepare_lod_def(scale_clusters[step - 1].certainty);
			scale_clusters[step].aweight = prepare_lod_def(scale_clusters[step - 1].aweight);
			scale_clusters[step].ax_parameters 
				= copy_ax_parameters(*local_ax_count, scale_clusters[step - 1].ax_parameters);

			double sf = scale_clusters[step - 1].input_scale / 2;
			scale_clusters[step].input_scale = sf;

			if (sf >= 1.0 || interpolant != NULL) {
				scale_clusters[step].input = scale_clusters[step - 1].input;
				scale_clusters[step].input_certainty = scale_clusters[step - 1].input_certainty;
				scale_ax_parameters(*local_ax_count, scale_clusters[step].ax_parameters, 0.5, 1);
			} else if (sf > 0.5) {
				scale_clusters[step].input = scale_clusters[step - 1].input->scale(sf, "alignment");
				scale_clusters[step].input_certainty = scale_clusters[step - 1].input->scale(sf, "alignment", 1);
				scale_ax_parameters(*local_ax_count, scale_clusters[step].ax_parameters, 0.5, sf);
			} else {
				scale_clusters[step].input = scale_clusters[step - 1].input->scale(0.5, "alignment");
				scale_clusters[step].input_certainty = scale_clusters[step - 1].input_certainty->scale(0.5, "alignment", 1);
				scale_ax_parameters(*local_ax_count, scale_clusters[step].ax_parameters, 0.5, 0.5);
			}

			init_nl_cluster(&(scale_clusters[step]));
		}

		return scale_clusters;
	}

	/*
	 * Destroy the first element in the scale cluster data structure.
	 */
	static void final_clusters(struct scale_cluster *scale_clusters, ale_pos scale_factor,
			unsigned int steps) {

		if (scale_clusters[0].input_scale < 1.0) {
			delete scale_clusters[0].input;
			delete scale_clusters[0].input_certainty;
		}

		free((void *)scale_clusters[0].ax_parameters);

		for (unsigned int step = 1; step < steps; step++) {
			delete scale_clusters[step].accum;
			delete scale_clusters[step].certainty;
			delete scale_clusters[step].aweight;

			if (scale_clusters[step].input_scale < 1.0) {
				delete scale_clusters[step].input;
				delete scale_clusters[step].input_certainty;
			}

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
	static ale_pos cp_rms_error(unsigned int m, transformation t) {
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

		return (ale_pos) sqrt(err / divisor);
	}

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
	static transformation follow(element_t *element, transformation offset, int lod) {

		transformation new_offset = offset;

		/*
		 * Criteria for using following.
		 */

		if (!element->old_is_default && !element->is_default && 
				default_initial_alignment_type == 1) {
			/*
			 * Ensure that the lod for the old initial and final
			 * alignments are equal to the lod for the new initial
			 * alignment.
			 */

			ui::get()->following();

			for (offset.set_current_index(0),
			     element->old_initial_alignment.set_current_index(0),
			     element->old_final_alignment.set_current_index(0),
			     new_offset.set_current_index(0);

			     offset.get_current_index() < _ma_card;

			     offset.push_element(),
			     new_offset.push_element()) {

				if (alignment_class == 0) {
					/*
					 * Translational transformations
					 */
		
					ale_pos t0 = -element->old_initial_alignment.eu_get(0) 
						   +  element->old_final_alignment.eu_get(0);
					ale_pos t1 = -element->old_initial_alignment.eu_get(1) 
						   +  element->old_final_alignment.eu_get(1);

					new_offset.eu_modify(0, t0);
					new_offset.eu_modify(1, t1);

				} else if (alignment_class == 1) {
					/*
					 * Euclidean transformations
					 */

					ale_pos t2 = -element->old_initial_alignment.eu_get(2) 
						   +  element->old_final_alignment.eu_get(2);

					new_offset.eu_modify(2, t2);

					point p( offset.scaled_height()/2 + offset.eu_get(0) - element->old_initial_alignment.eu_get(0),
						 offset.scaled_width()/2 + offset.eu_get(1) - element->old_initial_alignment.eu_get(1) );

					p = element->old_final_alignment.transform_scaled(p);

					new_offset.eu_modify(0, p[0] - offset.scaled_height()/2 - offset.eu_get(0));
					new_offset.eu_modify(1, p[1] - offset.scaled_width()/2 - offset.eu_get(1));
					
				} else if (alignment_class == 2) {
					/*
					 * Projective transformations
					 */

					point p[4];

					p[0] = element->old_final_alignment.transform_scaled(element->old_initial_alignment
					     . scaled_inverse_transform(offset.get_current_element().transform_scaled(point(      0        ,       0       ))));
					p[1] = element->old_final_alignment.transform_scaled(element->old_initial_alignment
					     . scaled_inverse_transform(offset.get_current_element().transform_scaled(point(offset.scaled_height(),       0       ))));
					p[2] = element->old_final_alignment.transform_scaled(element->old_initial_alignment
					     . scaled_inverse_transform(offset.get_current_element().transform_scaled(point(offset.scaled_height(), offset.scaled_width()))));
					p[3] = element->old_final_alignment.transform_scaled(element->old_initial_alignment
					     . scaled_inverse_transform(offset.get_current_element().transform_scaled(point(      0        , offset.scaled_width()))));

					new_offset.gpt_set(p);
				}
			}

			ui::get()->set_offset(offset);
		}

		return new_offset;
	}

	static void test_global(diff_stat_t *here, scale_cluster si, transformation t, 
			int local_ax_count, int m, ale_accum local_gs_mo, ale_pos perturb) {

		diff_stat_t test(*here);

		test.diff(si, perturb, t, local_ax_count, m);

		unsigned int ovl = overlap(si, t, local_ax_count);

		if (ovl >= local_gs_mo && test.better()) {
			test.confirm();
			*here = test;
			ui::get()->set_match(here->get_error());
			ui::get()->set_offset(here->get_offset());
		} else {
			test.discard();
		}

	}

	/*
	 * Get the set of global transformations for a given density
	 */
	static void test_globals(diff_stat_t *here, 
			scale_cluster si, transformation t, int local_gs, ale_pos adj_p,
			int local_ax_count, int m, ale_accum local_gs_mo, ale_pos perturb) {

		transformation offset = t;

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

		if (local_gs == 1 || local_gs == 3 || local_gs == 4 || local_gs == 6) {

			/*
			 * Inner
			 */

			for (ale_pos i = inner_min_t[0]; i <= inner_max_t[0]; i += adj_p)
			for (ale_pos j = inner_min_t[1]; j <= inner_max_t[1]; j += adj_p) {
				transformation test_t = offset;
				test_t.translate(point(i, j));
				test_global(here, si, test_t, local_ax_count, m, local_gs_mo, perturb);
			}
		} 
		
		if (local_gs == 2 || local_gs == 3 || local_gs == -1 || local_gs == 6) {

			/*
			 * Outer
			 */

			for (ale_pos i = outer_min_t[0]; i <= outer_max_t[0]; i += adj_p)
			for (ale_pos j = outer_min_t[1]; j <  inner_min_t[1]; j += adj_p) {
				transformation test_t = offset;
				test_t.translate(point(i, j));
				test_global(here, si, test_t, local_ax_count, m, local_gs_mo, perturb);
			}
			for (ale_pos i = outer_min_t[0]; i <= outer_max_t[0]; i += adj_p)
			for (ale_pos j = outer_max_t[1]; j >  inner_max_t[1]; j -= adj_p) {
				transformation test_t = offset;
				test_t.translate(point(i, j));
				test_global(here, si, test_t, local_ax_count, m, local_gs_mo, perturb);
			}
			for (ale_pos i = outer_min_t[0]; i <  inner_min_t[0]; i += adj_p)
			for (ale_pos j = outer_min_t[1]; j <= outer_max_t[1]; j += adj_p) {
				transformation test_t = offset;
				test_t.translate(point(i, j));
				test_global(here, si, test_t, local_ax_count, m, local_gs_mo, perturb);
			}
			for (ale_pos i = outer_max_t[0]; i >  inner_max_t[0]; i -= adj_p)
			for (ale_pos j = outer_min_t[1]; j <= outer_max_t[1]; j += adj_p) {
				transformation test_t = offset;
				test_t.translate(point(i, j));
				test_global(here, si, test_t, local_ax_count, m, local_gs_mo, perturb);
			}
		}
	}

	static void get_translational_set(std::vector<transformation> *set,
			transformation t, ale_pos adj_p) {

		ale_pos adj_s;

		transformation offset = t;
		transformation test_t;

		for (int i = 0; i < 2; i++)
		for (adj_s = -adj_p; adj_s <= adj_p; adj_s += 2 * adj_p) {

			test_t = offset;

			test_t.translate(i ? point(adj_s, 0) : point(0, adj_s));

			set->push_back(test_t);
		}
	}

	static int threshold_ok(ale_accum error) {
		if ((1 - error) * (ale_accum) 100 >= match_threshold)
			return 1;

		if (!(match_threshold >= 0))
			return 1;

		return 0;
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
	static diff_stat_t _align(int m, int local_gs, element_t *element) {

		const image *input_frame = element->input_frame;

		/*
		 * Local upper/lower data, possibly dependent on image
		 * dimensions.
		 */

		ale_pos local_lower, local_upper;
		ale_accum local_gs_mo;

		/*
		 * Select the minimum dimension as the reference.
		 */

		ale_pos reference_size = input_frame->height();
		if (input_frame->width() < reference_size)
			reference_size = input_frame->width();
		ale_accum reference_area = input_frame->height()
			                 * input_frame->width();

		if (perturb_lower_percent)
			local_lower = (double) perturb_lower
				    * (double) reference_size
				    * (double) 0.01
				    * (double) scale_factor;
		else
			local_lower = perturb_lower;

		if (perturb_upper_percent)
			local_upper = (double) perturb_upper
				    * (double) reference_size
				    * (double) 0.01
				    * (double) scale_factor;
		else
			local_upper = perturb_upper;

		local_upper = pow(2, floor(log(local_upper) / log(2)));

		if (gs_mo_percent)
			local_gs_mo = (double) _gs_mo
				    * (double) reference_area
				    * (double) 0.01
				    * (double) scale_factor;
		else
			local_gs_mo = _gs_mo;

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

		/*
		 * Plain (unsigned int) casting seems to be broken in some cases.
		 */

		unsigned int steps = (perturb > pow(2, lod_max)) 
			           ? (unsigned int) lrint(log(perturb) / log(2)) - lod_max + 1 : 1;

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

			// element->default_initial_alignment = orig_t;
			element->default_initial_alignment.set_current_element(orig_t.get_element(0));
			element->default_initial_alignment.set_dimensions(input_frame);

		} else if (default_initial_alignment_type == 1)

			/*
			 * Follow previous transformation, setting new image
			 * dimensions.
			 */

			element->default_initial_alignment.set_dimensions(input_frame);

		else
			assert(0);

		element->old_is_default = element->is_default;

		/*
		 * Set the default transformation.
		 */

		transformation offset = element->default_initial_alignment;

		/*
		 * Load any file-specified transformations
		 */

		for (offset.set_current_index(0); 
		     offset.get_current_index() < _ma_card; 
		     offset.push_element()) {

			offset = tload_next(tload, alignment_class == 2, 
					offset, 
					&element->is_default, 
					offset.get_current_index() == 0);

		}

		offset.set_current_index(0);

		if (perturb > 0) {

			/*
			 * Apply following logic
			 */

			transformation new_offset = follow(element, offset, lod);

			new_offset.set_current_index(0);
			
			element->old_initial_alignment = offset;
			offset = new_offset;

		} else {
			element->old_initial_alignment = offset;
		}

		/*
		 * Control point alignment
		 */

		if (local_gs == 5) {

			transformation o = offset;

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

			ale_pos lowest_error = cp_rms_error(m, o);

			ale_pos rot_lower = 2 * (double) local_lower
					  / sqrt(pow(scale_clusters[0].input->height(), 2)
					       + pow(scale_clusters[0].input->width(),  2))
					  * 180
					  / M_PI;

			if  (alignment_class > 0)
			for (double rot = 30; rot > rot_lower; rot /= 2) 
			for (double srot = -rot; srot < rot * 1.5; srot += rot * 2) {
				int is_improved = 1;
				while (is_improved) {
					is_improved = 0;
					transformation test_t = o;
					/*
					 * XXX: is this right?
					 */
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
				ale_pos adj_p = lowest_error;

				if (adj_p < local_lower)
					adj_p = local_lower;

				while (adj_p >= local_lower) {
					transformation test_t = o;
					int is_improved = 1;
					ale_pos test_v;
					ale_pos adj_s;

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
		}

		/*
		 * Pre-alignment exposure adjustment
		 */

		if (_exp_register) {
			ui::get()->exposure_1();
			set_exposure_ratio(m, scale_clusters[0], offset, local_ax_count, 0);
		}

		/*
		 * Scale transform for lod
		 */

		for (int lod_ = 0; lod_ < lod; lod_++) {
			transformation s = offset;
			transformation t = offset;

			t.rescale(1 / (double) 2);

			if (!(t.scaled_height() > 0 && t.scaled_height() < s.scaled_height())
			 || !(t.scaled_width() > 0  && t.scaled_width() < s.scaled_width())) {
				perturb /= pow(2, lod - lod_);
				lod = lod_;
				break;
			} else {
				offset = t;
			}
		}

		ui::get()->set_offset(offset);

		struct scale_cluster si = scale_clusters[lod];

		/*
		 * Projective adjustment value
		 */

		ale_pos adj_p = (perturb >= pow(2, lod_diff))
			     ? pow(2, lod_diff) : (double) perturb;

		/*
		 * Orientational adjustment value in degrees.
		 *
		 * Since rotational perturbation is now specified as an
		 * arclength, we have to convert.
		 */

		ale_pos adj_o = (double) 2 * (double) perturb 
				  / sqrt(pow((double) scale_clusters[0].input->height(), (double) 2)
				       + pow((double) scale_clusters[0].input->width(),  (double) 2))
				  * (double) 180
				  / M_PI;

		/*
		 * Barrel distortion adjustment value
		 */

		ale_pos adj_b = perturb * bda_mult;

		/*
		 * Global search overlap requirements.
		 */

		local_gs_mo = (double) local_gs_mo / pow(pow(2, lod), 2);

		/*
		 * Alignment statistics.
		 */

		diff_stat_t here;

		/*
		 * Current difference (error) value
		 */

		ui::get()->prematching();
		here.diff(si, perturb, offset, local_ax_count, m);
		ui::get()->set_match(here.get_error());

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

		if (perturb >= local_lower && local_gs != 0 && local_gs != 5
		 && (local_gs != 6 || element->is_default)) {
			
			ui::get()->global_alignment(perturb, lod);
			ui::get()->gs_mo(local_gs_mo);

			test_globals(&here, si, here.get_offset(), local_gs, adj_p,
					local_ax_count, m, local_gs_mo, perturb);

			ui::get()->set_match(here.get_error());
			ui::get()->set_offset(here.get_offset());
		}

		/*
		 * Announce perturbation size
		 */

		ui::get()->aligning(perturb, lod);

		/*
		 * Run initial tests to get perturbation multipliers and error
		 * centroids.
		 */

		std::vector<transformation> t_set;

		here.get_perturb_set(&t_set, adj_p, adj_o, adj_b, current_bd, modified_bd);

		/*
		 * Perturbation adjustment loop.  
		 */

		int stable_count = 0;

		while (perturb >= local_lower) {

			/*
			 * Orientational adjustment value in degrees.
			 *
			 * Since rotational perturbation is now specified as an
			 * arclength, we have to convert.
			 */

			ale_pos adj_o = 2 * (double) perturb 
				          / sqrt(pow(scale_clusters[0].input->height(), 2)
					       + pow(scale_clusters[0].input->width(),  2))
					  * 180
					  / M_PI;

			/*
			 * Barrel distortion adjustment value
			 */

			ale_pos adj_b = perturb * bda_mult;

			diff_stat_t old_here = here;

			here.perturb_test(perturb, adj_p, adj_o, adj_b, current_bd, modified_bd,
				stable_count);

			if (here.get_offset() == old_here.get_offset())
				stable_count++;
			else
				stable_count = 0;

			if (stable_count == 3) {

				stable_count = 0;

				here.calculate_element_region();

				if (here.get_current_index() + 1 < _ma_card) {
					here.push_element();
					here.make_element_nontrivial(adj_p, adj_o);
					element->is_primary = 0;
				} else {

					here.set_current_index(0);

					element->is_primary = 1;

					perturb *= 0.5;

					if (lod > 0) {

						/*
						 * Work with images twice as large
						 */

						lod--;
						si = scale_clusters[lod];

						/* 
						 * Rescale the transforms.
						 */

						ale_pos rescale_factor = (double) scale_factor
						                       / (double) pow(2, lod)
								       / (double) here.get_offset().scale();

						here.rescale(rescale_factor, si);

					} else {
						adj_p = perturb;
					}

					/*
					 * Announce changes
					 */

					ui::get()->alignment_perturbation_level(perturb, lod);

				}
			}

			ui::get()->set_match(here.get_error());
			ui::get()->set_offset(here.get_offset());
		}

		here.set_current_index(0);

		if (lod > 0) {
			ale_pos rescale_factor = (double) scale_factor
					       / (double) pow(2, lod)
					       / (double) here.get_offset().scale();

			here.rescale(rescale_factor, scale_clusters[0]);
		}

		offset = here.get_offset();

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
		offset.use_full_support();
		here.diff(scale_clusters[0], perturb, offset, local_ax_count, m);
		here.confirm();
		offset.use_restricted_support();
		ui::get()->set_match(here.get_error());

		/*
		 * Free the level-of-detail structures
		 */

		final_clusters(scale_clusters, scale_factor, steps);

		/*
		 * Ensure that the match meets the threshold.
		 */

		if (threshold_ok(here.get_error())) {
			/*
			 * Update alignment variables
			 */
			latest_ok = 1;
			element->default_initial_alignment = offset;
			element->old_final_alignment = offset;
			ui::get()->alignment_match_ok();
		} else if (local_gs == 4) {

			/*
			 * Align with outer starting points.
			 */

			/*
			 * XXX: This probably isn't exactly the right thing to do,
			 * since variables like old_initial_value have been overwritten.
			 */

			diff_stat_t nested_result = _align(m, -1, element);

			if (threshold_ok(nested_result.get_error())) {
				return nested_result;
			} else if (nested_result.get_error() < here.get_error()) {
				here = nested_result;
			}

			if (is_fail_default)
				offset = element->default_initial_alignment;

			ui::get()->set_match(here.get_error());
			ui::get()->alignment_no_match();
		
		} else if (local_gs == -1) {
			
			latest_ok = 0;
			latest_t = offset;
			return here;

		} else {
			if (is_fail_default)
				offset = element->default_initial_alignment;
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

		for (offset.set_current_index(0); 
		     offset.get_current_index() < _ma_card; 
		     offset.push_element()) {

			tsave_next(tsave, offset, alignment_class == 2, 
					  offset.get_current_index() == 0);
		}

		offset.set_current_index(0);

		/*
		 * Update match statistics.
		 */

		match_sum += (1 - here.get_error()) * (ale_accum) 100;
		match_count++;
		latest = m;
		latest_t = offset;

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

		alignment_weights = new_image_ale_real(rows, cols,
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
				alignment_weights->set_pixel(i, j, 
					alignment_weights->get_pixel(i, j) 
				      * weight_map->get_bl(map_weight_position));
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
				alignment_weights->set_pixel(i, j,
					(pixel) alignment_weights->get_pixel(i, j) 
				      * (pixel) wmx_weights->get_pixel(i, j));
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

		alignment_weights = new_image_ale_real(rows, cols,
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
				alignment_weights->set_chan(i / cols, i % cols, k,
					sqrt(pow(inout[i][0] / (rows * cols), 2)
					   + pow(inout[i][1] / (rows * cols), 2)));
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
		assert (n >= 0);

		static std::vector<element_t> elements;

		if (latest < 0) {

			elements.resize(1);

			/*
			 * Handle the initial frame
			 */

			elements[0].input_frame = image_rw::open(n);

			const image *i = elements[0].input_frame;
			int is_default;
			transformation result = alignment_class == 2
				? transformation::gpt_identity(i, scale_factor)
				: transformation::eu_identity(i, scale_factor);
			result = tload_first(tload, alignment_class == 2, result, &is_default);
			tsave_first(tsave, result, alignment_class == 2);

			if (_keep > 0) {
				kept_t = new transformation[image_rw::count()];
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

			elements[0].default_initial_alignment = result;
			orig_t = result;

			image_rw::close(n);
		}

		for (int i = latest + 1; i <= n; i++) {
			int j = 0;

			/*
			 * Handle supplemental frames.
			 */

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

			assert (reference_image != NULL);
			assert (reference_defined != NULL);

			elements[j].input_frame = image_rw::open(i);
			elements[j].is_primary = 1;

			_align(i, _gs, &elements[j]);

			image_rw::close(n);
		}

		if (elements.size() > _ma_card)
			elements.resize(_ma_card);
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
		} else if (!strcmp(type, "defaults")) {
			_gs = 6;
		} else if (!strcmp(type, "points")) {
			_gs = 5;
			keep();
		} else {
			ui::get()->error("bad global search type");
		}
	}

	/*
	 * Multi-alignment contiguity
	 */
	static void ma_cont(double value) {
		_ma_cont = value;
	}

	/*
	 * Multi-alignment subdivision limit
	 */
	static void ma_decomp(unsigned int value) {
		assert(0);
	}

	/*
	 * Multi-alignment cardinality
	 */
	static void ma_card(unsigned int value) {
		assert (value >= 1);
		_ma_card = value;
	}

	/*
	 * Set the minimum overlap for global searching
	 */
	static void gs_mo(ale_pos value, int _gs_mo_percent) {
		_gs_mo = value;
		gs_mo_percent = _gs_mo_percent;
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
		return match_sum / (ale_accum) match_count;
	}
};

#endif
