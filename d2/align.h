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
	 * Reference rendering to align against
	 */

	static render *reference;
	static filter::scaled_filter *interpolant;
	static const image *reference_image;
	static const image *reference_defined;

 	/*
	 * Per-pixel alignment weights
	 */

	static const image *alignment_weights;

	/*
	 * Frequency-dependent alignment weights
	 */

	static double horiz_freq_cut;
	static double vert_freq_cut;
	static double avg_freq_cut;
	static image *frequency_weights;
	static const char *fw_output;

	/*
	 * Algorithmic alignment weighting
	 */

	static const char *wmx_exec;
	static const char *wmx_file;
	static const char *wmx_defs;
	static image *wmx_weights;

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
	 * Helper variables for new --follow semantics (as of 0.5.0).
	 *
	 * These variables implement delta following.  The change between the
	 * non-default old initial alignment and old final alignment is used to
	 * adjust the non-default current initial alignment.  If either the old
	 * or new initial alignment is a default alignment, the old --follow
	 * semantics are preserved.
	 */

	static int is_default, old_is_default;
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
	static ale_pos perturb_upper;

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

	/*
	 * Exclusion regions
	 */

	static int *ax_parameters;
	static int ax_count;

	/*
	 * Not-quite-symmetric difference function.  Determines the difference in areas
	 * where the arrays overlap.  Uses the first array's notion of pixel positions.
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
	static ale_accum diff(const image *a, const image *b, const image *definition, transformation t, ale_pos _mc_arg,
			int ax_count, const int *ax_parameters) {
		assert (reference_image);
		ale_accum result = 0;
		ale_accum divisor = 0;

		point offset = a->offset();

		int i, j, k;

		if (interpolant != NULL) 
			interpolant->set_parameters(t, b, offset);

		/*
		 * We always the same code for exhaustive and Monte Carlo pixel
		 * sampling, setting _mc_arg = 1 when all pixels are to be
		 * sampled.
		 */

		if (_mc_arg <= 0 || _mc_arg >= 1)
			_mc_arg = 1;

		int index;
		int index_max = a->height() * a->width();
		
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
		 */

		srand(1);

		for(index = -1 + (int) ceil((mc_max+1) 
					  * ( (1 + ((ale_pos) rand()) ) 
					    / (1 + ((ale_pos) RAND_MAX)) ));
		    index < index_max;
		    index += (int) ceil((mc_max+1) 
				      * ( (1 + ((ale_pos) rand()) ) 
					/ (1 + ((ale_pos) RAND_MAX)) ))){

			i = index / a->width();
			j = index % a->width();

			int ax_ok = 1;

			for (int idx = 0; idx < ax_count; idx++)
				if (i + offset[0] >= ax_parameters[idx * 4 + 0]
				 && i + offset[0] <= ax_parameters[idx * 4 + 1]
				 && j + offset[1] >= ax_parameters[idx * 4 + 2]
				 && j + offset[1] <= ax_parameters[idx * 4 + 3])
					ax_ok = 0;

			if (ax_ok == 0) {
				// ((image *)a)->set_pixel(i, j, pixel(0, 0, 0));
				continue;
			}

			/*
			 * Transform
			 */

			struct point q;

			q = t.scaled_inverse_transform(
				point(i + offset[0], j + offset[1]));

			ale_pos ti = q[0];
			ale_pos tj = q[1];

			/*
			 * Check that the transformed coordinates are within
			 * the boundaries of array b.  
			 *
			 * Also, check that the weight value in the accumulated array
			 * is nonzero, unless we know it is nonzero by virtue of the
			 * fact that it falls within the region of the original frame
			 * (e.g. when we're not increasing image extents).
			 */

			if (ti >= 0
			 && ti <= b->height() - 1
			 && tj >= 0
			 && tj <= b->width() - 1
			 && definition->get_pixel(i, j)[0] != 0) { 



				pixel pa = a->get_pixel(i, j);
				pixel pb;
				pixel weight;

				if (interpolant != NULL)
					interpolant->filtered(i, j, &pb, &weight);
				else
					pb = b->get_bl(point(ti, tj));

				/*
				 * Calculate per-channel alignment weight
				 */

				pixel aweight = pixel(1, 1, 1);
				if (alignment_weights != NULL) {
					point aweight_offset = a->offset() - alignment_weights->offset();
					point aweight_position = aweight_offset + point(i, j);
					if (aweight_position[0] >= 0
					 && aweight_position[1] >= 0
					 && aweight_position[0] < alignment_weights->height()
					 && aweight_position[1] < alignment_weights->width())
						aweight = alignment_weights->get_pixel((unsigned int) aweight_position[0],
										       (unsigned int) aweight_position[1]);
				}

				/*
				 * Adjust per-channel alignment weight by
				 * frequency weighting.
				 */

				if (frequency_weights != NULL)
					aweight *= frequency_weights->get_pixel(i, j);

				/*
				 * Adjust per-channel alignment weight by
				 * algorithmic weighting.
				 */

				if (wmx_weights != NULL)
					aweight *= wmx_weights->get_pixel(i, j);

				/*
				 * Determine alignment type.
				 */

				if (channel_alignment_type == 0) {
					/*
					 * Align based on all channels.
					 */


					for (k = 0; k < 3; k++) {
						ale_real achan = pa[k];
						ale_real bchan = pb[k];

						result += aweight[k] * pow(fabs(achan - bchan), metric_exponent);
						divisor += aweight[k] * pow(achan > bchan ? achan : bchan, metric_exponent);
					}
				} else if (channel_alignment_type == 1) {
					/*
					 * Align based on the green channel.
					 */

					ale_real achan = pa[1];
					ale_real bchan = pb[1];

					result += aweight[1] * pow(fabs(achan - bchan), metric_exponent);
					divisor += aweight[1] * pow(achan > bchan ? achan : bchan, metric_exponent);
				} else if (channel_alignment_type == 2) {
					/*
					 * Align based on the sum of all channels.
					 */

					ale_real asum = 0;
					ale_real bsum = 0;
					ale_real wsum = 0;

					for (k = 0; k < 3; k++) {
						asum += pa[k];
						bsum += pb[k];
						wsum += aweight[k] / 3;
					}

					result += wsum * pow(fabs(asum - bsum), metric_exponent);
					divisor += wsum * pow(asum > bsum ? asum : bsum, metric_exponent);
				}
			}
		}
		return pow(result / divisor, 1/metric_exponent);
	}

	/*
	 * Adjust exposure for an aligned frame B against reference A.
	 */
	static void set_exposure_ratio(unsigned int m, const image *a, const image *b, 
			const image *definition, transformation t, int ax_count, const int *ax_parameters) {

		pixel_accum asum(0, 0, 0), bsum(0, 0, 0);

		point offset = a->offset();

		for (unsigned int i = 0; i < a->height(); i++)
		for (unsigned int j = 0; j < a->width(); j++) {

			int ax_ok = 1;

			for (int idx = 0; idx < ax_count; idx++)
				if (i + offset[0] >= ax_parameters[idx * 4 + 0]
				 && i + offset[0] <= ax_parameters[idx * 4 + 1]
				 && j + offset[1] >= ax_parameters[idx * 4 + 2]
				 && j + offset[1] <= ax_parameters[idx * 4 + 3])
					ax_ok = 0;

			if (ax_ok == 0) {
				// ((image *)a)->set_pixel(i, j, pixel(0, 0, 0));
				continue;
			}

			/*
			 * Transform
			 */

			struct point q;

			q = t.scaled_inverse_transform(
				point(i + offset[0], j + offset[1]));

			/*
			 * Check that the transformed coordinates are within
			 * the boundaries of array b.  
			 *
			 * Also, check that the weight value in the accumulated array
			 * is nonzero, unless we know it is nonzero by virtue of the
			 * fact that it falls within the region of the original frame
			 * (e.g. when we're not increasing image extents).
			 */

			if (q[0] >= 0
			 && q[0] <= b->height() - 1
			 && q[1] >= 0
			 && q[1] <= b->width() - 1
			 && definition->get_pixel(i, j).minabs_norm() != 0) { 
				asum += a->get_pixel(i, j);
				bsum += b->get_bl(q);
			}
		}

		// std::cerr << (asum / bsum) << " ";

		image_rw::exp(m).set_multiplier((asum / bsum)
				* image_rw::exp(m).get_multiplier());
	}

	static void init_ax_parameter_scales(int frame, int *local_ax_count, const int **ax_parameter_scales) {
		int *ax_parameter_scales_0 = (int *) malloc(ax_count * 4 * sizeof(int));

		assert (ax_parameter_scales_0);

		if (!ax_parameter_scales_0)  {
			fprintf(stderr, "Unable to allocate memory for exclusion regions.\n");
			exit(1);
		}

		*local_ax_count = 0;

		for (int idx = 0; idx < ax_count; idx++) {
			if (ax_parameters[6 * idx + 4] > frame 
			 || ax_parameters[6 * idx + 5] < frame)
				continue;

			ax_parameter_scales_0[4 * *local_ax_count + 0]
				= (int) floor(scale_factor * ax_parameters[6 * idx + 0]);
			ax_parameter_scales_0[4 * *local_ax_count + 1]
				= (int) ceil (scale_factor * ax_parameters[6 * idx + 1]);
			ax_parameter_scales_0[4 * *local_ax_count + 2]
				= (int) floor(scale_factor * ax_parameters[6 * idx + 2]);
			ax_parameter_scales_0[4 * *local_ax_count + 3]
				= (int) ceil (scale_factor * ax_parameters[6 * idx + 3]);

			(*local_ax_count)++;
		}

		// ax_parameter_scales_0 = (int *) realloc(ax_parameter_scales_0, *local_ax_count * 4 * sizeof(int));

		ax_parameter_scales[0] = ax_parameter_scales_0;
	}

	static void halve_ax_parameter_scales(int local_ax_count, const int **ax_parameter_scales) {
		int *ax_parameter_scales_new = (int *) malloc(local_ax_count * 4 * sizeof(int));

		assert (ax_parameter_scales_new);

		if (!ax_parameter_scales_new)  {
			fprintf(stderr, "Unable to allocate memory for exclusion regions.\n");
			exit(1);
		}

		for (int idx = 0; idx < local_ax_count; idx++) {
			ax_parameter_scales_new[idx * 4 + 0] = ax_parameter_scales[0][idx * 4 + 0] / 2;
			ax_parameter_scales_new[idx * 4 + 1] = (ax_parameter_scales[0][idx * 4 + 1] + 1) / 2;
			ax_parameter_scales_new[idx * 4 + 2] = ax_parameter_scales[0][idx * 4 + 2] / 2;
			ax_parameter_scales_new[idx * 4 + 3] = (ax_parameter_scales[0][idx * 4 + 3] + 1) / 2;
		}

		ax_parameter_scales[1] = ax_parameter_scales_new;
	}

	/*
	 * Align frame m against the reference.
	 */
	static void _align(int m) {
		ale_pos perturb = pow(2, floor(log(perturb_upper) / log(2)));
		transformation offset;
		ale_accum here;
		int lod;
		const image *input_frame = image_rw::open(m);

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

		int steps = (perturb > pow(2, lod_max)) ? (int) (log(perturb) / log(2)) - lod_max + 1 : 1;
		int step;
		int local_ax_count;
		const image **accum_scales = (const image **) malloc(steps * sizeof(image *));
		const image **input_scales = (const image **) malloc(steps * sizeof(image *));
		const image **defined_scales = (const image **) malloc(steps * sizeof(image *));
		const int   **ax_parameter_scales = (const int **) malloc(steps * sizeof(int   *));

		assert (accum_scales);
		assert (input_scales);
		assert (defined_scales);
		assert (ax_parameter_scales);

		if (!accum_scales || !input_scales || !defined_scales || !ax_parameter_scales) {
			fprintf(stderr, "Couldn't allocate memory for alignment.\n");
			exit(1);
		}

		/*
		 * Prepare multiple levels of detail.
		 */

		/*
		 * First, we prepare the highest levels of detail.  Then, we
		 * use scale_by_half() and defined_scale_by_half() to create
		 * reduced levels of detail.
		 */

		accum_scales[0] = reference_image;

		defined_scales[0] = reference_defined;

		if (scale_factor > 1.0)
			input_scales[0] = input_frame->scale(scale_factor, "alignment");
		else 
			input_scales[0] = input_frame;

		init_ax_parameter_scales(m, &local_ax_count, ax_parameter_scales);
		
		for (step = 1; step < steps; step++) {
			accum_scales[step] = accum_scales[step - 1]->scale_by_half("accum_scales[step]");
			input_scales[step] = input_scales[step - 1]->scale_by_half("accum_scales[step]");
			defined_scales[step] = defined_scales[step - 1]->defined_scale_by_half("defined_scales[step]");
			halve_ax_parameter_scales(local_ax_count, ax_parameter_scales + step - 1);
		}

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
			 * XXX: we assume that the lod for the old initial
			 * alignment is equal to the lod for the new initial
			 * alignment, both being equal to the variable 'lod'.
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

			old_final_alignment.rescale (1 / pow(2, lod));

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
		offset = new_offset;

		ale_pos _mc_arg = _mc * pow(2, 2 * lod);
		const image *ai = accum_scales[lod];
		const image *ii = input_scales[lod];
		const image *di = defined_scales[lod];
		const int   *pi = ax_parameter_scales[lod];

		/*
		 * Positional adjustment value
		 */

		ale_pos adj_p = (perturb >= pow(2, lod_diff))
			     ? pow(2, lod_diff) : (double) perturb;

		/*
		 * Pre-alignment exposure adjustment
		 */

		if (_exp_register)
			set_exposure_ratio(m, accum_scales[0], input_scales[0],
					defined_scales[0], offset, local_ax_count, ax_parameter_scales[0]);

		/*
		 * Current difference (error) value
		 */

		here = diff(ai, ii, di, offset, _mc_arg, local_ax_count, pi);

		/*
		 * Current and modified barrel distortion parameters
		 */

		ale_pos current_bd[BARREL_DEGREE];
		ale_pos modified_bd[BARREL_DEGREE];
		offset.bd_get(current_bd);
		offset.bd_get(modified_bd);

		/*
		 * Perturbation adjustment loop.  
		 */

		while (perturb >= perturb_lower) {

			ale_pos adj_s;

			/*
			 * Orientational adjustment value
			 */

			ale_pos adj_o = perturb;

			/*
			 * Barrel distortion adjustment value
			 */

			ale_pos adj_b = perturb * bda_mult;

			transformation test_t;
			ale_accum test_d;
			ale_accum old_here = here;

			if (alignment_class < 2 && alignment_class >= 0) {

				/* 
				 * Translational or euclidean transformation
				 */

				for (int i = 0; i < 2; i++)
				for (adj_s = -adj_p; adj_s <= adj_p; adj_s += 2 * adj_p) {

					test_t = offset;

					test_t.eu_modify(i, adj_s);

					test_d = diff(ai, ii, di, test_t, _mc_arg, local_ax_count, pi);

					if (test_d < here) {
						here = test_d;
						offset = test_t;
						goto done;
					}
				}

				if (alignment_class == 1 && adj_o < rot_max)
				for (adj_s = -adj_o; adj_s <= adj_o; adj_s += 2 * adj_o) {

					test_t = offset;

					test_t.eu_modify(2, adj_s);

					test_d = diff(ai, ii, di, test_t, _mc_arg, local_ax_count, pi);

					if (test_d < here) {
						here = test_d;
						offset = test_t;
						goto done;
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

					test_t.gpt_modify(j, i, adj_s);

					test_d = diff(ai, ii, di, test_t, _mc_arg, local_ax_count, pi);

					if (test_d < here) {
						here = test_d;
						offset = test_t;
						adj_s += 3 * adj_p;
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

					test_d = diff(ai, ii, di, test_t, _mc_arg, local_ax_count, pi);

					if (test_d < here) {
						here = test_d;
						offset = test_t;
						modified_bd[rd] += adj_s;
						goto done;
					}
				}
			}
			
	done:
			if (!(here < old_here)) {
				perturb *= 0.5;

				if (lod > 0) {

					/* 
					 * We're about to work with images
					 * twice as large, so rescale the 
					 * transform.
					 */

					offset.rescale(2);

					/*
					 * Work with images twice as large
					 */

					lod--;
					ai = accum_scales[lod];
					ii = input_scales[lod];
					di = defined_scales[lod];
					pi = ax_parameter_scales[lod];
					_mc_arg /= 4;

					here = diff(ai, ii, di, offset, _mc_arg, local_ax_count, pi);

				} else {
					adj_p = perturb;
				}

				/*
				 * Announce that we've dropped a perturbation level.
				 */

				fprintf(stderr, ".");
			}
		}

		if (lod > 0) {
			offset.rescale(pow(2, lod));
		}

		/*
		 * Post-alignment exposure adjustment
		 */

		if (_exp_register)
			set_exposure_ratio(m, accum_scales[0], input_scales[0],
					defined_scales[0], offset, local_ax_count, ax_parameter_scales[0]);

		/*
		 * Recalculate error
		 */

		here = diff(accum_scales[0], input_scales[0],
				defined_scales[0], offset, _mc, local_ax_count,
				ax_parameter_scales[0]);

		/*
		 * Free the level-of-detail structures
		 */

		if (scale_factor > 1.0)
			delete input_scales[0];

		free((void *)ax_parameter_scales[0]);

		for (step = 1; step < steps; step++) {
			delete accum_scales[step];
			delete input_scales[step];
			delete defined_scales[step];
			free((void *)ax_parameter_scales[step]);
		}

		free(accum_scales);
		free(input_scales);
		free(defined_scales);
		free(ax_parameter_scales);

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
			fprintf(stderr, " okay (%f%% match)", (double) (1 - here) * 100);
		} else {
			if (is_fail_default)
				offset = default_initial_alignment;
			latest_ok = 0;
			fprintf(stderr, " no match (%f%% match)", (double) (1 - here) * 100);
		}

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

		image_rw::close(m);
		latest = m;
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
	 * Update algorithmic weights
	 */
	static void wmx_update() {
#ifdef USE_UNIX

		static exposure *exp_def = new exposure_default();
		static exposure *exp_bool = new exposure_boolean();

		if (wmx_weights != NULL) {
			delete wmx_weights;
			wmx_weights = NULL;
		}

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

			fprintf(stderr, "\n\n*** An error occurred while running `%s %s`. ***\n\n\n", wmx_exec, wmx_file);
			exit(1);
		}

		wait(&exit_status);

		if (exit_status) {
			fprintf(stderr, "\n\n*** Could not fork in d2::align.  ***\n\n\n");
			exit(1);
		}

		wmx_weights = image_rw::read_image(wmx_file, exp_def);

		if (wmx_weights->height() != rows || wmx_weights->width() != cols) {
			fprintf(stderr, "\n\n*** Error: algorithmic weighting must not change image size. ***\n\n\b");
			exit(1);
		}
#endif
	}

	/*
	 * Update frequency weights
	 */
	static void fw_update() {
#ifdef USE_FFTW
		if (frequency_weights != NULL) {
			delete frequency_weights;
			frequency_weights = NULL;
		}

		if (horiz_freq_cut == 0
		 && vert_freq_cut  == 0
		 && avg_freq_cut   == 0)
			return;

		int rows = reference_image->height();
		int cols = reference_image->width();

		frequency_weights = new image_ale_real(rows, cols,
				reference_image->depth(),
				"frequency weights");

		assert (frequency_weights);

		fftw_complex *inout = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * rows * cols);

		assert (inout);

		fftw_plan pf = fftw_plan_dft_2d(rows, cols,
						inout, inout, 
						FFTW_FORWARD, FFTW_ESTIMATE);

		fftw_plan pb = fftw_plan_dft_2d(rows, cols,
						inout, inout, 
						FFTW_BACKWARD, FFTW_ESTIMATE);

		for (int k = 0; k < 3; k++) {
			for (int i = 0; i < rows * cols; i++) {
				inout[i][0] = reference_image->get_pixel(i / cols, i % cols)[k];
				inout[i][1] = 0;
			}

			fftw_execute(pf);
			hipass(rows, cols, inout);
			fftw_execute(pb);

			for (int i = 0; i < rows * cols; i++) {
#if 0
				frequency_weights->pix(i / cols, i % cols)[k] = fabs(inout[i][0] / (rows * cols));
#else
				frequency_weights->pix(i / cols, i % cols)[k] = 
					sqrt(pow(inout[i][0] / (rows * cols), 2)
					   + pow(inout[i][1] / (rows * cols), 2));
#endif
			}
		}

		fftw_destroy_plan(pf);
		fftw_destroy_plan(pb);
		fftw_free(inout);

		if (fw_output != NULL)
			image_rw::write_image(fw_output, frequency_weights);
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

				if (!kept_t || !kept_ok) {
					fprintf(stderr, 
					   "Couldn't allocate memory for alignment\n");
					exit(1);
				}

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

			reference->sync(i - 1);
			reference_image = reference->get_image();
			reference_defined = reference->get_defined();

			fw_update();
			wmx_update();

			assert (reference_image != NULL);
			assert (reference_defined != NULL);

			_align(i);
		}
	}

public:

	/*
	 * Register exposure
	 */
	static void exp_register() {
		_exp_register = 1;
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

	static void set_perturb_lower(ale_pos pl) {
		perturb_lower = pl;
	}

	static void set_perturb_upper(ale_pos pu) {
		perturb_upper = pu;
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
	static void set_alignment_weights(const image *i) {
		alignment_weights = i;
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
	 * Don't use Monte Carlo alignment sampling.
	 */
	static void no_mc() {
		_mc = 0;
	}

	/*
	 * Set alignment exclusion regions
	 */
	static void set_exclusion(int *_ax_parameters, int _ax_count) {
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
