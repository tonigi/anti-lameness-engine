// Copyright 2002 David Hilvert <dhilvert@ugcs.caltech.edu>

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

#ifndef __align_h__
#define __align_h__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "gpt.h"
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

	static double scale_factor;

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
	static const image  *reference_image;
	static const image_weights *reference_defined;

	/*
	 * New image to align
	 */

	static const image *input_frame;

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
	 * We may need to know whether extension is occurring.
	 */
	static int extend;
	static int extend_orig_height;
	static int extend_orig_width;

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

	static double metric_exponent;

	/*
	 * Match threshold
	 */

	static double match_threshold;

	/*
	 * Perturbation lower and upper bounds.
	 */

	static double perturb_lower;
	static double perturb_upper;

	/*
	 * Maximum level-of-detail scale factor is 2^lod_max/perturb.
	 */

	static int lod_max;

	/*
	 * Maximum rotational perturbation
	 */

	static double rot_max;

	/*
	 * Alignment match sum
	 */

	static double match_sum;

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

	static double _mc;

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
	static double diff(const image *a, const image *b, const image_weights *d, transformation t, double _mc_arg) {
		assert (reference_image);
		float result = 0;
		float divisor = 0;

		float apm_a = (extend || default_initial_alignment_type == 1)
			    ? (float) 1 : (float) a->avg_pixel_magnitude();
		float apm_b = (extend || default_initial_alignment_type == 1)
			    ? (float) 1 : (float) b->avg_pixel_magnitude();

		int i, j, k;

		/*
		 * We cut-and-paste to handle the difference between
		 * monte-carlo and non-monte-carlo cases.  This is amenable to
		 * abstraction by function definition.  An investigation into
		 * execution speed for the various possibilities would be
		 * appropriate.
		 */

		if (_mc_arg <= 0 || _mc_arg >= 1) {

			/*
			 * No monte-carlo pixel sampling; sample all pixels
			 */

			for (i = 0; (unsigned int) i < a->height(); i++)
			for (j = 0; (unsigned int) j < a->width();  j++) {

				/*
				 * Transform
				 */

				struct point q;
				struct point offset = a->offset();

				q = t.inverse_transform(
					point(i + offset[0], j + offset[1]));

				my_real ti = q[0];
				my_real tj = q[1];

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
				 && (extend == 0
				  // || (i >= -offset[0]
				  //  && j >= -offset[1] 
				  //  && i <  -offset[0] + extend_orig_height
				  //  && j <  -offset[1] + extend_orig_width)
				  || d->get_pixel_component(i, j, 0) != 0)) { 

					if (channel_alignment_type == 0) {
						/*
						 * Align based on all channels.
						 */

						for (k = 0; k < 3; k++) {
							float achan = a->get_pixel_component(i, j, k) / apm_a;
							float bchan = b->get_bl_component(ti, tj, k) / apm_b;

							result += pow(fabs(achan - bchan), metric_exponent);
							divisor += pow(achan > bchan ? achan : bchan, metric_exponent);
						}
					} else if (channel_alignment_type == 1) {
						/*
						 * Align based on the green channel.  XXX: apm_* shouldn't be used here.
						 */

						float achan = a->get_pixel_component(i, j, 1) / apm_a;
						float bchan = b->get_bl_component(ti, tj, 1) / apm_b;

						result += pow(fabs(achan - bchan), metric_exponent);
						divisor += pow(achan > bchan ? achan : bchan, metric_exponent);
					} else if (channel_alignment_type == 2) {
						/*
						 * Align based on the sum of all channels.
						 */

						float asum = 0;
						float bsum = 0;

						for (k = 0; k < 3; k++) {
							asum += a->get_pixel_component(i, j, k) / apm_a;
							bsum += b->get_bl_component(ti, tj, k) / apm_b;
						}

						result += pow(fabs(asum - bsum), metric_exponent);
						divisor += pow(asum > bsum ? asum : bsum, metric_exponent);
					}
				}
			}
		} else {

			/*
			 * Monte Carlo pixel sampling
			 */

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

			double u = (1 - _mc_arg) / _mc_arg;

			double mc_max = (floor(2*u) * (1 + floor(2*u)) + 2*u)
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
						  * ( (1 + ((double) rand()) ) 
						    / (1 + ((double) RAND_MAX)) ));
			    index < index_max;
			    index += (int) ceil((mc_max+1) 
				              * ( (1 + ((double) rand()) ) 
					        / (1 + ((double) RAND_MAX)) ))){

				i = index / a->width();
				j = index % a->width();

				/*
				 * Transform
				 */

				struct point q;
				struct point offset = a->offset();

				q = t.inverse_transform(
					point(i + offset[0], j + offset[1]));

				my_real ti = q[0];
				my_real tj = q[1];

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
				 && (extend == 0
				  // || (i >= -offset[0]
				  //  && j >= -offset[1] 
				  //  && i <  -offset[0] + extend_orig_height
				  //  && j <  -offset[1] + extend_orig_width)
				  || d->get_pixel_component(i, j, 0) != 0)) { 

					if (channel_alignment_type == 0) {
						/*
						 * Align based on all channels.
						 */

						for (k = 0; k < 3; k++) {
							float achan = a->get_pixel_component(i, j, k) / apm_a;
							float bchan = b->get_bl_component(ti, tj, k) / apm_b;

							result += pow(fabs(achan - bchan), metric_exponent);
							divisor += pow(achan > bchan ? achan : bchan, metric_exponent);
						}
					} else if (channel_alignment_type == 1) {
						/*
						 * Align based on the green channel.  XXX: apm_* shouldn't be used here.
						 */

						float achan = a->get_pixel_component(i, j, 1) / apm_a;
						float bchan = b->get_bl_component(ti, tj, 1) / apm_b;

						result += pow(fabs(achan - bchan), metric_exponent);
						divisor += pow(achan > bchan ? achan : bchan, metric_exponent);
					} else if (channel_alignment_type == 2) {
						/*
						 * Align based on the sum of all channels.
						 */

						float asum = 0;
						float bsum = 0;

						for (k = 0; k < 3; k++) {
							asum += a->get_pixel_component(i, j, k) / apm_a;
							bsum += b->get_bl_component(ti, tj, k) / apm_b;
						}

						result += pow(fabs(asum - bsum), metric_exponent);
						divisor += pow(asum > bsum ? asum : bsum, metric_exponent);
					}
				}
			}
		}

		return pow(result / divisor, 1/metric_exponent);
	}

	/*
	 * Update an accumulated image with a new input image.
	 */
	static void update() {
		double perturb = pow(2, floor(log(perturb_upper) / log(2)));
		transformation offset;
		double here;
		int lod;

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
		const image **accum_scales = (const image **) malloc(steps * sizeof(image *));
		const image **input_scales = (const image **) malloc(steps * sizeof(image *));
		const image_weights **defined_scales = (const image_weights **) malloc(steps * sizeof(image_weights *));

		assert (accum_scales);
		assert (input_scales);
		assert (defined_scales);

		if (!accum_scales || !input_scales || !defined_scales) {
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

		if (extend)
			defined_scales[0] = reference_defined;

		if (scale_factor != 1.0) {
			
			/*
			 * In cursory tests, scaling is much more expensive
			 * than cloning, so we don't go out of our way to avoid
			 * this cloning step, although it could be avoided with
			 * certain (possibly ugly) source code modifications.
			 */

			image *input_frame_clone = input_frame->clone();

			input_frame_clone->scale(scale_factor);

			input_scales[0] = input_frame_clone;

		} else {

			input_scales[0] = input_frame;

		}

		for (step = 1; step < steps; step++) {
			accum_scales[step] = accum_scales[step - 1]->scale_by_half();
			input_scales[step] = input_scales[step - 1]->scale_by_half();
			if (extend) {
				defined_scales[step] = defined_scales[step - 1]->defined_scale_by_half();
			}
		}

		/*
		 * Initialize variables used in the main loop.
		 */

		lod = (steps - 1);

		/*
		 * Initialize the default initial transform
		 */

		if (default_initial_alignment_type == 0) 
			
			/*
			 * Identity transformation
			 */

			default_initial_alignment = (alignment_class == 2)
					  ?  transformation::gpt_identity(input_frame, scale_factor)
					  :  transformation:: eu_identity(input_frame, scale_factor);

		else if (default_initial_alignment_type == 1)

			/*
			 * Follow previous transformation
			 *
			 * NB: (input_frame, scale_factor) is used as an
			 * idiom to mean "the current input scaled by the given
			 * factor".  Hence, we are 'rescaling' according to the
			 * new input frame dimensions; the scale factor merely
			 * adjusts these dimensions.
			 */

			default_initial_alignment.rescale(input_frame, scale_factor);

		else
			assert(0);

		old_is_default = is_default;
		offset = tload_next(tload, lod, alignment_class == 2, default_initial_alignment, &is_default);

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
	
				double t0 = -old_initial_alignment.eu_get(0) + old_final_alignment.eu_get(0);
				double t1 = -old_initial_alignment.eu_get(1) + old_final_alignment.eu_get(1);

				new_offset.eu_modify(0, t0);
				new_offset.eu_modify(1, t1);

			} else if (alignment_class == 1) {
				/*
				 * Euclidean transformations
				 */

				double t2 = -old_initial_alignment.eu_get(2) + old_final_alignment.eu_get(2);

				new_offset.eu_modify(2, t2);

				point p( offset.height()/2 + offset.eu_get(0) - old_initial_alignment.eu_get(0),
					 offset.width()/2 + offset.eu_get(1) - old_initial_alignment.eu_get(1) );

				p = old_final_alignment(p);

				new_offset.eu_modify(0, p[0] - offset.height()/2 - offset.eu_get(0));
				new_offset.eu_modify(1, p[1] - offset.width()/2 - offset.eu_get(1));
				
			} else if (alignment_class == 2) {
				/*
				 * Projective transformations
				 */

				point p[4];

				p[0] = old_final_alignment(old_initial_alignment
				     . inverse_transform(offset(point(      0        ,       0       ))));
				p[1] = old_final_alignment(old_initial_alignment
				     . inverse_transform(offset(point(offset.height(),       0       ))));
				p[2] = old_final_alignment(old_initial_alignment
				     . inverse_transform(offset(point(offset.height(), offset.width()))));
				p[3] = old_final_alignment(old_initial_alignment
				     . inverse_transform(offset(point(      0        , offset.width()))));

				new_offset.gpt_set(p);
			}
		}

		old_initial_alignment = offset;
		offset = new_offset;

		double _mc_arg = _mc * pow(2, 2 * lod);
		const image *ai = accum_scales[lod];
		const image *ii = input_scales[lod];
		const image_weights *di = defined_scales[lod];

		/*
		 * Positional adjustment value
		 */

		double adj_p = (perturb >= pow(2, lod_diff))
			     ? pow(2, lod_diff) : perturb;

		/*
		 * Current difference (error) value
		 */

		here = diff(ai, ii, di, offset, _mc_arg);

		/*
		 * Perturbation adjustment loop.  
		 */

		while (perturb >= perturb_lower) {

			double adj_s;

			/*
			 * Orientational adjustment value
			 */

			double adj_o = perturb;

			transformation test_t;
			double test_d;
			double old_here = here;

			if (alignment_class < 2 && alignment_class >= 0) {

				/* 
				 * Translational or euclidean transformation
				 */

				for (int i = 0; i < 2; i++)
				for (adj_s = -adj_p; adj_s <= adj_p; adj_s += 2 * adj_p) {

					test_t = offset;

					test_t.eu_modify(i, adj_s);

					test_d = diff(ai, ii, di, test_t, _mc_arg);

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

					test_d = diff(ai, ii, di, test_t, _mc_arg);

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

					test_d = diff(ai, ii, di, test_t, _mc_arg);

					if (test_d < here) {
						here = test_d;
						offset = test_t;
						adj_s += 3 * adj_p;
					}
				}

			} else assert(0);

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
					_mc_arg /= 4;

					here = diff(ai, ii, di, offset, _mc_arg);

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
			here = diff(accum_scales[0], input_scales[0], defined_scales[0], offset, _mc);
		}

		/*
		 * Free the level-of-detail structures
		 */

		if (scale_factor != 1.0)
			delete input_scales[0];

		for (step = 1; step < steps; step++) {
			delete accum_scales[step];
			delete input_scales[step];
			if (extend)
				delete defined_scales[step];
		}

		free(accum_scales);
		free(input_scales);
		free(defined_scales);

		/*
		 * Save the transformation information
		 */

		tsave_next(tsave, offset, alignment_class == 2);

		latest_t = offset;

		/*
		 * Ensure that the match meets the threshold.
		 */

		if ((1 - here) * 100 > match_threshold) {
			latest_ok = 1;
			default_initial_alignment = offset;
			old_final_alignment = offset;
			fprintf(stderr, " okay (%f%% match)", (1 - here) * 100);
		} else {
			latest_ok = 0;
			fprintf(stderr, " no match (%f%% match)", (1 - here) * 100);
		}

		match_sum += (1 - here) * 100;
		match_count++;
	}

	/*
	 * Update alignment to frame N.
	 */
	static void update(int n) {
		assert (n <= latest + 1);

		for (int i = latest + 1; i <= n; i++) {
			if (i > 0)  {
				assert (reference != NULL);

				reference->sync(i - 1);
				reference_image = reference->get_image();
				reference_defined = reference->get_defined();

				assert (reference_image != NULL);
				assert (reference_defined != NULL);
			}
			if (i == 1) {
				{
					const image *im = image_rw::open(0);

					extend_orig_height = im->height();
					extend_orig_width  = im->width();

					image_rw::close(0);
				}

				default_initial_alignment = (alignment_class == 2)
						  ? transformation
							:: gpt_identity(reference_image, scale_factor)
						  : transformation
							::  eu_identity(reference_image, scale_factor);

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

				}
			}
			if (i > latest) {
				if (_keep) {
					kept_t[latest] = latest_t;
					kept_ok[latest] = latest_ok;
				}
				input_frame = image_rw::open(i);
				update();
				image_rw::close(i);
				latest = i;
			}
		}
	}

public:

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
	 * Message that we are extending images.
	 */
	static void set_extend() {
		extend = 1;
	}

	/*
	 * Message that we are not extending images.
	 */
	static void no_extend() {
		extend = 0;
	}

	/*
	 * Error metric exponent
	 */

	static void set_metric_exponent(double me) {
		metric_exponent = me;
	}

	/*
	 * Match threshold
	 */

	static void set_match_threshold(double mt) {
		match_threshold = mt;
	}

	/*
	 * Perturbation lower and upper bounds.
	 */

	static void set_perturb_lower(double pl) {
		perturb_lower = pl;
	}

	static void set_perturb_upper(double pu) {
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
	 * Level-of-detail
	 */

	static void set_lod_max(int lm) {
		lod_max = lm;
	}

	/*
	 * Set the scale factor
	 */
	static void set_scale(double s) {
		scale_factor = s;
	}

	/*
	 * Set reference rendering to align against
	 */
	static void set_reference(render *r) {
		reference = r;
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
	 *
	 * Frame 0 always matches.
	 */
	static int match(int n) {
		update(n);

		if (n == 0)
			return 1;
		else if (n == latest)
			return latest_ok;
		else if (_keep)
			return kept_ok[n];
		else
			assert(0);

		return 0;
	}

	/*
	 * Message that old alignment data should be kept.
	 */
	static void keep() {
		assert (latest == 0);
		_keep = 1;
	}

	/*
	 * Get alignment for frame N.
	 *
	 * Frame 0 is always aligned with the identity transformation.
	 */
	static transformation of(int n) {
		update(n);
		if (n == 0) {
			const image *i = image_rw::open(n);
			transformation result = 
				transformation::eu_identity(i, scale_factor);
			image_rw::close(n);
			return result;
		} else if (n == latest)
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
	static void mc(double n) {
		_mc = n;
	}

	/*
	 * Don't use Monte Carlo alignment sampling.
	 */
	static void no_mc() {
		_mc = 0;
	}

	/*
	 * Get match summary statistics.
	 */
	static double match_summary() {
		return match_sum / match_count;
	}
};

#endif
