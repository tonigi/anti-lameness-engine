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

/*
 * align.h: Handle alignment of frames.
 */

#ifndef __d2align_h__
#define __d2align_h__

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

	static ale_pos scale_factor;

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
	static const image *reference_image;
	static const image *reference_defined;

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
	static ale_accum diff(const image *a, const image *b, const image *definition, transformation t, ale_pos _mc_arg) {
		assert (reference_image);
		ale_accum result = 0;
		ale_accum divisor = 0;

		int i, j, k;

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

			/*
			 * Transform
			 */

			struct point q;
			struct point offset = a->offset();

			q = t.inverse_transform(
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
			 && (extend == 0
			  || definition->get_pixel(i, j)[0] != 0)) { 

				pixel pa = a->get_pixel(i, j);
				pixel pb = b->get_bl(point(ti, tj));

				if (channel_alignment_type == 0) {
					/*
					 * Align based on all channels.
					 */


					for (k = 0; k < 3; k++) {
						ale_real achan = pa[k];
						ale_real bchan = pb[k];

						result += pow(fabs(achan - bchan), metric_exponent);
						divisor += pow(achan > bchan ? achan : bchan, metric_exponent);
					}
				} else if (channel_alignment_type == 1) {
					/*
					 * Align based on the green channel.
					 */

					ale_real achan = pa[1];
					ale_real bchan = pb[1];

					result += pow(fabs(achan - bchan), metric_exponent);
					divisor += pow(achan > bchan ? achan : bchan, metric_exponent);
				} else if (channel_alignment_type == 2) {
					/*
					 * Align based on the sum of all channels.
					 */

					ale_real asum = 0;
					ale_real bsum = 0;

					for (k = 0; k < 3; k++) {
						asum += pa[k];
						bsum += pb[k];
					}

					result += pow(fabs(asum - bsum), metric_exponent);
					divisor += pow(asum > bsum ? asum : bsum, metric_exponent);
				}
			}
		}
		return pow(result / divisor, 1/metric_exponent);
	}

	/*
	 * Adjust exposure for an aligned frame B against reference A.
	 */
	static void set_exposure_ratio(unsigned int m, const image *a, const image *b, 
			const image *definition, transformation t) {

		pixel_accum asum(0, 0, 0), bsum(0, 0, 0);

		for (unsigned int i = 0; i < a->height(); i++)
		for (unsigned int j = 0; j < a->width(); j++) {

			/*
			 * Transform
			 */

			struct point q;
			struct point offset = a->offset();

			q = t.inverse_transform(
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
			 && (extend == 0
			  || definition->get_pixel(i, j)[0] != 0)) { 
				asum += a->get_pixel(i, j);
				bsum += b->get_bl(q);
			}
		}

		// std::cerr << (asum / bsum) << " ";

		image_rw::exp(m).set_multiplier((asum / bsum)
				* image_rw::exp(m).get_multiplier());
	}


	/*
	 * Update an accumulated image with a new input image.
	 */
	static void update(int m) {
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
		const image **accum_scales = (const image **) malloc(steps * sizeof(image *));
		const image **input_scales = (const image **) malloc(steps * sizeof(image *));
		const image **defined_scales = (const image **) malloc(steps * sizeof(image *));

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

		if (scale_factor != 1.0)
			input_scales[0] = input_frame->scale(scale_factor, "input_scales[0]");
		else
			input_scales[0] = input_frame;

		for (step = 1; step < steps; step++) {
			accum_scales[step] = accum_scales[step - 1]->scale_by_half("accum_scales[step]");
			input_scales[step] = input_scales[step - 1]->scale_by_half("accum_scales[step]");
			if (extend) {
				defined_scales[step] = defined_scales[step - 1]->defined_scale_by_half("defined_scales[step]");
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

		ale_pos _mc_arg = _mc * pow(2, 2 * lod);
		const image *ai = accum_scales[lod];
		const image *ii = input_scales[lod];
		const image *di = defined_scales[lod];

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
					defined_scales[0], offset);

		/*
		 * Current difference (error) value
		 */

		here = diff(ai, ii, di, offset, _mc_arg);

		/*
		 * Perturbation adjustment loop.  
		 */

		while (perturb >= perturb_lower) {

			ale_pos adj_s;

			/*
			 * Orientational adjustment value
			 */

			ale_pos adj_o = perturb;

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
		}

		/*
		 * Post-alignment exposure adjustment
		 */

		if (_exp_register)
			set_exposure_ratio(m, accum_scales[0], input_scales[0],
					defined_scales[0], offset);

		/*
		 * Recalculate error
		 */

		here = diff(accum_scales[0], input_scales[0], defined_scales[0], offset, _mc);

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
			/*
			 * Update alignment variables
			 */
			latest_ok = 1;
			default_initial_alignment = offset;
			old_final_alignment = offset;
			fprintf(stderr, " okay (%f%% match)", (double) (1 - here) * 100);
		} else {
			latest_ok = 0;
			fprintf(stderr, " no match (%f%% match)", (double) (1 - here) * 100);
		}

		match_sum += (1 - here) * 100;
		match_count++;

		image_rw::close(m);
		latest = m;
	}

	/*
	 * Update alignment to frame N.
	 */
	static void update_to(int n) {
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
				update(i);
			}
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
		update_to(n);

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
		update_to(n);
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
	 * Get match summary statistics.
	 */
	static ale_accum match_summary() {
		return match_sum / match_count;
	}
};

#endif
