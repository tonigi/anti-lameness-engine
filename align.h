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
	 * Alignment match sum
	 */

	static double match_sum;

	/*
	 * Alignment match count.
	 */

	static int match_count;

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
	 */
	static double diff(image *a, image *b, transformation t) {
		assert (reference_image);
		float result = 0;
		float divisor = 0;

		float apm_a = (extend || default_initial_alignment_type == 1)
			    ? (float) 1 : (float) a->avg_pixel_magnitude();
		float apm_b = (extend || default_initial_alignment_type == 1)
			    ? (float) 1 : (float) b->avg_pixel_magnitude();

		int i, j, k;

		for (i = 0; (unsigned int) i < a->height(); i++)
		for (j = 0; (unsigned int) j < a->width();  j++) {

			/*
			 * Transform
			 */

			struct point q;
			struct point offset = reference_image->offset();

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
			  || (i >= -offset[0]
			   && j >= -offset[1] 
			   && i <  -offset[0] + extend_orig_height
			   && j <  -offset[1] + extend_orig_width)
			  || reference_defined->get_pixel_component(i, j, 0) != 0)) { 

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
		 * Determine how many whole-pixel adjustment steps will occur
		 */

		int steps = (perturb >= 1) ? (int) (log(perturb) / log(2)) + 1 : 1;
		int step;
		image **accum_scales = (image **) malloc(steps * sizeof(image *));
		image **input_scales = (image **) malloc(steps * sizeof(image *));

		assert (accum_scales);
		assert (input_scales);

		if (!accum_scales || !input_scales) {
			fprintf(stderr, "Couldn't allocate memory for alignment.\n");
			exit(1);
		}

		/*
		 * Level-of-detail difference.  Use a level of detail 2^lod_diff finer
		 * than the adjustment resolution.  Empirically, it's okay to use a
		 * level-of-detail equal to twice the resolution of the perturbation,
		 * so we set lod_diff to 1, as 2^1==2.
		 */

		const int lod_diff = 1;

		/*
		 * Prepare multiple levels of detail.
		 */

		accum_scales[0] = reference_image->clone();
		input_scales[0] = input_frame->clone();
		
		for (step = 1; step < steps; step++) {
			accum_scales[step] = accum_scales[step - 1]->clone();
			accum_scales[step]->scale(0.5);
			input_scales[step] = input_scales[step - 1]->clone();
			input_scales[step]->scale(0.5);
		}

		/*
		 * Initialize variables used in the main loop.
		 */

		lod = (steps - 1) - lod_diff;
		lod = (lod < 0) ? 0 : lod;

		/*
		 * Initialize the default initial transform
		 */

		if (default_initial_alignment_type == 0) 
			
			/*
			 * Identity transformation
			 */

			default_initial_alignment = (alignment_class == 2)
					  ?  transformation::gpt_identity(input_frame)
					  :  transformation:: eu_identity(input_frame);

		else if (default_initial_alignment_type == 1)

			/*
			 * Follow previous transformation
			 */

			default_initial_alignment.rescale(input_frame);

		else
			assert(0);

		offset = tload_next(tload, lod, alignment_class == 2, default_initial_alignment);

		here = diff(accum_scales[lod], input_scales[lod], offset);

		/*
		 * Simulated annealing perturbation adjustment loop.  
		 */

		while (perturb >= perturb_lower) {

			image *ai = accum_scales[lod];
			image *ii = input_scales[lod];

			double adj_p = (perturb >= pow(2, lod_diff)) 
				     ? pow(2, lod_diff) : perturb;
			double adj_s;

			transformation test_t;
			double test_d;
			double old_here = here;

			int i, j;

			if (alignment_class < 2 && alignment_class >= 0) {

				/* 
				 * Translational or euclidean transformation
				 */

				for (i = 0; i < 2 + alignment_class; i++)
				for (adj_s = -adj_p; adj_s <= adj_p; adj_s += 2 * adj_p) {

					test_t = offset;

					test_t.eu_modify(i, adj_s);

					test_d = diff(ai, ii, test_t);

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

				for (i = 0; i < 4; i++)
				for (j = 0; j < 2; j++)
				for (adj_s = -adj_p; adj_s <= adj_p; adj_s += 2 * adj_p) {

					test_t = offset;

					test_t.gpt_modify(j, i, adj_s);

					test_d = diff(ai, ii, test_t);

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

					if (perturb >= perturb_lower)
						here = diff(accum_scales[lod], 
							input_scales[lod], offset);
				}

				/*
				 * Announce that we've dropped a perturbation level.
				 */

				fprintf(stderr, ".");
			}
		}

		/*
		 * Free the level-of-detail structures
		 */
		for (step = 0; step < steps; step++) {
			delete accum_scales[step];
			delete input_scales[step];
		}
		free(accum_scales);
		free(input_scales);

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

				reference->operator()(i - 1);
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
							:: gpt_identity(reference_image)
						  : transformation
							::  eu_identity(reference_image);

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
				transformation::eu_identity(i);
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
	 * Get match summary statistics.
	 */
	static double match_summary() {
		return match_sum / match_count;
	}
};

#endif
