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

#include "image.h"
#include "hf-filter.h"
#include "image_rw.h"
#include "gpt.h"
#include "my_real.h"
#include "tfile.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

/*
 * Version Information
 */

char *version = "0.4.1"
#ifdef USE_MAGICK
		" (File handler: ImageMagick)";
#else
		" (File handler: PPM binary)";
#endif


/* 
 * Global image variables 
 *
 * XXX: these are confusing and should go away.
 *
 */

/* 
 * Accumulated merged image and associated weight array,
 */

image *accum_image;
image_weights *accum_weight;

/*
 * Drizzled image and associated weight array
 */

image *drizzle_image;
image_weights *drizzle_weight;

/*
 * Current supplemental image to be accumulated.
 */

image *input_image;

/*
 * Output image.
 */

image *output_image;

/*
 * Global flags
 */

/*
 * Transformation code.
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

static int transform_code = 1;

/*
 * Default offset type.
 *
 * 0. Start perturbation with an identity transformation by default.
 *
 * 1. Start perturbation with the most recently accepted frame's 
 * final transformation by default.
 */

static int default_transform_type = 0;
static transformation default_transform;

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

static int align_code = 2;

static double scale_factor = 1;
static double metric = 2;
static double match_threshold = 0;
static double perturb_lower = 0.125;
static double perturb_upper = 32;
static double hf_enhance = 0.0;
static double match_sum = 0;
static double drizzle_radius = -1;
static int drizzle_only = 0;
static int match_count = 0;
static int extend = 0;
static int extend_offset_i = 0;
static int extend_offset_j = 0;
static int extend_orig_height = 0;
static int extend_orig_width = 0;
static struct tload_t *tload = NULL;
static struct tsave_t *tsave = NULL;

/*
 * Output current image to a specified file.
 */
inline void output(char *filename) {
	if (drizzle_radius >= 0 && !drizzle_only) {

		/*
		 * Combine drizzled and normally rendered data if necessary.
		 */

		unsigned int i, j, k;

		assert (accum_image->width()  == drizzle_image->width());
		assert (accum_image->height() == drizzle_image->height());
		assert (accum_image->width()  == output_image->width());
		assert (accum_image->height() == output_image->height());

		for (i = 0; i < accum_image->height(); i++)
		for (j = 0; j < accum_image->width();  j++)
		for (k = 0; k < 3; k++) {
			output_image->set_pixel_component(i, j, k, 
				(drizzle_weight->get_pixel_component(i, j, 0) > 0)
				? drizzle_image->get_pixel_component(i, j, k)
				: accum_image->  get_pixel_component(i, j, k));
		}	
	}

	write_image(filename, output_image);
}

/*
 * Increase extents of accumulated, drizzled, and output
 * images according to a new image to be merged.
 */

void increase_extents(transformation t) {
	int extend_top = 0;
	int extend_bottom = 0;
	int extend_left = 0;
	int extend_right = 0;

	my_real height = t.height(), width = t.width();

	point p[4];

	p[0] = t(point(   0  ,   0  ));
	p[1] = t(point(height,   0  ));
	p[2] = t(point(height, width));
	p[3] = t(point(   0  , width));

	for (int n = 0; n < 4; n++) {

		if (p[n][0] < -extend_offset_i - extend_top) {
			extend_top = (int) ceil(-p[n][0] - extend_offset_i);
		}
		if (p[n][0] > accum_image->height() - extend_offset_i + extend_bottom) {
			extend_bottom = (int) ceil(p[n][0] - accum_image->height() +
					extend_offset_i);
		}
		if (p[n][1] < -extend_offset_j - extend_left) {
			extend_left = (int) ceil(-p[n][1] - extend_offset_j);
		}
		if (p[n][1] > accum_image->width() - extend_offset_j + extend_right) {
			extend_right = (int) ceil(p[n][1] - accum_image->width() + 
					extend_offset_j);
		}
	}

	extend_offset_i += extend_top;
	extend_offset_j += extend_left;

	accum_image->extend(extend_top, extend_bottom, 
		extend_left, extend_right);
	accum_weight->extend(extend_top, extend_bottom, 
		extend_left, extend_right);

	if (drizzle_radius >= 0) {
		drizzle_image->extend(extend_top, extend_bottom, 
			extend_left, extend_right);
		drizzle_weight->extend(extend_top, extend_bottom, 
			extend_left, extend_right);
		if (!drizzle_only)
			output_image->extend(extend_top, 
				extend_bottom, extend_left, extend_right);
	}
}	

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
inline double diff(image *a, image *b, transformation t) {
	float result = 0;
	float divisor = 0;

	float apm_a = (extend || default_transform_type == 1)
	       	    ? (float) 1 : (float) a->avg_pixel_magnitude();
	float apm_b = (extend || default_transform_type == 1)
		    ? (float) 1 : (float) b->avg_pixel_magnitude();

	int i, j, k;

	for (i = 0; (unsigned int) i < a->height(); i++)
	for (j = 0; (unsigned int) j < a->width();  j++) {

		/*
		 * Transform
		 */

		struct point q;

		q = t.inverse_transform(
			point(i - extend_offset_i, j - extend_offset_j));

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
		  || (i >= extend_offset_i 
		   && j >= extend_offset_j 
		   && i <  extend_offset_i + extend_orig_height
		   && j <  extend_offset_j + extend_orig_width)
		  || accum_weight->get_pixel_component(i, j, 0) > 0)) { 

			if (align_code == 0) {
				/*
				 * Align based on all channels.
				 */

				for (k = 0; k < 3; k++) {
					float achan = a->get_pixel_component(i, j, k) / apm_a;
					float bchan = b->get_bl_component(ti, tj, k) / apm_b;

					result += pow(fabs(achan - bchan), metric);
					divisor += pow(achan > bchan ? achan : bchan, metric);
				}
			} else if (align_code == 1) {
				/*
				 * Align based on the green channel.  XXX: apm_* shouldn't be used here.
				 */

				float achan = a->get_pixel_component(i, j, 1) / apm_a;
				float bchan = b->get_bl_component(ti, tj, 1) / apm_b;

				result += pow(fabs(achan - bchan), metric);
				divisor += pow(achan > bchan ? achan : bchan, metric);
			} else if (align_code == 2) {
				/*
				 * Align based on the sum of all channels.
				 */

				float asum = 0;
				float bsum = 0;

				for (k = 0; k < 3; k++) {
					asum += a->get_pixel_component(i, j, k) / apm_a;
					bsum += b->get_bl_component(ti, tj, k) / apm_b;
				}

				result += pow(fabs(asum - bsum), metric);
				divisor += pow(asum > bsum ? asum : bsum, metric);
			}
		}
	}

	return pow(result / divisor, 1/metric);
}

/*
 * Drizzle part of a delta frame with part of a target frame using the specified
 * transformation.
 */
inline void
drizzle(image *target, image *delta, transformation t) {
	int i, j, k;

	for (i = 0; (unsigned int) i < target->height(); i++)
	for (j = 0; (unsigned int) j < target->width();  j++) {
		int ii, jj;

		/*
		 * Transform
		 */

		point q;

		q = t.inverse_transform(
			point(i - extend_offset_i, j - extend_offset_j));

		my_real ti = q[0];
		my_real tj = q[1];

		q = t.inverse_transform(
			point(i - extend_offset_i + 1, j - extend_offset_j));

		my_real ui = fabs(q[0] - ti);
		my_real uj = fabs(q[1] - tj);

		q = t.inverse_transform(
			point(i - extend_offset_i, j - extend_offset_j + 1));

		my_real vi = fabs(q[0] - ti);
		my_real vj = fabs(q[1] - tj);

		/*
		 * We map the area of the target pixel onto the delta
		 * frame as a rectangular area oriented on the delta
		 * frame's axes.  Note that this results in an area
		 * that may be the wrong shape or orientation.
		 *
		 * We define two estimates of the rectangle's
		 * dimensions below.  For rotations of 0, 90, 180, or
		 * 270 degrees, max and sum are identical.  For
		 * other orientations, sum is too large and max is
		 * too small.  We use the mean of max and sum, which we
		 * then divide by two to obtain the distance between
		 * the center and the edge.
		 */

		my_real maxi = (ui > vi) ? ui : vi;
		my_real maxj = (uj > vj) ? uj : vj;
		my_real sumi = ui + vi;
		my_real sumj = uj + vj;

		my_real di = (maxi + sumi) / 4;
		my_real dj = (maxj + sumj) / 4;

		ti /= scale_factor;
		tj /= scale_factor;
		di /= scale_factor;
		dj /= scale_factor;

		for (ii = (int) floor(ti - di - drizzle_radius); 
			ii <= ceil(ti + di + drizzle_radius); ii++)
		for (jj = (int) floor(tj - dj - drizzle_radius); 
			jj <= ceil(tj + dj + drizzle_radius); jj++) {
	
			my_real top = (ti - di < ii - drizzle_radius)
				    ? (ii - drizzle_radius)
				    : (ti - di);
			my_real bot = (ti + di > ii + drizzle_radius)
				    ? (ii + drizzle_radius)
				    : (ti + di);
			my_real lef = (tj - dj < jj - drizzle_radius)
				    ? (jj - drizzle_radius)
				    : (tj - dj);
			my_real rig = (tj + dj > jj + drizzle_radius)
				    ? (jj + drizzle_radius)
				    : (tj + dj);

			if (top < bot 
			 && lef < rig 
			 && ii >= 0
			 && ii * scale_factor < delta->height()
			 && jj >= 0
			 && jj * scale_factor < delta->width()) {

				double weight = drizzle_weight->get_pixel_component(i, j, 0);
				double thisw  = (bot - top) * (rig - lef) / (di * dj);

				drizzle_weight->set_pixel_component(i, j, 0, weight + thisw);

				for (k = 0; k < 3; k++)
					target->set_pixel_component(i, j, k, 
						(int) ((weight * target->get_pixel_component(
							i, j, k) + thisw
							 * delta->get_pixel_component(
								 (int) (ii * scale_factor), 
								 (int) (jj * scale_factor), 
								 k))
						/ (weight + thisw)));
			}
		}
	}
}

/*
 * Merge part of a delta frame with part of a target frame using the specified
 * transformation.
 */
inline void
merge(image *target, image *delta, transformation t) {
	int i, j, k;

	for (i = 0; (unsigned int) i < target->height(); i++)
		for (j = 0; (unsigned int) j < target->width(); j++) {

			/*
			 * Transform
			 */

			struct point q;

			q = t.inverse_transform(
				point(i - extend_offset_i, j - extend_offset_j));	

			my_real ti = q[0];
			my_real tj = q[1];

			/*
			 * Check that the transformed coordinates are within
			 * the boundaries of the delta frame.
			 */

			if (ti >= 0
			 && ti <= delta->height() - 1
			 && tj >= 0
			 && tj <= delta->width() - 1){ 

				/*
				 * Determine and update merging weight at this pixel
				 */

				double weight = accum_weight->get_pixel_component(i, j, 0);
				accum_weight->set_pixel_component(i, j, 0, weight + 1);

				/*
				 * Update each channel
				 */
				
				for (k = 0; k < 3; k++)
					target->set_pixel_component(i, j, k,
						(weight * target->get_pixel_component(i, j, k)
					     +   delta->get_bl_component(ti, tj, k)) / (weight + 1));
			}
		}
}

void filter(image *target, image *source, double scale_factor, double hf_enhance) {
	unsigned int i, j, k;

	assert (target->height() == source->height());
	assert (target->width() == source->width());
	assert (target->depth() == source->depth());

	for (i = 0; i < target->height(); i++)
		for (j = 0; j < target->width(); j++) 
			for (k = 0; k < 3; k++) {
				int result = (int) (hf_enhance 
						  * hf_filter(scale_factor, i, j, k, source)
					          + source->get_pixel_component(i, j, k));

				if (result < 0)
					result = 0;
				if (result > 255)
					result = 255;

				target->set_pixel_component(i, j, k, result);
			}
}


/*
 * Update an accumulated image with a new input image.
 */
inline void update() {
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

	accum_scales[0] = accum_image->clone();
	input_scales[0] = input_image->clone();
	
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

	if (default_transform_type == 0) 
		
		/*
		 * Identity transformation
		 */

		default_transform = (transform_code == 2)
				  ?  transformation::gpt_identity(input_image)
				  :  transformation:: eu_identity(input_image);

	else if (default_transform_type == 1)

		/*
		 * Follow previous transformation
		 */

		default_transform.rescale(input_image);

	else
		assert(0);

	offset = tload_next(tload, lod, transform_code == 2, default_transform);

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

		if (transform_code < 2 && transform_code >= 0) {

			/* 
			 * Translational or euclidean transformation
			 */

			for (i = 0; i < 2 + transform_code; i++)
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
		
		} else if (transform_code == 2) {

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

	tsave_next(tsave, offset, transform_code == 2);

	/*
	 * Ensure that the match meets the threshold.
	 */

	if ((1 - here) * 100 > match_threshold) {
		if (extend == 1)
			increase_extents(offset);
		merge (accum_image, input_image, offset);
		if (drizzle_radius >= 0)
			drizzle (drizzle_image, input_image, offset);
		default_transform = offset;
		fprintf(stderr, " okay (%f%% match)", (1 - here) * 100);
	} else {
		fprintf(stderr, " no match (%f%% match)", (1 - here) * 100);
	}

	match_sum += (1 - here) * 100;
	match_count++;
}

int main(int argc, char *argv[]){
	int i;

	if (argc == 2 && !strcmp(argv[1], "--version")) {

		/*
		 * Output the version
		 */

		fprintf(stderr, "%s\n", version);

		return 0;
	}

	/* 
	 * Iterate through all arguments but the last argument; the last
	 * argument is the output file, which we update incrementally.
	 */

	for (i = 1; i < argc - 1; i++) {


		if (!strcmp(argv[i], "--scale2")) {
			scale_factor = 2;
			fprintf(stderr, "\n\n*** Warning: --scale2 is deprecated.  "
					"Use --scale=2 instead. ***\n\n\n");
		} else if (!strcmp(argv[i], "--scale4")) {
			scale_factor = 4;
			fprintf(stderr, "\n\n*** Warning: --scale4 is deprecated.  "
					"Use --scale=4 instead. ***\n\n\n");
		} else if (!strcmp(argv[i], "--scale8")) {
			scale_factor = 8;
			fprintf(stderr, "\n\n*** Warning: --scale8 is deprecated.  "
					"Use --scale=8 instead. ***\n\n\n");
		} else if (!strcmp(argv[i], "--align-all")) {
			align_code = 0;
		} else if (!strcmp(argv[i], "--align-green")) {
			align_code = 1;
		} else if (!strcmp(argv[i], "--align-sum")) {
			align_code = 2;
		} else if (!strcmp(argv[i], "--translation")) {
			transform_code = 0;
		} else if (!strcmp(argv[i], "--euclidean")) {
			transform_code = 1;
		} else if (!strcmp(argv[i], "--projective")) {
			transform_code = 2;
		} else if (!strcmp(argv[i], "--identity")) {
			default_transform_type = 0;
		} else if (!strcmp(argv[i], "--follow")) {
			default_transform_type = 1;
		} else if (!strcmp(argv[i], "--no-extend")) {
			extend = 0;
		} else if (!strcmp(argv[i], "--extend")) {
			extend = 1;
		} else if (!strcmp(argv[i], "--drizzle-only")) {
			drizzle_only = 1;
		} else if (!strncmp(argv[i], "--scale=", strlen("--scale="))) {
			sscanf(argv[i] + strlen("--scale="), "%lf", &scale_factor);
			if (scale_factor < 1.0) {
				fprintf(stderr, "\n\n*** Warning: Ignoring scale value "
						"smaller than 1.0 ***\n\n\n");
			}
		} else if (!strncmp(argv[i], "--metric=", strlen("--metric="))) {
			sscanf(argv[i] + strlen("--metric="), "%lf", &metric);
		} else if (!strncmp(argv[i], "--threshold=", strlen("--threshold="))) {
			sscanf(argv[i] + strlen("--threshold="), "%lf", &match_threshold);
		} else if (!strncmp(argv[i], "--drizzle-diam=", strlen("--drizzle-diam="))) {
			sscanf(argv[i] + strlen("--drizzle-diam="), "%lf", &drizzle_radius);
			drizzle_radius /= 2;
		} else if (!strncmp(argv[i], "--perturb-upper=", strlen("--perturb-upper="))) {
			sscanf(argv[i] + strlen("--perturb-upper="), "%lf", &perturb_upper);
		} else if (!strncmp(argv[i], "--stepsize=", strlen("--stepsize="))) {
			fprintf(stderr, "\n\n*** Warning: --stepsize is deprecated.  "
					"Use --perturb-lower instead. ***\n\n\n");
			sscanf(argv[i] + strlen("--stepsize="), "%lf", &perturb_lower);
		} else if (!strncmp(argv[i], "--hf-enhance=", strlen("--hf-enhance="))) {
			sscanf(argv[i] + strlen("--hf-enhance="), "%lf", &hf_enhance);
		} else if (!strncmp(argv[i], "--perturb-lower=", strlen("--perturb-lower="))) {
			sscanf(argv[i] + strlen("--perturb-lower="), "%lf", &perturb_lower);
		} else if (!strncmp(argv[i], "--trans-load=", strlen("--trans-load="))) {
			tload_delete(tload);
			tload = tload_new(argv[i] + strlen("--trans-load="));
		} else if (!strncmp(argv[i], "--trans-save=", strlen("--trans-save="))) {
			tsave_delete(tsave);
			tsave = tsave_new(argv[i] + strlen("--trans-save="));
		} else if (accum_image == NULL) {

			/* 
			 * First file argument.  Print general file information as well
			 * as information specific to this argument.  Initialize image
			 * file handler.
			 */

			init_image();

			tsave_orig(tsave, argv[i]);
			tsave_target(tsave, argv[argc - 1]);

			fprintf(stderr, "Output file will be '%s'.\n", 
					argv[argc - 1]);

			fprintf(stderr, "Reading original frame '%s'", argv[i]);

			accum_image = read_image(argv[i]);

			if (scale_factor != 1)
				accum_image->scale(scale_factor);

			extend_orig_height = accum_image->height();
			extend_orig_width  = accum_image->width();

			accum_weight = new_image_weights(accum_image->height(), 
					                 accum_image->width(), 1);

			merge(accum_image, accum_image, 
					transformation::eu_identity(accum_image));

			if (drizzle_radius >= 0) {

				drizzle_image  = new_image(accum_image->height(), 
							   accum_image->width(), 3);

				drizzle_weight = new_image_weights(accum_image->height(),
								   accum_image->width(), 1);

				drizzle(drizzle_image, accum_image, 
					transformation::eu_identity(accum_image));

			}

			if (drizzle_radius >= 0 && drizzle_only == 0)
				output_image = new_image(accum_image->height(),
							 accum_image->width(), 3);
			else if (drizzle_only)
				output_image = drizzle_image;
			else
				output_image = accum_image;

			default_transform = (transform_code == 2)
				          ? transformation
					  	:: gpt_identity(accum_image)
				          : transformation
						::  eu_identity(accum_image);

			output(argv[argc - 1]);

			fprintf(stderr, ".\n");

		} else {

			tsave_info (tsave, argv[i]);

			fprintf(stderr, "Merging supplemental frame '%s'", 
					argv[i]);

			input_image = read_image(argv[i]);
			if (scale_factor != 1)
				input_image->scale(scale_factor);
			update();

			delete input_image;

			output(argv[argc - 1]);

			fprintf(stderr, ".\n");
			
		}

	}

	if (accum_image) {

		if (hf_enhance != 0) {

			image *ppsource = output_image->clone();
						    
			fprintf(stderr, "Enhancing high frequencies");

			filter(output_image, ppsource, scale_factor, hf_enhance);

			write_image(argv[argc - 1], output_image);

			fprintf(stderr, ".\n");
		}

		/*
		 * Destroy the image file handler
		 */

		destroy_image();

		/*
		 * Delete the transformation file structures
		 */

		tsave_delete(tsave);
		tload_delete(tload);

		/*
		 * Output an summary match statistic.
		 */

		fprintf(stderr, "Done (%f%% average match).\n", match_sum / match_count);
		
	} else {

		/*
		 * If there was no output, the user might need more information.
		 */

		fprintf(stderr, 
			"\n"
			"Usage: %s [<options>] <original-frame> [<supplemental-frame> ...] <output-file>\n"
			"   or: %s --version\n"
			"\n\n"
			"Alignment channel options:\n\n"
			"--align-all       Align images using all color channels\n"
			"--align-green     Align images using the green channel\n"
			"--align-sum       Align images using a sum of channels [default]\n"
			"\n\n"
			"Transformation options:\n\n"
			"--translation     Only adjust the position of images\n"
			"--euclidean       Adjust the position and orientation of images [default]\n"
			"--projective      Use projective transformations.  Best quality, but slow.\n"
			"\n\n"
			"Image extents:\n\n"
			"--extend          Increase image extents to accommodate all pixel data.\n"
			"--no-extend       Don't increase extents; crop to original frame. [default]\n"
			"\n\n"
			"Transformation defaults:\n\n"
			"--identity        Default alignment begins with identity transform. [default]\n"
			"--follow          Default alignment begins with the previous frame's transform.\n"
			"\n\n"
			"Transformation file operations:\n\n"
			"--trans-load=x    Load initial transformation settings from file x\n"
			"--trans-save=x    Save final transformation data in file x\n"
			"\n\n"
			"Tunable parameters:\n\n"
			"--scale=x         Scale images by the factor x (where x is at least 1.0)\n"
			"--hf-enhance=x    Enhance high frequency details by factor x. (0.0 is default)\n"
			"--metric=x        Set the error metric exponent (2 is default)\n"
			"--threshold=x     Min. match threshold; a perfect match is 100.  (0 is default)\n"
			"--perturb-upper=x Max. correction step, in pixels or degrees (32.0 is default)\n"
			"--perturb-lower=x Min. correction step, in pixels or degrees (0.125 is default)\n"
			"\n\n"
			"Drizzling [very experimental]:\n\n"
			"--drizzle-diam=x  Drizzle with input pixel diameter x (between 0 and 1).\n"
			"--drizzle-only    If drizzling, output black for pixels with no drizzle data.\n"
			"\n",
			argv[0], argv[0]);
		return 1;
	}

	return 0;
}

