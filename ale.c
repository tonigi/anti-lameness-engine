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
#include "ppm.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

char *version = "0.0.0";

image *display_image, *input_image, *weight_image;

/*
 * Not-quite-symmetric difference function.  Determines the difference in areas
 * where the arrays overlap.  Uses the first array's notion of pixel positions.
 */
inline double diff(image *a, image *b, double xoff, double yoff) {
	double result = 0;
	double apm_a = avg_pixel_magnitude(a);
	double apm_b = avg_pixel_magnitude(b);

	unsigned int ha = height(a);
	unsigned int wa = width(a);
	unsigned int hb = height(b);
	unsigned int wb = width(b);

	/* Determine the start of the common area of the images */

	unsigned int ystart = (unsigned int)
		((yoff > 0) ? ceil(yoff) : 0);
	unsigned int xstart = (unsigned int)
		((xoff > 0) ? ceil(xoff) : 0);

	/* Determine the overlapping width and height of the images */

	unsigned int hc = (unsigned int)
		((yoff + hb > ha) ? ha : floor(yoff + hb));
	unsigned int wc = (unsigned int)
		((xoff + wb > wa) ? wa : floor(xoff + wb));

	int i, j, k;

	// assert(ha == hb);
	// assert(wa == wb);

	for (i = ystart; i < hc; i++)
		for (j = xstart; j < wc; j++)
			for (k = 0; k < 3; k++)
				result += 
					fabs(apm_b * get_pixel_component(a, i, j, k) - 
					     apm_a * get_bl_component(b, i - yoff, j - xoff, k));

	return (result / hc) / wc;
}

/*
 * Merge part of a delta frame with part of a target frame using the specified
 * offset.
 */
inline void merge(image *target, image *delta, double xoff, double yoff) {

	unsigned int ht = height(target);
	unsigned int wt = width(target);
	unsigned int hd = height(delta);
	unsigned int wd = width(delta);

	/* Determine the start of the common area of the images */

	unsigned int ystart = (unsigned int)
		((yoff > 0) ? ceil(yoff) : 0);
	unsigned int xstart = (unsigned int)
		((xoff > 0) ? ceil(xoff) : 0);	

	/* Determine the overlapping width and height of the images */

	unsigned int hc = (unsigned int)
		((yoff + hd > ht) ? ht : floor(yoff + hd));
	unsigned int wc = (unsigned int)
		((xoff + wd > wt) ? wt : floor(xoff + wd));

	int i, j, k;

	for (i = ystart; i < hc; i++)
		for (j = xstart; j < wc; j++) {

			int weight = get_pixel_component(weight_image, i, j, 0);
			set_pixel_component(weight_image, i, j, 0, ++weight);
			
			for (k = 0; k < 3; k++)
				set_pixel_component(target, i, j, k,
					(weight * get_pixel_component(target, i, j, k)
				     +   get_bl_component(delta, i - yoff, j - xoff, k)) / (weight + 1));
		}

}


/*
 * Update an accumulated image with a new input image.
 */
inline void update() {
	double perturb = 16;
	double xoff = 0, yoff = 0;
	double here = diff(display_image, input_image, 0, 0);

	double u = diff(display_image, input_image, xoff, yoff - perturb);
	double d = diff(display_image, input_image, xoff, yoff + perturb);
	double l = diff(display_image, input_image, xoff - perturb, yoff);
	double r = diff(display_image, input_image, xoff + perturb, yoff);

	while (perturb > .1) {
		// fprintf(stderr, "(%f, %f)\n", xoff, yoff);
		// fprintf(stderr, "diff: (h, u, d, l, r) %f, %f, %f, %f, %f\n", here, u, d, l, r);

		if (u < here) {
			yoff -= perturb;
			d = here;
			here = u;
			l = diff(display_image, input_image, xoff - perturb, yoff);
			r = diff(display_image, input_image, xoff + perturb, yoff);
			u = diff(display_image, input_image, xoff, yoff - perturb);
		} else if (d < here) {
			yoff += perturb;
			u = here;
			here = d;
			l = diff(display_image, input_image, xoff - perturb, yoff);
			r = diff(display_image, input_image, xoff + perturb, yoff);
			d = diff(display_image, input_image, xoff, yoff + perturb);
		} else if (l < here) {
			xoff -=perturb;
			r = here;
			here = l;
			u = diff(display_image, input_image, xoff, yoff - perturb);
			d = diff(display_image, input_image, xoff, yoff + perturb);
			l = diff(display_image, input_image, xoff - perturb, yoff);
		} else if (r < here) {
			xoff += perturb;
			l = here;
			here = r;
			u = diff(display_image, input_image, xoff, yoff - perturb);
			d = diff(display_image, input_image, xoff, yoff + perturb);
			r = diff(display_image, input_image, xoff + perturb, yoff);
		} else {
			perturb *= 0.5;
			u = diff(display_image, input_image, xoff, yoff - perturb);
			d = diff(display_image, input_image, xoff, yoff + perturb);
			l = diff(display_image, input_image, xoff - perturb, yoff);
			r = diff(display_image, input_image, xoff + perturb, yoff);

			/*
			 * Announce that we've dropped a perturbation level.
			 */

			fprintf(stderr, ".");

		}
	}


	merge (display_image, input_image, xoff, yoff);
}

int main(int argc, char *argv[]){
	double res_scale = 1;
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
			res_scale = 2;
		} else if (!strcmp(argv[i], "--scale4")) {
			res_scale = 4;
		} else if (!strcmp(argv[i], "--scale8")) {
			res_scale = 8;
		} else if (display_image == NULL) {

			/* 
			 * First file argument.  Print general file information as well
			 * as information specific to this argument.
			 */

			fprintf(stderr, "Output file will be '%s'.\n", 
					argv[argc - 1]);

			fprintf(stderr, "Reading original frame '%s'", argv[i]);

			display_image = read_ppm(argv[i]);

			if (res_scale != 1)
				scale(display_image, res_scale);
			weight_image = new_image(height(display_image), 
						 width(display_image), 1);

			write_ppm(argv[argc - 1], display_image);

			fprintf(stderr, ".\n");

		} else {

			fprintf(stderr, "Merging supplemental frame '%s'", 
					argv[i]);

			input_image = read_ppm(argv[i]);
			if (res_scale != 1)
				scale(input_image, res_scale);
			update();

			write_ppm(argv[argc - 1], display_image);

			fprintf(stderr, ".\n");
			
		}

	}

	/*
	 * If there was no output, the user might need more information.
	 */

	if (display_image == NULL) {
		fprintf(stderr, 
			"Usage: %s [--scale2] [--scale4] [--scale8] <input-files> ... <output-file>\n"
			"   or: %s --version\n",
			argv[0], argv[0]);
		return 1;
	} else {
		fprintf(stderr, "Done.\n");
	}

	return 0;
}

