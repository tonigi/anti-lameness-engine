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

#ifndef __ppm_h__
#define __ppm_h__

#include <stdio.h>
#include <assert.h>
#include "image.h"

static inline void error_ppm(char *filename) {
	fprintf(stderr, 
		"\n'%s' doesn't look like a binary PPM file.\n", 
		filename);
	exit(1);
}

static inline image *read_ppm(char *filename) {
	unsigned int i, j, k;
	image *im;
	char m1, m2, val;
	int w, h, mcv;
	int n;
	FILE *f = fopen(filename, "r");
	assert(f);

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

	/* Magic */

	n = fscanf(f, "%c%c", &m1, &m2);
	assert(n == 2 && m1 == 'P' && m2 == '6');

	if (n != 2 || m1 != 'P' || m2 != '6')
		error_ppm(filename);

	/* Width */

	n = fscanf(f, "%d", &w);
	assert(n == 1);

	if (n != 1)
		error_ppm(filename);

	/* Height */

	n = fscanf(f, "%d", &h);
	assert(n == 1);

	if (n != 1)
		error_ppm(filename);

	/* Maximum component value */

	n = fscanf(f, "%d", &mcv);
	assert(n == 1);
	assert(mcv == 255);

	if (n != 1)
		error_ppm(filename);

	/* Trailing whitespace */

	fgetc(f);

	/* Make a new image */

	im = new_image(h, w, 3);

	/* Pixels */

	for (i = 0; i < height(im); i++)
		for (j = 0; j < width(im); j++)
			for (k = 0; k < 3; k++) {
				fscanf(f, "%c", &val);
				set_pixel_component(im, i, j, k, val);
			}

	/* Done */

	fclose(f);

	return im;
}

static inline void write_ppm(char *filename, image *im) {
	unsigned int i, j, k;
	FILE *f = fopen(filename, "w");
	assert(f);

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

	/* Magic */

	fprintf(f, "P6 ");

	/* Width */

	fprintf(f, "%d ", width(im));

	/* Height */

	fprintf(f, "%d ", height(im));

	/* Maximum component value */

	fprintf(f, "255\n");

	/* Pixels */

	for (i = 0; i < height(im); i++)
		for (j = 0; j < width(im); j++)
			for (k = 0; k < 3; k++)
				fprintf(f, "%c", get_pixel_component(im, i, j, k));

	/* Done */

	fclose(f);
}

#endif
