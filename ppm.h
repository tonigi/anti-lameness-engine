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

static inline void error_ppm(const char *filename) {
	fprintf(stderr, 
		"\n'%s' doesn't look like a binary PPM file.\n", 
		filename);
	exit(1);
}

static inline void eat_comments(FILE *f, const char *filename) {
#if 0

	/*
	 * %[] is apparently not compatible with WINE 20021007.
	 */

	char comment[2];

	while (fscanf(f, " %1[#]", comment)) {
		int a = 0;
		while (a != '\n') {
			if (feof(f))
				error_ppm(filename);
			a = fgetc(f);
		}
	}
#else

	/*
	 * ... so we try a different approach.
	 */

	int next = ' ';

	while (next == ' ' || next == '\n' || next == '\t' || next == '#') {
		next = fgetc(f);
		if (next == '#')
			while (next != '\n' && next != EOF)
				next = fgetc(f);
		if (feof(f))
			error_ppm(filename);
	}

	if (ungetc(next, f) == EOF) {
		assert(0);
		fprintf(stderr, "Unable to ungetc().");
		exit(1);
	}
#endif
}

static inline image *read_ppm(const char *filename) {
	unsigned int i, j, k;
	image *im;
	char m1, m2, val;
	int w, h, mcv;
	int n;
	FILE *f = fopen(filename, "rb");
	assert(f);

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

	/* Magic */

	eat_comments(f, filename);	/* XXX - should we eat comments here? */
	n = fscanf(f, "%c%c", &m1, &m2);
	assert(n == 2 && m1 == 'P' && m2 == '6');

	if (n != 2 || m1 != 'P' || m2 != '6')
		error_ppm(filename);

	/* Width */

	eat_comments(f, filename);
	n = fscanf(f, " %d", &w);
	assert(n == 1);

	if (n != 1)
		error_ppm(filename);

	/* Height */

	eat_comments(f, filename);
	n = fscanf(f, "%d", &h);
	assert(n == 1);

	if (n != 1)
		error_ppm(filename);

	/* Maximum component value */

	eat_comments(f, filename);
	n = fscanf(f, "%d", &mcv);
	assert(n == 1);
	assert(mcv == 255);

	if (n != 1)
		error_ppm(filename);

	/* Trailing whitespace */

	fgetc(f);

	/* Make a new image */

	im = new image(h, w, 3);
	assert (im);

	/* Pixels */

#if 1
	for (i = 0; i < im->height(); i++)
	for (j = 0; j < im->width();  j++)
	for (k = 0; k < im->depth();  k++) {
		fscanf(f, "%c", &val);
		im->set_pixel_component(i, j, k, val);
	}
#else
	/* 15% improvement in speed; tested for ALE 0.4.5  */

	fread(im->get_pixel_array(), 1, im->height() * im->width() * im->depth(), f);
#endif

	/* Done */

	fclose(f);

	return im;
}

static inline void write_ppm(const char *filename, const image *im) {
	unsigned int i, j, k;
	FILE *f = fopen(filename, "wb");
	assert(f);

	if (f == NULL) {
		fprintf(stderr, "\n\nUnable to open '%s'.\n\n", filename);
		exit(1);
	}

	/* Magic */

	fprintf(f, "P6 ");

	/* Width */

	fprintf(f, "%d ", im->width());

	/* Height */

	fprintf(f, "%d ", im->height());

	/* Maximum component value */

	fprintf(f, "255\n");

	/* Pixels */

	for (i = 0; i < im->height(); i++)
	for (j = 0; j < im->width();  j++)
	for (k = 0; k < im->depth();  k++)
		fprintf(f, "%c", im->get_pixel_component(i, j, k));

	/* Done */

	fclose(f);
}

#endif
