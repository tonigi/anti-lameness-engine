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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#ifndef __image_h__
#define __image_h__

typedef struct {
	unsigned int dimx, dimy, depth;
	unsigned char *p;
	double apm_memo;
} image;

static inline image *new_image(unsigned int dimy, unsigned int dimx, unsigned int depth) {
	image *im = (image *) malloc(sizeof(image));

	assert(im);

	im->dimx = dimx;
	im->dimy = dimy;
	im->depth = depth;
	im->p = (unsigned char *) calloc(dimx * dimy * depth, sizeof(unsigned char));
	im->apm_memo = -1;

	assert (im->p);

	return im;
}

static inline void image_final(image *im) {
	free(im->p);
}

static inline unsigned int width(image *im) {
	return im->dimx;
}

static inline unsigned int height(image *im) {
	return im->dimy;
}

static inline unsigned char get_pixel_component(image *im, unsigned int y, unsigned int x, unsigned int color) {
	assert (x < im->dimx);
	assert (y < im->dimy);
	assert (color < im->depth);
	return im->p[y * im->dimx * im->depth + x * im->depth + color];
}

static inline void set_pixel_component(image *im, unsigned int y, unsigned int x, unsigned int color, unsigned char value) {
	assert (x < im->dimx);
	assert (y < im->dimy);
	assert (color < im->depth);

	/* XXX: We could do much better than this. */
	im->apm_memo = -1;

	im->p[y * im->dimx * im->depth + x * im->depth + color] = value;
}

/*
 * Get a color value at a given position using bilinear interpolation between the
 * four nearest pixels.
 */
static inline unsigned char get_bl_component(image *im, double y, double x, unsigned int color) {
	int lx = (int) floor(x);
	int hx = (int) floor(x) + 1;
	int ly = (int) floor(y);
	int hy = (int) floor(y) + 1;

	assert (y >= 0);
	assert (x >= 0);
	assert (y <= im->dimy - 1);
	assert (x <= im->dimx - 1);
	assert (color < im->depth);

	return (int) ((hx - x) * (hy - y)
			       * get_pixel_component(im, ly, lx, color)
		    + (hx - x) * (y - ly)
			       * get_pixel_component(im, hy % im->dimy, lx, color)
		    + (x - lx) * (y - ly)
			       * get_pixel_component(im, hy % im->dimy, hx % im->dimx, color)
		    + (x - lx) * (hy - y)
			       * get_pixel_component(im, ly, hx % im->dimx, color));
}

/*
 * Scale an image by some factor
 */
static inline void scale(image *im, double f) {
	image *is = new_image((int) floor(height(im) * f), 
			      (int) floor(width(im) * f), im->depth);

	unsigned int i, j, k;

	for (i = 0; i < height(is); i++)
		for (j = 0; j < width(is); j++) 
			for (k = 0; k < is->depth; k++)
				set_pixel_component(is, i, j, k,
					get_bl_component(im, (i/f <= height(im) - 1) ? (i/f) : (height(im) - 1), 
							     (j/f <= width(im) - 1)  ? (j/f) : (width(im) - 1), k));

	free(im->p);

	im->p = is->p;
	im->dimx = is->dimx;
	im->dimy = is->dimy;
	im->depth = is->depth;
	im->apm_memo = -1;

	free(is);
}

/*
 * Clone an image
 */
static inline image *clone(image *im) {
	image *ic = new_image(height(im), width(im), im->depth);

	unsigned int i, j, k;

	for (i = 0; i < height(ic); i++)
		for (j = 0; j < width(ic); j++)
			for (k = 0; k < ic->depth; k++)
				set_pixel_component(ic, i, j, k,
					get_pixel_component(im, i, j, k));

	ic->apm_memo = im->apm_memo;

	return ic;
}

static inline double avg_pixel_magnitude(image *im) {
	unsigned int i, j, k;

	if (im->apm_memo >= 0)
		return im->apm_memo;

	im->apm_memo = 0;

	for (i = 0; i < im->dimy; i++)
		for (j = 0; j < im->dimx; j++)
			for (k = 0; k < im->depth; k++)
				im->apm_memo += get_pixel_component(im, i, j, k);

	im->apm_memo /= (im->dimy * im->dimx * im->depth);

	return im->apm_memo;
}

#endif
