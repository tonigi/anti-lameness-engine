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
 * image_rw.h: Read and write images.
 */


#ifndef __image_rw_h__
#define __image_rw_h__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "image.h"
#include "ppm.h"

#ifdef USE_MAGICK
#include <magick/api.h>
#endif

static inline image *read_image(char *filename)
#ifdef USE_MAGICK
{
	/*
	 * Patterned after http://www.imagemagick.org/www/api.html
	 * and http://www.imagemagick.org/www/smile.c
	 */

	ExceptionInfo exception;
	Image *mi;
	ImageInfo *image_info;
	image *im;
	const PixelPacket *p;

	int i, j;

	InitializeMagick("ale");
	GetExceptionInfo(&exception);
	image_info = CloneImageInfo((ImageInfo *) NULL);

	strncpy(image_info->filename, filename, MaxTextExtent);
	mi = ReadImage(image_info, &exception);
	if (exception.severity != UndefinedException) {
		fprintf(stderr, "\n\n");
		CatchException(&exception);
		fprintf(stderr, "\n");
	}
	if (mi == (Image *) NULL)
		exit(1);

	im = new_image(mi->rows, mi->columns, 3);

	for (i = 0; i < mi->rows; i++) {
		p = AcquireImagePixels(mi, 0, i, mi->columns, 1, &exception);

		if (exception.severity != UndefinedException)
			CatchException(&exception);
		if (p == NULL)
			exit(1);
	
		for (j = 0; j < mi->columns; j++) {
			set_pixel_component(im, i, j, 0, p->red);
			set_pixel_component(im, i, j, 1, p->green);
			set_pixel_component(im, i, j, 2, p->blue);
			p++;
		}
	}

	DestroyImage(mi);
	DestroyImageInfo(image_info);
	DestroyMagick();

	return im;
}	
#else
{
	return read_ppm(filename);
}
#endif

static inline void write_image(char *filename, image *im) {
#ifdef USE_MAGICK
	
	/*
	 * Patterned after http://www.imagemagick.org/www/api.html
	 * and http://www.imagemagick.org/www/smile.c
	 */

	ExceptionInfo exception;
	Image *mi;
	ImageInfo *image_info;
	PixelPacket *p;

	int i, j;

	InitializeMagick("ale");
	GetExceptionInfo(&exception);
	image_info = CloneImageInfo((ImageInfo *) NULL);
	strncpy(image_info->filename, filename, MaxTextExtent);

	mi = AllocateImage(image_info);
	if (mi == (Image *) NULL) 
		MagickError(ResourceLimitError,
			"Unable to display image", "MemoryAllocationFailed");

	mi->columns = width(im);
	mi->rows = height(im);
	mi->depth = 8;

	for (i = 0; i < mi->rows; i++) {
		p = SetImagePixels(mi, 0, i, mi->columns, 1);
		if (p == NULL)
			break;
		for (j = 0; j < mi->columns; j++) {
			p->red = get_pixel_component(im, i, j, 0);
			p->green = get_pixel_component(im, i, j, 1);
			p->blue = get_pixel_component(im, i, j, 2);
			p->red |= p->red << 8;
			p->green |= p->green << 8;
			p->blue |= p->blue << 8;
			p++;
		}
		if (!SyncImagePixels(mi))
			break;
	}

	if (!WriteImage(image_info, mi)) {
		fprintf(stderr, "\n\nImageMagick Error.  Perhaps you need to indicate the file type of the\n");
		fprintf(stderr, "output file by providing an appropriate extension (e.g. .png, .jpg, .tiff).\n\n");
		CatchException(&mi->exception);
		exit(1);
	}

	DestroyImage(mi);
	DestroyImageInfo(image_info);
	DestroyMagick();
#else
	write_ppm(filename, im);
#endif
}

#endif
