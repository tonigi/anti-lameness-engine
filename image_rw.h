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
#include <unistd.h>
#include <sys/stat.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "image.h"
#include "ppm.h"

#ifdef USE_MAGICK
#include <magick/api.h>
#endif

class image_rw {

private:

	/*
	 * Private methods to init and shut down the file reader.
	 */
	
	/*
	 * Initialize the image file handler
	 */
	static void init_image() {
#ifdef USE_MAGICK
		InitializeMagick("ale");
#endif
	}

	/*
	 * Destroy the image file handler
	 */
	static void destroy_image() {
#ifdef USE_MAGICK
		DestroyMagick();
#endif
	}


	/*
	 * Private methods to read and write image files
	 */

	/*
	 * Read an image from a file
	 */
	static image *read_image(const char *filename)
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

		unsigned int i, j;

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
				im->set_pixel_component(i, j, 0, p->red);
				im->set_pixel_component(i, j, 1, p->green);
				im->set_pixel_component(i, j, 2, p->blue);
				p++;
			}
		}

		DestroyImage(mi);
		DestroyImageInfo(image_info);

		return im;
	}	
#else
	{
		return read_ppm(filename);
	}
#endif

	/*
	 * Write an image to a file
	 */
	static void write_image(const char *filename, const image *im) {
#ifdef USE_MAGICK
		
		/*
		 * Patterned after http://www.imagemagick.org/www/api.html
		 * and http://www.imagemagick.org/www/smile.c
		 */

		ExceptionInfo exception;
		Image *mi;
		ImageInfo *image_info;
		PixelPacket *p;

		unsigned int i, j;

		GetExceptionInfo(&exception);
		image_info = CloneImageInfo((ImageInfo *) NULL);
		strncpy(image_info->filename, filename, MaxTextExtent);

		mi = AllocateImage(image_info);
		if (mi == (Image *) NULL) 
			MagickError(ResourceLimitError,
				"Unable to display image", "MemoryAllocationFailed");

		mi->columns = im->width();
		mi->rows = im->height();
		mi->depth = 8;

		for (i = 0; i < mi->rows; i++) {
			p = SetImagePixels(mi, 0, i, mi->columns, 1);
			if (p == NULL)
				break;
			for (j = 0; j < mi->columns; j++) {
				p->red = im->get_pixel_component(i, j, 0);
				p->green = im->get_pixel_component(i, j, 1);
				p->blue = im->get_pixel_component(i, j, 2);
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
#else
		write_ppm(filename, im);
#endif
	}

	/*
	 * See if a file works with stat().
	 */
	static void try_stat(const char *filename) {
		struct stat s;
		if (stat(filename, &s)) {
			fprintf(stderr, "Error trying to stat file '%s'.\n", 
					filename);
			exit(1);
		}
	}

	/*
	 * Private data members
	 */

	static const char **filenames;
	static int filename_count;
	static const char *output_filename;
	static const image **images;
	static int *files_open;
	static double scale_factor;

	static int latest_close_num;
	static const image *latest_close;

public:

	/*
	 * Initializer.
	 *
	 * Handle COUNT input files with names in array FILENAMES and output
	 * file OUTPUT_FILENAME.  FILENAMES should be an array of char * that
	 * is never freed.  OUTPUT_FILENAME should be a char * that is never
	 * freed.
	 */
	static void init(int _filename_count, const char **_filenames, 
			const char *_output_filename, double _scale_factor){
		assert (_filename_count > 0);

		init_image();

		filenames = _filenames;
		filename_count = _filename_count;
		output_filename = _output_filename;
		scale_factor = _scale_factor;

		images = (const image **)malloc(filename_count * sizeof(image *));
		files_open = (int *)calloc(filename_count, sizeof(int));

		assert (images);
		assert (files_open);

		if (!images || !files_open) {
			fprintf(stderr, "Unable to allocate memory for images.\n");
			exit(1);
		}

		for (int i = 0; i < filename_count; i++) {
			try_stat(filenames[i]);
		}

		fprintf(stderr, "Output file will be '%s'.\n", 
				output_filename);

	}

	static void destroy() {
		assert (filename_count > 0);
		destroy_image();
	}

	static int count() {
		assert (filename_count > 0);
		return filename_count;
	}

	static const char *name(int image) {
		assert (image >= 0);
		assert (image <  filename_count);

		return filenames[image];
	}

	static const char *output_name() {
		assert (filename_count > 0);
		return output_filename;
	}

	static void output(const image *i) {
		assert (filename_count > 0);
		write_image(output_filename, i);
	}

	static const image *open(int n) {
		assert (n >= 0);
		assert (n <  filename_count);
		assert (!files_open[n]);

		if (n == latest_close_num) {
			images[n] = latest_close;
			latest_close_num = -1;
			files_open[n] = 1;
		} else {
			image *i = read_image(filenames[n]);

			i->scale(scale_factor);

			images[n] = i;
			files_open[n] = 1;
		}
		return images[n];
	}

	static image *copy(int n) {
		assert (n >= 0);
		assert (n <  filename_count);

		if (files_open[n])
			return images[n]->clone();
		else {
			image *i = read_image(filenames[n]);
			i->scale(scale_factor);
			return i;
		}
	}

	static void close(int image) {
		assert (image >= 0);
		assert (image <  filename_count);
		assert (files_open[image]);

		if (latest_close_num >= 0)
			delete latest_close;

		latest_close = images[image];
		latest_close_num = image;

		files_open[image] = 0;
	}

};

#endif
