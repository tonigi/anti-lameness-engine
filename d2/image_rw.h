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
 * image_rw.h: Read and write images.
 */


#ifndef __image_rw_h__
#define __image_rw_h__

#include "image.h"
#include "image_ale_real.h"
#include "image_bayer_ale_real.h"
#include "ppm.h"
#include "exposure/exposure.h"
#include "exposure/exposure_default.h"

class image_rw {

	/*
	 * Private data members
	 */

	/*
	 * PPM type
	 *
	 * 0 = No type selected
	 * 1 = PPM Raw
	 * 2 = PPM Plain
	 */
	static int ppm_type;

	/*
	 * Bit depth
	 */
	static unsigned int num_bits;
	static unsigned int mcv;

	/*
	 * Input and output exposure models
	 */
	static exposure *input_exposure;
	static exposure *output_exposure;
	static int exposure_scale;

	/*
	 * Default bayer pattern
	 */
	static unsigned int bayer;

	/*
	 * Pointer to the output filename
	 */
	static const char *output_filename;
	
	/*
	 * Variables relating to input image files and image data structures.
	 */
	static const char **filenames;
	static unsigned int file_count;
	static const image **images;
	static int *files_open;

	/*
	 * The most recently closed image number.
	 */
	static int latest_close_num;

	/*
	 * A cache of the data associated with the most recently closed image.
	 */
	static const image *latest_close;

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

public:

	/*
	 * Methods to read and write image files
	 */

	/*
	 * Read an image from a file
	 */
	static image *read_image(const char *filename, exposure *exp, char *name = "file") {
#ifdef USE_MAGICK
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

		if (bayer == IMAGE_BAYER_NONE)
			im = new image_ale_real(mi->rows, mi->columns, 3, name, exp);
		else
			im = new image_bayer_ale_real(mi->rows, mi->columns, 3, bayer, name, exp);

		for (i = 0; i < mi->rows; i++) {
			p = AcquireImagePixels(mi, 0, i, mi->columns, 1, &exception);

			if (exception.severity != UndefinedException)
				CatchException(&exception);
			if (p == NULL)
				exit(1);

			for (j = 0; j < mi->columns; j++) {
				pixel input( (ale_real) p->red / MaxRGB,
					     (ale_real) p->green / MaxRGB,
					     (ale_real) p->blue / MaxRGB );

				im->set_pixel(i, j, exp->linearize(input));

				p++;
			}
		}

		DestroyImage(mi);
		DestroyImageInfo(image_info);

		return im;
#else
		return read_ppm(filename, exp, bayer);
#endif
	}

	/*
	 * Initializer.
	 *
	 * Handle FILE_COUNT input files with names in array FILENAMES and
	 * output file OUTPUT_FILENAME.  FILENAMES should be an array of char *
	 * that is never freed.  OUTPUT_FILENAME should be a char * that is
	 * never freed.  
	 *
	 * INPUT_EXPOSURE should be an array of FILE_COUNT exposure objects
	 * that is never freed.  OUTPUT_EXPOSURE should be an exposure * that
	 * is never freed.
	 */
	static void init(unsigned int _file_count, const char **_filenames, 
			const char *_output_filename, exposure *_input_exposure,
			exposure *_output_exposure){
		assert (_file_count > 0);

		init_image();

		filenames = _filenames;
		file_count = _file_count;
		output_filename = _output_filename;
		input_exposure = _input_exposure;
		output_exposure = _output_exposure;

		images = (const image **)malloc(file_count * sizeof(image *));
		files_open = (int *)calloc(file_count, sizeof(int));

		assert (images);
		assert (files_open);

		if (!images || !files_open) {
			fprintf(stderr, "Unable to allocate memory for images.\n");
			exit(1);
		}

		fprintf(stderr, "Output file will be '%s'.\n", 
				output_filename);

	}

	static void ppm_plain() {
		ppm_type = 2;
	}

	static void ppm_raw() {
		ppm_type = 1;
	}

	static void ppm_auto() {
#ifdef USE_MAGICK
		ppm_type = 0;
#else
		fprintf(stderr, "\n\n*** Error: --auto flag not supported on this build. ***\n"
				    "*** (Hint: Rebuild with IMAGEMAGICK=1)              ***\n\n");
		exit(1);
#endif
	}

	static void set_default_bayer(unsigned int b) {
		bayer = b;
	}

	static void depth16() {
		num_bits = 16;
		mcv = 65535;
	}

	static void depth8() {
		num_bits = 8;
		mcv = 255;
	}

	static void destroy() {
		assert (file_count > 0);
		destroy_image();
	}

	static unsigned int count() {
		assert (file_count > 0);
		return file_count;
	}

	static const char *name(unsigned int image) {
		assert (image <  file_count);

		return filenames[image];
	}

	static const char *output_name() {
		assert (file_count > 0);
		return output_filename;
	}

	/*
	 * Write an image to a file
	 */
	static void write_image(const char *filename, const image *im, exposure *exp = output_exposure, int rezero = 0, int exp_scale_override = 0) {
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

		/*
		 * Set the output image depth
		 */

		if (MaxRGB < 65535 || mcv < 65535)
			mi->depth = 8;
		else
			mi->depth = 16;

		if (MaxRGB < 65535 && mcv == 65535) {
			fprintf(stderr, "\n\n*** Warning: 16 bit-per-channel file output was specified,\n");
			fprintf(stderr, "*** but ImageMagick has not been compiled with support for this.\n");
			fprintf(stderr, "*** Writing output using 8 bits per channel.\n");
		}

		/*
		 * Set compression type
		 */

		if (ppm_type == 2) {
			mi->compression = NoCompression;
			image_info->compression = NoCompression;
			strncpy(mi->magick, "PNM", MaxTextExtent);
			strncpy(image_info->magick, "PNM", MaxTextExtent);
		} else if (ppm_type == 1) {
			strncpy(mi->magick, "PNM", MaxTextExtent);
			strncpy(image_info->magick, "PNM", MaxTextExtent);
		}

		/*
		 * Automatic exposure adjustment (don't blow out highlights)
		 */
		ale_real maxval = 1;
		ale_real minval = (rezero ? im->minval() : 0);
		pixel minval_pixel(minval, minval, minval);

		if (exposure_scale || exp_scale_override) {
			ale_real new_maxval = im->maxval();

			if (new_maxval > maxval)
				maxval = new_maxval;
		}

		/*
		 * Write the image
		 */

		for (i = 0; i < mi->rows; i++) {
			p = SetImagePixels(mi, 0, i, mi->columns, 1);
			if (p == NULL)
				break;

			for (j = 0; j < mi->columns; j++) {
				pixel unlinearized(exp->unlinearize((im->get_pixel(i, j) - minval_pixel) 
							          / (maxval - minval)));

				unlinearized = unlinearized.clamp();

				p->red = (Quantum) round(unlinearized[0] * MaxRGB);
				p->green = (Quantum) round(unlinearized[1] * MaxRGB);
				p->blue = (Quantum) round(unlinearized[2] * MaxRGB);
				p++;
			}

			if (!SyncImagePixels(mi))
				break;
		}

		if (!WriteImage(image_info, mi)) {

			/*
			 * Perhaps file type was unknown?  Set to PNM by default.
			 */

			strncpy(mi->magick, "PNM", MaxTextExtent);
			strncpy(image_info->magick, "PNM", MaxTextExtent);

			if (!WriteImage(image_info, mi)) {
				fprintf(stderr, "\n\n");
				CatchException(&mi->exception);
				fprintf(stderr, "\n");
				exit(1);
			}
		}

		DestroyImage(mi);
		DestroyImageInfo(image_info);
#else
		write_ppm(filename, im, exp, mcv, ppm_type == 2, rezero, exposure_scale || exp_scale_override);
#endif
	}

	static void output(const image *i) {
		assert (file_count > 0);
		write_image(output_name(), i, output_exposure);
	}

	static void vise_write(const char *p, const char *s, const image *i) {
		static int count = 0;
		int length = strlen(p) + strlen(s) + 8;
		char *output_string = (char *) malloc(length * sizeof(char));

		snprintf(output_string, length, "%s%08d%s", p, count, s);

		write_image(output_string, i, output_exposure);

		count++;
	}

	static exposure &exp(int n) {
		return input_exposure[n];
	}

	static const exposure &const_exp(int n) {
		return input_exposure[n];
	}

	static exposure &exp() {
		return *output_exposure;
	}

	static void exp_scale() {
		exposure_scale = 1;
	}

	static void exp_noscale() {
		exposure_scale = 0;
	}

	static const exposure &const_exp() {
		return *output_exposure;
	}

	static const image *open(unsigned int n) {
		assert (n <  file_count);
		assert (!files_open[n]);

		if (latest_close_num >= 0 && n == (unsigned int) latest_close_num) {
			images[n] = latest_close;
			latest_close_num = -1;
			files_open[n] = 1;
		} else {

			image *i = read_image(filenames[n], input_exposure + n);

			images[n] = i;
			files_open[n] = 1;
		}
		return images[n];
	}

	static image *copy(unsigned int n, char *name) {
		assert (n <  file_count);

		if (files_open[n])
			return images[n]->clone(name);
		else {
			image *i = read_image(filenames[n], input_exposure + n, name);
			return i;
		}
	}

	static void close(unsigned int image) {
		assert (image <  file_count);
		assert (files_open[image]);

		if (latest_close_num >= 0)
			delete latest_close;
		latest_close = images[image];
		latest_close_num = image;

		files_open[image] = 0;
	}

};

#endif
