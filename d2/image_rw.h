// Copyright 2002 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
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
#include "image_accel.h"
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
	 * Nearest-neighbor defined value radius.
	 */
	static double nn_defined_radius;

	/*
	 * Input and output exposure models
	 */
	static exposure **input_exposure;
	static exposure *output_exposure;
	static int exposure_scale;

	/*
	 * Default bayer pattern
	 */
	static unsigned int bayer_default;

	/*
	 * Image-specific bayer patterns.
	 */
	static unsigned int *bayer_specific;

	/*
	 * Pointer to the output filename
	 */
	static const char *output_filename;
	
	/*
	 * Variables relating to input image files and image data structures.
	 */
	static const char **filenames;
	static unsigned int file_count;
	static ale_image *images;
	static int *files_open;

	/*
	 * The most recently closed image number.
	 */
	static int latest_close_num;

	/*
	 * Actual cache size.
	 */
	static double cache_size;

	/*
	 * Number of cached files.
	 */
	static unsigned int cache_count;

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

	static ale_image read_image_im(const char *filename, exposure *exp, const char *name, 
			unsigned int bayer, int init_reference_gain) {
		static int warned = 0;
#ifdef USE_MAGICK

                if (MaxRGB < 65535 && mcv == 65535 && !warned) {
                        fprintf(stderr, "\n\n*** Warning: " MagickPackageName " has not been compiled with 16 bit support.\n");
                        fprintf(stderr, "*** Reading input using 8 bits per channel.\n");
                        fprintf(stderr, "*** \n"); 
                        fprintf(stderr, "*** (To silence this warning, specify option --8bpc)\n\n\n");

			warned = 1;
                }

		/*
		 * Patterned after http://www.imagemagick.org/www/api.html
		 * and http://www.imagemagick.org/www/smile.c
		 */

		ExceptionInfo exception;
		Image *mi;
		ImageInfo *image_info;
		ale_image im;
		FILE *converted_f = tmpfile();
		const PixelPacket *p;

		unsigned int i, j;

		ale_real black_level = exp->get_black_level();

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

		im = ale_new_image(accel::context(),
			(bayer == IMAGE_BAYER_NONE) ? ALE_IMAGE_RGB : ALE_IMAGE_Y,
			(image_info->depth <= 8) ? ALE_TYPE_UINT_8 :
			((image_info->depth <= 16) ? ALE_TYPE_UINT_16 :
			((image_info->depth <= 32) ? ALE_TYPE_UINT_32 : ALE_TYPE_UINT_64)));

		for (i = 0; i < mi->rows; i++) {
			p = AcquireImagePixels(mi, 0, i, mi->columns, 1, &exception);

			if (exception.severity != UndefinedException)
				CatchException(&exception);
			if (p == NULL)
				exit(1);

			for (j = 0; j < mi->columns; j++) {
				
				Quantum ival[3] = { p->red, p->green, p->blue };

				for (int k = 0; k < 3; k++) {

					if (!ale_has_channel(i, j, k, bayer))
						continue;

					if (image_info->depth <= 8) {
						cl_uchar c = ival[k];
						fwrite(&c, sizeof(cl_uchar), 1, converted_f);
					} else if (image_info->depth <= 16) {
						cl_ushort c = ival[k];
						fwrite(&c, sizeof(cl_ushort), 1, converted_f);
					} else if (image_info->depth <= 32) {
						cl_uint c = ival[k];
						fwrite(&c, sizeof(cl_uint), 1, converted_f);
					} else if (image_info->depth <= 64) {
						cl_ulong c = ival[k];
						fwrite(&c, sizeof(cl_ulong), 1, converted_f);
					} else {
						fprintf(stderr, "error: unable to handle image depth %u\n", (unsigned int) image_info->depth);
						exit(1);
					}
				}

				p++;
			}
		}

		ale_image_set_file_static(im, mi->columns, mi->rows, converted_f, 0, ppm_void_file_close, converted_f);

		DestroyImage(mi);
		DestroyImageInfo(image_info);

		return im;
#else
		return NULL;
#endif
	}

public:

	/*
	 * Read an image from a file
	 */
	static ale_image read_image(const char *filename, exposure *exp, const char *name = "file", 
			unsigned int bayer = IMAGE_BAYER_DEFAULT, int init_reference_gain = 0) {
		ale_image result;

		if (bayer == IMAGE_BAYER_DEFAULT)
			bayer = bayer_default;

		if (is_eppm(filename)) {
			result = read_ppm(filename, exp, bayer, init_reference_gain);
		}

#ifdef USE_MAGICK
		result = read_image_im(filename, exp, name, bayer, init_reference_gain);
#else
		result = read_ppm(filename, exp, bayer);
#endif

		return result;
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
			const char *_output_filename, exposure **_input_exposure,
			exposure *_output_exposure){
		assert (_file_count > 0);

		init_image();

		filenames = _filenames;
		file_count = _file_count;
		output_filename = _output_filename;
		input_exposure = _input_exposure;
		output_exposure = _output_exposure;

		images = (ale_image *)malloc(file_count * sizeof(ale_image));
		bayer_specific = (unsigned int *)malloc(file_count * sizeof(unsigned int));
		files_open = (int *)calloc(file_count, sizeof(int));

		assert (images);
		assert (bayer_specific);
		assert (files_open);

		if (!images || !files_open || !bayer_specific) {
			fprintf(stderr, "Unable to allocate memory for images.\n");
			exit(1);
		}

		for (unsigned int i = 0; i < file_count; i++)
			bayer_specific[i] = IMAGE_BAYER_DEFAULT;

		ui::get()->identify_output(output_filename);
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
		bayer_default = b;
	}

	static void set_specific_bayer(unsigned int index, unsigned int b) {
		assert (bayer_specific);
		bayer_specific[index] = b;
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

	static void def_nn(double _nn) {
		nn_defined_radius = _nn;
	}

	static const char *output_name() {
		assert (file_count > 0);
		return output_filename;
	}

	/*
	 * Write an image to a file
	 */
	static void write_image(const char *filename, ale_image im, double gamma = 0.45, int rezero = 0, int exp_scale_override = 0, int use_weights = 0) {
		static int warned = 0;

		/*
		 * Handle ALE-specific magical filenames.
		 */

		if (!strcmp(filename, "dump:")) {
			FILE *image_data = ale_image_retain_file(im);

			int format = ale_image_get_format(im);
			int type = ale_image_get_type(im);

			fprintf(stderr, "Image dump: ");
			for (unsigned int i = 0; i < ale_image_get_height(im); i++)
			for (unsigned int j = 0; j < ale_image_get_width(im); j++) {
				fprintf(stderr, "(%d, %d): ", i, j);

				fprintf(stderr, "[");

				for (unsigned int k = 0; k < ale_image_get_depth(im); k++) {

					if (k != 0)
						fprintf(stderr, " ");

					switch (type) {
					ALE_TYPE_UINT_8: {
						cl_uchar data;
						assert(fread(&data, sizeof(cl_uchar), 1, image_data));
						fprintf(stderr, "%u", (unsigned int) data);
						break; }
					ALE_TYPE_UINT_16: {
						cl_ushort data;
						assert(fread(&data, sizeof(cl_ushort), 1, image_data));
						fprintf(stderr, "%u", (unsigned int) data);
						break; }
					ALE_TYPE_UINT_32: {
						cl_uint data;
						assert(fread(&data, sizeof(cl_uint), 1, image_data));
						fprintf(stderr, "%u", (unsigned int) data);
						break; }
					ALE_TYPE_UINT_64: {
						cl_ulong data;
						assert(fread(&data, sizeof(cl_ulong), 1, image_data));
						fprintf(stderr, "%llu", data);
						break; }
					ALE_TYPE_FLOAT_32: {
						cl_float data;
						assert(fread(&data, sizeof(cl_float), 1, image_data));
						fprintf(stderr, "%f", (double) data);
						break; }
					ALE_TYPE_FLOAT_64: {
						double data;
						assert(fread(&data, sizeof(double), 1, image_data));
						fprintf(stderr, "%f", data);
						break; }
					}
				}
				fprintf(stderr, "] ");
			}
			fprintf(stderr, "\n");

			ale_image_release_file(im, image_data);

			return;
		}

#ifdef USE_MAGICK
		/*
		 * If necessary, adjust mcv to match library limits.
		 */

		if (MaxRGB < 65535 && mcv == 65535) {
			fprintf(stderr, "\n\n*** Warning: " MagickPackageName " has not been compiled with 16 bit support.\n");
			fprintf(stderr, "*** Writing output using 8 bits per channel.\n");
                        fprintf(stderr, "*** \n"); 
                        fprintf(stderr, "*** (To silence this warning, specify option --8bpc)\n\n\n");

			mcv = 255;
		}

		/*
		 * We currently can't handle mcv greater than 65535 here.
		 */

		assert(mcv <= 65535);

		if (mcv > 65535) {
			fprintf(stderr, "error: I don't know how to produce greater than 16-bit output.\n");
			exit(1);
		} 
#endif

		/*
		 * Automatic exposure adjustment (don't blow out highlights)
		 */
		ale_real maxval = 1;
		ale_real minval = (rezero ? (ale_real) ale_image_minval(im) : (ale_real) 0);
		if (minval > 0)
			minval = 0;
		pixel minval_pixel(minval, minval, minval);


		if (exposure_scale || exp_scale_override) {
			ale_real new_maxval = ale_image_maxval(im);

			if (new_maxval > maxval)
				maxval = new_maxval;
		}

		/*
		 * Nearest-neighbor fill.
		 */

		 ale_image temp_image = use_weights ? ale_image_nn_fill(im, nn_defined_radius)
		 				    : ale_image_get_weights(im);

		/*
		 * Unlinearize and quantize
		 */

		ale_image quantized_image = ale_new_image(accel::context(),
				ale_image_get_format(temp_image), (mcv > 255) ?
				ALE_TYPE_UINT_16 : ALE_TYPE_UINT_8);

		ale_eval("MAP_PIXEL(%0I, p, CLAMP(pow((GET_PIXEL(%1i, p) - (PIXEL(1, 1, 1) * %2f)) / (%3f - %2f), %4f)) * %5f)",
			quantized_image, temp_image, minval, maxval, gamma, (double) ((mcv > 255) ? 65535 : 255));

		ale_image_release(temp_image);

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

		mi->columns = ale_image_get_width(im);
		mi->rows = ale_image_get_height(im);

		/*
		 * Set the output image depth
		 */

		if (mcv < 65535)
			mi->depth = 8;
		else
			mi->depth = 16;

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
		 * Write the image
		 */

		FILE *image_data = ale_image_retain_file(quantized_image);

		for (i = 0; i < mi->rows; i++) {
			p = SetImagePixels(mi, 0, i, mi->columns, 1);
			if (p == NULL)
				break;

			for (j = 0; j < mi->columns; j++) {

				cl_ushort val[3];

				for (int k = 0; k < 3; k++) {
					if (mcv <= 255) {
						val[k] = fgetc(image_data);
					} else {
						assert(fread(&(val[k]), sizeof(cl_ushort), 1, image_data));
					}
				}

				p->red =   (Quantum) val[0];
				p->green = (Quantum) val[1];
				p->blue =  (Quantum) val[2];
				p++;
			}

			if (!SyncImagePixels(mi))
				break;
		}

		ale_image_release_file(quantized_image, image_data);

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
		write_ppm(filename, quantized_image, mcv, ppm_type == 2);
#endif

		ale_image_release(quantized_image);

	
	}

	static void output(ale_image i) {
		assert (file_count > 0);
		write_image(output_name(), i, output_exposure->get_gamma());
	}

	static void vise_write(const char *p, const char *s, ale_image i) {
		static int count = 0;
		int length = strlen(p) + strlen(s) + 8;
		char *output_string = (char *) malloc(length * sizeof(char));

		snprintf(output_string, length, "%s%08d%s", p, count, s);

		write_image(output_string, i, output_exposure->get_gamma());

		count++;
	}

	static exposure &exp(int n) {
		return *input_exposure[n];
	}

	static const exposure &const_exp(int n) {
		return *input_exposure[n];
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

	static const unsigned int bayer(unsigned int n) {
		if (bayer_specific[n] == IMAGE_BAYER_DEFAULT)
			return bayer_default;
		else
			return bayer_specific[n];
	}

	static const ale_image open_simple(unsigned int n) {
		assert (n <  file_count);

		return read_image(filenames[n], input_exposure[n], "file", bayer(n), (n == 0));
	}

	static const ale_image open(unsigned int n) {
		assert (n <  file_count);
		assert (!files_open[n]);

		files_open[n] = 1;

		if (latest_close_num >= 0 && n == (unsigned int) latest_close_num) {
			latest_close_num = -1;
			return images[n];
		} 
		
		if (n < cache_count)
			return images[n];

		ui::get()->loading_file();
		ale_image i = read_image(filenames[n], input_exposure[n], "file", bayer(n), (n == 0));

		images[n] = i;

		return images[n];
	}

	static void open_all() {
		for (unsigned int n = 0; n < file_count; n++) 
			open(n);
	}

	static ale_image get_open(unsigned int n) {
		assert (files_open[n]);
		return images[n];
	}

#if 0
	/*
	 * XXX: as far as I can tell, 'copy' is currently always used
	 * as a kind of 'retain', so that calls to this can probably be
	 * replaced by calls to an appropriate retain method.
	 */

	static ale_image copy(unsigned int n, const char *name) {
		assert (n <  file_count);

		if (files_open[n])
			return images[n]->clone(name);
		else {
			image *i = read_image(filenames[n], input_exposure[n], name, bayer(n), (n == 0));
			return i;
		}
	}
#endif

#if 0
	/*
	 * This method should be replaced by ale_sequence_release_image.
	 */

	static void close(unsigned int image) {
		assert (image <  file_count);
		assert (files_open[image]);

		files_open[image] = 0;

		if (image < cache_count)
			return;

		if (image == cache_count) {
			double image_size = ((double) images[image]->storage_size()) / pow(2, 20);

			if (image_size + cache_size < cache_size_max) {
				cache_size += image_size;
				cache_count++;
				ui::get()->cache(cache_size, cache_size_max);
				return;
			} else {
				ui::get()->cache_status(0);
			}
		}

		if (latest_close_num >= 0)
			delete images[latest_close_num];

		latest_close_num = image;
	}

	static void close_all() {
		for (unsigned int n = 0; n < file_count; n++) 
			close(n);
	}
#endif

};

#endif
