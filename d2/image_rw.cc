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

#include "image_rw.h"

/*
 * See image_rw.h for details on these variables.
 */

int image_rw::ppm_type = 0;

unsigned int image_rw::num_bits = 8;
unsigned int image_rw::mcv = 255;

unsigned int image_rw::file_count = 0;
const char *image_rw::output_filename = NULL;
const char **image_rw::filenames = NULL;
const image **image_rw::images = NULL;
int *image_rw::files_open;

int image_rw::latest_close_num = -1;
const image *image_rw::latest_close = NULL;

double image_rw::nn_defined_radius = 0;

exposure **image_rw::input_exposure = NULL;
exposure *image_rw::output_exposure = NULL;
unsigned int image_rw::bayer_default = 0;
unsigned int *image_rw::bayer_specific = NULL;
int image_rw::exposure_scale = 1;
