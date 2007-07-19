// Copyright 2002, 2003, 2004, 2005, 2006 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                                      <dhilvert@ugcs.caltech.edu>

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
 * ale.cc: The main module of ale.  This calls the input handler.
 */

/*
 * Configuration
 */

#include <config.h>

/*
 * Types
 */

#include "ale_pos.h"
#include "ale_real.h"

/*
 * Version Information
 */

char *package_name = PACKAGE_NAME;

char *short_version = VERSION;

char *version = PACKAGE_NAME " Version:      " VERSION "\n"
#ifdef USE_MAGICK
		"File handler:     ImageMagick\n"
#else
		"File handler:     PPM\n"
#endif
		"Color data:       " ALE_REAL_PRECISION_STRING "\n"
		"Coordinate data:  " ALE_POS_PRECISION_STRING "\n"
#ifdef USE_FFTW
		"DFT:              FFTW3\n"
#else
		"DFT:              Built-in\n"
#endif
#if defined USE_PTHREAD 
		"Threads:          POSIX\n"
#else
		"Threads:          Disabled\n"
#endif
#if defined NDEBUG && !defined DEBUG
                "Assertions:       Disabled\n"
#elif defined DEBUG && !defined NDEBUG
                "Assertions:       Enabled\n"
#elif defined NDEBUG
                "Assertions:       Probably disabled\n"
#else
                "Assertions:       Probably enabled\n"
#endif
#if OPTIMIZATIONS == 1
		"Optimizations:    Enabled\n"
#else
		"Optimizations:    Disabled\n"
#endif
;

/*
 * User interface includes
 */

#include "ui/ui.h"

/*
 * main() calls the input handler.
 */

int main(int argc, const char *argv[]){

	/* 
	 * Call UI routine to handle options and other interface input.
	 * Returning from this function indicates program success.
	 */

	ui::handle_input(argc, argv, package_name, short_version, version);

	return 0;

}
