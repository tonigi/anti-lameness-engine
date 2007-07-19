// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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
 * vise.h: A superclass for all video stabilization engine classes.
 */

#ifndef __vise_h__
#define __vise_h__

#include "transformation.h"
#include "image.h"
#include "point.h"

/*
 * Class vise accepts images from vise_core.  It is assumed, when an image is
 * received, that all transformations required for stabilization are available.
 * This class is abstract, and must be subclassed to be instantiated.
 */

class vise {
	const char *prefix;
	const char *suffix;
protected:
	ale_real scale_factor;

	render *r;

	/*
	 * Write a frame.  This function does not support sequences having
	 * more than 100,000,000 frames (starting count with frame zero).
	 */
	void write_frame(const image *im, unsigned int frame_number) {
		if (frame_number > 99999999) {
			fprintf(stderr, "\n\n *** Frame count too high for d2::vise::write_frame() ***\n\n\n");
			exit(1);
		}

		int length = strlen(prefix) + strlen(suffix) + 8 + 1;
		char *filename_string = (char *) malloc(length * sizeof(char));

		snprintf(filename_string, length, "%s%08d%s", prefix, frame_number, suffix);

		image_rw::write_image(filename_string, im, &image_rw::exp());

		free(filename_string);
	}
public:
	vise(render *r, const char *prefix, const char *suffix, ale_real scale_factor) {
		this->prefix = prefix;
		this->suffix = suffix;
		this->r = r;
		this->scale_factor = scale_factor;
	}

	/*
	 * Accept an image for rendering. 
	 */

	virtual void render_frame(unsigned int frame_number) = 0;

	/* 
	 * Report the frame lag for this stabilizer.
	 */

	virtual unsigned int lag() = 0;

	virtual ~vise() {
	}

};

#endif
