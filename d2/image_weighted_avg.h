// Copyright 2002, 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                          <dhilvert@ugcs.caltech.edu>

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
 * image_weighted_avg.h: Image representing a weighted average of inputs.
 */

#ifndef __image_weighted_avg_h__
#define __image_weighted_avg_h__

#include "image_ale_real.h"
#include "exposure/exposure.h"
#include "point.h"
#include "image.h"

#define ALE_GLSL_IMAGE_WEIGHTED_AVG_INCLUDE \
ALE_GLSL_IMAGE_INCLUDE \
"struct image_weighted_avg {\n"\
"	int type;\n"\
"	vec2 offset;\n"\
"};\n"\
"bool image_weighted_avg_accumulate_norender(inout image_weighted_avg _this, vec4 pos);\n"\
"void image_weighted_avg_accumulate(inout image_weighted_avg _this, vec4 pos, int frame, vec3 value, vec3 confidence);\n"\
"image image_weighted_avg_get_weights(inout image_weighted_avg _this);\n"\
"vec3 image_weighted_avg_get_pixel(inout image_weighted_avg _this, vec4 pos);\n"

class image_weighted_avg : public image, public gpu::program::library {
private:
	void trigger(pixel multiplier) {
		assert(0);
	}

public:
	image_weighted_avg (unsigned int dimy, unsigned int dimx, unsigned int
			depth, const char *name = "anonymous") 
			: image(dimy, dimx, depth, name, NULL) {
	}

	virtual ~image_weighted_avg() {
	}

	void set_pixel(unsigned int y, unsigned int x, spixel p) {
		assert(0);
	}

	spixel get_pixel(unsigned int y, unsigned int x) const {
		assert(0);

		return spixel(0, 0, 0);
	}

	void set_chan(unsigned int y, unsigned int x, unsigned int k, ale_sreal c) {
		assert(0);
	}

	ale_sreal get_chan(unsigned int y, unsigned int x, unsigned int k) const {
		assert(0);

		return 0;
	}

	/*
	 * Make a new image suitable for receiving scaled values.
	 */
	virtual image *scale_generator(int height, int width, int depth, const char *name) const {
		return new_image_ale_real(height, width, depth, name, _exp);
	}

	/*
	 * Pre-transformation check for whether an area should be skipped.
	 * Takes image weights as an argument.
	 */
	virtual int accumulate_norender(int i, int j) = 0;

	/*
	 * Accumulate pixels 
	 */
	virtual void accumulate(int i, int j, int f, pixel new_value, pixel new_weight) = 0;

	virtual void accumulate_accel(const gpu::program *p) {
		assert(0);
	}

	/*
	 * Get color map
	 */
	virtual image *get_colors() = 0;

	/*
	 * Get weight map
	 */
	virtual image *get_weights() = 0;
};

#endif
