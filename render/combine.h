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
 * combine.h: A renderer that combines two renderings.
 */

#ifndef __combine_h__
#define __combine_h__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../gpt.h"
#include "../image.h"
#include "../point.h"

/*
 * Combine two renderings. 
 *
 * Available data is taken from the PARTIAL rendering.  When no data from
 * the PARTIAL rendering is available, data from the DEFAULT rendering
 * is substituted.
 */

class combine : public render {
private:
	render *_default;
	render *partial;
	image *output_image;
public:

	/*
	 * Constructor
	 */
	combine(render *_default, render *partial) {
		this->_default = _default;
		this->partial = partial;
		this->output_image = NULL;
	}

	virtual ~combine() {
		if (output_image)
			delete output_image;
	}

	/*
	 * Result of rendering.
	 */

	virtual const image *get_image() {
		unsigned int i, j, k;

		const image *default_image = _default->get_image();
		const image *partial_image = partial->get_image();

		assert (default_image->width()  == partial_image->width());
		assert (default_image->height() == partial_image->height());

		if (output_image)
			return output_image;

		output_image = new image(default_image->height(),
				default_image->width(), 3);

		const image_weights *partial_weight = partial->get_defined();

		for (i = 0; i < default_image->height(); i++)
		for (j = 0; j < default_image->width();  j++)
		for (k = 0; k < 3; k++) {
			output_image->set_pixel_component(i, j, k, 
				(partial_weight->get_pixel_component(i, j, 0) != 0)
				? partial_image->get_pixel_component(i, j, k)
				: default_image->  get_pixel_component(i, j, k));
		}	
		return output_image;
	}

	/*
	 * Definition map.  Unit-depth image whose pixels are nonzero where
	 * the image is defined.
	 */

	virtual const image_weights *get_defined() {
		assert (0);
		return NULL;
	}

	/*
	 * Perform rendering steps requiring no frames beyond frame N.
	 */

	virtual void sync(int n) {
		if (output_image) {
			delete output_image;
			output_image = NULL;
		}
		_default->sync(n);
		partial->sync(n);
	}

};

#endif
