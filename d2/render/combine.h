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
 * combine.h: A renderer that combines two renderings.
 */

#ifndef __combine_h__
#define __combine_h__

#include "../transformation.h"
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
	image *defined_image;
public:

	/*
	 * Constructor
	 */
	combine(render *_default, render *partial) {
		this->_default = _default;
		this->partial = partial;
		this->output_image = NULL;
		this->defined_image = NULL;
	}

	virtual ~combine() {
		if (output_image)
			delete output_image;
		if (defined_image)
			delete defined_image;
	}

	/*
	 * Result of rendering.
	 */

	virtual const image *get_image() {
		const image *default_image = _default->get_image();
		const image *partial_image = partial->get_image();

		assert (default_image->width()  == partial_image->width());
		assert (default_image->height() == partial_image->height());

		if (output_image)
			return output_image;

		output_image = new image_ale_real(default_image->height(),
				default_image->width(), 3, NULL);

		output_image->set_offset(default_image->offset());

		const image *partial_weight = partial->get_defined();

		for (unsigned int i = 0; i < default_image->height(); i++)
		for (unsigned int j = 0; j < default_image->width();  j++)
			output_image->set_pixel(i, j, 
				(partial_weight->get_pixel(i, j).min_norm() >= render::get_wt())
				? partial_image->get_pixel(i, j)
				: default_image->get_pixel(i, j));

		return output_image;
	}

	/*
	 * Definition map.  Unit-depth image whose pixels are nonzero where
	 * the image is defined.
	 */

	virtual const image *get_defined() {
		unsigned int i, j, k;

		if (defined_image)
			return defined_image;

		const image *partial_weight = partial->get_defined();
		const image *default_weight = _default->get_defined();

		assert (default_weight->width()  == partial_weight->width());
		assert (default_weight->height() == partial_weight->height());
		
		defined_image = new image_ale_real(default_weight->height(),
				default_weight->width(), 3, NULL);

		defined_image->set_offset(default_weight->offset());

		for (i = 0; i < default_weight->height(); i++)
		for (j = 0; j < default_weight->width();  j++)
		for (k = 0; k < default_weight->depth();  k++)
			defined_image->set_pixel(i, j, 
					(partial_weight->get_pixel(i, j).min_norm() >= render::get_wt())
					? partial_weight->get_pixel(i, j)
					: default_weight->get_pixel(i, j));

		return defined_image;
	}

	/*
	 * Perform rendering steps requiring no frames beyond frame N.
	 */

	virtual void sync(int n) {
		render::sync(n);
		if (output_image) {
			delete output_image;
			output_image = NULL;
		}
		if (defined_image) {
			delete defined_image;
			defined_image = NULL;
		}
		_default->sync(n);
		partial->sync(n);
	}

	virtual void step() {
	}

	virtual void init_point_renderer(unsigned int h, unsigned int w, unsigned int d) {
		_default->init_point_renderer(h, w, d);
		partial->init_point_renderer(h, w, d);
		output_image = new image_zero(h, w, d);
		defined_image = new image_zero(h, w, d);
	}

	virtual void point_render(unsigned int i, unsigned int j, unsigned int f, transformation t) {
		_default->point_render(i, j, f, t);
		partial->point_render(i, j, f, t);
	}

	virtual void finish_point_rendering() {
		_default->finish_point_rendering();
		partial->finish_point_rendering();
		delete defined_image;
		delete output_image;

		/*
		 * These will be generated upon a call to get_image() or
		 * get_defined().
		 */

		defined_image = NULL;
		output_image = NULL;
	}

	const render *get_default() const {
		return _default;
	}

	const render *get_partial() const {
		return partial;
	}

	void free_memory() {
		delete output_image;
		delete defined_image;
	}
};

#endif
