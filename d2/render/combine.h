// Copyright 2002, 2007 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

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
 * combine.h: A renderer that combines two renderings.
 */

#ifndef __combine_h__
#define __combine_h__

#include "../transformation.h"
#include "../image.h"
#include "../point.h"
#include "incremental.h"
#include "../filter/filter.h"

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
	mutable image *output_image;
	mutable image *defined_image;
	int synced;

	class refilter : public thread::decompose_domain {
		combine *c;
		const render *fine;
		const render *coarse;
		const filter::filter *f;
		const image *fine_weight;
		const image *fine_image;
		const image *coarse_defined;
		image *output_image;

	protected:
		void subdomain_algorithm(unsigned int thread,
				int i_min, int i_max, int j_min, int j_max) {

			for (int i = i_min; i < i_max; i++)
			for (int j = j_min; j < j_max; j++) 
			for (unsigned int k = 0; k < 3; k++){

				if (!(coarse_defined->chan(i, j, k) > 0))
					continue;

				ale_real filter_scale = 1;
				ale_real filtered_weight;
				ale_real filtered_value;

				do {

					filtered_weight = 0;
					filtered_value = 0;

					/* 
					 * lrintf() may be faster than ceil/floor() on some architectures.
					 * See render/psf/raster.h for more details.
					 */

					for (int ii = (int) lrintf(-f->support() * filter_scale); 
						 ii < (int) lrintf( f->support() * filter_scale); ii++)
					for (int jj = (int) lrintf(-f->support() * filter_scale);
						 jj < (int) lrintf( f->support() * filter_scale); jj++) {

						if (ii + i < 0
						 || jj + j < 0
						 || ii + i >= (int) fine_weight->height()
						 || jj + j >= (int) fine_weight->width())
							continue;

						ale_real pw = fine_weight->get_pixel(i + ii, j + jj)[k];

						if (!(pw > 0))
							continue;

						ale_real w = pw * f->response(point(ii / filter_scale, 
										    jj / filter_scale));
					
						ale_real v = fine_image->get_pixel(i + ii, j + jj)[k];
						
						if (!finite(w) || !finite(v))
							continue;

						filtered_weight += w;
						filtered_value  += w * v;
					}

					if (filtered_weight < render::get_wt())
						/* filter_scale += 1; */
						filter_scale *= 2;

				} while (filtered_weight < render::get_wt()
				    && filter_scale < coarse_defined->width() 
						    + coarse_defined->height());

				output_image->chan(i, j, k) = filtered_value / filtered_weight;
			}
		}
	public:
		refilter(combine *_c,
		         const render *_fine,
			 const render *_coarse,
			 const filter::filter *_f,
			 const image *_fine_weight,
			 const image *_fine_image,
			 const image *_coarse_defined,
			 image *_output_image) : decompose_domain(0, _coarse_defined->height(),
			                                          0, _coarse_defined->width()) {

			c = _c;
			fine = _fine;
			coarse = _coarse;
			f = _f;
			fine_weight = _fine_weight;
			fine_image = _fine_image;
			coarse_defined = _coarse_defined;
			output_image = _output_image;
		}
	};

	const image *get_image_dynamic() const {
		assert(typeid(*partial) == typeid(incremental));

		if (typeid(*_default) != typeid(combine) || !synced) {
			/*
			 * Degenerate case.
			 */
			output_image = _default->get_image()->clone("degenerate dynamic filter");
			return output_image;
		}

		combine *c = (combine *)_default;
		const render *fine = c->get_partial();
		const render *coarse = c->get_default();
		const filter::filter *f = ((incremental *)partial)->get_invariant()->ssfe()->
				get_scaled_filter()->get_filter();
		const image *fine_weight = fine->get_defined();
		const image *fine_image = fine->get_image();
		const image *coarse_defined = coarse->get_defined();

		output_image = new image_ale_real(coarse_defined->height(),
				coarse_defined->width(), 3, NULL);

		output_image->set_offset(coarse_defined->offset());

		assert (coarse_defined->width()  == fine_image->width());
		assert (coarse_defined->height() == fine_image->height());
		assert (coarse_defined->width()  == fine_weight->width());
		assert (coarse_defined->height() == fine_weight->height());

		ui::get()->refilter_start();

		refilter r(c, fine, coarse, f, fine_weight, fine_image, 
		           coarse_defined, output_image);
		r.run();
		
		ui::get()->refilter_done();

		return output_image;
	}
public:

	/*
	 * Constructor
	 */
	combine(render *_default, render *partial) {
		this->_default = _default;
		this->partial = partial;
		this->output_image = NULL;
		this->defined_image = NULL;
		this->synced = 0;
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

	virtual const image *get_image() const {

		if (output_image)
			return output_image;

		assert(typeid(*partial) != typeid(combine));

		/*
		 * Dynamic filtering is handled separately.
		 */
		if (typeid(*partial) == typeid(incremental)
		 && (((incremental *)partial)->get_invariant()->
		 		ssfe()->get_scaled_filter()->is_dynamic()))
			return get_image_dynamic();

		const image *default_image = _default->get_image();

		output_image = new image_ale_real(default_image->height(),
				default_image->width(), 3, NULL);

		output_image->set_offset(default_image->offset());

		const image *partial_image = partial->get_image();
		const image *partial_weight = partial->get_defined();

		assert (default_image->width()  == partial_image->width());
		assert (default_image->height() == partial_image->height());

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

	virtual const image *get_defined() const {
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

	virtual int sync() {
		if (output_image) {
			delete output_image;
			output_image = NULL;
		}
		if (defined_image) {
			delete defined_image;
			defined_image = NULL;
		}
		_default->sync();
		partial->sync();
		synced = 1;

		return 1;
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
		output_image = NULL;
		defined_image = NULL;
	}
};

#endif
