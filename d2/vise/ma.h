// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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
 * ma.h: A video stabilizer that uses moving averages to calculate
 * transformations.
 */

#ifndef __ma_h__
#define __ma_h__

#include "../vise.h"
#include "../image.h"
#include "../point.h"

/*
 * Stabilize using moving averages.
 *
 * For a given frame x, the moving averages are calculated over frames
 * ranging from x - r to x + r, where r is the specified RANGE size.
 */

class ma : public vise {
	unsigned int range;
public:
	ma(render *r, unsigned int range, const char *prefix, const char *suffix, 
			ale_real scale_factor) : vise(r, prefix, suffix, scale_factor) {
		r->extend_queue(range);
		this->range = range;
	}

	/*
	 * Accept an image for rendering. 
	 */

	void render_frame(unsigned int frame_number) {

		const image *im = r->get_image(frame_number);
		int replace = 0;
		int replace_ex = 0;
		const filter::scaled_filter *scf = NULL;
		unsigned int rx_count = render::get_rx_count();
		const ale_pos *rx_parameters = render::get_rx_parameters();
		int rx_show = render::is_rx_show();

		/*
		 * Determine, for single-invariant chains,  whether replacement 
		 * is occurring, and, if so, determine whether we are honoring
		 * exclusion regions.
		 */
		if (typeid(*r) == typeid(incremental)
		 && ((incremental *)r)->get_invariant()->is_last()) {
			scf = ((incremental *)r)->get_invariant()->ssfe()->get_scaled_filter();
			if (scf->is_coarse() && scf->get_filter()->support() >= 1) {
				replace = 1;
				replace_ex = ((incremental *)r)->get_invariant()->ssfe()->ex_is_honored();
			}
		}

		/*
		 * Calculate the parameters for the desired 
		 * transformation.
		 */

		point p[4] = {
			point(0, 0),
			point(0, 0),
			point(0, 0),
			point(0, 0)
		};

		ale_pos bd[BARREL_DEGREE] = {0, /* ... */};
		unsigned int bd_count = BARREL_DEGREE;
		int frame_count = 0;

		for (int f = frame_number - range; f <= (int) (frame_number + range); f++) {
			if (f < 0
			 || f >= (int) image_rw::count())
				continue;

			frame_count++;

			transformation t = align::of(f);

			for (unsigned int i = 0; i < 4; i++)
				p[i] = p[i] + t.transform_scaled(point((i == 1 || i == 2) ? t.scaled_height() : 0,
						      (i > 1)            ? t.scaled_width()  : 0));

			if (t.bd_count() < bd_count)
				bd_count = t.bd_count();

			for (unsigned int i = 0; i < bd_count; i++)
				bd[i] += t.bd_get(i);
		}

		for (unsigned int i = 0; i < 4; i++)
			p[i] = p[i] / (ale_real) frame_count;

		for (unsigned int i = 0; i < bd_count; i++)
			bd[i] /= frame_count;

		/*
		 * Generate the desired transformation.
		 */
		transformation t = align::of(frame_number);
		transformation s = t;
		unsigned int new_height = (unsigned int) 
			floor(s.unscaled_height() * scale_factor);
		unsigned int new_width = (unsigned int) 
			floor(s.unscaled_width() * scale_factor);
		s.set_domain(new_height, new_width);
		s.gpt_set(p);
		s.bd_set(bd_count, bd);

		image *rendered = new image_ale_real(new_height, new_width, 3);

		if (replace) {
			const image *replace_image = NULL;
			replace_image = image_rw::open(frame_number);
			scf->set_parameters(t, s, replace_image);
		}

		for (unsigned int i = 0; i < rendered->height(); i++)
		for (unsigned int j = 0; j < rendered->width();  j++) {
			point unoffset_p = s.transform_scaled(point(i, j));
			point p = unoffset_p - im->offset();
			// point p_replace = t.inverse_transform(s(point(i, j)));
			double shading = 1;


			if (rx_show || (replace && replace_ex))
			for (unsigned int param = 0; param < rx_count; param++)
				if (unoffset_p[0] >= rx_parameters[6 * param + 0]
				 && unoffset_p[0] <= rx_parameters[6 * param + 1]
				 && unoffset_p[1] >= rx_parameters[6 * param + 2]
				 && unoffset_p[1] <= rx_parameters[6 * param + 3]
				 && frame_number >= (unsigned) rx_parameters[6 * param + 4]
				 && frame_number <= (unsigned) rx_parameters[6 * param + 5])
					shading *= 0.5;

			if (shading < 1 && !rx_show && replace && replace_ex)
				continue;

			if (replace) {
				pixel value, weight;

				scf->filtered(i, j, &value, &weight);

				if (weight.min_norm() > 0) {
					rendered->set_pixel(i, j, shading * value);
					continue;
				}
			}


			if (p[0] >= 0 
			 && p[0] <= im->height() - 1
			 && p[1] >= 0
			 && p[1] <= im->width() - 1)
				rendered->set_pixel(i, j, shading * im->get_bl(p));
			else
				rendered->set_pixel(i, j, pixel(0, 0, 0));
		}

		if (replace)
			image_rw::close(frame_number);

		write_frame(rendered, frame_number);

		delete rendered;
	}

	/* 
	 * Report the frame lag for this stabilizer.
	 */

	unsigned int lag() {
		return range;
	}
};

#endif
