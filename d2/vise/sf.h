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
 * sf.h: A video stabilizer that uses a single frame's transformation.
 */

#ifndef __sf_h__
#define __sf_h__

#include "../vise.h"
#include "../image.h"
#include "../point.h"

/*
 * Stabilize to the viewpoint of a single frame.
 */

class sf : public vise {
	unsigned int frame;
public:
	sf(render *r, unsigned int frame, const char *prefix, const char *suffix,
			ale_real scale_factor) : vise(r, prefix, suffix, scale_factor) {
		if (frame > 10) {
			fprintf(stderr, "\n\n*** Warning: large values <x> for VISE sf:<x> are not recommended ***\n\n");
		}
		r->extend_queue(frame);
		this->frame = frame;
	}

	/*
	 * Accept an image for rendering. 
	 */

	void render_frame(unsigned int frame_number) {

		const image *im = r->get_image(frame_number);
		int replace = 0;     // Are image regions being replaced?
		int replace_ex = 0;  // If image regions are being replaced, are we honoring exclusion regions?
		unsigned int rx_count = render::get_rx_count();
		const filter::scaled_filter *scf = NULL;
		const int *rx_parameters = render::get_rx_parameters();
		int rx_show = render::is_rx_show();

		/*
		 * Determine, for single-invariant chains, whether replacement
		 * is occurring, and, if so, determine whether we are honoring
		 * exclusion regions.
		 */
		if (typeid(*r) == typeid(incremental)
		 && ((incremental *)r)->get_invariant()->is_last()) {
			scf = ((incremental *)r)->get_invariant()->ssfe()->get_scaled_filter();
			assert(scf);
			if (scf->is_coarse() && scf->get_filter()->support() >= 1) {
				replace = 1;
				replace_ex = ((incremental *)r)->get_invariant()->ssfe()->ex_is_honored();
			}
		}

		/*
		 * Generate the desired transformation.
		 */
		transformation t = align::of(frame_number);
		transformation s = align::of(frame);
		unsigned int new_height = (unsigned int) 
			floor(s.unscaled_height() * scale_factor);
		unsigned int new_width = (unsigned int) 
			floor(s.unscaled_width() * scale_factor);
		s.set_domain(new_height, new_width);

		image *rendered = new image_ale_real(new_height, new_width, 3);

		const image *replace_image = NULL;
		
		if (replace) {
			replace_image = image_rw::open(frame_number);
			scf->set_parameters(t, s, replace_image);
		}

		for (unsigned int i = 0; i < rendered->height(); i++)
		for (unsigned int j = 0; j < rendered->width();  j++) {
			point unoffset_p = s.transform_scaled(point(i, j));
			point p = unoffset_p - im->offset();
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
		return frame;
	}
};

#endif
