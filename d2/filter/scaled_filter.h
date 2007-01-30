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

#ifndef __scaled_filter_h__
#define __scaled_filter_h__

#include "filter.h"
#include "mult.h"
#include "sinc.h"
#include "triangle.h"
#include "box.h"
#include "gauss.h"
#include "lanczos.h"

/*
 * Useful constants.
 */

static const double sqrt2          = 1.41421356237309504880;
static const double one_over_sqrt2 = 0.70710678118654752440;

/*
 * Scaled filter class.
 */

class scaled_filter {
private:
	/*
	 * Frequency limit:
	 *
	 *    0: indicates the higher limit (identical to output frequency)
	 *
	 *    1: indicates the safer limit (minimum of input and output
	 *       frequencies)
	 */
	int frequency_limit;
	filter *f;

	/*
	 * filter function parameters
	 *
	 * Parameters include data regarding the current image and
	 * transformation.  All parameters are mutable.
	 */

	/*
	 * A pointer to the current image is stored, as well as the image
	 * offset and bayer type.
	 */

	mutable const image *im;
	mutable unsigned int bayer;
	mutable point offset;

	/* We are either using one transformation or a combination of two
	 * transformations.  T_TWO indicates whether two transformations are
	 * being used.  T[] stores the transformations.  When only one
	 * transformation is used, OFFSET stores the image offset.
	 */

	mutable unsigned int t_two;
	mutable transformation t[2];
	mutable int _is_projective;

	/*
	 * Transform a point using the current transformation.
	 */
	point transform(point p) const {
		if (t_two)
			return t[1].unscaled_inverse_transform(t[0].transform_unscaled(p));

		return t[0].transform_unscaled(p);
	}

	/*
	 * Inverse of the above.
	 */
	point transform_inverse(point p) const {
		if (t_two)
			return t[0].unscaled_inverse_transform(t[1].transform_unscaled(p));

		return t[0].unscaled_inverse_transform(p);
	}

	/*
	 * Returns non-zero if the transformation might be non-Euclidean.
	 */
	int is_projective() const {
		return _is_projective;
	}

	/*
	 * If we are limiting to source frequencies, then scale a filter to
	 * accept only frequencies we know to be expressible in the source.
	 * (Or do this approximately.)
	 */
	void freq_limit(point p, point mapped_p, ale_pos *hscale_g, 
			ale_pos *hscale_rb, ale_pos *wscale_g, ale_pos *wscale_rb) const {

		if (frequency_limit == 0)
			return;

		ale_real hnorm, wnorm;

		point dh = transform_inverse(p + point(1, 0));
		point dw = transform_inverse(p + point(0, 1));

		hnorm = (mapped_p - dh).norm();
		wnorm = (mapped_p - dw).norm();

		if (bayer == IMAGE_BAYER_NONE) {
			if (hnorm < 1) {
				*hscale_g = hnorm;
				*hscale_rb = hnorm;
			}
			if (wnorm < 1) {
				*wscale_g = wnorm;
				*wscale_rb = wnorm;
			}
		} else {
			if (hnorm < sqrt2) {
				*hscale_g = hnorm / sqrt2;
				*hscale_rb = hnorm / 2;
			} else if (hnorm < 2) {
				*hscale_rb = hnorm / 2;
			}
			if (wnorm < sqrt2) {
				*wscale_g = wnorm / sqrt2;
				*wscale_rb = wnorm / 2;
			} else if (wnorm < 2) {
				*wscale_rb = wnorm / 2;
			}
		}
	}

	void filter_channel(point p, point mapped_p, unsigned int k, ale_pos hscale,
			ale_pos wscale, pixel *result, pixel *weight, int honor_exclusion, int frame) const {

#if 1
		/*
		 * This test matches the technical description.
		 */
		if (hscale < 1 && wscale < 1) {
#else
		/*
		 * This approach is ~33% faster for Euclidean transformations,
		 * but is likely to produce different results in some cases.
		 */
		if (hscale <= 1 && wscale <= 1) {
#endif

			/*
			 * Handle the especially coarse case.
			 */

			ale_pos fscale;

			if (frequency_limit) {
				fscale  = (bayer == IMAGE_BAYER_NONE)
					? 1
					: (k == 1) ? sqrt2 : 2;
			} else {
				fscale = 1;
			}
			
			ale_pos support = f->support() * fscale;

			point min = mapped_p - point(support, support);
			point max = mapped_p + point(support, support);

			if (min[0] < 0 || 1 / support == 0)
				min[0] = 0;
			if (max[0] > im->height() - 1 || 1 / support == 0)
				max[0] = im->height() - 1;
			if (min[1] < 0 || 1 / support == 0)
				min[1] = 0;
			if (max[1] > im->width() - 1 || 1 / support == 0)
				max[1] = im->width() - 1;
			/*
			 * Iterate over the source pixels.
			 */

			for (int i = (int) ceil (min[0]); 
				 i<= (int) floor(max[0]); i++)
			for (int j = (int) ceil (min[1]); 
				 j<= (int) floor(max[1]); j++) {

				if (render::is_excluded_f(i, j, frame))
					continue;

				if (bayer != IMAGE_BAYER_NONE
				 && ((image_bayer_ale_real *) im)->bayer_color(i, j) != k)
					continue;

				point a = point(i, j);

				ale_real v = im->chan(i, j, k);
				ale_real w = im->exp().confidence(k, v) * f->response((a - mapped_p) / fscale);

				(*weight)[k] += w;
				(*result)[k] += w * v;
			}
		} else {
			/*
			 * Establish the boundaries of filtering in the source.
			 */

			point min = mapped_p;
			point max = mapped_p;

			for (int ii = 0; ii < 2; ii++)
			for (int jj = 0; jj < 2; jj++) {

				point b = transform_inverse(point((f->support() / hscale) * (ii ? -1 : 1),
									     (f->support() / wscale) * (jj ? -1 : 1))
								     + p);

				for (int d = 0; d < 2; d++) {
					if (b[d] < min[d])
						min[d] = b[d];
					if (b[d] > max[d])
						max[d] = b[d];
				}

			}

			if (min[0] < 0 || 1 / f->support() == 0)
				min[0] = 0;
			if (max[0] > im->height() - 1 || 1 / f->support() == 0)
				max[0] = im->height() - 1;
			if (min[1] < 0 || 1 / f->support() == 0)
				min[1] = 0;
			if (max[1] > im->width() - 1 || 1 / f->support() == 0)
				max[1] = im->width() - 1;

			/*
			 * Iterate over the source pixels.
			 */

			for (int i = (int) ceil (min[0]); 
				 i<= (int) floor(max[0]); i++)
			for (int j = (int) ceil (min[1]); 
				 j<= (int) floor(max[1]); j++) {

				if (render::is_excluded_f(i, j, frame))
					continue;

				if (bayer != IMAGE_BAYER_NONE
				 && ((image_bayer_ale_real *) im)->bayer_color(i, j) != k)
					continue;

				point a = transform(point(i, j));

				ale_real v = im->chan(i, j, k);
				ale_real w = im->exp().confidence(k, v) * f->response((a - p) * point(hscale, wscale));

				(*weight)[k] += w;
				(*result)[k] += w * v;
			}

		}
	}

public:

	scaled_filter(filter *f, int frequency_limit) {
		this->frequency_limit = frequency_limit;
		this->f = f;
	}

	const filter *get_filter() const {
		return f;
	}

	int equals(scaled_filter *s) {
		return (frequency_limit == s->frequency_limit
		     && f->equals(s->f));
	}

	int is_coarse() const {
		return frequency_limit == 1;
	}

	int is_fine() const {
		return frequency_limit == 0;
	}

	/*
	 * Set the transformation T and image IM.
	 */
	void set_parameters(transformation _t, const image *_im, point _offset) const {
		t_two = 0;
		t[0] = _t;
		im = _im;

		bayer = im->get_bayer();
		offset = _offset;

		_is_projective = _t.is_projective();
	}

	/*
	 * Set the transformations T and S, and image IM.
	 */
	void set_parameters(transformation _t, transformation _s, const image *_im) const {
		t_two = 1;
		t[0] = _t;
		t[1] = _s;
		im = _im;

		bayer = im->get_bayer();
		offset = point(0, 0);

		_is_projective = t[0].is_projective() || t[1].is_projective();
	}

	

	/*
	 * Return filtered RESULT and WEIGHT at point P in a coordinate system
	 * specified by the inverse of transformation T based on data taken
	 * from image IM.
	 */
	void filtered(int i, int j, pixel *result, pixel *weight, int honor_exclusion, int frame) const {

		point p = point(i, j) + offset;

		*result = pixel(0, 0, 0);
		*weight = pixel(0, 0, 0);
		
		ale_pos hscale_g = 1;
		ale_pos hscale_rb = 1;
		ale_pos wscale_g = 1;
		ale_pos wscale_rb = 1;

		point mapped_p = transform_inverse(p);

		/*
		 * Allowing points such as these results in problems that are
		 * annoying to resolve.  Solving such problems is not obviously
		 * worth the effort, given that they can be avoided entirely by
		 * disallowing the points.
		 */

		if (mapped_p[0] < 0 || mapped_p[0] > im->height() - 1
		 || mapped_p[1] < 0 || mapped_p[1] > im->width() - 1)
			return;

		freq_limit(p, mapped_p, &hscale_g, &hscale_rb, &wscale_g, &wscale_rb);

		filter_channel(p, mapped_p, 0, hscale_rb, wscale_rb, result, weight, honor_exclusion, frame);
		filter_channel(p, mapped_p, 2, hscale_rb, hscale_rb, result, weight, honor_exclusion, frame);
		filter_channel(p, mapped_p, 1, hscale_g , hscale_g , result, weight, honor_exclusion, frame);

		for (unsigned int k = 0; k < 3; k++) {
			if (fabs((*weight)[k]) < 0.0001) {
				(*weight)[k] = 0;
				(*result)[k] = 0;
			} else {
				(*result)[k] /= (*weight)[k];
			}
		}
	}
};
#endif
