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

// static const ale_pos sqrt2          = (ale_pos) 1.41421356237309504880;
// static const ale_pos one_over_sqrt2 = (ale_pos) 0.70710678118654752440;

static const ale_pos sqrt2          = sqrt((ale_pos) 2);
static const ale_pos one_over_sqrt2 = 1 / sqrt2;

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
	 *
	 *    2: indicates a limit dynamically derived from subsequent chain
	 *       elements.
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
	mutable transformation t0, t1;
	mutable trans_single ts0, ts1;
	mutable int _is_projective;

	/*
	 * Transform a point using the current transformation.
	 */
	point transform(point p) const {
		if (t_two)
			return ts1.unscaled_inverse_transform(ts0.transform_unscaled(p));

		return ts0.transform_unscaled(p);
	}

	/*
	 * Inverse of the above.
	 */
	point transform_inverse(point p) const {
		if (t_two)
			return ts0.unscaled_inverse_transform(ts1.transform_unscaled(p));

		return ts0.unscaled_inverse_transform(p);
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

		ale_pos hnorm, wnorm;

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
			ale_pos wscale, pixel *result, pixel *weight, int honor_exclusion, 
			int frame, ale_real prev_value = 0, ale_real prev_weight = 0) const {

		ale_real temp_result = (*result)[k], temp_weight = (*weight)[k];
		ale_real certainty;

		if (prev_weight > 0)
			certainty = im->exp().confidence(k, prev_value);
		else
			certainty = 1;   /* We calculate certainty later */

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
					? (ale_pos) 1
					: (k == 1) ? sqrt2 : (ale_pos) 2;
			} else {
				fscale = 1;
			}
			
			int min_i, max_i, min_j, max_j;

			ale_pos support = f->support() * fscale;

			if (1 / support == 0) {
				min_i = 0;
				max_i = im->height() - 1;
				min_j = 0;
				max_j = im->width() - 1;
			} else {

				point min = mapped_p - point(support, support);
				point max = mapped_p + point(support, support);

				/*
				 * lrintf() may be faster than ceil/floor() on some architectures.
				 * See render/psf/raster.h for more details.
				 */

				min_i = (int) lrintf(min[0]);
				max_i = (int) lrintf(max[0]);
				min_j = (int) lrintf(min[1]);
				max_j = (int) lrintf(max[1]);

				if (min_i < 0)
					min_i = 0;
				if (max_i > (int) im->height() - 1)
					max_i = (int) im->height() - 1;
				if (min_j < 0)
					min_j = 0;
				if (max_j > (int) im->width() - 1)
					max_j = (int) im->width() - 1;
			}

			/*
			 * Iterate over the source pixels.
			 */

			for (int i = min_i; i <= max_i; i++)
			for (int j = min_j; j <= max_j; j++) {

				if (honor_exclusion && render::is_excluded_f(i, j, frame))
					continue;

				if (bayer != IMAGE_BAYER_NONE
				 && (im->get_channels(i, j) & (1 << k)) == 0)
					continue;

				point a = point(i, j);

				ale_real v = im->get_chan(i, j, k);

				ale_real response = f->response((a - mapped_p) / fscale);

				ale_real w = certainty * response;

				temp_weight += w;
				temp_result += w * v;
			}
		} else {
			/*
			 * Establish the boundaries of filtering in the source.
			 */

			point min = mapped_p;
			point max = mapped_p;

			int imin[2];
			int imax[2];

			ale_pos sup = f->support();

			if (1 / sup == 0) {
				imin[0] = 0;
				imax[0] = im->height() - 1;
				imin[1] = 0;
				imax[1] = im->width() - 1;
			} else {

				ale_pos hsup = sup / hscale;
				ale_pos wsup = sup / wscale;

				for (int ii = -1; ii <= +1; ii+=2)
				for (int jj = -1; jj <= +1; jj+=2) {

					point b = transform_inverse(point(hsup * ii, wsup * jj) + p);
					for (int d = 0; d < 2; d++) {
						if (b[d] < min[d])
							min[d] = b[d];
						if (b[d] > max[d])
							max[d] = b[d];
					}

				}

				/*
				 * lrintf() may be faster than ceil/floor() on some architectures.
				 * See render/psf/raster.h for more details.
				 */

				imin[0] = lrintf(min[0]);
				imax[0] = lrintf(max[0]);
				imin[1] = lrintf(min[1]);
				imax[1] = lrintf(max[1]);

				if (imin[0] < 0)
					imin[0] = 0;
				if (imax[0] > (int) im->height() - 1)
					imax[0] = (int) im->height() - 1;
				if (imin[1] < 0)
					imin[1] = 0;
				if (imax[1] > (int) im->width() - 1)
					imax[1] = (int) im->width() - 1;
			}

			/*
			 * Iterate over the source pixels.
			 */

			for (int i = imin[0]; i <= imax[0]; i++)
			for (int j = imin[1]; j <= imax[1]; j++) {

				if (honor_exclusion && render::is_excluded_f(i, j, frame))
					continue;

				if (bayer != IMAGE_BAYER_NONE
				 && (im->get_channels(i, j) & (1 << k)) == 0)
					continue;

				point a = transform(point(i, j));

				ale_real v = im->get_chan(i, j, k);

				ale_real response = f->response((a - p) * point(hscale, wscale));

				ale_real w = certainty * response;

				temp_weight += w;
				temp_result += w * v;
			}
		}

		if (!(prev_weight > 0) 
		 && d2::exposure::get_confidence() != 0
		 && temp_weight > 0) {

			/*
			 * Calculate certainty for the first pass.
			 */

			certainty = im->exp().confidence(k, temp_result / temp_weight);
			temp_weight *= certainty;
			temp_result *= certainty;
		}

		(*result)[k] = temp_result;
		(*weight)[k] = temp_weight;
	}

public:

	scaled_filter(filter *f, int frequency_limit) :
			t0(transformation::eu_identity()),
			t1(transformation::eu_identity()) {

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

	int is_dynamic() const {
		return frequency_limit == 2;
	}

	/*
	 * Set the transformation T and image IM.
	 */
	void set_parameters(transformation _t, const image *_im, point _offset) const {
		t_two = 0;
		t0 = _t;
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
		t0 = _t;
		t1 = _s;
		im = _im;

		bayer = im->get_bayer();
		offset = point(0, 0);

		_is_projective = t0.is_projective() || t1.is_projective();
	}

	

	/*
	 * Return filtered RESULT and WEIGHT at point P in a coordinate system
	 * specified by the inverse of transformation T based on data taken
	 * from image IM.
	 */
	void filtered(int i, int j, pixel *result, pixel *weight, 
			int honor_exclusion, int frame, 
			pixel prev_value = pixel(0, 0, 0), 
			pixel prev_weight = pixel(0, 0, 0)) const {

		point p = point(i, j) + offset;

		*result = pixel(0, 0, 0);
		*weight = pixel(0, 0, 0);
		
		ale_pos hscale_g = 1;
		ale_pos hscale_rb = 1;
		ale_pos wscale_g = 1;
		ale_pos wscale_rb = 1;

		if (t_two) {
			ts1 = t_at_point(p);
			ts0 = t_at_inv_point(ts1.transform_unscaled(p));
		} else {
			ts0 = t_at_inv_point(p);
		}

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

		filter_channel(p, mapped_p, 0, hscale_rb, wscale_rb, result, weight, honor_exclusion, frame, prev_value[0], prev_weight[0]);
		filter_channel(p, mapped_p, 2, hscale_rb, hscale_rb, result, weight, honor_exclusion, frame, prev_value[2], prev_weight[2]);
		filter_channel(p, mapped_p, 1, hscale_g , hscale_g , result, weight, honor_exclusion, frame, prev_value[1], prev_weight[1]);
	}
};
#endif
