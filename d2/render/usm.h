// Copyright 2002, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

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
 * usm.h: A render subclass that implements an unsharp mask
 * postprocessing algorithm.
 */

#ifndef __usm_h__
#define __usm_h__

#include "../image.h"
#include "../render.h"
#include "psf/psf.h"

class usm : public render {

	int done;
	int inc;
	image *done_image;
	render *input;
	const image *input_image;
	const image *input_defined;
	ale_real scale_factor;
	ale_real usm_multiplier;
	psf *lresponse, *nlresponse;
	const exposure *exp;

	/* 
	 * USM value for point (i, j).
	 */
	pixel _usm(int i, int j, const image *im) {

		pixel result;

		ale_real d = scale_factor / 2;

		pixel weight;

		/*
		 * Convolve with the linear filter, iterating over pixels 
		 * according to the filter support, and tracking contribution
		 * weights in the variable WEIGHT.
		 */

		for (int ii = (int) floor(i - d + lresponse->min_i());
			ii <= ceil(i + d + lresponse->max_i()); ii++)
		for (int jj = (int) floor(j - d + lresponse->min_j());
			jj <= ceil(j + d + lresponse->max_j()); jj++) {

			ale_real top = ii - d;
			ale_real bot = ii + d;
			ale_real lef = jj - d;
			ale_real rig = jj + d;

			if (ii >= (int) 0
			 && ii < (int) im->height()
			 && jj >= (int) 0
			 && jj < (int) im->width()
			 && input_defined->get_pixel(ii, jj)[0]) {

				class psf::psf_result r = (*lresponse)((top - i) / scale_factor,
								(bot - i) / scale_factor,
								(lef - j) / scale_factor,
								(rig - j) / scale_factor);

				if (nlresponse) {

					/*
					 * Convolve with the non-linear filter,
					 * iterating over pixels according to the
					 * filter support, and tracking contribution
					 * weights in the variable WWEIGHT.
					 *
					 * Note:  This approach is efficient
					 * space-wise, but inefficient timewise.  There
					 * is probably a better approach to this.
					 */

					pixel rresult(0, 0, 0), wweight(0, 0, 0);

					for (int iii = (int) floor(ii - d + nlresponse->min_i());
						iii <= ceil(ii + d + lresponse->max_i()); iii++)
					for (int jjj = (int) floor(jj - d + nlresponse->min_j());
						jjj <= ceil(jj + d + lresponse->max_j()); jjj++) {

						ale_real top = iii - d;
						ale_real bot = iii + d;
						ale_real lef = jjj - d;
						ale_real rig = jjj + d;

						if (iii >= (int) 0
						 && iii < (int) im->height()
						 && jjj >= (int) 0
						 && jjj < (int) im->width()
						 && input_defined->get_pixel(iii, jjj)[0]) {

							class psf::psf_result r = (*nlresponse)((top - ii) / scale_factor,
											(bot - ii) / scale_factor,
											(lef - jj) / scale_factor,
											(rig - jj) / scale_factor);

							wweight += r.weight();
							rresult += r(exp->unlinearize(im->get_pixel(iii, jjj)));
						}
					}

					result += r(exp->linearize(rresult / wweight));
				} else {
					result += r(im->get_pixel(ii, jj));
				}

				weight += r.weight();

			}
		}

		result /= weight;
		result = im->get_pixel(i, j) - result;

		if (finite(result[0])
		 && finite(result[1])
		 && finite(result[2]))
			return result;
		else
			return pixel(0, 0, 0);
	}

	void _filter() {
		assert (done_image->height() == input_image->height());
		assert (done_image->width() == input_image->width());
		assert (done_image->depth() == input_image->depth());

		for (unsigned int i = 0; i < done_image->height(); i++)
		for (unsigned int j = 0; j < done_image->width(); j++) {

			if (!input_defined->get_pixel(i, j)[0])
				continue;

			done_image->set_pixel(i, j, 
				input_image->get_pixel(i, j)
			      + usm_multiplier
			      * _usm(i, j, input_image));
		}
	}

public:

	usm(render *input, ale_real scale_factor, ale_real usm_multiplier, int _inc, 
			psf *lresponse, psf *nlresponse, exposure *exp) {
		this->input = input;
		done = 0;
		inc = _inc;
		this->scale_factor = scale_factor;
		this->usm_multiplier = usm_multiplier;
		this->lresponse = lresponse;
		this->nlresponse = nlresponse;
		this->exp = exp;
	}

	const image *get_image() {
		if (done)
			return done_image;
		else
			return input->get_image();
	}

	const image *get_defined() {
		return input->get_defined();
	}

	void sync(int n) {
		input->sync(n);
	}

	int sync() {
		input->sync();
		fprintf(stderr, "Applying USM");
		done = 1;
		done_image = input->get_image()->clone("USM done_image");
		input_image = input->get_image();
		input_defined = input->get_defined();
		_filter();

		if (inc)
			image_rw::output(done_image);

		fprintf(stderr, ".\n");

		return 0;
	}

	virtual ~usm() {
	}
};

#endif
