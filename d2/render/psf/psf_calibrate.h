// Copyright 2003 David Hilvert <dhilvert@auricle.dyndns.org>,
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

#ifndef __psf_calibrate_h__
#define __psf_calibrate_h__

#include "../../image.h"
#include "../../render.h"
#include "../ipc.h"

class psf_calibrate : public ipc {
public:
	psf_calibrate(render *input, unsigned int iterations, int _inc, psf *lresponse, psf *nlresponse) 
			: ipc(input, iterations, _inc, lresponse, nlresponse, 1) {
		fprintf(stderr, "\nIPC Calibration module.\n\n");
		fprintf(stderr, "This module is designed for use with a calibration script.\n\n");
	}

        virtual void _ip_frame(ale_accum *diff, int m) {

		/*
		 * Get alignment information for frame m.
		 */

                transformation t = align::of(m);

		/*
		 * We create real and simulated input-frame data structures
		 * REAL and SIMULATED, as well as simulated input-frame weights
		 * SIM_WEIGHTS, used to track the weights of contributions to each
		 * simulated input-frame pixel component.
		 */

                const image *real = image_rw::open(m);
                image *simulated = new image_ale_real(
				real->height(),
				real->width(), 3);

		/*
		 * Calculate the simulated input frame SIMULATED from the image
		 * approximation APPROXIMATION, iterating over image
		 * approximation pixels and tracking contributions to simulated
		 * frame pixels in the data structure SIM_WEIGHTS.
		 */

                image *sim_weights = new image_ale_real(
                                simulated->height(),
                                simulated->width(), 3);

                for (unsigned int i = 0; i < approximation->height(); i++)
                for (unsigned int j = 0; j < approximation->width();  j++) {

			/*
			 * Obtain the position Q and dimensions D of
			 * image approximation pixel (i, j) in the coordinate
			 * system of the simulated frame.
			 */

                        point p = point(i + approximation->offset()[0], j + approximation->offset()[1]);
			point q;
			ale_pos d[2];

			t.unscaled_map_area_inverse(p, &q, d);
			
			/*
			 * Iterate over all simulated frame pixels influenced
			 * by the scene pixel (i, j), as determined by the 
			 * response function.
			 */

                        for (int ii = (int) floor(q[0] - d[0] + lresponse->min_i());
                                ii <= ceil(q[0] + d[0] + lresponse->max_i()); ii++)
                        for (int jj = (int) floor(q[1] - d[1] + lresponse->min_j());
                                jj <= ceil(q[1] + d[1] + lresponse->max_j()); jj++) {

                                ale_pos top = q[0] - d[0];
                                ale_pos bot = q[0] + d[0];
                                ale_pos lef = q[1] - d[1];
                                ale_pos rig = q[1] + d[1];

                                if (ii >= (int) 0
                                 && ii < (int) real->height()
                                 && jj >= (int) 0
                                 && jj < (int) real->width()) {

					psf::psf_result r =
						(*lresponse)(top - ii, bot - ii,
							    lef - jj, rig - jj,
							    lresponse->select(ii, jj));

					sim_weights->pix(ii, jj) += r.weight();
					simulated->pix(ii, jj) += r(approximation->get_pixel(i, j));
                                }
                        }
                }

		/*
		 * Normalize SIMULATED by SIM_WEIGHTS
		 */

		for (unsigned int i = 0; i < simulated->height(); i++)
		for (unsigned int j = 0; j < simulated->width();  j++)
			simulated->pix(i, j) 
				/= sim_weights->get_pixel(i, j);

		delete sim_weights;

		/*
		 * If NLRESPONSE is defined, then redefine SIMULATED to account
		 * for this.
		 */

		if (nlresponse != NULL) {
			image *nlsimulated = new image_ale_real(
					simulated->height(),
					simulated->width(), 3);

			image *nlsim_weights = new image_ale_real(
					simulated->height(),
					simulated->width(), 3);

			for (unsigned int i = 0; i < simulated->height(); i++)
			for (unsigned int j = 0; j < simulated->width();  j++) {

				for (int ii = (int) floor(i - 0.5 + nlresponse->min_i());
					ii <= ceil(i + 0.5 + nlresponse->max_i()); ii++)
				for (int jj = (int) floor(i - 0.5 + nlresponse->min_j());
					jj <= ceil(i + 0.5 + nlresponse->max_j()); jj++) {

					ale_pos top = i - 0.5;
					ale_pos bot = i + 0.5;
					ale_pos lef = j - 0.5;
					ale_pos rig = j + 0.5;

					if (ii >= (int) 0
					 && ii < (int) nlsimulated->height()
					 && jj >= (int) 0
					 && jj < (int) nlsimulated->width()) {

						psf::psf_result r =
							(*nlresponse)(top - ii, bot - ii,
								    lef - jj, rig - jj,
								    nlresponse->select(i, j));

						nlsim_weights->pix(ii, jj) += r.weight();
						nlsimulated->pix(ii, jj) += r(real->exp().unlinearize(simulated->get_pixel(i, j)));
					}
				}
			}

			/*
			 * Normalize nlsimulated.
			 */

			for (unsigned int i = 0; i < simulated->height(); i++)
			for (unsigned int j = 0; j < simulated->width();  j++)
				nlsimulated->pix(i, j) 
					/= nlsim_weights->get_pixel(i, j);

			/*
			 * Linearize nlsimulated
			 */

			for (unsigned int i = 0; i < simulated->height(); i++)
			for (unsigned int j = 0; j < simulated->width();  j++)
				nlsimulated->pix(i, j) =
					real->exp().linearize(nlsimulated->get_pixel(i, j));

			delete simulated;
			delete nlsim_weights;

			simulated = nlsimulated;
		}

		/*
		 * For each SIMULATED pixel, calculate the difference from
		 * the corresponding REAL pixel, and update the sum of squares
		 * of differences.
		 */

		ale_pos margin_i1 = lresponse->min_i() + (nlresponse ? nlresponse->min_i() : 0);
		ale_pos margin_i2 = lresponse->max_i() + (nlresponse ? nlresponse->max_i() : 0);
		ale_pos margin_j1 = lresponse->min_j() + (nlresponse ? nlresponse->min_j() : 0);
		ale_pos margin_j2 = lresponse->max_j() + (nlresponse ? nlresponse->max_j() : 0);

		for (unsigned int i = 0; i < simulated->height(); i++)
		for (unsigned int j = 0; j < simulated->width();  j++) {

			/*
			 * Establish margins.  This is designed to reduce the
			 * influence of boundary conditions.
			 */

			point p;

			p = t.transform_unscaled(point(i + margin_i1, j + margin_j1));
			if (p[0] < 0 || p[0] > approximation->height()
			 || p[1] < 0 || p[1] > approximation->width())
				continue;

			p = t.transform_unscaled(point(i + margin_i1, j + margin_j2));
			if (p[0] < 0 || p[0] > approximation->height()
			 || p[1] < 0 || p[1] > approximation->width())
				continue;

			p = t.transform_unscaled(point(i + margin_i2, j + margin_j1));
			if (p[0] < 0 || p[0] > approximation->height()
			 || p[1] < 0 || p[1] > approximation->width())
				continue;

			p = t.transform_unscaled(point(i + margin_i2, j + margin_j2));
			if (p[0] < 0 || p[0] > approximation->height()
			 || p[1] < 0 || p[1] > approximation->width())
				continue;

			/*
			 * Real and simulated responses
			 */

			pixel comp_real = real->get_pixel(i, j);
			pixel comp_simu = simulated->get_pixel(i, j);

			for (unsigned int k = 0; k < simulated->depth();  k++) {

				if (!finite(comp_simu[k]))
					continue;


				/*
				 * Error calculation
				 */

				if ((comp_real[k] < 1.0 || comp_simu[k] < 1.0 )
				 && (comp_real[k] > 0 || comp_simu[k] > 0))

					/*
					 * real and simulated are distinguishable
					 * within the dynamic range of the program
					 * inputs, so calculate the error for this
					 * channel.
					 */

					(*diff) += pow(comp_simu[k] - comp_real[k], 2);

			}
		}

		image_rw::close(m); 
		delete simulated;

	}

	virtual void _ip() {

		/*
		 * Input images 0 through count()-2 are frames captured with
		 * the device to be calibrated, so we combine the difference
		 * values for all of these frames against the calibration image
		 * count()-1.
		 */

		ale_accum diff = 0;

		approximation = image_rw::copy(image_rw::count() - 1, "PSF_CALIBRATE reference");

		for (int m = 0; m < image_rw::count() - 1; m++) {
			_ip_frame(&diff, m);
		}

		diff = pow(diff, 0.5);

		fprintf(stderr, "\n\nPSF Error:: %e\n\n", (double) diff);

		delete approximation;
        }
};

#endif
