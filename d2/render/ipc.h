// Copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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
 * ipc.h: A render subclass implementing an iterative image reconstruction
 * algorithm based on the work of Michal Irani and Shmuel Peleg.  The
 * projection and backprojection functions used are user-configurable.
 * For more information on the theory behind the algorithm, see the last
 * section and appendix of:
 *
 * 	http://www.wisdom.weizmann.ac.il/~irani/abstracts/superResolution.html
 *
 * The algorithm in the source paper looks something like this (PSF' is the
 * backprojection kernel, and corresponds to what the authors of the paper call
 * AUX):
 *
 * ===============================================================
 *    Forward         Backward           Points of comparison
 * ---------------------------------------------------------------
 *
 *    scene(n)        scene(n+1)
 *
 *      |                 ^
 *      |                 |
 *     PSF               PSF'
 *      |                 |
 *      |        ---------+              <--- difference
 *      V       /         |
 *
 *   simulated(n)       real
 *
 * ===============================================================
 *
 * This assumes a single colorspace representation.  However, consumer cameras
 * sometimes perform sharpening in non-linear colorspace, whereas lens and
 * sensor blurring occurs in linear colorspace.  Hence, there can be two
 * colorspaces involved; ALE accounts for this with linear and non-linear
 * colorspace PSFs.  Hence, the algorithm we use looks something like:
 *
 * ===============================================================
 *    Forward         Backward           Points of comparison
 * ---------------------------------------------------------------
 *
 *    scene(n)        scene(n+1)
 *
 *      |                 ^
 *      |                 |
 *    LPSF              LPSF'
 *      |                 |
 *      |       ----------+               <--- difference,
 *      V      /          |                    exposure
 *                                             re-estimation
 *  lsimulated(n)       lreal         
 *                                   
 *      |                 ^         
 *      |                 |
 *   unlinearize       linearize
 *      |                 |
 *      V                 |
 *
 *  lsim_nl(n)        lreal_nl(n)
 *   
 *      |                 ^
 *      |                 |
 *    NLPSF             NLPSF'
 *      |                 |
 *      |       ----------+               <--- difference
 *      V      /          |
 *
 *  nlsimulated(n)     real_nl
 *
 *                        ^
 *                        |
 *                    unlinearize
 *                        |
 *                        |
 *
 *                      real
 *
 * ===============================================================
 */

#ifndef __ipc_h__
#define __ipc_h__

#include "../image.h"
#include "../render.h"
#include "psf/rasterizer.h"
#include "psf/normalizer.h"
#include "psf/backprojector.h"

class ipc : public render {
protected:
        int done;
	int inc;
        image *approximation;
        render *input;
        unsigned int iterations;
	psf *lresponse, *nlresponse;
	int exposure_register;

	/*
	 * Calculate the simulated input frame SIMULATED from the image
	 * approximation APPROXIMATION, iterating over image approximation
	 * pixels.  
	 */

	void _ip_frame_simulate(int frame_num, image *approximation, 
			image *lsimulated, image *nlsimulated, 
			transformation t, const raster *lresponse, 
			const raster *nlresponse, const exposure &exp) {

		/*
		 * Initializations for linear filtering
		 */

		image_ale_real *lsim_weights = new image_ale_real(
				lsimulated->height(),
				lsimulated->width(),
				lsimulated->depth());

		/*
		 * Linear filtering, iterating over approximation pixels
		 */

                for (unsigned int i = 0; i < approximation->height(); i++)
                for (unsigned int j = 0; j < approximation->width();  j++) {

			if (is_excluded(approximation->offset(), i, j, frame_num))
				continue;

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
			 * Convenient variables for expressing the boundaries
			 * of the mapped area.
			 */
			
			ale_pos top = q[0] - d[0];
			ale_pos bot = q[0] + d[0];
			ale_pos lef = q[1] - d[1];
			ale_pos rig = q[1] + d[1];

			/*
			 * Iterate over all simulated frame pixels influenced
			 * by the scene pixel (i, j), as determined by the 
			 * response function.
			 */

                        for (int ii = (int) floor(top + lresponse->min_i());
                                ii <= ceil(bot + lresponse->max_i()); ii++)
                        for (int jj = (int) floor(lef + lresponse->min_j());
                                jj <= ceil(rig + lresponse->max_j()); jj++) {

				if (ii < (int) 0
				 || ii >= (int) lsimulated->height()
				 || jj < (int) 0
				 || jj >= (int) lsimulated->width())
					continue;

				class rasterizer::psf_result r =
					(*lresponse)(top - ii, bot - ii,
						    lef - jj, rig - jj,
						    lresponse->select(ii, jj));

				lsimulated->pix(ii, jj) +=
					r(approximation->get_pixel(i, j));
				lsim_weights->pix(ii, jj) +=
					r.weight();
                        }
                }

		/*
		 * Normalize linear
		 */

		for (unsigned int ii = 0; ii < lsimulated->height(); ii++)
		for (unsigned int jj = 0; jj < lsimulated->width(); jj++) {
			pixel weight = lsim_weights->get_pixel(ii, jj);
			const ale_real weight_floor = 0.00001;

			if (weight[0] > weight_floor
			 && weight[1] > weight_floor
			 && weight[2] > weight_floor)
				lsimulated->pix(ii, jj)
					/= weight;
			else
				lsimulated->pix(ii, jj)
					/= 0;  /* Generate a non-finite value */
		}

		/*
		 * Finalize linear
		 */

		delete lsim_weights;

		/*
		 * Return if there is no non-linear step.
		 */

		if (nlsimulated == NULL)
			return;

		/*
		 * Initialize non-linear
		 */

		image_ale_real *nlsim_weights = new image_ale_real(
				nlsimulated->height(),
				nlsimulated->width(),
				nlsimulated->depth());

		/*
		 * Iterate non-linear
		 */

                for (unsigned int i = 0; i < lsimulated->height(); i++)
                for (unsigned int j = 0; j < lsimulated->width();  j++) {

			/*
			 * Convenient variables for expressing the boundaries
			 * of the mapped area.
			 */
			
			ale_pos top = i - 0.5;
			ale_pos bot = i + 0.5;
			ale_pos lef = j - 0.5;
			ale_pos rig = j + 0.5;

			/*
			 * Iterate over all simulated frame pixels influenced
			 * by the scene pixel (i, j), as determined by the 
			 * response function.
			 */

                        for (int ii = (int) floor(top + nlresponse->min_i());
                                ii <= ceil(bot + nlresponse->max_i()); ii++)
                        for (int jj = (int) floor(lef + nlresponse->min_j());
                                jj <= ceil(rig + nlresponse->max_j()); jj++) {

				if (ii < (int) 0
				 || ii >= (int) nlsimulated->height()
				 || jj < (int) 0
				 || jj >= (int) nlsimulated->width())
					continue;

				class rasterizer::psf_result r =
					(*nlresponse)(top - ii, bot - ii,
						    lef - jj, rig - jj,
						    nlresponse->select(ii, jj));

				nlsimulated->pix(ii, jj) +=
					r(exp.unlinearize(lsimulated->get_pixel(i, j)));
				nlsim_weights->pix(ii, jj) +=
					r.weight();
                        }
                }

		/*
		 * Normalize non-linear
		 */

		for (unsigned int ii = 0; ii < nlsimulated->height(); ii++)
		for (unsigned int jj = 0; jj < nlsimulated->width(); jj++) {
			pixel weight = nlsim_weights->get_pixel(ii, jj);
			ale_real weight_floor = 0.00001;

			if (weight[0] > weight_floor
			 && weight[1] > weight_floor
			 && weight[2] > weight_floor)
				nlsimulated->pix(ii, jj)
					/= weight;
			else
				nlsimulated->pix(ii, jj)
					/= 0;  /* Generate a non-finite value */
		}

		/*
		 * Finalize non-linear
		 */

		delete nlsim_weights;
	}



	/*
	 * For each pixel in APPROXIMATION, calculate the differences
	 * between SIMULATED and REAL pixels influenced by this pixel,
	 * and back-project the differences onto the correction array
	 * C.  The number of frames backprojected to each pixel in C is
	 * counted in CC.  
	 *
	 * Since APPROXIMATION can always, given sufficient computational
	 * resources, be configured to be of finer resolution than SIMULATED,
	 * we map APPROXIMATION pixels onto the SIMULATED grid rather than vice
	 * versa.  This should reduce the error incurred by approximations in
	 * d2::transformation::unscaled_map_area*().  
	 *
	 * This approach requires multiplying the resultant integral by a
	 * corrective factor, since the grids are generally of different
	 * densities.
	 */

	void _ip_frame_correct(int frame_num, image *approximation, image *c, 
			image *cc, const image *real, image *lsimulated, 
			image *nlsimulated, transformation t, 
			const backprojector *lresponse, 
			const backprojector *nlresponse) {

		/*
		 * Generate the image to compare lsimulated with.
		 */

		const image *lreal;

		if (nlsimulated == NULL)
			lreal = real;
		else {

			image *new_lreal = new image_ale_real(
					real->height(),
					real->width(),
					real->depth(),
					"IPC lreal", 
					&real->exp());

			/*
			 * Unlinearize values from lsimulated.
			 */

			for (unsigned int i = 0; i < lsimulated->height(); i++)
			for (unsigned int j = 0; j < lsimulated->width(); j++)
				new_lreal->set_pixel(i, j, real->exp().unlinearize(
							lsimulated->get_pixel(i, j)));

			/*
			 * Backproject from real to lreal, iterating over all pixels
			 * in lreal.
			 */

			for (unsigned int i = 0; i < new_lreal->height(); i++)
			for (unsigned int j = 0; j < new_lreal->width();  j++) {

				if (is_excluded(approximation->offset(), i, j, frame_num))
					continue;


				/*
				 * Convenient variables for expressing the boundaries
				 * of the mapped area.
				 */
				
				ale_pos top = i - 0.5;
				ale_pos bot = i + 0.5;
				ale_pos lef = j - 0.5;
				ale_pos rig = j + 0.5;

				/*
				 * Iterate over non-linear pixels influenced by linear
				 * pixels.
				 */

				for (int ii = (int) floor(top + nlresponse->min_i());
					ii <= ceil(bot + nlresponse->max_i()); ii++)
				for (int jj = (int) floor(lef + nlresponse->min_j());
					jj <= ceil(rig + nlresponse->max_j()); jj++) {

					if (ii < (int) 0
					 || ii >= (int) nlsimulated->height()
					 || jj < (int) 0
					 || jj >= (int) nlsimulated->width())
						continue;

					class backprojector::psf_result r =
						(*nlresponse)(top - ii, bot - ii,
							    lef - jj, rig - jj,
							    nlresponse->select(ii, jj));


					pixel comp_real =
						real->exp().unlinearize(real->get_pixel(ii, jj));

					pixel comp_simu =
						nlsimulated->get_pixel(ii, jj);

					if (!finite(comp_simu[0])
					 || !finite(comp_simu[1])
					 || !finite(comp_simu[2]))
						continue;

					/*
					 * Backprojection value.
					 */

					pixel bpv = r(comp_real - comp_simu);

					/*
					 * Error calculation
					 */

					new_lreal->pix(i, j) += bpv;
				}
			}

			/*
			 * Linearize lreal.
			 */

			for (unsigned int i = 0; i < new_lreal->height(); i++)
			for (unsigned int j = 0; j < new_lreal->width(); j++)
				new_lreal->set_pixel(i, j, real->exp().linearize(
							new_lreal->get_pixel(i, j)));

			lreal = new_lreal;
		}

		/*
		 * Perform exposure adjustment.
		 *
		 * XXX: it would be cleaner to remove the 'nlsimulated' term
		 * from the following test, but, empirically, this causes
		 * problems with raw Bayer pattern data, so we leave this as-is
		 * for now.
		 */

		if (exposure_register && nlsimulated) {

			pixel ec;

#if 0
			ec = lsimulated->avg_channel_magnitude()
			   / lreal->avg_channel_magnitude();
#elsif 0
			pixel_accum ratio_sum;
			pixel_accum weight_sum;

			for (unsigned int i = 0; i < lreal->height(); i++)
			for (unsigned int j = 0; j < lreal->width(); j++) {
				pixel s = lsimulated->get_pixel(i, j);
				pixel r = lreal->get_pixel(i, j);
				pixel confidence = real->exp().confidence(r);

				if (s[0] > 0.001 
				 && s[1] > 0.001
				 && s[2] > 0.001
				 && r[0] > 0.001
				 && r[1] > 0.001
				 && r[2] > 0.001) {
					ratio_sum += confidence * s / r;
					weight_sum += confidence;
				}
			}

			ec = ratio_sum / weight_sum;
#else
			pixel_accum ssum, rsum;

			for (unsigned int i = 0; i < lreal->height(); i++)
			for (unsigned int j = 0; j < lreal->width(); j++) {
				pixel s = lsimulated->get_pixel(i, j);
				pixel r = lreal->get_pixel(i, j);
				pixel confidence = real->exp().confidence(r)
					         * real->exp().confidence(s);

				ssum += confidence * s;
				rsum += confidence * r;
			}

			ec = ssum / rsum;

#endif

			real->exp().set_multiplier(
				real->exp().get_multiplier() * ec);
		}


		/*
		 * Iterate over all pixels in the approximation.
		 */

                for (unsigned int i = 0; i < approximation->height(); i++)
                for (unsigned int j = 0; j < approximation->width();  j++) {

			/*
			 * Obtain the position Q and dimensions D of image
			 * approximation pixel (i, j) in the coordinate system
			 * of the simulated (and real) frame.
			 */

                        point p = point(i + approximation->offset()[0], j + approximation->offset()[1]);
			point q;
			ale_pos d[2];

			t.unscaled_map_area_inverse(p, &q, d);
			
			/*
			 * Convenient variables for expressing the boundaries
			 * of the mapped area.
			 */
			
			ale_pos top = q[0] - d[0];
			ale_pos bot = q[0] + d[0];
			ale_pos lef = q[1] - d[1];
			ale_pos rig = q[1] + d[1];

			/*
			 * Iterate over frame pixels influenced by the scene
			 * pixel.
			 */

                        for (int ii = (int) floor(top + lresponse->min_i());
                                ii <= ceil(bot + lresponse->max_i()); ii++)
                        for (int jj = (int) floor(lef + lresponse->min_j());
                                jj <= ceil(rig + lresponse->max_j()); jj++) {

				if (ii < (int) 0
				 || ii >= (int) lreal->height()
				 || jj < (int) 0
				 || jj >= (int) lreal->width())
					continue;

				unsigned int selection = lresponse->select(ii, jj);

				class backprojector::psf_result r =
					(*lresponse)(top - ii, bot - ii,
						    lef - jj, rig - jj,
						    selection);


				/*
				 * R is the result of integration in the
				 * coordinate system of the simulated frame.
				 * We must rescale to get the result of
				 * integration in the coordinate system of the
				 * approximation image.
				 */

				r *= (1 / (bot - top) / (rig - lef));

				pixel comp_lreal =
					lreal->get_pixel(ii, jj);

				// pixel comp_real = 
				//	real->get_pixel(ii, jj);

				pixel comp_simu =
					lsimulated->get_pixel(ii, jj);

				if (!finite(comp_simu[0])
				 || !finite(comp_simu[1])
				 || !finite(comp_simu[2])
				 || !finite(comp_lreal[0])
				 || !finite(comp_lreal[1])
				 || !finite(comp_lreal[2]))
					continue;

				/*
				 * Backprojection value unadjusted 
				 * for confidence.
				 */

				pixel bpv = r(comp_lreal - comp_simu);

				/*
				 * Confidence [equal to (1, 1, 1) when
				 * confidence is uniform].
				 */

				// Ordinary certainty
				// pixel conf = real->exp().confidence(comp_lreal);

				// One-sided certainty
				// pixel conf = real->exp().one_sided_confidence(comp_lreal, bpv);
				      // conf = real->exp().one_sided_confidence(comp_real, bpv);
				      
				// Estimate-based certainty
				// pixel conf = real->exp().confidence(comp_simu);
					
				// One-sided estimate-based certainty
				pixel conf = real->exp().one_sided_confidence(comp_simu, bpv);
					
				/*
				 * If a color is bayer-interpolated, then we have no confidence in its
				 * value.
				 */

				if (real->get_bayer() != IMAGE_BAYER_NONE) {
					int color = ((image_bayer_ale_real *)real)->bayer_color(ii, jj);
					conf[(color + 1) % 3] = 0;
					conf[(color + 2) % 3] = 0;
				}

				/*
				 * Error calculation
				 */

				c->pix(i, j) += bpv * conf;

				/*
				 * Increment the backprojection weight.  When
				 * confidence is uniform, this should weight
				 * each frame's correction equally.
				 */

				cc->pix(i, j) += conf * r.weight() 
					       / lresponse->integral(selection);
			}
                }

		if (nlsimulated)
			delete lreal;
	}

	/* 
	 * Adjust correction array C based on the difference between the
	 * simulated projected frame and real frame M.  Update the correction
	 * count CC for affected pixels in C.
	 */

	virtual void _ip_frame(int frame_num, image *c, image *cc, const image *real, 
			transformation t, const raster *f, const backprojector *b,
			const raster *nlf, const backprojector *nlb) {

		/*
		 * Initialize simulated data structures
		 */

                image *lsimulated = new image_ale_real(
				real->height(),
				real->width(), 
				real->depth());

                image *nlsimulated = NULL;

		if (nlf)
			nlsimulated = new image_ale_real(
				real->height(),
				real->width(), 
				real->depth());

		/*
		 * Create simulated frames with forward projection.
		 */

		_ip_frame_simulate(frame_num, approximation, lsimulated, nlsimulated, t, f, nlf, real->exp());

		/*
		 * Update the correction array using backprojection.
		 */

		_ip_frame_correct(frame_num, approximation, c, cc, real, lsimulated, nlsimulated, t, b, nlb);

		/*
		 * Finalize data structures.
		 */

                delete lsimulated;
		delete nlsimulated;
        }

	/*
	 * Iterate _ip_frame() over all frames, and update APPROXIMATION after
	 * corrections from all frames have been summed.  Repeat for the number
	 * of iterations specified by the user.
	 */

        virtual void _ip() {

		/*
		 * Create rasterized PSF and backprojection kernel AUX.
		 */

		raster **f = (raster **) malloc(image_rw::count() * sizeof(raster *));
		backprojector **b = (backprojector **) malloc(image_rw::count() * sizeof(backprojector *));

		for (unsigned int m = 0; m < image_rw::count(); m++) {

			if (!align::match(m))
				continue;

			transformation t = align::of(m);

			f[m] = new normalizer(new rasterizer(lresponse, t));
			b[m] = new backprojector(f[m]);
		}

		raster *nlf = NULL;
		backprojector *nlb = NULL;

		if (nlresponse) {
			nlf = new normalizer(new rasterizer(nlresponse, transformation::eu_identity()));
			nlb = new backprojector(nlf);
		}

                for (unsigned int n = 0; n < iterations; n++) {

                        image *correction = new image_ale_real (
                                        approximation->height(),
                                        approximation->width(),
                                        approximation->depth());

			image *correction_count = new image_ale_real(
					approximation->height(),
					approximation->width(),
					approximation->depth());

			/*
			 * Iterate over all frames
			 */

                        for (unsigned int m = 0; m < image_rw::count(); m++) {

				if (!align::match(m))
					continue;

				ui::get()->ip_frame_start(m);

				transformation t = align::of(m);
				const image *real = image_rw::open(m);

                                _ip_frame(m, correction, correction_count, real,
					t, f[m], b[m], nlf, nlb);

				image_rw::close(m);
			}

			/*
			 * Update the approximation.
			 */

			ui::get()->ip_update();

			for (unsigned int i = 0; i < approximation->height(); i++)
			for (unsigned int j = 0; j < approximation->width();  j++) {

				pixel  cpix = correction       ->get_pixel(i, j);
				pixel ccpix = correction_count ->get_pixel(i, j);
				pixel  apix = approximation    ->get_pixel(i, j);

				for (unsigned int k = 0; k < 3; k++) {

					const ale_real cc_floor = 0.00001;

					if (ccpix[k] < cc_floor)
						continue;

					ale_real new_value = cpix[k] / ccpix[k] + apix[k];

					assert (finite(apix[k]));
					assert (finite(ccpix[k]));
					assert (finite(cpix[k]));
					assert (finite(new_value));

					/*
					 * Negative light doesn't make sense.
					 */
					if (new_value < 0)
						new_value = 0;

					approximation->chan(i, j, k) = new_value;
				}
			}

                        delete correction;
			delete correction_count;

			if (inc) {
				ui::get()->ip_write();
				image_rw::output(approximation);
			}

			ui::get()->ip_step_done();

                }

		for (unsigned int m = 0; m < image_rw::count(); m++) {

			if (!align::match(m))
				continue;

			delete f[m];
			delete b[m];
		}

		free(f);
		free(b);

		delete nlf;
		delete nlb;
        }

public:

        ipc(render *input, unsigned int iterations, int _inc, psf *lresponse, psf *nlresponse, int exposure_register) {
                this->input = input;
                done = 0;
		inc = _inc;
                this->iterations = iterations;
		this->lresponse = lresponse;
		this->nlresponse = nlresponse;
		this->exposure_register = exposure_register;
        }

        const image *get_image() {
                if (done)
                        return approximation;
                else
                        return input->get_image();
        }

        const image *get_defined() {
                return input->get_defined();
        }

        void sync(int n) {
		render::sync(n);
                input->sync(n);
        }

	void step() {
		return;
	}

        virtual int sync() {
		input->sync();
		ui::get()->ip_start();
                done = 1;
                approximation = optimizations::get_ip_working_image(input->get_image());
                _ip();
		ui::get()->ip_done();

                return 0;
        }

	virtual ~ipc() {
	}

	void free_memory() {
	}
};

#endif