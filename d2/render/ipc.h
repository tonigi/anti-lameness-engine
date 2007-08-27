// Copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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
        int synced;
	int inc;
        image *approximation;
        render *input;
        unsigned int iterations;
	psf *lresponse, *nlresponse;
	int exposure_register;
	int use_weighted_median;
	double weight_limit;

	/*
	 * Calculate the simulated input frame SIMULATED from the image
	 * approximation APPROXIMATION, iterating over image approximation
	 * pixels.  
	 */

	struct sim_args {
		int frame_num;
		image *approximation;
		image *lsimulated;
		image *nlsimulated;
		image *lsim_weights;
		image *nlsim_weights;
		transformation t;
		const raster *lresponse;
		const raster *nlresponse;
		const exposure *exp;
		point *extents;
	};

	class simulate_linear : public thread::decompose_domain, 
	                        private sim_args {
		point *subdomain_extents;

	protected:
		void prepare_subdomains(unsigned int N) {
			subdomain_extents = new point[N * 2];

			for (unsigned int n = 0; n < N; n++) {
				point *se = subdomain_extents + 2 * n;

				for (int d = 0; d < 2; d++) {
					se[0][d] = extents[0][d];
					se[1][d] = extents[1][d];
				}
			}
		}

		void subdomain_algorithm(unsigned int thread,
				int i_min, int i_max, int j_min, int j_max) {

			point *extents = subdomain_extents + thread * 2;

			/*
			 * Linear filtering, iterating over approximation pixels
			 */

			for (int i = i_min; i < i_max; i++)
			for (int j = j_min; j < j_max; j++) {

				if (is_excluded_r(approximation->offset(), i, j, frame_num))
					continue;

				/*
				 * Obtain the position Q and dimensions D of
				 * image approximation pixel (i, j) in the coordinate
				 * system of the simulated frame.
				 */

				point p = point(i + approximation->offset()[0], j + approximation->offset()[1]);
				point q;
				ale_pos d[2];

				/*
				 * XXX: This appears to calculate the wrong thing.
				 */

				if (is_excluded_f(p, frame_num))
					continue;

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

					if (ii < extents[0][0])
						extents[0][0] = ii;
					if (jj < extents[0][1])
						extents[0][1] = jj;
					if (ii > extents[1][0])
						extents[1][0] = ii;
					if (jj > extents[1][1])
						extents[1][1] = jj;

					class rasterizer::psf_result r =
						(*lresponse)(top - ii, bot - ii,
							    lef - jj, rig - jj,
							    lresponse->select(ii, jj),
							    lsimulated->get_channels(ii, jj));

					lock();

					if (lsimulated->get_bayer() == IMAGE_BAYER_NONE) {
						lsimulated->pix(ii, jj) +=
							r(approximation->get_pixel(i, j));
						lsim_weights->pix(ii, jj) +=
							r.weight();
					} else {
						int k = ((image_bayer_ale_real *)lsimulated)->bayer_color(ii, jj);
						lsimulated->chan(ii, jj, k) +=
							r(approximation->get_pixel(i, j))[k];
						lsim_weights->chan(ii, jj, k) +=
							r.weight()[k];
					}

					unlock();
				}
			}
		}

		void finish_subdomains(unsigned int N) {
			/*
			 * Determine extents
			 */

			for (unsigned int n = 0; n < N; n++) {
				point *se = subdomain_extents + 2 * n;

				for (int d = 0; d < 2; d++) {
					if (se[0][d] < extents[0][d])
						extents[0][d] = se[0][d];
					if (se[1][d] > extents[1][d])
						extents[1][d] = se[1][d];
				}
			}

			delete[] subdomain_extents;
		}
	public:

		simulate_linear(sim_args s) : decompose_domain(0, s.approximation->height(),
				                               0, s.approximation->width()),
					      sim_args(s) {
		}
	};


	class simulate_nonlinear : public thread::decompose_domain, 
	                           private sim_args {
		point *subdomain_extents;

	protected:
		void prepare_subdomains(unsigned int N) {
			subdomain_extents = new point[N * 2];

			for (unsigned int n = 0; n < N; n++) {
				point *se = subdomain_extents + 2 * n;

				for (int d = 0; d < 2; d++) {
					se[0][d] = extents[0][d];
					se[1][d] = extents[1][d];
				}
			}
		}

		void subdomain_algorithm(unsigned int thread,
				int i_min, int i_max, int j_min, int j_max) {

			point *extents = subdomain_extents + thread * 2;

			/*
			 * Iterate non-linear
			 */

			for (int i = i_min; i < i_max; i++)
			for (int j = j_min; j < j_max; j++) {

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

					if (ii < extents[0][0])
						extents[0][0] = ii;
					if (jj < extents[0][1])
						extents[0][1] = jj;
					if (ii > extents[1][0])
						extents[1][0] = ii;
					if (jj > extents[1][1])
						extents[1][1] = jj;

					class rasterizer::psf_result r =
						(*nlresponse)(top - ii, bot - ii,
							    lef - jj, rig - jj,
							    nlresponse->select(ii, jj));

					lock();

					nlsimulated->pix(ii, jj) +=
						r(exp->unlinearize(lsimulated->get_pixel(i, j)));
					nlsim_weights->pix(ii, jj) +=
						r.weight();

					unlock();
				}
			}
		}

		void finish_subdomains(unsigned int N) {
			/*
			 * Determine extents
			 */

			for (unsigned int n = 0; n < N; n++) {
				point *se = subdomain_extents + 2 * n;

				for (int d = 0; d < 2; d++) {
					if (se[0][d] < extents[0][d])
						extents[0][d] = se[0][d];
					if (se[1][d] > extents[1][d])
						extents[1][d] = se[1][d];
				}
			}

			delete[] subdomain_extents;
		}
	public:

		simulate_nonlinear(sim_args s) : decompose_domain(0, s.lsimulated->height(),
						                  0, s.lsimulated->width()),
					         sim_args(s) {
		}
	};

	void _ip_frame_simulate(int frame_num, image *approximation, 
			image *lsimulated, image *nlsimulated, 
			transformation t, const raster *lresponse, 
			const raster *nlresponse, const exposure &exp,
			point *extents) {

		sim_args args;

		/*
		 * Initializations for linear filtering
		 */

		image *lsim_weights = NULL;

		if (lsimulated->get_bayer() == IMAGE_BAYER_NONE) {
			lsim_weights = new image_ale_real(
				lsimulated->height(),
				lsimulated->width(),
				lsimulated->depth());
		} else {
			lsim_weights = new image_bayer_ale_real(
				lsimulated->height(),
				lsimulated->width(),
				lsimulated->depth(), 
				lsimulated->get_bayer());
		}
		assert (lsim_weights);

		args.frame_num = frame_num;
		args.approximation = approximation;
		args.lsimulated = lsimulated;
		args.nlsimulated = nlsimulated;
		args.lsim_weights = lsim_weights;
		args.t = t;
		args.lresponse = lresponse;
		args.nlresponse = nlresponse;
		args.exp = &exp;
		args.extents = extents;

		/*
		 * Simulate linear
		 */

		simulate_linear sl(args);
		sl.run();

		/*
		 * Normalize linear
		 */

		for (unsigned int ii = 0; ii < lsimulated->height(); ii++)
		for (unsigned int jj = 0; jj < lsimulated->width(); jj++) {
			const ale_real weight_floor = 1e-10;
			ale_accum zero = 0;

			for (int k = 0; k < 3; k++) {
				if ((lsimulated->get_channels(ii, jj) & (1 << k)) == 0)
					continue;

				ale_real weight = lsim_weights->chan(ii, jj, k);

				if (!(weight > weight_floor))
					lsimulated->chan(ii, jj, k)
						/= zero;  /* Generate a non-finite value */

				lsimulated->chan(ii, jj, k) /= weight;
			}
		}

		/*
		 * Finalize linear
		 */

		delete lsim_weights;

		/*
		 * Return if there is no non-linear step.
		 */

		if (nlsimulated == NULL) {
			return;
		}

		/*
		 * Initialize non-linear
		 */

		image_ale_real *nlsim_weights = new image_ale_real(
				nlsimulated->height(),
				nlsimulated->width(),
				nlsimulated->depth());

		args.nlsim_weights = nlsim_weights;

		/*
		 * Simulate non-linear
		 */

		simulate_nonlinear snl(args);
		snl.run();

		/*
		 * Normalize non-linear
		 */

		for (unsigned int ii = 0; ii < nlsimulated->height(); ii++)
		for (unsigned int jj = 0; jj < nlsimulated->width(); jj++) {
			pixel weight = nlsim_weights->get_pixel(ii, jj);
			ale_real weight_floor = 1e-10;
			ale_accum zero = 0;

			for (int k = 0; k < 3; k++)
				if (!(weight[k] > weight_floor))
					nlsimulated->pix(ii, jj)[k]
						/= zero;  /* Generate a non-finite value */

			nlsimulated->pix(ii, jj)
				/= weight;
		}

		/*
		 * Finalize non-linear
		 */

		delete nlsim_weights;
	}

	struct correction_t {
		/*
		 * Type
		 */
		int use_weighted_median;

		/*
		 * Weighted Median
		 */
		image_weighted_median *iwm;

		/*
		 * Weight limit
		 */
		double weight_limit;

		/*
		 * Common
		 */
		image *c;
		image *cc;

		/*
		 * Create correction structures.
		 */
		correction_t(int use_weighted_median, double weight_limit, unsigned int h, unsigned int w, unsigned int d) {
			this->use_weighted_median = use_weighted_median;
			if (use_weighted_median)
				iwm = new image_weighted_median(h, w, d);
			this->weight_limit = weight_limit;
			c = new image_ale_real(h, w, d);
			cc = new image_ale_real(h, w, d);
		}

		/*
		 * Destroy correction structures
		 */
		~correction_t() {
			if (use_weighted_median)
				delete iwm;
			delete c;
			delete cc;
		}

		/*
		 * Correction count
		 */
		pixel get_count(int i, int j) {
			if (use_weighted_median)
				return iwm->get_weights()->get_pixel(i, j);
			else
				return cc->get_pixel(i, j);
		}

		/*
		 * Check for weight limit
		 */
		int weight_limited(int i, int j) {
			if (weight_limit
			 && get_count(i, j)[0] >= weight_limit
			 && get_count(i, j)[1] >= weight_limit
			 && get_count(i, j)[2] >= weight_limit)
				return 1;

			return 0;
		}

		/*
		 * Correction value
		 */
		pixel get_correction(int i, int j) {
			if (use_weighted_median)
				return iwm->get_colors()->get_pixel(i, j);
			else
				return c->get_pixel(i, j)
				     / cc->get_pixel(i, j);
		}

		/*
		 * Frame end
		 */
		void frame_end(int frame_num) {
			if (use_weighted_median) {

				for (unsigned int i = 0; i < c->height(); i++)
				for (unsigned int j = 0; j < c->width(); j++) {

					/*
					 * Update the median calculator
					 */
					pixel cval = c->get_pixel(i, j);
					pixel ccval = cc->get_pixel(i, j);
					iwm->accumulate(i, j, frame_num, cval / ccval, ccval);

					/*
					 * Reset the counters
					 */
					c->pix(i, j) = pixel::zero();
					cc->pix(i, j) = pixel::zero();
				}
			}
		}

		/*
		 * Update correction structures, using either a weighted mean update or
		 * a weighted median update.
		 */
		void update(int i, int j, pixel value_times_weight, pixel weight) {
			c->pix(i, j) += value_times_weight;
			cc->pix(i, j) += weight;
		}
	};

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

	struct correct_args {
		int frame_num;
		image *approximation;
		correction_t *cu;
		const image *real;
		image *lreal;
		image *lsimulated;
		image *nlsimulated;
		transformation t;
		const backprojector *lresponse;
		const backprojector *nlresponse;
		double weight_limit;
	};

	class correct_nonlinear : public thread::decompose_domain, 
	                          private correct_args {

	protected:
		void subdomain_algorithm(unsigned int thread,
				int i_min, int i_max, int j_min, int j_max) {

			/*
			 * Unlinearize values from lsimulated.
			 */

			for (int i = i_min; i < i_max; i++)
			for (int j = j_min; j < j_max; j++)
				lreal->set_pixel(i, j, real->exp().unlinearize(
							lsimulated->get_pixel(i, j)));


			/*
			 * Backproject from real to lreal, iterating over all pixels
			 * in lreal.
			 */

			for (int i = i_min; i < i_max; i++)
			for (int j = j_min; j < j_max; j++) {

				/*
				 * XXX: Is this right?
				 */

				if (is_excluded_r(approximation->offset(), i, j, frame_num))
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

					lreal->pix(i, j) += bpv;
				}
			}

			/*
			 * Linearize lreal.
			 */

			for (int i = i_min; i < i_max; i++)
			for (int j = j_min; j < j_max; j++)
				lreal->set_pixel(i, j, real->exp().linearize(
							lreal->get_pixel(i, j)));
		}

	public:

		correct_nonlinear(correct_args c) : decompose_domain(0, c.lreal->height(),
				                                     0, c.lreal->width()),
					            correct_args(c) {
		}
	};

	class correct_linear : public thread::decompose_domain, 
	                       private correct_args {

	protected:
		void subdomain_algorithm(unsigned int thread,
				int i_min, int i_max, int j_min, int j_max) {

			/*
			 * Iterate over all pixels in the approximation.
			 */

			for (int i = i_min; i < i_max; i++)
			for (int j = j_min; j < j_max; j++) {

				/*
				 * Check correction count against any weight limit.
				 */

				if (cu->weight_limited(i, j))
					continue;

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

					if (is_excluded_f(ii, jj, frame_num))
						continue;

					unsigned int selection = lresponse->select(ii, jj);

					char channels = lreal->get_channels(ii, jj);

					class backprojector::psf_result r =
						(*lresponse)(top - ii, bot - ii,
							    lef - jj, rig - jj,
							    selection, channels);


					/*
					 * R is the result of integration in the
					 * coordinate system of the simulated frame.
					 * We must rescale to get the result of
					 * integration in the coordinate system of the
					 * approximation image.
					 */

					r *= (1 / (bot - top) / (rig - lef));

					pixel comp_lreal =
						lreal->get_raw_pixel(ii, jj);

					// pixel comp_real = 
					//	real->get_pixel(ii, jj);

					pixel comp_simu =
						lsimulated->get_raw_pixel(ii, jj);

#if 1

					/*
					 * Under the assumption that finite() testing
					 * may be expensive, limit such tests to active
					 * channels.
					 */
					int found_nonfinite = 0;
					for (int k = 0; k < 3; k++) {
						if (((1 << k) & channels)
						 && (!finite(comp_simu[k])
						  || !finite(comp_lreal[k]))) {
							found_nonfinite = 1;
							break;
						}
					}

					if (found_nonfinite)
						continue;
#else 

					if (!finite(comp_simu[0])
					 || !finite(comp_simu[1])
					 || !finite(comp_simu[2])
					 || !finite(comp_lreal[0])
					 || !finite(comp_lreal[1])
					 || !finite(comp_lreal[2]))
						continue;
#endif

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

					// One-sided approximation-based certainty
					// pixel conf = real->exp().one_sided_confidence(approximation->pix(i, j), bpv);
						
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

					// c->pix(i, j) += bpv * conf;

					/*
					 * Increment the backprojection weight.  When
					 * confidence is uniform, this should weight
					 * each frame's correction equally.
					 */

					// cc->pix(i, j) += conf * r.weight() 
					//	       / lresponse->integral(selection);

					cu->update(i, j, bpv * conf, conf * r.weight() / lresponse->integral(selection));
				}
			}
		}

	public:

		correct_linear(correct_args c) : decompose_domain(0, c.approximation->height(),
				                                  0, c.approximation->width()),
					         correct_args(c) {
		}
	};

	void _ip_frame_correct(int frame_num, image *approximation,
			correction_t *cu, const image *real, image *lsimulated, 
			image *nlsimulated, transformation t, 
			const backprojector *lresponse, 
			const backprojector *nlresponse, point *extents) {

		correct_args args;

		args.frame_num = frame_num;
		args.approximation = approximation;
		args.cu = cu;
		args.real = real;
		args.lsimulated = lsimulated;
		args.nlsimulated = nlsimulated;
		args.t = t;
		args.lresponse = lresponse;
		args.nlresponse = nlresponse;

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

			args.lreal = new_lreal;

			correct_nonlinear cn(args);
			cn.run();

			lreal = new_lreal;
		}

		/*
		 * Perform exposure adjustment.
		 */

		if (exposure_register) {

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

//			for (unsigned int i = 0; i < lreal->height(); i++)
//			for (unsigned int j = 0; j < lreal->width(); j++) {
			for (unsigned int i = (unsigned int) floor(extents[0][0]); 
			                  i < (unsigned int)  ceil(extents[1][0]); i++)
			for (unsigned int j = (unsigned int) floor(extents[0][1]); 
			                  j < (unsigned int)  ceil(extents[1][1]); j++) {

				pixel s = lsimulated->get_pixel(i, j);
				pixel r = lreal->get_pixel(i, j);

				pixel confidence = real->exp().confidence(s);

				if (!s.finite()
				 || !r.finite()
				 || !confidence.finite())
					continue;

				if (real->get_bayer() != IMAGE_BAYER_NONE) {
					int color = ((image_bayer_ale_real *)real)->bayer_color(i, j);
					ssum[color] += confidence[color] * s[color];
					rsum[color] += confidence[color] * r[color];
				} else {
					ssum += confidence * s;
					rsum += confidence * r;
				}

			}

			ec = ssum / rsum;

#endif

			if (ec.finite() 
			 && rsum[0] > 0.001
			 && rsum[1] > 0.001
			 && rsum[2] > 0.001)
				real->exp().set_multiplier(
					real->exp().get_multiplier() * ec);
		}

		args.lreal = (d2::image *) lreal;

		correct_linear cl(args);
		cl.run();

		if (nlsimulated)
			delete lreal;
	}

	/* 
	 * Adjust correction array C based on the difference between the
	 * simulated projected frame and real frame M.  Update the correction
	 * count CC for affected pixels in C.
	 */

	virtual void _ip_frame(int frame_num, correction_t *cu, const image *real, 
			transformation t, const raster *f, const backprojector *b,
			const raster *nlf, const backprojector *nlb) {

		ui::get()->d2_irani_peleg_start();

		/*
		 * Initialize simulated data structures
		 */

                image *lsimulated;
		
		if (real->get_bayer() == IMAGE_BAYER_NONE) {
			lsimulated = new image_ale_real(
				real->height(),
				real->width(), 
				real->depth());
		} else {
			lsimulated = new image_bayer_ale_real(
				real->height(),
				real->width(),
				real->depth(),
				real->get_bayer());
		}

                image *nlsimulated = NULL;

		if (nlf)
			nlsimulated = new image_ale_real(
				real->height(),
				real->width(), 
				real->depth());

		/*
		 * Create simulated frames with forward projection.
		 */

		ui::get()->ip_frame_simulate_start();

		point extents[2] = { point::posinf(), point::neginf() };

		_ip_frame_simulate(frame_num, approximation, lsimulated, nlsimulated, t, f, nlf, real->exp(), extents);

		/*
		 * Update the correction array using backprojection.
		 */

		ui::get()->ip_frame_correct_start();
		_ip_frame_correct(frame_num, approximation, cu, real, lsimulated, nlsimulated, t, b, nlb, extents);

		/*
		 * Finalize data structures.
		 */

                delete lsimulated;
		delete nlsimulated;

		ui::get()->d2_irani_peleg_stop();
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

			correction_t *correction = new correction_t(
					use_weighted_median,
					weight_limit,
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

                                _ip_frame(m, correction, real,
					t, f[m], b[m], nlf, nlb);

				image_rw::close(m);

				correction->frame_end(m);
			}

			/*
			 * Update the approximation.
			 */

			ui::get()->ip_update();

			for (unsigned int i = 0; i < approximation->height(); i++)
			for (unsigned int j = 0; j < approximation->width();  j++) {

				pixel  cpix = correction->get_correction(i, j);
				pixel ccpix = correction->get_count(i, j);
				pixel  apix = approximation->get_pixel(i, j);

				for (unsigned int k = 0; k < 3; k++) {

					const ale_real cc_floor = 1e-20;

					if (ccpix[k] < cc_floor)
						continue;

					if (!finite(cpix[k]))
						continue;

					if (!finite(apix[k]))
						continue;

					ale_real new_value = cpix[k] + apix[k];

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

        ipc(render *input, unsigned int iterations, int _inc, psf *lresponse, 
			psf *nlresponse, int exposure_register,
			int use_weighted_median, double ipwl) {
                this->input = input;
                synced = 0;
		inc = _inc;
                this->iterations = iterations;
		this->lresponse = lresponse;
		this->nlresponse = nlresponse;
		this->exposure_register = exposure_register;
		this->use_weighted_median = use_weighted_median;
		this->weight_limit = ipwl;
        }

        const image *get_image() const {
                if (synced)
                        return approximation;
                else
                        return input->get_image();
        }

        const image *get_defined() const {
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
                synced = 1;
                approximation = optimizations::get_ip_working_image(input->get_image());
		ui::get()->ip_start();
                _ip();
		ui::get()->ip_done();

		/*
		 * Since we write output internally for --inc, no external
		 * update is necessary in this case.
		 */

                return 0;
        }

	virtual ~ipc() {
	}

	void free_memory() {
	}
};

#endif
