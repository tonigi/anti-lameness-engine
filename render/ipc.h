// Copyright 2003 David Hilvert <dhilvert@ugcs.caltech.edu>

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
 * projection and backprojection function used is configurable by template.
 * For more information on the algorithm used, see:
 *
 * 	http://www.wisdom.weizmann.ac.il/~irani/abstracts/superResolution.html
 */

#ifndef __ipc_h__
#define __ipc_h__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "../image.h"
#include "../render.h"

template <class Response>
class ipc : public render {

        int done;
	int inc;
        image *done_image;
        render *input;
        double scale_factor;
        unsigned int iterations;

	Response response;

	/* 
	 * Adjust correction array C based on the difference between the
	 * simulated projected image and actual image M. 
	 */

        void _ip_frame(image_weights *c, int m) {

                transformation t = align::of(m);

                const image *im = image_rw::open(m);
                image *forward = image_rw::copy(m);
                image_weights *iw = new image_weights(
                                im->height(),
                                im->width(), 3);

		/*
		 * First, calculate the simulated projected image.
		 */

                for (unsigned int i = 0; i < done_image->height(); i++)
                for (unsigned int j = 0; j < done_image->width();  j++) {

			// fprintf(stderr, "(%d, %d)\n", i, j);
			
                        /*
                         * Here we more-or-less cut-and-paste-and-modify the
                         * drizzling code from drizzle.h.  Surely there's a more
                         * elegant way to do this.
                         */

                        point p  = point(i + im->offset()[0], j + im->offset()[1]);
                        point p0 = p, p1 = p;

                        p0[0]++;
                        p1[1]++;

                        point q  = t.inverse_transform(p);
                        point q0 = t.inverse_transform(p0);
                        point q1 = t.inverse_transform(p1);

                        double ui = fabs(q0[0] - q[0]);
                        double uj = fabs(q0[1] - q[1]);
                        double vi = fabs(q1[0] - q[0]);
                        double vj = fabs(q1[1] - q[1]);

                        double maxi = (ui > vi) ? ui : vi;
                        double maxj = (uj > vj) ? uj : vj;
                        double sumi = ui + vi;
                        double sumj = uj + vj;

                        double di = (maxi + sumi) / 4 / scale_factor;
                        double dj = (maxj + sumj) / 4 / scale_factor;

                        double ti = q[0] / scale_factor;
                        double tj = q[1] / scale_factor;

                        for (int ii = (int) floor(ti - di + response.min_i());
                                ii <= ceil(ti + di + response.max_i()); ii++)
                        for (int jj = (int) floor(tj - dj + response.min_j());
                                jj <= ceil(tj + dj + response.max_j()); jj++) {

                                my_real top = ti - di;
                                my_real bot = ti + di;
                                my_real lef = tj - dj;
                                my_real rig = tj + dj;

                                if (ii >= (int) 0
                                 && ii < (int) im->height()
                                 && jj >= (int) 0
                                 && jj < (int) im->width()) {

					class Response::ipc_result r = response(top - ii,
							                        bot - ii,
										lef - jj,
										rig - jj, i, j);

                                        for (int k_in = 0; k_in < 3; k_in++) {
					int k_out = k_in;

						if (r(k_in, k_out) != 0) {

							double w_old = iw->get_pixel_component(ii, jj, k_out);
							double w_new = w_old + fabs(r(k_in, k_out));

							iw->set_pixel_component(ii, jj, k_out, w_new);

							// if (k_out == 0)
								// fprintf(stderr, "(%d, %d, %.3lf, %.3lf)\n", 
								// 	ii, jj, (double) r(k_in, k_out), w_new);
							double new_value = 
								(w_old * forward->get_pixel_component(
										ii, jj, k_out)
							       + (r(k_in, k_out) 
								      * done_image->get_pixel_component(
										i, j, k_in)))
							     / w_new;

							// fprintf(stderr, "%lf\n", new_value);

							if (new_value > 255)
								new_value = 255;
							if (new_value < 0)
								new_value = 0;

							forward->set_pixel_component(ii, jj, k_out,
								(int) new_value);
						}
					}
                                }
                        }
                }

		// write_ppm("ale-internal-debug-ipc.ppm", forward);

		/*
		 * Now calculate the differences between the simulated
		 * image and the actual image, and add this difference, 
		 * for each pixel, to the corresponding pixel in the
		 * correction array C.
		 */

                for (unsigned int i = 0; i < done_image->height(); i++)
                for (unsigned int j = 0; j < done_image->width();  j++) {

			// fprintf(stderr, "[%d, %d]\n", i, j);

                        /*
                         * Here we more-or-less cut-and-paste-and-modify the
                         * code from above.  Surely there's a more elegant way
                         * to do this.
                         */

                        point p  = point(i + im->offset()[0], j + im->offset()[1]);
                        point p0 = p, p1 = p;

                        p0[0]++;
                        p1[1]++;

                        point q  = t.inverse_transform(p);
                        point q0 = t.inverse_transform(p0);
                        point q1 = t.inverse_transform(p1);

                        double ui = fabs(q0[0] - q[0]);
                        double uj = fabs(q0[1] - q[1]);
                        double vi = fabs(q1[0] - q[0]);
                        double vj = fabs(q1[1] - q[1]);

                        double maxi = (ui > vi) ? ui : vi;
                        double maxj = (uj > vj) ? uj : vj;
                        double sumi = ui + vi;
                        double sumj = uj + vj;

                        double di = (maxi + sumi) / 4 / scale_factor;
                        double dj = (maxj + sumj) / 4 / scale_factor;

                        double ti = q[0] / scale_factor;
                        double tj = q[1] / scale_factor;

                        for (int ii = (int) floor(ti - di - response.max_i());
                                ii <= ceil(ti + di - response.min_i()); ii++)
                        for (int jj = (int) floor(tj - dj - response.max_j());
                                jj <= ceil(tj + dj - response.min_j()); jj++) {

                                my_real top = ti - di;
                                my_real bot = ti + di;
                                my_real lef = tj - dj;
                                my_real rig = tj + dj;

                                if (ii >= (int) 0
                                 && ii < (int) im->height()
                                 && jj >= (int) 0
                                 && jj < (int) im->width()) {

					struct Response::ipc_result r = response(top - ii,
							                         bot - ii,
										 lef - jj,
										 rig - jj, i, j);

                                        for (int k_in = 0; k_in < 3; k_in++) {
					int k_out = k_in;

						if (r(k_in, k_out) != 0) {

							double w_old = c->get_pixel_component(i, j, k_in + 3);
							double w_new = w_old + fabs(r(k_in, k_out));

							// fprintf(stderr, "[%d, %d, %d, %lf]\n", 
							//         ii, jj, r(k_in, k_out), w_new);

							c->set_pixel_component(i, j, k_in + 3, w_new);

							c->set_pixel_component(i, j, k_in, 
								((w_old 
								 * c->get_pixel_component(i, j, k_in)
								+ r(k_in, k_out)
								  * (im->get_pixel_component(ii, jj, k_out)
								   - forward->get_pixel_component(ii, jj, k_out)))
								 / w_new));

							// fprintf(stderr, "[%lf]\n", c->get_pixel_component(
							//			i, j, k_in));
						}
					}
                                }
                        }
                }


                image_rw::close(m);
                delete iw;
                delete forward;
        }

	/*
	 * Iterate _ip_frame() over all frames, and update DONE_IMAGE after
	 * corrections from all frames have been summed.  Repeat for the number
	 * of iterations specified by the user.
	 */

        void _ip() {
                for (unsigned int n = 0; n < iterations; n++) {

                        image_weights *correction = new image_weights(
                                        done_image->height(),
                                        done_image->width(),
                                        done_image->depth() * 2);

                        for (int m = 0; m < image_rw::count(); m++)
                                _ip_frame(correction, m);


			for (unsigned int i = 0; i < done_image->height(); i++)
			for (unsigned int j = 0; j < done_image->width();  j++) 
			for (unsigned int k = 0; k < done_image->depth();  k++) {

				// if (k == 0)
					// fprintf(stderr, "(%d, %d, %lf)\n", i, j, 
					//		correction->get_pixel_component(
					//			i, j, 3));

				// fprintf(stderr, "[%lf]\n", correction->get_pixel_component(i, j, k));
				
				double new_value = 
					correction->get_pixel_component(
						i, j, k)
				      + done_image->get_pixel_component(
						i, j, k);

				if (new_value < 0)
					new_value = 0;
				else if (new_value > 255)
					new_value = 255;

				done_image->set_pixel_component(i, j, k, 
						(unsigned char) new_value);
			}

                        delete correction;

			if (inc)
				image_rw::output(done_image);

                        fprintf(stderr, ".");

                }
        }

public:

        ipc(render *input, double scale_factor, unsigned int iterations, int _inc) {
                this->input = input;
                done = 0;
		inc = _inc;
                this->scale_factor = scale_factor;
                this->iterations = iterations;
        }

        const image *get_image() {
                if (done)
                        return done_image;
                else
                        return input->get_image();
        }

        const image_weights *get_defined() {
                return input->get_defined();
        }

        void sync(int n) {
                input->sync(n);
        }

        int sync() {
		input->sync();
                fprintf(stderr, "Iterating Irani-Peleg");
                done = 1;
                done_image = input->get_image()->clone();
                _ip();

		fprintf(stderr, "\n");

                return 0;
        }

	Response *module() {
		return &response;
	}

	virtual ~ipc() {
	}

};

#endif
