// Copyright 2002 David Hilvert <dhilvert@ugcs.caltech.edu>

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

#ifndef __gpt_h__
#define __gpt_h__

#include "my_real.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Structure to describe a general projective transformation
 *
 * Member names roughly correspond to a typical treatment of projective
 * transformations from:
 *
 *	Heckbert, Paul.  "Projective Mappings for Image Warping."  Excerpted 
 *		from his Master's Thesis (UC Berkeley, 1989).  1995.
 *
 * http://www.cs.cmu.edu/afs/cs/project/classes-ph/862.95/www/notes/proj.ps
 *
 * For convenience, Heckbert's 'x' and 'y' are noted here numerically by '0'
 * and '1', respectively.  Thus, 'x0' is denoted 'x[0][0]'; 'y0' is 'x[1][0]'.
 *
 * eu[i] are the parameters for euclidean transformations.
 *
 */
struct gpt {
	my_real input_width, input_height;

	my_real x[2][4];

	my_real eu[3];

	my_real a, b, c, d, e, f, g, h;
};

/*
 * Calculate resultant values for a general projective transformation given
 * that we are transforming from a unit square to a specified arbitrary
 * quadrilateral.  Follow the calculations outlined in the document by Paul
 * Heckbert cited above.
 */
inline struct gpt gpt_resultant(struct gpt io) {
	my_real delta_01 = io.x[0][1] - io.x[0][2];
	my_real delta_02 = io.x[0][3] - io.x[0][2];
	my_real sigma_0  = io.x[0][0] - io.x[0][1] + io.x[0][2] - io.x[0][3];
	my_real delta_11 = io.x[1][1] - io.x[1][2];
	my_real delta_12 = io.x[1][3] - io.x[1][2];
	my_real sigma_1  = io.x[1][0] - io.x[1][1] + io.x[1][2] - io.x[1][3];

	io.g = (sigma_0  * delta_12 - sigma_1  * delta_02)
	     / (delta_01 * delta_12 - delta_11 * delta_02)
	     / io.input_width;

	io.h = (delta_01 * sigma_1  - delta_11 * sigma_0 )
	     / (delta_01 * delta_12 - delta_11 * delta_02)
	     / io.input_height;

	io.a = (io.x[0][1] - io.x[0][0] + io.g * io.x[0][1])
	     / io.input_width;
	io.b = (io.x[0][3] - io.x[0][0] + io.h * io.x[0][3])
	     / io.input_height;
	io.c = io.x[0][0];

	io.d = (io.x[1][1] - io.x[1][0] + io.g * io.x[1][1])
	     / io.input_width;
	io.e = (io.x[1][3] - io.x[1][0] + io.h * io.x[1][3])
	     / io.input_height;
	io.f = io.x[1][0];

	return io;
}

/*
 * Calculate resultant values for a euclidean transformation.
 */
inline struct gpt eu_resultant(struct gpt io) {
	int i;
	
	io.x[0][0] = 0;              io.x[1][0] = 0;
	io.x[0][1] = io.input_width; io.x[1][1] = 0;
	io.x[0][2] = io.input_width; io.x[1][2] = io.input_height;
	io.x[0][3] = 0;              io.x[1][3] = io.input_height;

	/*
	 * Rotate
	 */

	{
		my_real theta = io.eu[2] * M_PI / 180;

		for (i = 0; i < 4; i++) {
			int x[2];

			x[0] = (io.x[0][i] - io.input_width/2)  * cos(theta)
			     + (io.x[1][i] - io.input_height/2) * sin(theta)
			     + io.input_width/2;
			x[1] = (io.x[1][i] - io.input_height/2) * cos(theta)
			     - (io.x[0][i] - io.input_width/2)  * sin(theta)
			     + io.input_height/2;
			
			io.x[0][i] = x[0];
			io.x[1][i] = x[1];
		}
	}

	/*
	 * Translate
	 */

	for (i = 0; i < 4; i++) {
		io.x[0][i] += io.eu[0];
		io.x[1][i] += io.eu[1];
	}

	return gpt_resultant(io);
}

/*
 * Calculate a euclidean transform modified in the indicated manner.
 */
inline struct gpt eu_modify(struct gpt io, int i1, my_real diff) {
	io.eu[i1] += diff;
	return io;
}

/*
 * Calculate a general projective transform modified in the indicated manner.
 */
inline struct gpt gpt_modify (struct gpt io, int i1, int i2, my_real diff) {
	io.x[i1][i2] += diff;
	return io;
}


#endif
