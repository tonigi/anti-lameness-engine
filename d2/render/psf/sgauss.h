// code by HJ Hornbeck, based on code copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//				    <dhilvert@ugcs.caltech.edu>

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

#ifndef __psf_sgauss_h__
#define __psf_sgauss_h__

#include "../../point.h"
#include "psf.h"

/*
 * A "Gaussian" point-spread function. It's aimed at duplicating some
 *  of the blurring inherent in a poor quality lens. "Gaussian" is in
 *  quotes because a true Gaussian curve is only at maximum at one
 *  point, whereas this filter can expand that to a disc; handy for
 *  duplicating some forms of blur. Thus the name: "sgauss".
 */

class sgauss : public psf {
       ale_pos sigma;		  // in pixels per standard deviation
       ale_pos centre_radius;	  // in pixels
       static const ale_pos cutoff = 2;	// standard deviations before we cut off
public:

       /*
	* The following four functions indicate filter boundaries.  Filter
	* support may include everything up to and including the boundaries
	* specified here.
	*/
       inline float max_i() const { return sigma * cutoff + centre_radius; }
       inline float min_i() const { return -max_i(); } // too lazy to cut/paste the above ;)
       inline float min_j() const { return -max_i(); }
       inline float max_j() const { return max_i(); }

       /*
	* Response function
	*
	* Get the response to the rectangle bounded by (top, bot, lef, rig).
	* This function must correctly handle points which fall outside of the
	* filter support.  The variety of the responding pixel is provided, in
	* case response is not uniform for all pixels (e.g. some sensor arrays
	* stagger red, green, and blue sensors).
	*/
       psf_result operator()(float top, float bot, float lef, float rig,
		       unsigned int variety) const {

	       psf_result result;

	       // pre-fill the colour result matrix
	       for (int k = 0; k < 3; k++)
		       result.matrix(k, k) = 0;

	       // calculate some needed values
	       ale_pos radius = sigma * cutoff + centre_radius;
	       ale_pos area_premult = (bot - top) * (rig - lef) / 25;
	       float vert_step = (bot - top) / 4;
	       float horiz_step = (rig - lef) / 4;
	       float total = 0;


	       // determine the final value by simple sampling:
	       for (float i = top; i <= bot; i += vert_step)
	       for (float j = lef; j <= rig; j += horiz_step) {

		       // calculate radius for given sample
		       ale_pos r = sqrt( i*i + j*j );
		       if ( r <= centre_radius )
			       total += 1;
		       else if ( r < radius ) // calculate gaussian falloff
			       total += exp( -(r-centre_radius)*(r-centre_radius) / (sigma*sigma) ) ;
		       // outside our radius? must be 0...
		       }

	       // adjust for point sampling and area
	       total *= area_premult;

	       // fill in the results
	       for (int k = 0; k < 3; k++)
		       result.matrix(k, k) = total;
	      
	       return result;
       }

       sgauss(ale_pos sig, ale_pos rad) {
	       sigma = sig;
	       centre_radius = rad;
       }
};

#endif
