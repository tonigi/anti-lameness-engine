// Copyright 2002 David Hilvert <dhilvert@auricle.dyndns.org>, 
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
 * gpt.h: Represent transformations of projective or Euclidean type.
 */

#ifndef __gpt_h__
#define __gpt_h__

#include "image.h"
#include "point.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Structure to describe a transformation of euclidean or projective kind.
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
 * and '1', respectively.  'x0' is denoted 'x[0][0]'; 'y0' is 'x[0][1]'.
 *
 * eu[i] are the parameters for euclidean transformations.
 *
 * We consider points to be transformed as homogeneous coordinate vectors
 * multiplied on the right of the transformation matrix, and so we consider the
 * transformation matrix as
 *
 *	-       -
 * 	| a b c |
 *	| d e f |
 * 	| g h i |
 *	-       -
 *
 * where element i is always equal to 1.
 *
 */
struct transformation {
private:
	ale_pos input_height, input_width;
	point x[4];
	ale_pos eu[3];
	ale_pos a, b, c, d, e, f, g, h;			// matrix
	ale_pos _a, _b, _c, _d, _e, _f, _g, _h;		// matrix inverse
	int is_projective;
	int resultant_memo;
	int resultant_inverse_memo;
	ale_pos scale_factor;

	/*
	 * Calculate resultant matrix values.
	 */
	void resultant() {

		/*
		 * If we already know the answers, don't bother calculating 
		 * them again.
		 */

		if (resultant_memo)
			return;

		if (is_projective) {

			/*
			 * Calculate resultant matrix values for a general
			 * projective transformation given that we are mapping
			 * from the source domain of dimension (input_height *
			 * input_width) to a specified arbitrary quadrilateral.
			 * Follow the calculations outlined in the document by
			 * Paul Heckbert cited above for the case in which the
			 * source domain is a unit square and then divide to
			 * correct for the scale factor in each dimension.
			 */

			/*
			 * First, perform calculations as outlined in Heckbert.
			 */

			ale_pos delta_01 = x[1][0] - x[2][0];
			ale_pos delta_02 = x[3][0] - x[2][0];
			ale_pos sigma_0  = x[0][0] - x[1][0] + x[2][0] - x[3][0];
			ale_pos delta_11 = x[1][1] - x[2][1];                   
			ale_pos delta_12 = x[3][1] - x[2][1];                   
			ale_pos sigma_1  = x[0][1] - x[1][1] + x[2][1] - x[3][1];

			g = (sigma_0  * delta_12 - sigma_1  * delta_02)
			  / (delta_01 * delta_12 - delta_11 * delta_02);

			h = (delta_01 * sigma_1  - delta_11 * sigma_0 )
			  / (delta_01 * delta_12 - delta_11 * delta_02);

			a = (x[1][0] - x[0][0] + g * x[1][0]);
			b = (x[3][0] - x[0][0] + h * x[3][0]);
			c =  x[0][0];                       
							    
			d = (x[1][1] - x[0][1] + g * x[1][1]);
			e = (x[3][1] - x[0][1] + h * x[3][1]);
			f =  x[0][1];

			/* 
			 * Finish by scaling so that our transformation maps
			 * from a rectangle of width and height matching the
			 * width and height of the input image.
			 */

			a /= input_height;
			b /= input_width;
			d /= input_height;
			e /= input_width;
			g /= input_height;
			h /= input_width;

		} else {

			/*
			 * Calculate matrix values for a euclidean
			 * transformation.  
			 *
			 * We want to translate the image center by (eu[0],
			 * eu[1]) and rotate the image about the center by
			 * eu[2] degrees.  This is equivalent to the following
			 * sequence of affine transformations applied to the
			 * point to be transformed:
			 *
			 * translate by (-h/2, -w/2)
			 * rotate    by eu[2] degrees about the origin
			 * translate by (h/2, w/2)
			 * translate by (eu[0], eu[1])
			 *
			 * The matrix assigned below represents the result of
			 * combining all of these transformations.  Matrix
			 * elements g and h are always zero in an affine
			 * transformation.
			 */ 
			
			ale_pos theta = eu[2] * M_PI / 180;

			a = cos(theta);
			b = sin(theta);
			c = 0.5 * (input_height * (1 - a) - input_width * b) + eu[0];
			d = -b;
			e = a;
			f = 0.5 * (input_height * b + input_width * (1 - a)) + eu[1];
			g = 0;
			h = 0;
		} 

		resultant_memo = 1;
	}

	/*
	 * Calculate the inverse transform matrix values.
	 */
	void resultant_inverse () {

		/*
		 * If we already know the answers, don't bother calculating 
		 * them again.
		 */

		if (resultant_inverse_memo)
			return;

		if (is_projective) {

			resultant();

			/*
			 * For projective transformations, we calculate
			 * the inverse of the forward transformation 
			 * matrix.
			 */

			ale_pos scale = a * e - b * d;

			_a = (e * 1 - f * h) / scale;
			_b = (h * c - 1 * b) / scale;
			_c = (b * f - c * e) / scale;
			_d = (f * g - d * 1) / scale;
			_e = (1 * a - g * c) / scale;
			_f = (c * d - a * f) / scale;
			_g = (d * h - e * g) / scale;
			_h = (g * b - h * a) / scale;

		} else {

			/*
			 * In cursory trials, using an explicit matrix inverse
			 * with Euclidean transformations can result in
			 * different (and at least in some cases worse)
			 * alignment characteristics than the method employed
			 * below, which doesn't appear to yield appreciably
			 * worse performance than an explicit matrix inverse.
			 *
			 * For Euclidean transformations, we invert the
			 * constituent transformations of the forward
			 * transformation and apply them in reverse order:
			 *
			 * translate by (-eu[0], -eu[1]) 
			 * translate by (-h/2, -w/2) 
			 * rotate    by -eu[2] degrees about the origin
			 * translate by (h/2, w/2)
			 *
			 * The matrix assigned below represents the result of
			 * combining all of these transformations.  Matrix
			 * elements g and h are always zero in an affine
			 * transformation.
			 */ 

			if (!resultant_memo) {
				ale_pos theta = eu[2] * M_PI / 180;

				a = cos(theta);
				b = sin(theta);
			}

			_a = a;
			_b = -b;
			_c = (-eu[0] - input_height/2) * a 
			   - (-eu[1] - input_width/2) * b
			   + input_height/2;
			_d = b;
			_e = a;
			_f = (-eu[0] - input_height/2) * b 
			   + (-eu[1] - input_width/2) * a
			   + input_width/2;
			_g = 0;
			_h = 0;

		}

		resultant_inverse_memo = 1;
	}

public:	
	/*
	 * Get width of input image.
	 */
	ale_pos width() {
		return input_width;
	}

	/*
	 * Get unscaled width of input image.
	 */
	unsigned int unscaled_width() {
		return (unsigned int) (input_width / scale_factor);
	}

	/*
	 * Get height of input image;
	 */
	ale_pos height() {
	       return input_height;
	}

	/*
	 * Get unscaled height of input image.
	 */
	unsigned int unscaled_height() {
		return (unsigned int) (input_height / scale_factor);
	}

	/*
	 * Transform point p.
	 */
	struct point transform(struct point p) {
		struct point result;

		resultant();

		result[0] = (a * p[0] + b * p[1] + c)
			  / (g * p[0] + h * p[1] + 1);
		result[1] = (d * p[0] + e * p[1] + f)
			  / (g * p[0] + h * p[1] + 1);

		return result;
	}

	/*
	 * operator() is the transformation operator.
	 */
	struct point operator()(struct point p) {
		return transform(p);
	}

	/*
	 * Map point p using the inverse of the transform.
	 */
	struct point inverse_transform(struct point p) {
		struct point result;

		resultant_inverse();

		result[0] = (_a * p[0] + _b * p[1] + _c)
			  / (_g * p[0] + _h * p[1] +  1);
		result[1] = (_d * p[0] + _e * p[1] + _f)
			  / (_g * p[0] + _h * p[1] +  1);

		return result;
	}

	/*
	 * Map unscaled point p.
	 */
	struct point transform_unscaled(struct point p) {
		p[0] *= scale_factor;
		p[1] *= scale_factor;

		return transform(p);
	}

	/*
	 * Map point p using the inverse of the transform into
	 * the unscaled image space.
	 */
	struct point unscaled_inverse_transform(struct point p) {
		point q = inverse_transform(p);

		q[0] /= scale_factor;
		q[1] /= scale_factor;

		return q;
	}
	
	/*
	 * Calculate projective transformation parameters from a euclidean
	 * transformation.
	 */
	void eu_to_gpt() {

		assert(!is_projective);

		x[0] = transform( point(      0      ,      0      ) );
		x[1] = transform( point( input_height,      0      ) );
		x[2] = transform( point( input_height, input_width ) );
		x[3] = transform( point(      0      , input_width ) );

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		is_projective = 1;
	}

	/*
	 * Calculate euclidean identity transform for a given image.
	 */
	static struct transformation eu_identity(const image *i = NULL, ale_pos scale_factor = 1) {
		struct transformation r;

		r.resultant_memo = 0;
		r.resultant_inverse_memo = 0;

		r.eu[0] = 0;
		r.eu[1] = 0;
		r.eu[2] = 0;
		r.input_width = i ? i->width() * scale_factor : 2;
		r.input_height = i ? i->height() * scale_factor : 2;
		r.scale_factor = scale_factor;

		r.is_projective = 0;

		return r;
	}

	/*
	 * Calculate projective identity transform for a given image.
	 */
	static transformation gpt_identity(const image *i, ale_pos scale_factor) {
		struct transformation r = eu_identity(i, scale_factor);

		r.eu_to_gpt();

		return r;
	}

	/*
	 * Modify a euclidean transform in the indicated manner.
	 */
	void eu_modify(int i1, ale_pos diff) {

		assert(!is_projective);

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		eu[i1] += diff;

	}

	/*
	 * Modify all euclidean parameters at once.
	 */
	void eu_set(ale_pos eu[3]) {
		
		resultant_memo = 0;
		resultant_inverse_memo = 0;

		for (int i = 0; i < 3; i++)
			this->eu[i] = eu[i];

		if (is_projective) {
			is_projective = 0;
			eu_to_gpt();
		}

	}

	/*
	 * Get the specified euclidean parameter
	 */
	ale_pos eu_get(int param) {
		assert (!is_projective);
		assert (param >= 0);
		assert (param < 3);
		return eu[param];
	}

	/*
	 * Modify a projective transform in the indicated manner.
	 */
	void gpt_modify(int i1, int i2, ale_pos diff) {

		assert (is_projective);

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		x[i2][i1] += diff;

	}

	/*
	 * Modify all projective parameters at once.
	 */
	void gpt_set(point x[4]) {

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		is_projective = 1;

		for (int i = 0; i < 4; i++)
			this->x[i] = x[i];

	}

	/*
	 * Get the specified projective parameter
	 */
	point gpt_get(int point) {
		assert (is_projective);
		assert (point >= 0);
		assert (point <  4);

		return x[point];
	}

	/*
	 * Get the specified projective parameter
	 */
	ale_pos gpt_get(int point, int dim) {
		assert (is_projective);
		assert (dim   >= 0);
		assert (dim   <  2);

		return gpt_get(point)[dim];
	}

	/*
	 * Rescale a transform with a given factor.
	 */
	void rescale(ale_pos factor) {

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		if (is_projective) {

			for (int i = 0; i < 4; i++)
			for (int j = 0; j < 2; j++)
				x[i][j] *= factor;

		} else {

			for (int i = 0; i < 2; i++)
				eu[i] *= factor;

		}

		input_width  *= factor;
		input_height *= factor;
		scale_factor *= factor;
	}

	/*
	 * Set the dimensions of the image.
	 */
	void set_dimensions(const image *im) {

		int new_height = (int) (im->height() * scale_factor);
		int new_width  = (int) (im->width() * scale_factor);

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		if (is_projective) {

			/*
			 * If P(w, x, y, z) is a projective transform mapping
			 * the corners of the unit square to points w, x, y, z,
			 * and Q(w, x, y, z)(i, j) == P(w, x, y, z)(ai, bj),
			 * then we have:
			 *
			 * Q(w, x, y, z) == P(  P(w, x, y, z)(0, 0),
			 *			P(w, x, y, z)(a, 0),
			 *			P(w, x, y, z)(a, b),
			 *			P(w, x, y, z)(0, b)  )
			 *
			 * If we note that P(w, x, y, z)(0, 0) == w, we can
			 * omit a calculation.  
			 *
			 * We take 'a' as the ratio (new_height /
			 * old_height) and 'b' as the ratio (new_width /
			 * old_width) if we want the common upper left-hand
			 * region of both new and old images to map to the same
			 * area.
			 *
			 * Since we're not mapping from the unit square, we
			 * take 'a' as new_height and 'b' as new_width to
			 * accommodate the existing scale factor.
			 */

			point _x, _y, _z;

			_x = transform(point(new_height,     0    ));
			_y = transform(point(new_height, new_width));
			_z = transform(point(    0     , new_width));

			x[1] = _x;
			x[2] = _y;
			x[3] = _z;

		}

		input_height = new_height;
		input_width  = new_width;
	} 

	/*
	 * Get the position and dimensions of a pixel P mapped from one
	 * coordinate system to another, using the forward transformation.  
	 * This function uses scaled input coordinates.
	 */
	void map_area(point p, point *q, ale_pos d[2]) {

		/*
		 * Determine the coordinates in the target frame for the source
		 * image pixel P and two adjacent source pixels.
		 */

		    (*q) = transform(p);
		point q0 = transform(point(p[0] + 1, p[1]));
		point q1 = transform(point(p[0], p[1] + 1));

		/*
		 * Calculate the distance between source image pixel and
		 * adjacent source pixels, measured in the coordinate system of
		 * the target frame.
		 */

		ale_pos ui = fabs(q0[0] - (*q)[0]);
		ale_pos uj = fabs(q0[1] - (*q)[1]);
		ale_pos vi = fabs(q1[0] - (*q)[0]);
		ale_pos vj = fabs(q1[1] - (*q)[1]);

		/*
		 * We map the area of the source image pixel P onto the target
		 * frame as a rectangular area oriented on the target frame's
		 * axes.  Note that this results in an area that may be the
		 * wrong shape or orientation.
		 *
		 * We define two estimates of the rectangle's dimensions below.
		 * For rotations of 0, 90, 180, or 270 degrees, max and sum are
		 * identical.  For other orientations, sum is too large and max
		 * is too small.  We use the mean of max and sum, which we then
		 * divide by two to obtain the distance between the center and
		 * the edge.
		 */

		ale_pos maxi = (ui > vi) ? ui : vi;
		ale_pos maxj = (uj > vj) ? uj : vj;
		ale_pos sumi = ui + vi;
		ale_pos sumj = uj + vj;

		d[0] = (maxi + sumi) / 4;
		d[1] = (maxj + sumj) / 4;
	}

	/*
	 * Get the position and dimensions of a pixel P mapped from one
	 * coordinate system to another, using the forward transformation.  If
	 * SCALE_FACTOR is not equal to one, multiply by the scale factor to
	 * obtain scaled coordinates.
	 */
	void map_area_unscaled(point p, point *q, ale_pos d[2]) {

		/*
		 * Apply any necessary scaling
		 */

		if (scale_factor != 1) {
			p[0] *= scale_factor;
			p[1] *= scale_factor;
		}

		map_area(p, q, d);
	}

	/*
	 * Get the position and dimensions of a pixel P mapped from one
	 * coordinate system to another, using the inverse transformation.  If
	 * SCALE_FACTOR is not equal to one, divide out the scale factor to
	 * obtain unscaled coordinates.  This method is very similar to the 
	 * map_area method above.
	 */
	void unscaled_map_area_inverse(point p, point *q, ale_pos d[2]) {

		/*
		 * Determine the coordinates in the target frame for the source
		 * image pixel P and two adjacent source pixels.
		 */

		    (*q) = inverse_transform(p);
		point q0 = inverse_transform(point(p[0] + 1, p[1]));
		point q1 = inverse_transform(point(p[0], p[1] + 1));


		/*
		 * Calculate the distance between source image pixel and
		 * adjacent source pixels, measured in the coordinate system of
		 * the target frame.
		 */

		ale_pos ui = fabs(q0[0] - (*q)[0]);
		ale_pos uj = fabs(q0[1] - (*q)[1]);
		ale_pos vi = fabs(q1[0] - (*q)[0]);
		ale_pos vj = fabs(q1[1] - (*q)[1]);

		/*
		 * We map the area of the source image pixel P onto the target
		 * frame as a rectangular area oriented on the target frame's
		 * axes.  Note that this results in an area that may be the
		 * wrong shape or orientation.
		 *
		 * We define two estimates of the rectangle's dimensions below.
		 * For rotations of 0, 90, 180, or 270 degrees, max and sum are
		 * identical.  For other orientations, sum is too large and max
		 * is too small.  We use the mean of max and sum, which we then
		 * divide by two to obtain the distance between the center and
		 * the edge.
		 */

		ale_pos maxi = (ui > vi) ? ui : vi;
		ale_pos maxj = (uj > vj) ? uj : vj;
		ale_pos sumi = ui + vi;
		ale_pos sumj = uj + vj;

		d[0] = (maxi + sumi) / 4;
		d[1] = (maxj + sumj) / 4;

		if (scale_factor != 1) {
			d[0] /= scale_factor;
			d[1] /= scale_factor;
			(*q)[0] /= scale_factor;
			(*q)[1] /= scale_factor;
		}
	}

	/*
	 * Modify all projective parameters at once.  Accommodate bugs in the
	 * version 0 transformation file handler (ALE versions 0.4.0p1 and
	 * earlier).  This code is only called when using a transformation data
	 * file created with an old version of ALE.
	 */
	void gpt_v0_set(point x[4]) {

		is_projective = 1;

		/*
		 * This is slightly modified code from version
		 * 0.4.0p1.
		 */

		ale_pos delta_01 = x[1][0] - x[2][0];
		ale_pos delta_02 = x[3][0] - x[2][0];
		ale_pos sigma_0  = x[0][0] - x[1][0] + x[2][0] - x[3][0];
		ale_pos delta_11 = x[1][1] - x[2][1];                   
		ale_pos delta_12 = x[3][1] - x[2][1];                   
		ale_pos sigma_1  = x[0][1] - x[1][1] + x[2][1] - x[3][1];

		g = (sigma_0  * delta_12 - sigma_1  * delta_02)
		  / (delta_01 * delta_12 - delta_11 * delta_02)
		  / input_width;

		h = (delta_01 * sigma_1  - delta_11 * sigma_0 )
		  / (delta_01 * delta_12 - delta_11 * delta_02)
		  / input_height;

		a = (x[1][0] - x[0][0] + g * x[1][0])
		  / input_width;                       
		b = (x[3][0] - x[0][0] + h * x[3][0])
		  / input_height;
		c = x[0][0];

		d = (x[1][1] - x[0][1] + g * x[1][1])
		  /  input_width;                       
		e = (x[3][1] - x[0][1] + h * x[3][1])
		  /  input_height;
		f =  x[0][1];

		resultant_memo = 1;
		resultant_inverse_memo = 0;

		this->x[0] = inverse_transform( point(      0      ,      0      ) );
		this->x[1] = inverse_transform( point( input_height,      0      ) );
		this->x[2] = inverse_transform( point( input_height, input_width ) );
		this->x[3] = inverse_transform( point(      0      , input_width ) );

		resultant_memo = 0;
		resultant_inverse_memo = 0;
	}

	/*
	 * Modify all euclidean parameters at once.  Accommodate bugs in the
	 * version 0 transformation file handler (ALE versions 0.4.0p1 and 
	 * earlier).  This code is only called when using a transformation data
	 * file created with an old version of ALE.
	 */
	void eu_v0_set(ale_pos eu[3]) {

		/*
		 * This is slightly modified code from version
		 * 0.4.0p1.
		 */

		int i;
		
		x[0][0] = 0;              x[0][1] = 0;
		x[1][0] = input_width;    x[1][1] = 0;
		x[2][0] = input_width;    x[2][1] = input_height;
		x[3][0] = 0;              x[3][1] = input_height;

		/*
		 * Rotate
		 */

		ale_pos theta = eu[2] * M_PI / 180;

		for (i = 0; i < 4; i++) {
			ale_pos _x[2];

			_x[0] = (x[i][0] - input_width/2)  * cos(theta)
			      + (x[i][1] - input_height/2) * sin(theta)
			      +  input_width/2;
			_x[1] = (x[i][1] - input_height/2) * cos(theta)
			      - (x[i][0] - input_width/2)  * sin(theta)
			      +  input_height/2;
			
			x[i][0] = _x[0];
			x[i][1] = _x[1];
		}

		/*
		 * Translate
		 */

		for (i = 0; i < 4; i++) {
			x[i][0] += eu[0];
			x[i][1] += eu[1];
		}

		if (is_projective) {
			gpt_v0_set(x);
			return;
		}

		/*
		 * Reconstruct euclidean parameters
		 */
			
		gpt_v0_set(x);

		point center(input_height / 2, input_width / 2);
		point center_image = transform(center);

		this->eu[0] = center_image[0] - center[0];
		this->eu[1] = center_image[1] - center[1];

		point center_left(input_height / 2, 0);
		point center_left_image = transform(center_left);

		ale_pos displacement = center_image[0] - center_left_image[0];

		this->eu[2] = asin(2 * displacement / input_width) / M_PI * 180;

		if (center_left_image[1] > center_image[1])
			this->eu[2] = this->eu[2] + 180;

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		is_projective = 0;

	}

};

#endif
