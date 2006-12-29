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
 * transformation.h: Represent transformations of the kind q = c(b^-1(p)),
 * where p is a point in the source coordinate system, q is a point in the
 * target coordinate system, b^-1 is a transformation correcting barrel
 * distortion, and c is a transformation of projective or Euclidean type.
 * (Note that ^-1 in this context indicates the function inverse rather than
 * the exponential.)
 */

#ifndef __transformation_h__
#define __transformation_h__

#include "image.h"
#include "point.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Number of coefficients used in correcting barrel distortion.
 */

#define BARREL_DEGREE 5

/*
 * Acceptable error for inverse barrel distortion, measured in scaled output
 * pixels.
 */

#define BARREL_INV_ERROR 0.01

/*
 * transformation: a structure to describe a transformation of kind q =
 * c(b^-1(p)), where p is a point in the source coordinate system, q is a point
 * in the target coordinate system, b^-1 is a transformation correcting barrel
 * distortion, and c is a projective or Euclidean transformation.  (Note that
 * ^-1 in this case indicates a function inverse, not exponentiation.)  Data
 * elements are divided into those describing barrel distortion correction and
 * those describing projective/Euclidean transformations.
 *
 * Barrel distortion correction estimates barrel distortion using polynomial
 * functions of distance from the center of an image, following (roughly) the
 * example set by Helmut Dersch in his PanoTools software:
 *
 * 	http://www.path.unimelb.edu.au/~dersch/barrel/barrel.html
 *
 * Projective transformation data member names roughly correspond to a typical
 * treatment of projective transformations from:
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
	unsigned int input_height, input_width;
	point x[4];
	ale_pos eu[3];
	mutable ale_pos a, b, c, d, e, f, g, h;			// matrix
	mutable ale_pos _a, _b, _c, _d, _e, _f, _g, _h;		// matrix inverse
	ale_pos bdc[BARREL_DEGREE];				// barrel-dist. coeffs.
	unsigned int bdcnum;					// number of bdcs
	int _is_projective;
	mutable int resultant_memo;
	mutable int resultant_inverse_memo;
	ale_pos scale_factor;

	/*
	 * Calculate resultant matrix values.
	 */
	void resultant() const {

		/*
		 * If we already know the answers, don't bother calculating 
		 * them again.
		 */

		if (resultant_memo)
			return;

		if (_is_projective) {

			/*
			 * Calculate resultant matrix values for a general
			 * projective transformation given that we are mapping
			 * from the source domain of dimension input_height *
			 * input_width to a specified arbitrary quadrilateral.
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

			a = cos(theta) * scale_factor;
			b = sin(theta) * scale_factor;
			c = 0.5 * (input_height * (scale_factor - a) - input_width * b) + eu[0] * scale_factor;
			d = -b;
			e = a;
			f = 0.5 * (input_height * b + input_width * (scale_factor - a)) + eu[1] * scale_factor;
			g = 0;
			h = 0;
		} 

		resultant_memo = 1;
	}

	/*
	 * Calculate the inverse transform matrix values.
	 */
	void resultant_inverse () const {

		/*
		 * If we already know the answers, don't bother calculating 
		 * them again.
		 */

		if (resultant_inverse_memo)
			return;

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

		resultant_inverse_memo = 1;
	}

public:	

	/*
	 * Returns non-zero if the transformation might be non-Euclidean.
	 */
	int is_projective() {
		return _is_projective;
	}

	/*
	 * Get scale factor.
	 */

	ale_pos scale() const {
		return scale_factor;
	}

	/*
	 * Get width of input image.
	 */
	ale_pos scaled_width() const {
		return (input_width * scale_factor);
	}

	/*
	 * Get unscaled width of input image.
	 */
	unsigned int unscaled_width() const {
		return (unsigned int) input_width;
	}

	/*
	 * Get height of input image;
	 */
	ale_pos scaled_height() const {
	       return (input_height * scale_factor);
	}

	/*
	 * Get unscaled height of input image.
	 */
	unsigned int unscaled_height() const {
		return (unsigned int) input_height;
	}

	/*
	 * Barrel distortion radial component.
	 */
	ale_pos bdr(ale_pos r) const {
		ale_pos s = r;
		for (unsigned int d = 0; d < bdcnum; d++)
			s += bdc[d] * (pow(r, d + 2) - r);
		return s;
	}

	/*
	 * Derivative of the barrel distortion radial component.
	 */
	ale_pos bdrd(ale_pos r) const {
		ale_pos s = 1;
		for (unsigned int d = 0; d < bdcnum; d++)
			s += bdc[d] * (pow(r, d + 1) - 1);
		return s;
	}

	/*
	 * Barrel distortion.
	 */
	struct point bd(struct point p) const {
		if (bdcnum > 0) {
			point half_diag = point(unscaled_height(), unscaled_width()) / 2;

			p -= half_diag;

			ale_pos r = p.norm() / half_diag.norm();

			if (r > 0.00001)
				p *= bdr(r)/r;

			p += half_diag;
		}

		assert (!isnan(p[0]) && !isnan(p[1]));

		return p;
	}

	/*
	 * Barrel distortion inverse.
	 */
	struct point bdi(struct point p) const {
		if (bdcnum > 0) {
			point half_diag = point(unscaled_height(), unscaled_width()) / 2;

			p -= half_diag;

			ale_pos r = p.norm() / half_diag.norm();
			ale_pos s = r;

			while (fabs(r - bdr(s)) * half_diag.norm() > BARREL_INV_ERROR)
				s += (r - bdr(s)) / bdrd(s);

			if (r > 0.0001)
				p *= s / r;

			p += half_diag;
		}

		assert (!isnan(p[0]) && !isnan(p[1]));
		
		return p;
	}

	/*
	 * Projective/Euclidean transformation
	 */
	struct point pe(struct point p) const {
		struct point result;

		resultant();

		result[0] = (a * p[0] + b * p[1] + c)
			  / (g * p[0] + h * p[1] + 1);
		result[1] = (d * p[0] + e * p[1] + f)
			  / (g * p[0] + h * p[1] + 1);

		return result;
	}

	/*
	 * Projective/Euclidean inverse
	 */
	struct point pei(struct point p) const {
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
	struct point transform_unscaled(struct point p) const {
		return pe(bdi(p));
	}

	/*
	 * Transform point p.
	 *
	 * Barrel distortion correction followed by a projective/euclidean
	 * transformation.
	 */
	struct point transform_scaled(struct point p) const {
		return transform_unscaled(p / scale_factor);
	}

#if 0
	/*
	 * operator() is the transformation operator.
	 */
	struct point operator()(struct point p) {
		return transform(p);
	}
#endif

	/*
	 * Map point p using the inverse of the transform into
	 * the unscaled image space.
	 */
	struct point unscaled_inverse_transform(struct point p) const {
		return bd(pei(p));
	}
	
	/*
	 * Map point p using the inverse of the transform.
	 *
	 * Projective/euclidean inverse followed by barrel distortion.
	 */
	struct point scaled_inverse_transform(struct point p) const {
		point q = unscaled_inverse_transform(p);

		q[0] *= scale_factor;
		q[1] *= scale_factor;

		return q;
	}

	/*
	 * Calculate projective transformation parameters from a euclidean
	 * transformation.
	 */
	void eu_to_gpt() {

		assert(!_is_projective);

		x[0] = transform_unscaled(point(      0      ,      0      ) );
		x[1] = transform_unscaled(point( input_height,      0      ) );
		x[2] = transform_unscaled(point( input_height, input_width ) );
		x[3] = transform_unscaled(point(      0      , input_width ) );

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		_is_projective = 1;
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
		r.input_width = i ? i->width() : 2;
		r.input_height = i ? i->height() : 2;
		r.scale_factor = scale_factor;

		r._is_projective = 0;

		r.bdcnum = 0;

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

		assert(!_is_projective);

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		if (i1 < 2)
			eu[i1] += diff / scale_factor;
		else
			eu[i1] += diff;

	}

	/*
	 * Rotate about a given point in the original reference frame.
	 */
	void eu_rotate_about_scaled(point center, ale_pos diff) {
		point fixpoint = scaled_inverse_transform(center);
		eu_modify(2, diff);
		point offset = center - transform_scaled(fixpoint);

		eu_modify(0, offset[0]);
		eu_modify(1, offset[1]);
	}

	/*
	 * Modify all euclidean parameters at once.
	 */
	void eu_set(ale_pos eu[3]) {
		
		resultant_memo = 0;
		resultant_inverse_memo = 0;

		this->eu[0] = eu[0] / scale_factor;
		this->eu[1] = eu[1] / scale_factor;
		this->eu[2] = eu[2];

		if (_is_projective) {
			_is_projective = 0;
			eu_to_gpt();
		}

	}

	/*
	 * Get the specified euclidean parameter
	 */
	ale_pos eu_get(int param) {
		assert (!_is_projective);
		assert (param >= 0);
		assert (param < 3);

		if (param < 2)
			return eu[param] * scale_factor;
		else
			return eu[param];
	}

	/*
	 * Modify a projective transform in the indicated manner.
	 */
	void gpt_modify(int i1, int i2, ale_pos diff) {

		assert (_is_projective);

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		x[i2][i1] += diff;

	}

	/*
	 * Modify a projective transform according to the group operation.
	 */
	void gr_modify(int i1, int i2, ale_pos diff) {
		assert (_is_projective);
		assert (i1 == 0 || i1 == 1);

		point diff_vector = (i1 == 0) 
			          ? point(diff, 0)
				  : point(0, diff);

		transformation t = *this;

		t.resultant_memo = 0;
		t.resultant_inverse_memo = 0;
		t.input_height = (unsigned int) scaled_height();
		t.input_width  = (unsigned int) scaled_width();
		t.scale_factor = 1;
		t.bdcnum = 0;

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		x[i2] = t.transform_scaled(t.scaled_inverse_transform(x[i2]) + diff_vector);
	}

	/*
	 * Modify all projective parameters at once.
	 */
	void gpt_set(point x[4]) {

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		_is_projective = 1;

		for (int i = 0; i < 4; i++)
			this->x[i] = x[i];

	}

	void gpt_set(point x1, point x2, point x3, point x4) {
		point x[4] = {x1, x2, x3, x4};
		gpt_set(x);
	}

	/*
	 * Get the specified projective parameter
	 */
	point gpt_get(int point) {
		assert (_is_projective);
		assert (point >= 0);
		assert (point <  4);

		return x[point];
	}

	/*
	 * Get the specified projective parameter
	 */
	ale_pos gpt_get(int point, int dim) {
		assert (_is_projective);
		assert (dim   >= 0);
		assert (dim   <  2);

		return gpt_get(point)[dim];
	}

	/*
	 * Translate by a given amount
	 */
	void translate(point p) {

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		if  (_is_projective)
		for (int i = 0; i < 4; i++)
			x[i]+=p;
		else {
			eu[0] += p[0] / scale_factor;
			eu[1] += p[1] / scale_factor;
		}
	}

	/*
	 * Rotate by a given amount about a given point.
	 */
	void rotate(point p, ale_pos degrees) {

		ale_pos radians = degrees * M_PI / (double) 180;

		if (_is_projective)
		for (int i = 0; i <= 4; i++) {
			x[i] -= p;
			x[i] = point(
				x[i][0] * cos(radians) + x[i][1] * sin(radians),
				x[i][1] * cos(radians) - x[i][0] * sin(radians));
			x[i] += p;
		} else {
			point usi = unscaled_inverse_transform(p);

			eu[2] += degrees;
			resultant_memo = 0;
			
			point p_rot = transform_unscaled(usi);

			eu[0] -= (p_rot[0] - p[0]) / scale_factor;
			eu[1] -= (p_rot[1] - p[1]) / scale_factor;
		}

		resultant_memo = 0;
		resultant_inverse_memo = 0;
	}

	/*
	 * Set the specified barrel distortion parameter.
	 */
	void bd_set(unsigned int degree, ale_pos value) {
		assert (degree < bdcnum);
		bdc[degree] = value;
	}

	/*
	 * Set all barrel distortion parameters.
	 */
	void bd_set(unsigned int degree, ale_pos values[BARREL_DEGREE]) {
		assert (degree <= BARREL_DEGREE);
		bdcnum = degree;
		for (unsigned int d = 0; d < degree; d++)
			bdc[d] = values[d];
	}

	/*
	 * Get the all barrel distortion parameters.
	 */
	void bd_get(ale_pos result[BARREL_DEGREE]) {
		for (unsigned int d = 0; d < bdcnum; d++)
			result[d] = bdc[d];
	}

	/*
	 * Get the specified barrel distortion parameter.
	 */
	ale_pos bd_get(unsigned int degree) {
		assert (degree < bdcnum);
		return bdc[degree];
	}

	/*
	 * Get the number of barrel distortion parameters.
	 */
	unsigned int bd_count() {
		return bdcnum;
	}

	/*
	 * Get the maximum allowable number of barrel distortion parameters.
	 */
	unsigned int bd_max() {
		return BARREL_DEGREE;
	}

	/*
	 * Modify the specified barrel distortion parameter.
	 */
	void bd_modify(unsigned int degree, ale_pos diff) {
		assert (degree < bdcnum);
		bd_set(degree, bd_get(degree) + diff);
	}

	/*
	 * Rescale a transform with a given factor.
	 */
	void rescale(ale_pos factor) {

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		if (_is_projective) {

			for (int i = 0; i < 4; i++)
			for (int j = 0; j < 2; j++)
				x[i][j] *= factor;

		} else {
#if 0
			/*
			 * Euclidean scaling is handled in resultant().
			 */
			for (int i = 0; i < 2; i++)
				eu[i] *= factor;
#endif
		}

		scale_factor *= factor;
	}

	/*
	 * Set a new domain.
	 */

	void set_domain(unsigned int new_height, unsigned int new_width) {
		resultant_memo = 0;
		resultant_inverse_memo = 0;

		input_width = new_width;
		input_height = new_height;
	}

	/*
	 * Set the dimensions of the image.
	 */
	void set_dimensions(const image *im) {

		int new_height = (int) im->height();
		int new_width  = (int) im->width();

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		if (_is_projective) {

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

			_x = transform_unscaled(point(new_height,     0    ));
			_y = transform_unscaled(point(new_height, new_width));
			_z = transform_unscaled(point(    0     , new_width));

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

		    (*q) = transform_scaled(p);
		point q0 = transform_scaled(point(p[0] + 1, p[1]));
		point q1 = transform_scaled(point(p[0], p[1] + 1));

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
	 * coordinate system to another, using the forward transformation.  
	 * This function uses unscaled input coordinates.
	 */
	void map_area_unscaled(point p, point *q, ale_pos d[2]) {

		/*
		 * Determine the coordinates in the target frame for the source
		 * image pixel P and two adjacent source pixels.
		 */

		    (*q) = transform_unscaled(p);
		point q0 = transform_unscaled(point(p[0] + 1, p[1]));
		point q1 = transform_unscaled(point(p[0], p[1] + 1));

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

		    (*q) = scaled_inverse_transform(p);
		point q0 = scaled_inverse_transform(point(p[0] + 1, p[1]));
		point q1 = scaled_inverse_transform(point(p[0], p[1] + 1));


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

		_is_projective = 1;

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
		  / (input_width * scale_factor);

		h = (delta_01 * sigma_1  - delta_11 * sigma_0 )
		  / (delta_01 * delta_12 - delta_11 * delta_02)
		  / (input_height * scale_factor);

		a = (x[1][0] - x[0][0] + g * x[1][0])
		  / (input_width * scale_factor);                       
		b = (x[3][0] - x[0][0] + h * x[3][0])
		  / (input_height * scale_factor);
		c = x[0][0];

		d = (x[1][1] - x[0][1] + g * x[1][1])
		  /  (input_width * scale_factor);                       
		e = (x[3][1] - x[0][1] + h * x[3][1])
		  /  (input_height * scale_factor);
		f =  x[0][1];

		resultant_memo = 1;
		resultant_inverse_memo = 0;

		this->x[0] = scaled_inverse_transform( point(      0      ,      0      ) );
		this->x[1] = scaled_inverse_transform( point( (input_height * scale_factor),      0      ) );
		this->x[2] = scaled_inverse_transform( point( (input_height * scale_factor), (input_width * scale_factor) ) );
		this->x[3] = scaled_inverse_transform( point(      0      , (input_width * scale_factor) ) );

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
		x[1][0] = (input_width * scale_factor);    x[1][1] = 0;
		x[2][0] = (input_width * scale_factor);    x[2][1] = (input_height * scale_factor);
		x[3][0] = 0;              x[3][1] = (input_height * scale_factor);

		/*
		 * Rotate
		 */

		ale_pos theta = eu[2] * M_PI / 180;

		for (i = 0; i < 4; i++) {
			ale_pos _x[2];

			_x[0] = (x[i][0] - (input_width * scale_factor)/2)  * cos(theta)
			      + (x[i][1] - (input_height * scale_factor)/2) * sin(theta)
			      +  (input_width * scale_factor)/2;
			_x[1] = (x[i][1] - (input_height * scale_factor)/2) * cos(theta)
			      - (x[i][0] - (input_width * scale_factor)/2)  * sin(theta)
			      +  (input_height * scale_factor)/2;
			
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

		if (_is_projective) {
			gpt_v0_set(x);
			return;
		}

		/*
		 * Reconstruct euclidean parameters
		 */
			
		gpt_v0_set(x);

		point center((input_height * scale_factor) / 2, (input_width * scale_factor) / 2);
		point center_image = transform_scaled(center);

		this->eu[0] = (center_image[0] - center[0]) / scale_factor;
		this->eu[1] = (center_image[1] - center[1]) / scale_factor;

		point center_left((input_height * scale_factor) / 2, 0);
		point center_left_image = transform_scaled(center_left);

		ale_pos displacement = center_image[0] - center_left_image[0];

		this->eu[2] = asin(2 * displacement / (input_width * scale_factor)) / M_PI * 180;

		if (center_left_image[1] > center_image[1])
			this->eu[2] = this->eu[2] + 180;

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		_is_projective = 0;

	}

	void debug_output() {
		fprintf(stderr, "[t.do ih=%u, iw=%d x=[[%f %f] [%f %f] [%f %f] [%f %f]] eu=[%f %f %f]\n"
				"      a-f=[%f %f %f %f %f %f %f %f] _a-_f=[%f %f %f %f %f %f %f %f]\n"
				"      bdcnm=%d ip=%d rm=%d rim=%d sf=%f]\n", 
				input_height, input_width, 
				x[0][0], x[0][1], 
				x[1][0], x[1][1], 
				x[2][0], x[2][1], 
				x[3][0], x[3][1], 
				eu[0], eu[1], eu[2],
				a, b, c, d, e, f, g, h,
				_a, _b, _c, _d, _e, _f, _g, _h,
				bdcnum,
				_is_projective,
				resultant_memo,
				resultant_inverse_memo,
				scale_factor);
	}

};

#endif
