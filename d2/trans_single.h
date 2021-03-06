// Copyright 2002, 2004 David Hilvert <dhilvert@auricle.dyndns.org>, 
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
 * trans_single.h: Represent transformations of the kind q = c(b^-1(p)),
 * where p is a point in the source coordinate system, q is a point in the
 * target coordinate system, b^-1 is a transformation correcting barrel
 * distortion, and c is a transformation of projective or Euclidean type.
 * (Note that ^-1 in this context indicates the function inverse rather than
 * the exponential.)
 */

#ifndef __trans_single_h__
#define __trans_single_h__

#include "trans_abstract.h"

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
struct trans_single : public trans_abstract {
private:
	point x[4];
	ale_pos eu[3];
	mutable ale_pos a, b, c, d, e, f, g, h;			// matrix
	mutable ale_pos _a, _b, _c, _d, _e, _f, _g, _h;		// matrix inverse
	int _is_projective;
	pixel tonal_multiplier;
	mutable int resultant_memo;
	mutable int resultant_inverse_memo;

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

		int ale_pos_casting = ale_pos_casting_status();
		ale_pos_enable_casting();

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
			
			ale_pos theta = (double) eu[2] * M_PI / 180;

			a = (double) cos(theta) * (double) scale_factor;
			b = (double) sin(theta) * (double) scale_factor;
			c = 0.5 * ((double) input_height * ((double) scale_factor - (double) a) 
			  - (double) input_width * (double) b) + (double) eu[0] 
			  * (double) scale_factor;
			d = -b;
			e = a;
			f = 0.5 * ((double) input_height * (double) b 
			         + (double) input_width * ((double) scale_factor 
				                         - (double) a)) 
				 + (double) eu[1] * (double) scale_factor;
			g = 0;
			h = 0;
		} 

		resultant_memo = 1;

		if (!ale_pos_casting)
			ale_pos_disable_casting();
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

		int ale_pos_casting = ale_pos_casting_status();
		ale_pos_enable_casting();

		/*
		 * For projective transformations, we calculate
		 * the inverse of the forward transformation 
		 * matrix.
		 */

		double scale = (double) a * (double) e - (double) b * (double) d;

		_a = ((double) e * 1 - (double) f * (double) h) / scale;
		_b = ((double) h * (double) c - 1 * (double) b) / scale;
		_c = ((double) b * (double) f - (double) c * (double) e) / scale;
		_d = ((double) f * (double) g - (double) d * 1) / scale;
		_e = (1 * (double) a - (double) g * (double) c) / scale;
		_f = ((double) c * (double) d - (double) a * (double) f) / scale;
		_g = ((double) d * (double) h - (double) e * (double) g) / scale;
		_h = ((double) g * (double) b - (double) h * (double) a) / scale;

		resultant_inverse_memo = 1;

		if (!ale_pos_casting)
			ale_pos_disable_casting();
	}

public:	
	
	trans_single &operator=(const trans_single &ta) {

		this->trans_abstract::operator=(*((trans_abstract *) &ta));

		for (int i = 0; i < 4; i++) {
			x[i] = ta.x[i];
		}

		for (int i = 0; i < 3; i++) {
			eu[i] = ta.eu[i];
		}

		_is_projective = ta._is_projective;

		tonal_multiplier = ta.tonal_multiplier;

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		return *this;
	}

	trans_single(const trans_single &ta) {
		operator=(ta);
	}

	trans_single() {
	}


	/*
	 * Returns non-zero if the transformation might be non-Euclidean.
	 */
	int is_projective() const {
		return _is_projective;
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


#if 0
	/*
	 * operator() is the transformation operator.
	 */
	struct point operator()(struct point p) {
		return transform(p);
	}
#endif

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
	 * Identity transform for given dimensions
	 */
	static struct trans_single eu_identity(unsigned int height = 2, unsigned int width = 2, ale_pos scale_factor = 1) {
		struct trans_single r;

		r.resultant_memo = 0;
		r.resultant_inverse_memo = 0;

		r.eu[0] = 0;
		r.eu[1] = 0;
		r.eu[2] = 0;
		r.input_width = width;
		r.input_height = height;
		r.scale_factor = scale_factor;

		r._is_projective = 0;

		r.tonal_multiplier = pixel(1, 1, 1);

		r.bd_set(0, (ale_pos *) NULL);

		return r;
	}

	/*
	 * Calculate euclidean identity transform for a given image.
	 */
	static struct trans_single eu_identity(const image *i, ale_pos scale_factor = 1) {
		if (i == NULL)
			return eu_identity(2, 2, scale_factor);

		return eu_identity(i->height(), i->width(), scale_factor);
	}

	/*
	 * Projective identity transform for given dimensions.
	 */
	static trans_single gpt_identity(unsigned int height, unsigned int width, ale_pos scale_factor) {
		struct trans_single r = eu_identity(height, width, scale_factor);

		r.eu_to_gpt();

		return r;
	}

	/*
	 * Calculate projective identity transform for a given image.
	 */
	static trans_single gpt_identity(const image *i, ale_pos scale_factor) {
		struct trans_single r = eu_identity(i, scale_factor);

		r.eu_to_gpt();

		return r;
	}

	/*
	 * Set the tonal multiplier
	 */
	void set_tonal_multiplier(pixel p) {
		tonal_multiplier = p;
	}

	pixel get_tonal_multiplier(struct point p) const {
		return tonal_multiplier;
	}

	pixel get_inverse_tonal_multiplier(struct point p) const {
		return tonal_multiplier;
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
		assert(center.defined());
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
	ale_pos eu_get(int param) const {
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
	 * Modify projective transform with control points at the corners of a
	 * specified bounded region.
	 */

	void gpt_modify_bounded(int i1, int i2, ale_pos diff, elem_bounds_int_t eb) {

		point new_x[4];

		trans_single t = gpt_identity(eb.imax - eb.imin, eb.jmax - eb.jmin, 1);

		t.gpt_modify(i1, i2, diff);

		for (int p = 0; p < 4; p++) 
			new_x[p] = t.transform_scaled(x[p] - point(eb.imin, eb.jmin)) + point(eb.imin, eb.jmin);

		gpt_set(new_x);
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

		trans_single t = *this;

		t.resultant_memo = 0;
		t.resultant_inverse_memo = 0;
		t.input_height = (unsigned int) scaled_height();
		t.input_width  = (unsigned int) scaled_width();
		t.scale_factor = 1;
		t.bd_set(0, (ale_pos *) NULL);

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

	void snap(ale_pos interval) {
		for (int i = 0; i < 4; i++)
		for (int j = 0; j < 2; j++)
			x[i][j] = round(x[i][j] / interval) * interval;

		interval /= scale();

		for (int i = 0; i < 2; i++)
			eu[i] = round(eu[i] / interval) * interval;

		interval *= 2 * 180 / M_PI
		          / sqrt(pow(unscaled_height(), 2) 
			       + pow(unscaled_width(),  2));

		eu[2] = round(eu[2] / interval) * interval;

		resultant_memo = 0;
		resultant_inverse_memo = 0;
	}

	/*
	 * Get the specified projective parameter
	 */
	point gpt_get(int point) const {
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
			x[i] += p;
		else {
			eu[0] += p[0] / scale_factor;
			eu[1] += p[1] / scale_factor;
		}
	}

	/*
	 * Rotate by a given amount about a given point.
	 */
	void rotate(point center, ale_pos degrees) {
		if (_is_projective)
		for (int i = 0; i <= 4; i++) {
			ale_pos radians = (double) degrees * M_PI / (double) 180;

			x[i] -= center;
			x[i] = point(
				(double) x[i][0] * cos(radians) + (double) x[i][1] * sin(radians),
				(double) x[i][1] * cos(radians) - (double) x[i][0] * sin(radians));
			x[i] += center;

			resultant_memo = 0;
			resultant_inverse_memo = 0;
		} else {
			assert(center.defined());
			point fixpoint = scaled_inverse_transform(center);
			eu_modify(2, degrees);
			point offset = center - transform_scaled(fixpoint);

			eu_modify(0, offset[0]);
			eu_modify(1, offset[1]);
		}
	}

	void reset_memos() {
		resultant_memo = 0;
		resultant_inverse_memo = 0;

	}

	/*
	 * Rescale a transform with a given factor.
	 */
	void specific_rescale(ale_pos factor) {

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
	}

	/*
	 * Set the dimensions of the image.
	 */
	void specific_set_dimensions(const image *im) {

		int new_height = (int) im->height();
		int new_width  = (int) im->width();

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

		ale_pos theta = (double) eu[2] * M_PI / 180;

		for (i = 0; i < 4; i++) {
			ale_pos _x[2];

			_x[0] = ((double) x[i][0] - (double) (input_width * scale_factor)/2)  * cos(theta)
			      + ((double) x[i][1] - (double) (input_height * scale_factor)/2) * sin(theta)
			      +  (input_width * scale_factor)/2;
			_x[1] = ((double) x[i][1] - (double) (input_height * scale_factor)/2) * cos(theta)
			      - ((double) x[i][0] - (double) (input_width * scale_factor)/2)  * sin(theta)
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
				(double) x[0][0], (double) x[0][1], 
				(double) x[1][0], (double) x[1][1], 
				(double) x[2][0], (double) x[2][1], 
				(double) x[3][0], (double) x[3][1], 
				(double) eu[0], (double) eu[1], (double) eu[2],
				(double) a, (double) b, (double) c, (double) d, 
				(double) e, (double) f, (double) g, (double) h,
				(double) _a, (double) _b, (double) _c, (double) _d, 
				(double) _e, (double) _f, (double) _g, (double) _h,
				bd_count(),
				_is_projective,
				resultant_memo,
				resultant_inverse_memo,
				(double) scale_factor);
	}

};

#endif
