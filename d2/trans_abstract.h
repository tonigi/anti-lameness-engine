// Copyright 2002, 2004, 2007 David Hilvert <dhilvert@auricle.dyndns.org>, 
//                                          <dhilvert@ugcs.caltech.edu>

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
 * trans_abstract.h: Abstract transformation superclass.
 */

#ifndef __trans_abstract_h__
#define __trans_abstract_h__

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

struct trans_abstract {
private:
	ale_pos bdc[BARREL_DEGREE];				// barrel-dist. coeffs.
	unsigned int bdcnum;					// number of bdcs

protected:
	ale_pos scale_factor;
	unsigned int input_height, input_width;

	virtual void specific_rescale(ale_pos factor) = 0;
	virtual void reset_memos() = 0;
	virtual void specific_set_dimensions(image *im) = 0;

public:	

	/*
	 * Returns non-zero if the transformation might be non-Euclidean.
	 */
	virtual int is_projective() = 0;

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
	 * Transformation sans barrel distortion
	 */
	virtual struct point pe(struct point p) const = 0;

	/*
	 * Transformation inverse sans barrel distortion
	 */
	virtual struct point pei(struct point p) const = 0;

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
		assert (p.defined());
		point q = unscaled_inverse_transform(p);

		q[0] *= scale_factor;
		q[1] *= scale_factor;

		return q;
	}

	/*
	 * Calculate projective transformation parameters from a euclidean
	 * transformation.
	 */
	virtual void eu_to_gpt() = 0;

	/*
	 * Modify a euclidean transform in the indicated manner.
	 */
	virtual void eu_modify(int i1, ale_pos diff) = 0;

	/*
	 * Rotate about a given point in the original reference frame.
	 */
	virtual void eu_rotate_about_scaled(point center, ale_pos diff) = 0;

	/*
	 * Modify all euclidean parameters at once.
	 */
	virtual void eu_set(ale_pos eu[3]) = 0;

	/*
	 * Get the specified euclidean parameter
	 */
	virtual ale_pos eu_get(int param) = 0;

	/*
	 * Modify a projective transform in the indicated manner.
	 */
	virtual void gpt_modify(int i1, int i2, ale_pos diff) = 0;

	/*
	 * Modify a projective transform according to the group operation.
	 */
	virtual void gr_modify(int i1, int i2, ale_pos diff) = 0;

	/*
	 * Modify all projective parameters at once.
	 */
	virtual void gpt_set(point x[4]) = 0;

	virtual void gpt_set(point x1, point x2, point x3, point x4) = 0;

	/*
	 * Get the specified projective parameter
	 */
	virtual point gpt_get(int point) = 0;

	/*
	 * Get the specified projective parameter
	 */
	virtual ale_pos gpt_get(int point, int dim) = 0;

	/*
	 * Translate by a given amount
	 */
	virtual void translate(point p) = 0;

	/*
	 * Rotate by a given amount about a given point.
	 */
	virtual void rotate(point p, ale_pos degrees) = 0;

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
	 * Get all barrel distortion parameters.
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
		specific_rescale();
		scale_factor *= factor;
	}

	/*
	 * Set a new domain.
	 */

	void set_domain(unsigned int new_height, unsigned int new_width) {
		reset_memos();
		input_width = new_width;
		input_height = new_height;
	}

	/*
	 * Set the dimensions of the image.
	 */
	void set_dimensions(const image *im) {
		reset_memos();
		specific_set_dimensions(im);
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
	virtual void gpt_v0_set(point x[4]) = 0;

	/*
	 * Modify all euclidean parameters at once.  Accommodate bugs in the
	 * version 0 transformation file handler (ALE versions 0.4.0p1 and 
	 * earlier).  This code is only called when using a transformation data
	 * file created with an old version of ALE.
	 */
	virtual void eu_v0_set(ale_pos eu[3]) = 0;

	virtual void debug_output() = 0;
};

#endif
