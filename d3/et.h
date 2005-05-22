// Copyright 2003 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the License,
    or (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 * d3/et.h: Represent three-dimensional Euclidean transformations.
 */

#ifndef __et_h__
#define __et_h__

#include "point.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Structure to describe a three-dimensional euclidean transformation.
 *
 * We consider points to be transformed as homogeneous coordinate vectors
 * multiplied on the right of the transformation matrix.  The transformations
 * are assumed to be small, so xyz Euler angles are used to specify rotation.
 * The order of transformations is (starting with a point representation in the
 * original coordinate system):
 *
 * i)   translation of the point
 * ii)  x-rotation about the origin
 * iii) y-rotation about the origin
 * iv)  z-rotation about the origin
 *
 * For more information on Euler angles, see:
 *
 * http://mathworld.wolfram.com/EulerAngles.html
 *
 * For more information on rotation matrices, see:
 *
 * http://mathworld.wolfram.com/RotationMatrix.html
 *
 * XXX: inconsistently with the 2D class, this class uses radians to represent
 * rotation values.
 *
 */

struct et {
private:
	ale_pos translation[3];
	ale_pos rotation[3];
	mutable ale_pos matrix[4][4];
	mutable ale_pos matrix_inverse[4][4];
	mutable int resultant_memo;
	mutable int resultant_inverse_memo;

	
	/*
	 * Create an identity matrix
	 */
	void mident(ale_pos m1[4][4]) const {
		for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			m1[i][j] = (i == j) ? 1 : 0;
	}

	/*
	 * Create a rotation matrix
	 */
	void mrot(ale_pos m1[4][4], ale_pos angle, int index) const {
		mident(m1);
		m1[(index+1)%3][(index+1)%3] = cos(angle);
		m1[(index+2)%3][(index+2)%3] = cos(angle);
		m1[(index+1)%3][(index+2)%3] = sin(angle);
		m1[(index+2)%3][(index+1)%3] = -sin(angle);
	}

	/*
	 * Multiply matrices M1 and M2; return result M3.
	 */
	void mmult(ale_pos m3[4][4], ale_pos m1[4][4], ale_pos m2[4][4]) const {
		int k;

		for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		for (m3[i][j] = 0, k = 0; k < 4; k++) 
			m3[i][j] += m1[i][k] * m2[k][j];
	}

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

		/*
		 * Multiply matrices.
		 */

		ale_pos m1[4][4], m2[4][4], m3[4][4];

		/*
		 * Translation
		 */

		mident(m1);

		for (int i = 0; i < 4; i++)
			m1[i][3] = translation[i];

		/*
		 * Rotation about x
		 */

		mrot(m2, rotation[0], 0);
		mmult(m3, m2, m1);

		/*
		 * Rotation about y
		 */

		mrot(m1, rotation[1], 1);
		mmult(m2, m1, m3);

		/*
		 * Rotation about z
		 */

		mrot(m3, rotation[2], 2);
		mmult(matrix, m3, m2);

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

		/*
		 * The inverse can be explicitly calculated from the rotation
		 * and translation parameters.  We invert the individual
		 * matrices and reverse the order of multiplication.
		 */

		ale_pos m1[4][4], m2[4][4], m3[4][4];

		/*
		 * Translation
		 */

		mident(m1);

		for (int i = 0; i < 4; i++)
			m1[i][3] = -translation[i];

		/*
		 * Rotation about x
		 */

		mrot(m2, -rotation[0], 0);
		mmult(m3, m1, m2);

		/*
		 * Rotation about y
		 */

		mrot(m1, -rotation[1], 1);
		mmult(m2, m3, m1);

		/*
		 * Rotation about z
		 */

		mrot(m3, -rotation[2], 2);
		mmult(matrix_inverse, m2, m3);

		resultant_inverse_memo = 1;
	}

public:	
	/*
	 * Transform point p.
	 */
	struct point transform(struct point p) const {
		struct point result;

		resultant();

		for (int i = 0; i < 3; i++) {
			for (int k = 0; k < 3; k++)
				result[i] += matrix[i][k] * p[k];
			result[i] += matrix[i][3];
		}

		return result;
	}

	/*
	 * operator() is the transformation operator.
	 */
	struct point operator()(struct point p) const {
		return transform(p);
	}

	/*
	 * Map point p using the inverse of the transform.
	 */
	struct point inverse_transform(struct point p) const {
		struct point result;

		resultant_inverse();

		for (int i = 0; i < 3; i++) {
			for (int k = 0; k < 3; k++)
				result[i] += matrix_inverse[i][k] * p[k];
			result[i] += matrix_inverse[i][3];
		}

		return result;
	}

	
	/*
	 * Return identity transformation.
	 */
	static struct et identity() {
		struct et r;

		r.resultant_memo = 0;
		r.resultant_inverse_memo = 0;

		r.translation[0] = 0;
		r.translation[1] = 0;
		r.translation[2] = 0;

		r.rotation[0] = 0;
		r.rotation[1] = 0;
		r.rotation[2] = 0;

		return r;
	}

	/*
	 * Adjust translation in the indicated manner.
	 */
	void modify_translation(int i1, ale_pos diff) {
		assert (i1 >= 0);
		assert (i1 < 3);

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		translation[i1] += diff;
	}

	/*
	 * Adjust rotation in the indicated manner.
	 */
	void modify_rotation(int i1, ale_pos diff) {
		assert (i1 >= 0);
		assert (i1 < 3);

		resultant_memo = 0;
		resultant_inverse_memo = 0;

		rotation[i1] += diff;
	}

};

#endif
