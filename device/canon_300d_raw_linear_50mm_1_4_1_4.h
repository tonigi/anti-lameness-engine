// Copyright 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 3 of the License,
    or (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef __canon_300d_raw_linear_50mm_1_4_1_4_h__
#define __canon_300d_raw_linear_50mm_1_4_1_4_h__

#include "../d2.h"

/*
 * Extend the Canon 300D (Digital Rebel) device with an EF 50mm f/1.4 lens at 1.4.
 */

#define LPSF_ROWS 9
#define LPSF_COLS 9

/*
 * This filter (50-1.4-1.4-1) can produce haloing.
 */

static const ale_real canon_300d_raw_linear_50mm_1_4_1_4_lpsf_calibrated_response[LPSF_ROWS][LPSF_COLS][3] = {
	{
		{ 41.6, 32, 118.4, },
		{ 6.4, -3.2, 32, },
		{ -9.6, 35.2, 9.6, },
		{ 12.8, -25.6, 9.6, },
		{ -16, -3.2, 6.4, },
		{ 9.6, 3.2, 35.2, },
		{ 6.4, 8.88178419700125e-16, 25.6, },
		{ -3.2, -6.4, 16, },
		{ -3.2, 8.88178419700125e-16, 92.8000000000001, },
	}, {
		{ 22.4, 16, 32, },
		{ -41.6, -35.2, -35.2, },
		{ 12.8, -9.6, -19.2, },
		{ -38.4, 3.2, 35.2, },
		{ 25.6, -12.8, -48, },
		{ 19.2, 16, 35.2, },
		{ -38.4, 3.2, -35.2, },
		{ -6.4, 8.88178419700125e-16, -35.2, },
		{ 19.2, 9.6, 28.8, },
	}, {
		{ 16, 9.6, 25.6, },
		{ 8.88178419700125e-16, -3.2, -9.6, },
		{ 19.2, 25.6, 64, },
		{ 44.8, 73.6, 67.2, },
		{ 67.2, 38.4, 118.4, },
		{ 57.6, 44.8, 105.6, },
		{ 12.8, -22.4, 48, },
		{ -9.6, -9.6, -16, },
		{ -9.6, -16, -2.66453525910038e-15, },
	}, {
		{ 6.4, -12.8, 25.6, },
		{ 6.4, 28.8, -44.8, },
		{ 57.6, 54.4, 102.4, },
		{ 38.4, 44.8, 64, },
		{ 115.2, 92.8000000000001, 44.8, },
		{ 115.2, 108.8, 115.2, },
		{ 28.8, 44.8, 108.8, },
		{ 12.8, 12.8, 48, },
		{ -12.8, -16, 12.8, },
	}, {
		{ 9.6, -6.4, 41.6, },
		{ -3.2, 22.4, 6.4, },
		{ 76.8, 115.2, 115.2, },
		{ 9.59999999999999, -28.8, -48, },
		{ 884.799999999998, 948.799999999999, 884.799999999998, },
		{ 115.2, 115.2, 115.2, },
		{ 76.8, 83.2, 115.2, },
		{ 25.6, 25.6, 51.2, },
		{ -25.6, 3.2, 9.6, },
	}, {
		{ 3.2, 9.6, 60.8, },
		{ 9.6, 19.2, -48, },
		{ 32, 73.6, 115.2, },
		{ 19.2, 3.2, 9.6, },
		{ 64, 73.6, 25.6, },
		{ 112, 96, 115.2, },
		{ 57.6, 64, 92.8, },
		{ -2.66453525910038e-15, -35.2, -22.4, },
		{ 6.4, -25.6, 12.8, },
	}, {
		{ 22.4, 0, 54.4, },
		{ 9.6, -3.2, -19.2, },
		{ 28.8, 51.2, 76.8, },
		{ 38.4, 32, 67.2, },
		{ 73.6, 41.6, 112, },
		{ 48, 60.8, 80, },
		{ -9.6, -12.8, 60.8, },
		{ -16, -32, -28.8, },
		{ -6.4, 19.2, -2.66453525910038e-15, },
	}, {
		{ 6.4, 3.2, 22.4, },
		{ -28.8, -22.4, -41.6, },
		{ -12.8, 3.2, -16, },
		{ 12.8, -9.6, 3.2, },
		{ -22.4, 16, -41.6, },
		{ 9.6, -22.4, 6.4, },
		{ -12.8, -12.8, -22.4, },
		{ 12.8, 6.4, -44.8, },
		{ 6.4, 19.2, 60.8, },
	}, {
		{ 35.2, 28.8, 108.8, },
		{ 3.2, 8.88178419700125e-16, 28.8, },
		{ -9.6, -6.4, 22.4, },
		{ -9.6, -16, 9.6, },
		{ -12.8, -3.2, 16, },
		{ 12.8, 9.6, 32, },
		{ 9.6, 12.8, 28.8, },
		{ -3.2, -3.2, 19.2, },
		{ -9.6, -16, 76.8 }  
	}
};

class canon_300d_raw_linear_50mm_1_4_1_4 : public canon_300d_raw_linear_50mm_1_4 {
public:

	/*
	 * Linear colorspace PSF
	 */

	class lpsf : public d2::psf_template<LPSF_ROWS, LPSF_COLS> {
	public:
		lpsf() : d2::psf_template<LPSF_ROWS, LPSF_COLS> (3, 3, 
				canon_300d_raw_linear_50mm_1_4_1_4_lpsf_calibrated_response) {
		}
	};

};

#undef LPSF_ROWS
#undef LPSF_COLS
#undef NLPSF_ROWS
#undef NLPSF_COLS

#endif
