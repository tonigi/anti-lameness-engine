// Copyright 2002, 2003, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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

#ifndef __optimizations_h__
#define __optimizations_h__

/*
 * Class for implementing certain optimizations
 */

class optimizations {
public:

	/*
	 * Instead of allocating new memory for the Irani-Peleg approximation,
	 * use allocated memory from the default rendering chain.
	 *
	 * XXX: This could break if de-allocation of the IPC approximation
	 * actually occurred.
	 */
	static d2::image *get_ip_working_image(const d2::image *im) {
#if OPTIMIZATIONS == 1
		return (d2::image *) im;
#else
		return im->clone("IPC Approximation");
#endif
	}

	/*
	 * When starting work on the 3D scene, get rid of memory allocated
	 * for 2D rendering chains.
	 */
	static void begin_3d_work() {
#if OPTIMIZATIONS == 1
		d2::render::free_all_memory();
#endif
	}
};


#endif
