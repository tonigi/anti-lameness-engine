// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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

#ifndef __util_h__
#define __util_h__

#ifdef USE_IOCTL
#include <sys/ioctl.h>
#endif

/*
 * Returns the terminal width if possible.  Else, returns
 * -1.  Based on the function 'determine_screen_width'
 * in wget 1.9.1.
 */
static inline int get_terminal_width (FILE *tty_file) {
	int error_code = -1;

#if defined USE_IOCTL && defined TIOCGWINSZ
	struct winsize wsz;

	if(ioctl(fileno(tty_file), TIOCGWINSZ, &wsz) < 0)
		return error_code;

	return wsz.ws_col;
#else
	return error_code;
#endif
}
#endif
