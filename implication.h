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

#ifndef __implication_h__
#define __implication_h__

/*
 * Information about program implication logic.
 */
class implication {
public:
	/*
	 * Describe an aspect of implication logic that
	 * results in a program parameter being changed.
	 */
	static void changed(const char *description, const char *changes, const char *option = NULL) {
		fprintf(stderr, "\n\n");
		fprintf(stderr, "Program options have been automatically modified to satisfy the following:\n\n");
		fprintf(stderr, description);
		fprintf(stderr, "\n\n");

		fprintf(stderr, "Changes are as follows:\n\n");
		fprintf(stderr, changes);
		fprintf(stderr, "\n\n");

		if (option) {
			fprintf(stderr, "This is equivalent to manually setting the following options:\n\n");
			fprintf(stderr, option);
			fprintf(stderr, "\n\n");
		}
	}
};

#endif
