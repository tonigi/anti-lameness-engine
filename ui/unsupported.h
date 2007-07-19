// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>, 
//                              <dhilvert@ugcs.caltech.edu>

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

#ifndef __unsupported_h__
#define __unsupported_h__

/*
 * Information about unsupported features.
 */
class unsupported {
public:
	/*
	 * Describe a feature that is unsupported for now.
	 */
	static void fornow(const char *description) {
		fprintf(stderr, "\n\n");
		fprintf(stderr, "The following feature is currently unsupported:\n\n");
		fprintf(stderr, description);
		fprintf(stderr, "\n\n");
#if 0
		fprintf(stderr, "For more information, see http://auricle.dyndns.org/ALE/unsupported/currently/\n\n");
#endif

		exit(1);
	}

	/*
	 * Describe an option that is undocumented.
	 */
	static void undocumented(const char *description) {
		fprintf(stderr, "\n\n");
		fprintf(stderr, "Warning: %s is undocumented.\n", description);
	}

	/*
	 * Describe an option that is no longer supported.
	 */
	static void discontinued(const char *description, const char *alternative = NULL, const char *alternative2 = NULL) {
		fprintf(stderr, "\n\n");
		fprintf(stderr, "Error: %s is no longer supported.\n", description);

		if (alternative && alternative2) {
			fprintf(stderr, "  Use either: %s\n", alternative);
			fprintf(stderr, "          or: %s\n", alternative2);
		} else if (alternative)
			fprintf(stderr, "  Use: %s\n", alternative);

		fprintf(stderr, "\n");

		exit(1);
	}
};

#endif
