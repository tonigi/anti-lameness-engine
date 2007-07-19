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

/*
 * vise_core.h: Manages instances of vise. 
 */

#ifndef __vise_core_h__
#define __vise_core_h__

#include "vise.h"
#include "image.h"
#include "vise/ma.h"
#include "vise/sf.h"

/*
 * Vise_core initializes, and maintains shared variables for, all instances of
 * vise.
 */

class vise_core {
	static vise **active;
	static unsigned int active_size;
	static ale_real scale_factor;

public:

	/*
	 * Set the VISE scale factor.
	 */
	static void set_scale(ale_real factor) {
		scale_factor = factor;
	}

	/*
	 * Add a new video stabilization engine.
	 */
	static void add(render *chain, const char *type, const char *prefix, const char *suffix) {

		/*
		 * Instantiate an engine of the appropriate type.
		 */
		if (!strcmp(type, "identity")) {

			/*
			 * Identity is a moving average 0 frames to either side.
			 */
			active = (vise **) realloc(active, ++active_size * sizeof(vise *));
			assert(active);
			if (active == NULL) {
				fprintf(stderr, "\n\n*** VISE: Unable to allocate memory ***\n\n\n");
				exit(1);
			}
			active[active_size - 1] = new ma(chain, 0, prefix, suffix, scale_factor);

		} else if (!strncmp(type, "ma:", 3)) {

			/*
			 * Moving average with an unsigned range parameter.
			 */

			unsigned int range;

			if(sscanf(type + 3, "%u", &range) != 1) {
				fprintf(stderr, "\n\n*** VISE: 'ma:' type requires an unsigned argument. ***\n\n\n");
				exit(1);
			}

			active = (vise **) realloc(active, ++active_size * sizeof(vise *));
			assert(active);
			if (active == NULL) {
				fprintf(stderr, "\n\n*** VISE: Unable to allocate memory ***\n\n\n");
				exit(1);
			}
			active[active_size - 1] = new ma(chain, range, prefix, suffix, scale_factor);
		} else if (!strncmp(type, "sf:", 3)) {

			/*
			 * Single frame with an unsigned frame parameter.
			 */

			unsigned int frame;

			if(sscanf(type + 3, "%u", &frame) != 1) {
				fprintf(stderr, "\n\n*** VISE: 'sf:' type requires an unsigned argument. ***\n\n\n");
				exit(1);
			}

			active = (vise **) realloc(active, ++active_size * sizeof(vise *));
			assert(active);
			if (active == NULL) {
				fprintf(stderr, "\n\n*** VISE: Unable to allocate memory ***\n\n\n");
				exit(1);
			}
			active[active_size - 1] = new sf(chain, frame, prefix, suffix, scale_factor);

		} else {
			fprintf(stderr, "\n\n*** VISE: Unknown type '%s' ***\n\n\n", type);
			exit(1);
		}

	}

	/*
	 * Add a new image to the rendering queue.
	 */
	static void frame_queue_add(unsigned int frame_number) {

		/*
		 * Process the current queue.
		 */
		for (unsigned int i = 0; i < active_size; i++) {
			int lag = active[i]->lag();

			if ((int) frame_number - lag >= 0)
				active[i]->render_frame(frame_number - lag);
		}

		/*
		 * If this is the last frame, then complete all rendering.
		 */
		if (frame_number == image_rw::count() - 1) 
		for (unsigned int i = 0; i < active_size; i++)
		for (int j = active[i]->lag() - 1; j >= 0; j--)
			active[i]->render_frame(frame_number - j);
	}
};

#endif
