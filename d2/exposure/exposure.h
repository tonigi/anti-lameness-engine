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

/*
 * exposure.h: A superclass for all exposure classes.
 */

#ifndef __exposure_h__
#define __exposure_h__

#include "../pixel.h"

/*
 * This class models a non-linear response function.  More information can be
 * found here:
 *
 * http://wearcam.org/comparametrics.pdf
 */

class exposure {
private:
	static float confidence_exponent;
	pixel _multiplier;

public:
	/*
	 * confidence/uniform static mutators
	 */
	static void set_confidence(float exponent) {
		confidence_exponent = exponent;
	}

	/*
	 * Event listener interface
	 */

	class listener {
		friend class exposure;
	private:
		listener *next;
		char *name;
		const exposure *target;
	public:
		virtual void trigger(pixel multiplier) = 0;

		listener () {
			next = NULL;
			target = NULL;
		}

		virtual ~listener() {
			if (target) {
				target->remove_listener(this);
			}
		}
	};

private:
	mutable listener *listener_head;
public:

	void add_listener(listener *l, char *name) const {
	
		/*
		 * This is a metafunction, so we consider it
		 * to leave the object constant.
		 */
		assert (l->next == NULL);
		assert (l->target == NULL);
		l->next = listener_head;
		l->target = this;
		l->name = name;
		listener_head = l;
	}

	void remove_listener(listener *l) const {

		assert (listener_head != NULL);

		if (listener_head == l) {
			listener_head = listener_head->next;
		} else {
			assert (listener_head->next != NULL);

			listener *a = listener_head;
			listener *b = listener_head->next;

			while (b != l) {
				assert (b->next != NULL);
				a = b;
				b = b->next;
			}

			a->next = b->next;
		}

		l->target = NULL;
	}
	
	void set_multiplier(pixel _multiplier) {
		listener *cl = listener_head;

		while(cl != NULL) {
			// fprintf(stderr, "Triggering '%s'.\n", cl->name);
			cl->trigger(_multiplier / this->_multiplier);
			cl = cl->next;
		}

		this->_multiplier = _multiplier;
	}

	pixel get_multiplier() const {
		return _multiplier;
	}

	/*
	 * This is a very hackish confidence function.  It's zero at the
	 * extremes of camera response and maximal at the center.
	 */
	virtual pixel confidence(pixel input) const {

		if (confidence_exponent) {
			pixel result;
			for (unsigned int k = 0; k < 3; k++) {

				result[k] = (0.25 - pow(0.5 - input[k] / _multiplier[k], 2)) * 4;

				if (result[k] < 0) {
					result[k] = 0.001;
				} else {
					result[k] = pow(result[k], confidence_exponent);

					if (result[k] < 0.001 || !finite(result[k]))
						result[k] = 0.001;
				}
			}
			return result;
		} else {
			return pixel(1, 1, 1);
		}
	}

	/*
	 * Confidence that the real value is lower or higher than the given
	 * value.
	 */
	virtual pixel one_sided_confidence(pixel input, pixel sign) const {
		if (confidence_exponent) {
			pixel result = confidence(input);
			for (unsigned int k = 0; k < 3; k++) {
				if ((sign[k] > 0 && input[k] / _multiplier[k] > 0.5)
				 || (sign[k] < 0 && input[k] / _multiplier[k] < 0.5))
					result[k] = 1;
			}
			return result;
		} else {
			return pixel(1, 1, 1);
		}
	}

	virtual pixel linearize(pixel input) const = 0;
	virtual pixel unlinearize(pixel input) const = 0;

	exposure() {
		listener_head = NULL;
		_multiplier = pixel(1, 1, 1);
	}
};

#endif
