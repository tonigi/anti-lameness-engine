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
	static ale_real confidence_exponent;
	pixel _multiplier;
	ale_real _gain_multiplier;
	static ale_real _gain_reference;
	ale_real _black_level;

public:

	/*
	 * confidence/uniform static mutators
	 */
	static void set_confidence(ale_real exponent) {
		confidence_exponent = exponent;
	}

	/*
	 * confidence accessor
	 */
	static ale_real get_confidence() {
		return confidence_exponent;
	}

	/*
	 * Event listener interface
	 */

	class listener {
		friend class exposure;
	private:
		listener *next;
		const char *name;
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

	void add_listener(listener *l, const char *name) const {
	
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

		fprintf(stderr, "_m=%f %f %f t->_m=%f %f %f\n",
			(double) _multiplier[0],
			(double) _multiplier[1],
			(double) _multiplier[2],
			(double) this->_multiplier[0],
			(double) this->_multiplier[1],
			(double) this->_multiplier[2]);

		while(cl != NULL) {
			fprintf(stderr, "Triggering '%s' with %f %f %f.\n", cl->name, 
					(double) (_multiplier[0] / this->_multiplier[0]),
					(double) (_multiplier[1] / this->_multiplier[1]),
					(double) (_multiplier[2] / this->_multiplier[2]));
			cl->trigger(_multiplier / this->_multiplier);
			cl = cl->next;
		}

		this->_multiplier = _multiplier;
	}

	void set_gain_multiplier(ale_real g) {
		_gain_multiplier = g;
	}

	ale_real get_gain_multiplier() {
		return _gain_multiplier;
	}

	void set_black_level(ale_real b) {
		_black_level = b;
	}

	ale_real get_black_level() {
		return _black_level;
	}

	static void set_gain_reference(ale_real r) {
		_gain_reference = r;
	}

	static ale_real get_gain_reference() {
		return _gain_reference;
	}

	pixel get_multiplier() const {
		return _multiplier;
	}

	virtual ale_real confidence(unsigned int k, ale_real input, 
			ale_real confidence_floor = ale_real_confidence_floor) const {

		ale_real _0 = (ale_real) 0;
		ale_real _4 = (ale_real) 4;
		ale_real _2 = (ale_real) 2;
		ale_real _1 = (ale_real) 1;

		if (confidence_exponent == _0)
			return 1;

		ale_real input_scaled = input / _multiplier[k];

		ale_real unexponentiated = _4 * ((_1 / _4) - (ale_real) pow((_1 / _2) 
		                                                      - input_scaled, _2));
		// ale_real unexponentiated = 4 * input_scaled * (0.25 - pow(0.5 - input_scaled, 2));

		if (unexponentiated < _0) 
			return confidence_floor;

		ale_real exponentiated = pow(unexponentiated, confidence_exponent);

		if (exponentiated < confidence_floor || !finite(exponentiated))
			return confidence_floor;

		return exponentiated;
	}

	/*
	 * This is a very hackish confidence function.  It's zero at the
	 * extremes of camera response and maximal at the center.
	 */
	virtual pixel confidence(pixel input, 
			ale_real confidence_floor = ale_real_confidence_floor) const {

		if (confidence_exponent != 0) {
			return pixel(confidence(0, input[0], confidence_floor),
				     confidence(1, input[1], confidence_floor),
				     confidence(2, input[2], confidence_floor));
		} else {
			return pixel(1, 1, 1);
		}
	}

	/*
	 * Confidence that the real value is lower or higher than the given
	 * value.  
	 *
	 * XXX: This function now applies the one-sided condition only to
	 * responses greater than 50%.  Should it be called
	 * 'upper_one_sided_confidence' instead?
	 */
	virtual pixel one_sided_confidence(pixel input, pixel sign) const {
		if (confidence_exponent != 0) {
			pixel result = confidence(input);
			for (unsigned int k = 0; k < 3; k++) {
				if (sign[k] > 0 && input[k] / _multiplier[k] > 1 / (ale_real) 2)
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
		_gain_multiplier = 1;
		_black_level = 0;
	}

	virtual ~exposure() {
	}
};

#endif
