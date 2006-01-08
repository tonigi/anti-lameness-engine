// Copyright 2003, 2004, 2005 David Hilvert <dhilvert@auricle.dyndns.org>,
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
 * d3/scene.h: Representation of a 3D scene.
 */

#ifndef __scene_h__
#define __scene_h__

#include "point.h"

/*
 * View angle multiplier.  
 *
 * Setting this to a value larger than one can be useful for debugging.
 */

#define VIEW_ANGLE_MULTIPLIER 1

class scene {

	/*
	 * Clipping planes
	 */
	static ale_pos front_clip;
	static ale_pos rear_clip;

	/*
	 * Decimation exponents for geometry
	 */
	static ale_pos primary_decimation_upper;
	static ale_pos input_decimation_lower;
	static ale_pos output_decimation_preferred;

	/*
	 * Input resolution divisor
	 */
	static ale_pos input_resolution_divisor;

	/*
	 * Output clipping
	 */
	static int output_clip;

	/*
	 * Model files
	 */
	static const char *load_model_name;
	static const char *save_model_name;

	/*
	 * Occupancy attenuation
	 */

	static double occ_att;

	/*
	 * Normalization of output by weight
	 */

	static int normalize_weights;

	/*
	 * Falloff exponent
	 */

	static double falloff_exponent;

	/*
	 * Third-camera error multiplier
	 */
	static double tc_multiplier;

	/*
	 * Occupancy update iterations
	 */
	static unsigned int ou_iterations;

	/*
	 * Pairwise ambiguity
	 */
	static unsigned int pairwise_ambiguity;

	/*
	 * Pairwise comparisons
	 */
	static const char *pairwise_comparisons;

	/*
	 * 3D Post-exclusion
	 */
	static int d3px_count;
	static double *d3px_parameters;

	/*
	 * Nearness threshold
	 */
	static const ale_real nearness;

	/*
	 * Encounter threshold for defined pixels.
	 */
	static double encounter_threshold;

	/*
	 * Structure to hold input frame information at levels of 
	 * detail between full detail and full decimation.
	 */
	class lod_image {
		unsigned int f;
		unsigned int entries;
		std::vector<const d2::image *> im;
		std::vector<pt> transformation;

	public:
		/*
		 * Constructor
		 */
		lod_image(unsigned int _f) {

			pt _pt;
			
			f = _f;
			im.push_back(d2::image_rw::copy(f, "3D reference image"));
			assert(im.back());
			entries = 1;
			_pt = d3::align::projective(f);
			_pt.scale(1 / _pt.scale_2d());
			transformation.push_back(_pt);

			while (im.back()->height() > 4
			    && im.back()->width() > 4) {

				im.push_back(im.back()->scale_by_half("3D, reduced LOD"));
				assert(im.back());

				_pt.scale(1 / _pt.scale_2d() / pow(2, entries));
				transformation.push_back(_pt);

				entries += 1;
			}
		}

		/*
		 * Get the number of scales
		 */
		unsigned int count() {
			return entries;
		}

		/*
		 * Get an image.
		 */
		const d2::image *get_image(unsigned int i) {
			assert(i >= 0);
			assert(i < entries);
			return im[i];
		}

		/*
		 * Get the transformation
		 */
		pt get_t(unsigned int i) {
			assert(i >= 0);
			assert(i < entries);
			return transformation[i];
		}

		/*
		 * Destructor
		 */
		~lod_image() {
			for (unsigned int i = 0; i < entries; i++) {
				delete im[i];
			}
		}
	};

	/*
	 * Structure to hold input frame information at all levels of detail.
	 */
	class lod_images {

		/*
		 * All images.
		 */

		std::vector<lod_image *> images;

	public:

		lod_images() {
			images.resize(d2::image_rw::count(), NULL);
		}

		void open(unsigned int f) {
			assert (images[f] = NULL);

			if (images[f] == NULL)
				images[f] = new lod_image(f);
		}

		void open_all() {
			for (unsigned int f = 0; f < d2::image_rw::count(); f++)
				open(f);
		}

		lod_image *get(unsigned int f) {
			assert (images[f] != NULL);
			return images[f];
		}

		void close(unsigned int f) {
			assert (images[f] != NULL);
			delete images[f];
			images[f] = NULL;
		}

		void close_all() {
			for (unsigned int f = 0; f < d2::image_rw::count(); f++)
				close(f);
		}

		~lod_images() {
			close_all();
		}
	};

	/*
	 * All levels-of-detail
	 */

	static struct lod_images *al;

	/*
	 * List for calculating weighted median.
	 */
	class wml {
		ale_real *data;
		unsigned int size;
		unsigned int used;

		ale_real &_w(unsigned int i) {
			assert(i <= used);
			return data[i * 2];
		}

		ale_real &_d(unsigned int i) {
			assert(i <= used);
			return data[i * 2 + 1];
		}

		void increase_capacity() {

			if (size > 0)
				size *= 2;
			else
				size = 1;

			data = (ale_real *) realloc(data, sizeof(ale_real) * 2 * (size * 1));

			assert(data);
			assert (size > used);

			if (!data) {
				fprintf(stderr, "Unable to allocate %d bytes of memory\n",
						sizeof(ale_real) * 2 * (size * 1));
				exit(1);
			} 
		}

		void insert_weight(unsigned int i, ale_real v, ale_real w) {
			assert(used < size);
			assert(used >= i);
			for (unsigned int j = used; j > i; j--) {
				_w(j) = _w(j - 1);
				_d(j) = _d(j - 1);
			}

			_w(i) = w;
			_d(i) = v;

			used++;
		}

	public:

		unsigned int get_size() {
			return size;
		}

		unsigned int get_used() {
			return used;
		}

		void print_info() {
			fprintf(stderr, "[st %p sz %u el", this, size);
			for (unsigned int i = 0; i < used; i++)
				fprintf(stderr, " (%f, %f)", _d(i), _w(i));
			fprintf(stderr, "]\n");
		}

		void clear() {
			used = 0;
		}

		void insert_weight(ale_real v, ale_real w) {
			for (unsigned int i = 0; i < used; i++) {
				if (_d(i) == v) {
					_w(i) += w;
					return;
				}
				if (_d(i) > v) {
					if (used == size)
						increase_capacity();
					insert_weight(i, v, w);
					return;
				}
			}
			if (used == size)
				increase_capacity();
			insert_weight(used, v, w);
		}

		/*
		 * Finds the median at half-weight, or between half-weight
		 * and zero-weight, depending on the attenuation value.
		 */

		ale_real find_median(double attenuation = 0) {

			assert(attenuation >= 0);
			assert(attenuation <= 1);

			ale_real zero1 = 0;
			ale_real zero2 = 0;
			ale_real undefined = zero1 / zero2;

			ale_accum weight_sum = 0;

			for (unsigned int i = 0; i < used; i++)
				weight_sum += _w(i);

//			if (weight_sum == 0)
//				return undefined;

			if (used == 0 || used == 1)
				return undefined;

			if (weight_sum == 0) {
				ale_accum data_sum = 0;
				for (unsigned int i = 0; i < used; i++)
					data_sum += _d(i);
				return data_sum / used;
			}
					

			ale_accum midpoint = weight_sum * (0.5 - 0.5 * attenuation);

			ale_accum weight_sum_2 = 0;

			for (unsigned int i = 0; i < used && weight_sum_2 < midpoint; i++) {
				weight_sum_2 += _w(i);

				if (weight_sum_2 > midpoint) {
					return _d(i);
				} else if (weight_sum_2 == midpoint) {
					assert (i + 1 < used);
					return (_d(i) + _d(i + 1)) / 2;
				} 
			}

			return undefined;
		}

		wml(int initial_size = 0) {

//			if (initial_size == 0) {
//				initial_size = (int) (d2::image_rw::count() * 1.5);
//			}

			size = initial_size;
			used = 0;

			if (size > 0) {
				data = (ale_real *) malloc(size * sizeof(ale_real) * 2);
				assert(data);
			} else {
				data = NULL;
			}
		}

		/*
		 * copy constructor.  This is required to avoid undesired frees.
		 */

		wml(const wml &w) {
			size = w.size;
			used = w.used;
			data = (ale_real *) malloc(size * sizeof(ale_real) * 2);
			assert(data);

			memcpy(data, w.data, size * sizeof(ale_real) * 2);
		}

		~wml() {
			free(data);
		}
	};

	/*
	 * Class for information regarding spatial regions of interest.
	 *
	 * This class is configured for convenience in cases where sampling is
	 * performed using an approximation of the fine:box:1,triangle:2 chain.
	 * In this case, the *_1 variables would store the fine data and the
	 * *_2 variables would store the coarse data.  Other uses are also
	 * possible.
	 */

	class spatial_info {
		/*
		 * Map channel value --> weight.
		 */
		wml color_weights_1[3];
		wml color_weights_2[3];

		/*
		 * Current color.
		 */
		d2::pixel color;

		/*
		 * Map occupancy value --> weight.
		 */
		wml occupancy_weights_1;
		wml occupancy_weights_2;

		/*
		 * Current occupancy value.
		 */
		ale_real occupancy;

		/*
		 * pocc/socc density
		 */

		unsigned int pocc_density;
		unsigned int socc_density;

		/*
		 * Insert a weight into a list.
		 */
		void insert_weight(wml *m, ale_real v, ale_real w) {
			m->insert_weight(v, w);
		}

		/*
		 * Find the median of a weighted list.  Uses NaN for undefined.
		 */
		ale_real find_median(wml *m, double attenuation = 0) {
			return m->find_median(attenuation);
		}

	public:
		/*
		 * Constructor.
		 */
		spatial_info() {
			color = d2::pixel::zero();
			occupancy = 0;
			pocc_density = 0;
			socc_density = 0;
		}

		/*
		 * Accumulate color; primary data set.
		 */
		void accumulate_color_1(int f, int i, int j, d2::pixel color, d2::pixel weight) {
			for (int k = 0; k < 3; k++)
				insert_weight(&color_weights_1[k], color[k], weight[k]);
		}

		/*
		 * Accumulate color; secondary data set.
		 */
		void accumulate_color_2(d2::pixel color, d2::pixel weight) {
			for (int k = 0; k < 3; k++)
				insert_weight(&color_weights_2[k], color[k], weight[k]);
		}

		/*
		 * Accumulate occupancy; primary data set.
		 */
		void accumulate_occupancy_1(int f, int i, int j, ale_real occupancy, ale_real weight) {
			insert_weight(&occupancy_weights_1, occupancy, weight);
		}

		/*
		 * Accumulate occupancy; secondary data set.
		 */
		void accumulate_occupancy_2(ale_real occupancy, ale_real weight) {
			insert_weight(&occupancy_weights_2, occupancy, weight);
			
			if (occupancy == 0 || occupancy_weights_2.get_size() < 96)
				return;

			// fprintf(stderr, "%p updated socc with: %f %f\n", this, occupancy, weight);
			// occupancy_weights_2.print_info();
		}

		/*
		 * Update color (and clear accumulation structures).
		 */
		void update_color() {
			for (int d = 0; d < 3; d++) {
				ale_real c = find_median(&color_weights_1[d]);
				if (isnan(c))
					c = find_median(&color_weights_2[d]);
				if (isnan(c))
					c = 0;

				color[d] = c;

				color_weights_1[d].clear();
				color_weights_2[d].clear();
			}
		}

		/*
		 * Update occupancy (and clear accumulation structures).
		 */
		void update_occupancy() {
			ale_real o = find_median(&occupancy_weights_1, occ_att);
			if (isnan(o))
				o = find_median(&occupancy_weights_2, occ_att);
			if (isnan(o))
				o = 0;

			occupancy = o;

			pocc_density = occupancy_weights_1.get_used();
			socc_density = occupancy_weights_2.get_used();

			occupancy_weights_1.clear();
			occupancy_weights_2.clear();

		}

		/*
		 * Get current color.
		 */
		d2::pixel get_color() {
			return color;
		}

		/*
		 * Get current occupancy.
		 */
		ale_real get_occupancy() {
			return occupancy;
		}

		/*
		 * Get primary color density.
		 */

		unsigned int get_pocc_density() {
			return pocc_density;
		}

		unsigned int get_socc_density() {
			return socc_density;
		}
	};

	/*
	 * Map spatial regions of interest to spatial info structures.  XXX:
	 * This may get very poor cache behavior in comparison with, say, an
	 * array.  Unfortunately, there is no immediately obvious array
	 * representation.  If some kind of array representation were adopted,
	 * it would probably cluster regions of similar depth from the
	 * perspective of the typical camera.  In particular, for a
	 * stereoscopic view, depth ordering for two random points tends to be
	 * similar between cameras, I think.  Unfortunately, it is never
	 * identical for all points (unless cameras are co-located).  One
	 * possible approach would be to order based on, say, camera 0's idea
	 * of depth.
	 */

	typedef std::map<struct space::node *, spatial_info> spatial_info_map_t;

	static spatial_info_map_t spatial_info_map;

public:

	/*
	 * Debugging variables.
	 */

	static unsigned long total_ambiguity;
	static unsigned long total_pixels;
	static unsigned long total_divisions;
	static unsigned long total_tsteps;

	/*
	 * Member functions
	 */

	static void et(double et_parameter) {
		encounter_threshold = et_parameter;
	}

	static void load_model(const char *name) {
		load_model_name = name;
	}

	static void save_model(const char *name) {
		save_model_name = name;
	}

	static void fc(ale_pos fc) {
		front_clip = fc;
	}

	static void di_upper(ale_pos _dgi) {
		primary_decimation_upper = _dgi;
	}

	static void do_try(ale_pos _dgo) {
		output_decimation_preferred = _dgo;
	}

	static void di_lower(ale_pos _idiv) {
		input_decimation_lower = _idiv;
	}

	static void oc() {
		output_clip = 1;
	}

	static void no_oc() {
		output_clip = 0;
	}

	static void rc(ale_pos rc) {
		rear_clip = rc;
	}

	/*
	 * Initialize 3D scene from 2D scene, using 2D and 3D alignment
	 * information.
	 */
	static void init_from_d2() {

		/*
		 * Rear clip value of 0 is converted to infinity.
		 */

		if (rear_clip == 0) {
			ale_pos one = +1;
			ale_pos zero = +0;

			rear_clip = one / zero;
			assert(isinf(rear_clip) == +1);
		}

		/*
		 * Scale and translate clipping plane depths.
		 */

		ale_pos cp_scalar = d3::align::projective(0).wc(point(0, 0, 0))[2];

		front_clip = front_clip * cp_scalar - cp_scalar;
		rear_clip = rear_clip * cp_scalar - cp_scalar;


		/*
		 * Allocate image structures.
		 */

		al = new lod_images;

		if (tc_multiplier != 0) {
			al->open_all();
		}
	}

	/*
	 * Perform spatial_info updating on a given subspace, for given
	 * parameters.
	 */
	static void subspace_info_update(space::iterate si, int f, d2::image *weights, const d2::image *im, pt _pt) {
		while(!si.done()) {

			space::traverse st = si.get();

			/*
			 * XXX: This could be more efficient, perhaps.
			 */

			if (spatial_info_map.count(st.get_node()) == 0) {
				si.next();
				continue;
			}

			spatial_info *sn = &spatial_info_map[st.get_node()];

			/*
			 * Get information on the subspace.
			 */

			d2::pixel color = sn->get_color();
			ale_real occupancy = sn->get_occupancy();

			/*
			 * Determine the view-local bounding box for the
			 * subspace.
			 */

			point bb[2];

			st.get_view_local_bb(_pt, bb);

			point min = bb[0];
			point max = bb[1];

//				fprintf(stderr, "frame %d color update space pointer %p, bb (%f, %f) -> (%f, %f)\n", 
//						f, st.get_node(), min[0], min[1], max[0], max[1]);
//
//				fprintf(stderr, "space %p c=[%f %f %f]\n", st.get_node(), color[0], color[1], color[2]);
//				fprintf(stderr, "space %p occ=[%g]\n", st.get_node(), occupancy);

			/*
			 * Use the center of the bounding box to grab interpolation data.
			 */

			d2::point interp((min[0] + max[0]) / 2, (min[1] + max[1]) / 2);

//				fprintf(stderr, "interp=(%f, %f)\n", interp[0], interp[1]);

#if 0
			/*
			 * For interpolation points, ensure that the
			 * bounding box area is at least 0.25. XXX: Why?
			 * Remove this constraint.
			 *
			 * XXX: Is interpolation useful for anything, given
			 * that we're using spatial info registers at multiple
			 * resolutions?
			 */

			if (/* (max[0] - min[0]) * (max[1] - min[1]) > 0.25
			 && */ max[0] > min[0]
			 && max[1] > min[1]) {
				d2::pixel encounter = (d2::pixel(1, 1, 1) - weights->get_bl(interp));
				d2::pixel pcolor = im->get_bl(interp);
				d2::pixel colordiff = (color - pcolor) * (ale_real) 256;

				if (falloff_exponent != 0) {
					d2::pixel max_diff = im->get_max_diff(interp) * (ale_real) 256;

					for (int k = 0; k < 3; k++)
					if (max_diff[k] > 1)
						colordiff[k] /= pow(max_diff[k], falloff_exponent);
				}

//					fprintf(stderr, "color_interp=(%f, %f, %f)\n", pcolor[0], pcolor[1], pcolor[2]);

//					sn->accumulate_color_2(pcolor, encounter);
				d2::pixel channel_occ = pexp(-colordiff * colordiff);
//					fprintf(stderr, "color_diff=(%f, %f, %f)\n", colordiff[0], colordiff[1], colordiff[2]);
//					fprintf(stderr, "channel_occ=(%g, %g, %g)\n", channel_occ[0], channel_occ[1], channel_occ[2]);

				/*
				 * XXX: the best approach is probably to use 3 separate occupancy
				 * data sets, just as there are 3 separate color data sets.
				 */

				ale_accum occ = channel_occ[0];

				for (int k = 1; k < 3; k++)
					if (channel_occ[k] < occ)
						occ = channel_occ[k];

				sn->accumulate_occupancy_2(occ, encounter[0]);
#if 0
				for (int k = 0; k < 3; k++)
					sn->accumulate_occupancy_2(channel_occ[k], encounter[k]);
#endif
			}
#endif
			
			
			/*
			 * Data structure to check modification of weights by
			 * higher-resolution subspaces.
			 */

			std::queue<d2::pixel> weight_queue;

			/*
			 * Check for higher resolution subspaces, and
			 * update the space iterator.
			 */

			if (st.get_node()->positive
			 || st.get_node()->negative) {

				/*
				 * Store information about current weights,
				 * so we will know which areas have been
				 * covered by higher-resolution subspaces.
				 */

				for (int i = (int) ceil(min[0]); i <= (int) floor(max[0]); i++)
				for (int j = (int) ceil(min[1]); j <= (int) floor(max[1]); j++) {
					if (i < 0 || j < 0)
						continue;
					weight_queue.push(weights->get_pixel(i, j));
				}
				
				/*
				 * Cleave space for the higher-resolution pass,
				 * skipping the current space, since we will
				 * process that later.
				 */

				space::iterate cleaved_space = si.cleave();

				cleaved_space.next();

				subspace_info_update(cleaved_space, f, weights, im, _pt);
			} else {
				si.next();
			}
				

			/*
			 * Iterate over pixels in the bounding box,
			 * adding new data to the subspace.  XXX:
			 * assume for now that all pixels in the
			 * bounding box intersect the subspace.
			 */

			for (int i = (int) ceil(min[0]); i <= (int) floor(max[0]); i++)
			for (int j = (int) ceil(min[1]); j <= (int) floor(max[1]); j++) {

				if (i < 0 || j < 0)
					continue;

				d2::pixel pcolor = im->get_pixel(i, j);
				d2::pixel colordiff = (color - pcolor) * (ale_real) 256;

				if (falloff_exponent != 0) {
					d2::pixel max_diff = im->get_max_diff(interp) * (ale_real) 256;

					for (int k = 0; k < 3; k++)
					if (max_diff[k] > 1)
						colordiff[k] /= pow(max_diff[k], falloff_exponent);
				}

//					fprintf(stderr, "(i, j) == (%d, %d); c=[%f %f %f]\n",
//							i, j, pcolor[0], pcolor[1], pcolor[2]);

				/*
				 * Determine the probability of
				 * encounter, divided by the occupancy.
				 */

				d2::pixel encounter = (d2::pixel(1, 1, 1) - weights->get_pixel(i, j));

				/*
				 * Check for higher-resolution modifications.
				 */

				int high_res_mod = 0;

				if (weight_queue.size()) {
					if (weight_queue.front() != weights->get_pixel(i, j)) {
						high_res_mod = 1;
						encounter = d2::pixel(1, 1, 1) - weight_queue.front();
					}
					weight_queue.pop();
				}

				/*
				 * Update subspace.
				 */

				sn->accumulate_color_1(f, i, j, pcolor, encounter);
				d2::pixel channel_occ = pexp(-colordiff * colordiff);
//					fprintf(stderr, "encounter=(%f, %f, %f)\n", encounter[0], encounter[1], encounter[2]);
//					fprintf(stderr, "color_diff=(%f, %f, %f)\n", colordiff[0], colordiff[1], colordiff[2]);
//					fprintf(stderr, "channel_occ=(%g, %g, %g)\n", channel_occ[0], channel_occ[1], channel_occ[2]);

				ale_accum occ = channel_occ[0];

				for (int k = 1; k < 3; k++)
					if (channel_occ[k] < occ)
						occ = channel_occ[k];

				sn->accumulate_occupancy_1(f, i, j, occ, encounter[0]);

				/*
				 * If weights have not been updated by
				 * higher-resolution cells, then update
				 * weights at the current resolution.
				 */

				if (!high_res_mod)
					weights->pix(i, j) += encounter * occupancy;
			}

		}
	}

	/*
	 * Run a single iteration of the spatial_info update cycle.
	 */
	static void spatial_info_update() {
		/*
		 * Iterate through each frame.
		 */
		for (unsigned int f = 0; f < d2::image_rw::count(); f++) {
			/*
			 * Get transformation data.
			 */
			pt _pt = align::projective(f);
			_pt.scale(1 / _pt.scale_2d() / pow(2, ceil(input_decimation_exponent)));

			/*
			 * Get color data for the frames.
			 */

			d2::image *im = d2::image_rw::copy(f, "3D reference image");

			assert(im);

			ale_pos decimation_index = input_decimation_exponent;
			while (decimation_index > 0
			    && im->height() > 2
			    && im->width() > 2) {
				d2::image *iim = im->scale_by_half("3D, reduced LOD");
				assert(iim);
				delete im;
				im = iim;
				decimation_index -= 1;
			}

			assert(im);

			/*
			 * Allocate an image for storing encounter probabilities.
			 */
			d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()), 
					(int) floor(_pt.scaled_width()), 3);

			assert(weights);

			/*
			 * Call subspace_info_update for the root space.
			 */

			subspace_info_update(space::iterate(_pt), f, weights, im, _pt);

			delete im;
			delete weights;
		}

		/*
		 * Update all spatial_info structures.
		 */
		for (spatial_info_map_t::iterator i = spatial_info_map.begin(); i != spatial_info_map.end(); i++) {
			i->second.update_color();
			i->second.update_occupancy();

//			d2::pixel color = i->second.get_color();

//			fprintf(stderr, "space p=%p updated to c=[%f %f %f] o=%f\n",
//					i->first, color[0], color[1], color[2], 
//					i->second.get_occupancy());
		}
	}

	/*
	 * Support function for view() and depth().
	 */

	static const void view_recurse(int type, d2::image *im, d2::image *weights, space::iterate si, pt _pt) {
		while (!si.done()) {
			space::traverse st = si.get();

			/*
			 * XXX: This could be more efficient, perhaps.
			 */

			if (spatial_info_map.count(st.get_node()) == 0) {
				si.next();
				continue;
			}

			spatial_info sn = spatial_info_map[st.get_node()];

			/*
			 * Get information on the subspace.
			 */

			d2::pixel color = sn.get_color();
			// d2::pixel color = d2::pixel(1, 1, 1) * (double) (((unsigned int) (st.get_node()) / sizeof(space)) % 65535);
			ale_real occupancy = sn.get_occupancy();

			/*
			 * Determine the view-local bounding box for the
			 * subspace.
			 */

			point bb[2];

			st.get_view_local_bb(_pt, bb);

			point min = bb[0];
			point max = bb[1];

			/*
			 * Data structure to check modification of weights by
			 * higher-resolution subspaces.
			 */

			std::queue<d2::pixel> weight_queue;

			/*
			 * Check for higher resolution subspaces, and
			 * update the space iterator.
			 */

			if (st.get_node()->positive
			 || st.get_node()->negative) {

				/*
				 * Store information about current weights,
				 * so we will know which areas have been
				 * covered by higher-resolution subspaces.
				 */

				for (int i = (int) ceil(min[0]); i <= (int) floor(max[0]); i++)
				for (int j = (int) ceil(min[1]); j <= (int) floor(max[1]); j++)
					weight_queue.push(weights->get_pixel(i, j));
				
				/*
				 * Cleave space for the higher-resolution pass,
				 * skipping the current space, since we will
				 * process that afterward.
				 */

				space::iterate cleaved_space = si.cleave();

				cleaved_space.next();

				view_recurse(type, im, weights, cleaved_space, _pt);

			} else {
				si.next();
			}
				

			/*
			 * Iterate over pixels in the bounding box, finding
			 * pixels that intersect the subspace.  XXX: assume
			 * for now that all pixels in the bounding box
			 * intersect the subspace.
			 */

			for (int i = (int) ceil(min[0]); i <= (int) floor(max[0]); i++)
			for (int j = (int) ceil(min[1]); j <= (int) floor(max[1]); j++) {

				/*
				 * Check for higher-resolution updates.
				 */

				if (weight_queue.size()) {
					if (weight_queue.front() != weights->get_pixel(i, j)) {
						weight_queue.pop();
						continue;
					}
					weight_queue.pop();
				}

				/*
				 * Determine the probability of encounter.
				 */

				d2::pixel encounter = (d2::pixel(1, 1, 1) - weights->get_pixel(i, j)) * occupancy;

				/*
				 * Update images.
				 */

				if (type == 0) {

					/*
					 * Color view
					 */

					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * color;

				} else if (type == 1) {

					/*
					 * Weighted (transparent) depth display
					 */

					ale_pos depth_value = _pt.wp_scaled(st.get_min())[2];
					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * depth_value;

				} else if (type == 2) {

					/*
					 * Ambiguity (ambivalence) measure.
					 */

					weights->pix(i, j) = d2::pixel(1, 1, 1);
					im->pix(i, j) += 0.1 * d2::pixel(1, 1, 1);

				} else if (type == 3) {

					/*
					 * Closeness measure.
					 */

					ale_pos depth_value = _pt.wp_scaled(st.get_min())[2];
					if (weights->pix(i, j)[0] == 0) {
						weights->pix(i, j) = d2::pixel(1, 1, 1);
						im->pix(i, j) = d2::pixel(1, 1, 1) * depth_value;
					} else if (im->pix(i, j)[2] < depth_value) {
						im->pix(i, j) = d2::pixel(1, 1, 1) * depth_value;
					} else {
						continue;
					}

				} else if (type == 4) {

					/*
					 * Weighted (transparent) contribution display
					 */

					ale_pos contribution_value = sn.get_pocc_density() + sn.get_socc_density();
					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * contribution_value;

				} else if (type == 5) {

					/*
					 * Weighted (transparent) occupancy display
					 */

					ale_pos contribution_value = occupancy;
					weights->pix(i, j) += encounter;
					im->pix(i, j)      += encounter * contribution_value;

				} else if (type == 6) {
					
					/*
					 * (Depth, xres, yres) triple
					 */

					ale_pos depth_value = _pt.wp_scaled(st.get_min())[2];
					weights->pix(i, j)[0] += encounter[0];
					if (weights->pix(i, j)[1] < encounter[0]) {
						weights->pix(i, j)[1] = encounter[0];
						im->pix(i, j)[0] = weights->pix(i, j)[1] * depth_value;
						im->pix(i, j)[1] = max[0] - min[0];
						im->pix(i, j)[2] = max[1] - min[1];
					}

				} else if (type == 7) {
					
					/*
					 * (xoff, yoff, 0) triple
					 */

					weights->pix(i, j)[0] += encounter[0];
					if (weights->pix(i, j)[1] < encounter[0]) {
						weights->pix(i, j)[1] = encounter[0];
						im->pix(i, j)[0] = i - min[0];
						im->pix(i, j)[1] = j - min[1];
						im->pix(i, j)[2] = 0;
					}

				} else 
					assert(0);
			}
		}
	}

	/*
	 * Generate an depth image from a specified view.
	 */
	static const d2::image *depth(pt _pt, int n = -1) {
		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		if (n >= 0) {
			assert((int) floor(d2::align::of(n).scaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(n).scaled_width())
			     == (int) floor(_pt.scaled_width()));
		}

		d2::image *im1 = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		d2::image *im2 = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		d2::image *im3 = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

		/*
		 * Use adaptive subspace data.
		 */

		d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		/*
		 * Iterate through subspaces.
		 */

		space::iterate si(_pt);

		view_recurse(6, im1, weights, si, _pt);

		delete weights;
		weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		view_recurse(7, im2, weights, si, _pt);

		/*
		 * Normalize depths by weights
		 */

		if (normalize_weights)
		for (unsigned int i = 0; i < im1->height(); i++)
		for (unsigned int j = 0; j < im1->width();  j++)
			im1->pix(i, j)[0] /= weights->pix(i, j)[1];

	
		for (unsigned int i = 0; i < im1->height(); i++)
		for (unsigned int j = 0; j < im1->width();  j++) {

			/*
			 * Handle interpolation.
			 */

			d2::point x;
			d2::point blx;
			d2::point res(im1->pix(i, j)[1], im1->pix(i, j)[2]);

			for (int d = 0; d < 2; d++) {

				if (im2->pix(i, j)[d] < res[d] / 2)
					x[d] = (ale_pos) (d?j:i) - res[d] / 2 - im2->pix(i, j)[d];
				else
					x[d] = (ale_pos) (d?j:i) + res[d] / 2 - im2->pix(i, j)[d];

				blx[d] = 1 - ((d?j:i) - x[d]) / res[d];
			}

			ale_real depth_val = 0;
			ale_real depth_weight = 0;

			for (int ii = 0; ii < 2; ii++)
			for (int jj = 0; jj < 2; jj++) {
				d2::point p = x + d2::point(ii, jj) * res;
				if (im1->in_bounds(p)) {

					ale_real d = im1->get_bl(p)[0];

					if (isnan(d))
						continue;

					ale_real w = ((ii ? (1 - blx[0]) : blx[0]) * (jj ? (1 - blx[1]) : blx[1]));
					depth_weight += w;
					depth_val += w * d;
				}
			}

			ale_real depth = depth_val / depth_weight;

			/*
			 * Handle exclusions and encounter thresholds
			 */

			point w = _pt.pw_scaled(point(i, j, depth));

			if (weights->pix(i, j)[0] < encounter_threshold || excluded(w)) {
				im3->pix(i, j) = d2::pixel::zero() / d2::pixel::zero();
			} else {
				im3->pix(i, j) = d2::pixel(1, 1, 1) * depth;
			}
		}

		delete weights;
		delete im1;
		delete im2;

		return im3;
	}

	static const d2::image *depth(unsigned int n) {

		assert (n < d2::image_rw::count());

		pt _pt = align::projective(n);

		return depth(_pt, n);
	}

	/*
	 * Generate an image from a specified view.
	 */
	static const d2::image *view(pt _pt, int n = -1) {
		assert ((unsigned int) n < d2::image_rw::count() || n < 0);

		if (n >= 0) {
			assert((int) floor(d2::align::of(n).scaled_height())
			     == (int) floor(_pt.scaled_height()));
			assert((int) floor(d2::align::of(n).scaled_width())
			     == (int) floor(_pt.scaled_width()));
		}

		const d2::image *depths = NULL;

		if (d3px_count > 0)
			depths = depth(_pt, n);

		d2::image *im = new d2::image_ale_real((int) floor(_pt.scaled_height()),
				               (int) floor(_pt.scaled_width()), 3);

		_pt.view_angle(_pt.view_angle() * VIEW_ANGLE_MULTIPLIER);

		/*
		 * Use adaptive subspace data.
		 */

		d2::image *weights = new d2::image_ale_real((int) floor(_pt.scaled_height()),
						(int) floor(_pt.scaled_width()), 3);

		/*
		 * Iterate through subspaces.
		 */

		space::iterate si(_pt);

		view_recurse(0, im, weights, si, _pt);

		for (unsigned int i = 0; i < im->height(); i++)
		for (unsigned int j = 0; j < im->width();  j++) {
			if (weights->pix(i, j).min_norm() < encounter_threshold
			 || (d3px_count > 0 && isnan(depths->pix(i, j)[0]))) {
				im->pix(i, j) = d2::pixel::zero() / d2::pixel::zero();
				weights->pix(i, j) = d2::pixel::zero();
			} else if (normalize_weights)
				im->pix(i, j) /= weights->pix(i, j);
		}

		if (d3px_count > 0)
			delete depths;

		delete weights;

		return im;
	}

	static void tcem(double _tcem) {
		tc_multiplier = _tcem;
	}

	static void oui(unsigned int _oui) {
		ou_iterations = _oui;
	}

	static void pa(unsigned int _pa) {
		pairwise_ambiguity = _pa;
	}

	static void pc(const char *_pc) {
		pairwise_comparisons = _pc;
	}

	static void d3px(int _d3px_count, double *_d3px_parameters) {
		d3px_count = _d3px_count;
		d3px_parameters = _d3px_parameters;
	}

	static void fx(double _fx) {
		falloff_exponent = _fx;
	}

	static void nw() {
		normalize_weights = 1;
	}

	static void no_nw() {
		normalize_weights = 0;
	}

	static int excluded(point p) {
		for (int n = 0; n < d3px_count; n++) {
			double *region = d3px_parameters + (6 * n);
			if (p[0] >= region[0]
			 && p[0] <= region[1]
			 && p[1] >= region[2]
			 && p[1] <= region[3]
			 && p[2] >= region[4]
			 && p[2] <= region[5])
				return 1;
		}

		return 0;
	}

	static const d2::image *view(unsigned int n) {

		assert (n < d2::image_rw::count());

		pt _pt = align::projective(n);

		return view(_pt, n);
	}

	/*
	 * Add specified control points to the model
	 */
	static void add_control_points() {
	}

	typedef struct {point iw; point ip, is;} analytic;
	typedef std::multimap<ale_real,analytic> score_map;
	typedef std::pair<ale_real,analytic> score_map_element;

	static void solve_point_pair(unsigned int f1, unsigned int f2, int i, int j, ale_pos ii, ale_pos jj, 
			ale_real *normalized_score, analytic *_a, 
			const pt &_pt1, const pt &_pt2, const d2::image *if1, const d2::image *if2) {

		*normalized_score = 0;
		*normalized_score /= *normalized_score;
		assert(isnan(*normalized_score));

		/*
		 * Create an orthonormal basis to
		 * determine line intersection.
		 */

		point bp0 = _pt1.pw_scaled(point(i, j, 0));
		point bp1 = _pt1.pw_scaled(point(i, j, 10));
		point bp2 = _pt2.pw_scaled(point(ii, jj, 0));

		point b0  = (bp1 - bp0).normalize();
		point b1n = bp2 - bp0;
		point b1  = (b1n - b1n.dproduct(b0) * b0).normalize();
		point b2  = point(0, 0, 0).xproduct(b0, b1).normalize(); // Should already have norm=1
		
		/*
		 * Select a fourth point to define a second line.
		 */

		point p3  = _pt2.pw_scaled(point(ii, jj, 10));

		/*
		 * Representation in the new basis.
		 */

		d2::point nbp0 = d2::point(bp0.dproduct(b0), bp0.dproduct(b1));
		// d2::point nbp1 = d2::point(bp1.dproduct(b0), bp1.dproduct(b1));
		d2::point nbp2 = d2::point(bp2.dproduct(b0), bp2.dproduct(b1));
		d2::point np3  = d2::point( p3.dproduct(b0),  p3.dproduct(b1));

		/*
		 * Determine intersection of line
		 * (nbp0, nbp1), which is parallel to
		 * b0, with line (nbp2, np3).
		 */

		/*
		 * XXX: a stronger check would be
		 * better here, e.g., involving the
		 * ratio (np3[0] - nbp2[0]) / (np3[1] -
		 * nbp2[1]).  Also, acceptance of these
		 * cases is probably better than
		 * rejection.
		 */

		if (np3[1] - nbp2[1] == 0)
			return;

		d2::point intersection = d2::point(nbp2[0] 
				+ (nbp0[1] - nbp2[1]) * (np3[0] - nbp2[0]) / (np3[1] - nbp2[1]),
				nbp0[1]);

		ale_pos b2_offset = b2.dproduct(bp0);

		/*
		 * Map the intersection back to the world
		 * basis.
		 */

		point iw = intersection[0] * b0 + intersection[1] * b1 + b2_offset * b2;

		/*
		 * Reject intersection points behind a
		 * camera.
		 */

		point icp = _pt1.wc(iw);
		point ics = _pt2.wc(iw);


		if (icp[2] >= 0 || ics[2] >= 0)
			return;

		/*
		 * Reject clipping plane violations.
		 */

		if (iw[2] > front_clip
		 || iw[2] < rear_clip)
			return;

		/*
		 * Score the point.
		 */

		point ip = _pt1.wp_scaled(iw);

		point is = _pt2.wp_scaled(iw);

		_a->iw = iw;
		_a->ip = ip;
		_a->is = is;

		/*
		 * Calculate score from color match.  Assume for now
		 * that the transformation can be approximated locally
		 * with a translation.
		 */

		ale_pos score = 0;
		ale_pos divisor = 0;
		ale_pos l1_multiplier = 0.125;

		if (if1->in_bounds(ip.xy())
		 && if2->in_bounds(is.xy())) {
			divisor += 1 - l1_multiplier;
			score += (1 - l1_multiplier)
			       * (if1->get_bl(ip.xy()) - if2->get_bl(is.xy())).normsq();
		}

		for (int iii = -1; iii <= 1; iii++)
		for (int jjj = -1; jjj <= 1; jjj++) {
			d2::point t(iii, jjj);

			if (!if1->in_bounds(ip.xy() + t)
			 || !if2->in_bounds(is.xy() + t))
				continue;

			divisor += l1_multiplier;
			score   += l1_multiplier
				 * (if1->get_bl(ip.xy() + t) - if2->get_bl(is.xy() + t)).normsq();
				 
		}

		/*
		 * Include third-camera contributions in the score.
		 */

		if (tc_multiplier != 0)
		for (unsigned int f = 0; f < d2::image_rw::count(); f++) {
			if (f == f1 || f == f2)
				continue;

			assert(cl);
			assert(cl->reference);
			assert(cl->reference[f]);

			pt _ptf = align::projective(f);
			_ptf.scale(1 / _ptf.scale_2d() / pow(2, ceil(input_decimation_exponent)));

			point p = _ptf.wp_scaled(iw);

			if (!cl->reference[f]->in_bounds(p.xy())
			 || !if2->in_bounds(ip.xy()))
				continue;

			divisor += tc_multiplier;
			score   += tc_multiplier
				 * (if1->get_bl(ip.xy()) - cl->reference[f]->get_bl(p.xy())).normsq();
		}

		*normalized_score = score / divisor;
	}

	/*
	 * Generate a map from scores to 3D points for various depths at point (i, j) in f1.
	 */
	static score_map p2f_score_map(unsigned int i, unsigned int j, lod_image *if1, lod_image *if2,
			std::vector<pt> pt_outputs) {

		score_map result;

		/*
		 * Map two depths to the secondary frame.
		 */

		point p1 = _pt2.wp_scaled(_pt1.pw_scaled(point(i, j,  1000)));
		point p2 = _pt2.wp_scaled(_pt1.pw_scaled(point(i, j, -1000)));

		/*
		 * For cases where the mapped points define a
		 * line and where points on the line fall
		 * within the defined area of the frame,
		 * determine the starting point for inspection.
		 * In other cases, continue to the next pixel.
		 */

		ale_pos diff_i = p2[0] - p1[0];
		ale_pos diff_j = p2[1] - p1[1];
		ale_pos slope = diff_j / diff_i;

		if (isnan(slope)) {
			fprintf(stderr, "%d->%d (%d, %d) has undefined slope\n",
					f1, f2, i, j);
			return result;
		}

		/*
		 * Make absurdly large/small slopes either infinity, negative infinity, or zero.
		 */

		if (fabs(slope) > if2->width() * 100) {
			double zero = 0;
			double one  = 1;
			double inf  = one / zero;
			slope = inf;
		} else if (slope < 1 / (double) if2->height() / 100
		        && slope > -1/ (double) if2->height() / 100) {
			slope = 0;
		}

		ale_pos top_intersect = p1[1] - p1[0] * slope;
		ale_pos lef_intersect = p1[0] - p1[1] / slope;
		ale_pos rig_intersect = p1[0] - (p1[1] - if2->width() + 2) / slope;
		ale_pos sp_i, sp_j;

		if (slope == 0) {
			sp_i = lef_intersect;
			sp_j = 0;
		} else if (finite(slope) && top_intersect >= 0 && top_intersect < if2->width() - 1) {
			sp_i = 0;
			sp_j = top_intersect;
		} else if (slope > 0 && lef_intersect >= 0 && lef_intersect <= if2->height() - 1) {
			sp_i = lef_intersect;
			sp_j = 0;
		} else if (slope < 0 && rig_intersect >= 0 && rig_intersect <= if2->height() - 1) {
			sp_i = rig_intersect;
			sp_j = if2->width() - 2;
		} else {
			return result;
		}

		/*
		 * Determine increment values for examining
		 * point, ensuring that incr_i is always
		 * positive.
		 */

		ale_pos incr_i, incr_j;

		if (fabs(diff_i) > fabs(diff_j)) {
			incr_i = 1;
			incr_j = slope;
		} else if (slope > 0) {
			incr_i = 1 / slope;
			incr_j = 1;
		} else {
			incr_i = -1 / slope;
			incr_j = -1;
		}
		
		/*
		 * Examine regions near the projected line.
		 */

		for (ale_pos ii = sp_i, jj = sp_j; 
			ii <= if2->height() - 1 && jj <= if2->width() - 1 && ii >= 0 && jj >= 0; 
			ii += incr_i, jj += incr_j) {

			ale_real normalized_score;
			analytic _a;

			solve_point_pair(f1, f2, i, j, ii, jj, &normalized_score, &_a, _pt1, _pt2, if1, if2);

			/*
			 * Reject points with undefined score.
			 */

			if (!finite(normalized_score))
				continue;

			/*
			 * Add the point to the score map.
			 */

			result.insert(score_map_element(normalized_score, _a));
		}

		/*
		 * Iterate through regions and add new locations with sub-pixel resolution
		 */
		int accumulated_passes = 0;
		int max_acc_passes = pairwise_ambiguity;
		for (score_map::iterator smi = result.begin(); smi != result.end(); smi++) {
			point is = smi->second.is;

			if (accumulated_passes++ >= max_acc_passes)
				break;

			for (ale_pos epsilon = -0.5; epsilon <= 0.5; epsilon += 0.125) {

				if (fabs(epsilon) < 0.001)
					continue;

				ale_real normalized_score;
				analytic _a;

				solve_point_pair(f1, f2, i, j, is[0] + epsilon * incr_i, is[1] + epsilon * incr_j,
						&normalized_score, &_a, _pt1, _pt2, if1, if2);

				if (!finite(normalized_score))
					continue;

				result.insert(score_map_element(normalized_score, _a));
			}
		}

		return result;
	}

	/*
	 * Attempt to refine space around a point, to high and low resolutions
	 * for two cameras, resulting in four resolutions in total.
	 */

	static void refine_space(point iw, pt _pt1, pt _pt2, std::vector<pt> pt_outputs) {

		space::traverse st = space::traverse::root();

		if (!st.includes(iw)) {
			assert(0);
			return;
		}

		int camera_highres[2] = {0, 0};
		int camera_lowres[2] = {0, 0};

		std::vector<int> output_highres;
		std::vector<int> output_lowres;

		output_highres.resize(pt_outputs.size(), 0);
		output_lowres.resize(pt_outputs.size(), 0);

		/*
		 * Loop until all resolutions of interest have been generated.
		 */
		
		for(;;) {

			point frame_min[2] = { point::posinf(), point::posinf() },
			      frame_max[2] = { point::neginf(), point::neginf() };

			std::vector<point> output_max;
			std::vector<point> output_min;

			output_max.resize(pt_outputs.size(), point::posinf());
			output_min.resize(pt_outputs.size(), point::neginf());

			point p[2] = { st.get_min(), st.get_max() };

			/*
			 * Cycle through the corner points bounding the
			 * subspace to determine a bounding box.  
			 *
			 * NB: This code is not identical to
			 * get_view_local_bb(), as it does not clip the
			 * results.
			 */

			for (int ibit = 0; ibit < 2; ibit++)
			for (int jbit = 0; jbit < 2; jbit++)
			for (int kbit = 0; kbit < 2; kbit++) {
				point pp = point(p[ibit][0], p[jbit][1], p[kbit][2]);

				point ppp[2] = {_pt1.wp_scaled(pp), _pt2.wp_scaled(pp)};

				/*
				 * Inputs
				 */

				for (int f = 0; f < 2; f++)
				for (int d = 0; d < 3; d++) {
					if (ppp[f][d] < frame_min[f][d] || isnan(ppp[f][d]))
						frame_min[f][d] = ppp[f][d];
					if (ppp[f][d] > frame_max[f][d] || isnan(ppp[f][d]))
						frame_max[f][d] = ppp[f][d];
				}

				/*
				 * Outputs
				 */

				for (unsigned int n = 0; n < pt_outputs.size(); n++) {

					point ppp_pt = pt_outputs[n].wp_scaled(pp);

					for (int d = 0; d < 3; d++) {
						if (ppp_pt[d] < output_min[n][d] || isnan(ppp_pt[d]))
							output_min[n][d] = ppp_pt[d];
						if (ppp_pt[d] > output_max[n][d] || isnan(ppp_pt[d]))
							output_max[n][d] = ppp_pt[d];
					}
				}
			}

			/*
			 * Generate any new desired spatial registers.
			 */

			/*
			 * Inputs
			 */

			for (int f = 0; f < 2; f++) {

				/*
				 * Low resolution
				 */

				if (frame_max[f][0] - frame_min[f][0] < 2
				 && frame_max[f][1] - frame_min[f][1] < 2
				 && camera_lowres[f] == 0) {
					spatial_info_map[st.get_node()];
					camera_lowres[f] = 1;
				}

				/*
				 * High resolution.
				 */

				if (frame_max[f][0] - frame_min[f][0] < 1
				 && frame_max[f][1] - frame_min[f][1] < 1
				 && camera_highres[f] == 0) {
					spatial_info_map[st.get_node()];
					camera_highres[f] = 1;
				}
			}

			/*
			 * Outputs
			 */

			for (unsigned int n = 0; n < pt_outputs.size(); n++) {

				/*
				 * Low resolution
				 */

				if (output_max[n][0] - output_min[n][0] < 2
				 && output_max[n][1] - output_min[n][1] < 2
				 && output_lowres[n] == 0) {
					spatial_info_map[st.get_node()];
					output_lowres[n] = 1;
				}

				/*
				 * High resolution.
				 */

				if (output_max[n][0] - output_min[n][0] < 1
				 && output_max[n][1] - output_min[n][1] < 1
				 && camera_highres[n] == 0) {
					spatial_info_map[st.get_node()];
					output_highres[n] = 1;
				}
			}

			/*
			 * Check for completion
			 */

			if (camera_highres[0]
			 && camera_highres[1]
			 && camera_lowres[0]
			 && camera_lowres[1]
			 && !count(output_lowres.begin(), output_lowres.end(), 0)
			 && !count(output_highres.begin(), output_highres.end(), 0));
				return;

			/*
			 * Check precision before analyzing space further.
			 */

			if (st.precision_wall()) {
				fprintf(stderr, "\n\n*** Error: reached subspace precision wall ***\n\n");
				assert(0);
				return;
			}

			if (st.positive().includes(iw)) {
				st = st.positive();
				total_tsteps++;
			} else if (st.negative().includes(iw)) {
				st = st.negative();
				total_tsteps++;
			} else {
				fprintf(stderr, "failed iw = (%f, %f, %f)\n", 
						iw[0], iw[1], iw[2]);
				assert(0);
			}
		}
	}

	/*
	 * Analyze space in a manner dependent on the score map.
	 */

	static void analyze_space_from_map(unsigned int f1, unsigned int f2, unsigned int i, 
			unsigned int j, pt _pt1, pt _pt2, score_map _sm, std::vector<pt> pt_outputs) {

		int accumulated_ambiguity = 0;
		int max_acc_amb = pairwise_ambiguity;

		for(score_map::iterator smi = _sm.begin(); smi != _sm.end(); smi++) {

			point iw = smi->second.iw;
			point ip = smi->second.ip;
			// point is = smi->second.is;

			if (accumulated_ambiguity++ >= max_acc_amb)
				break;

			total_ambiguity++;

			/*
			 * Attempt to refine space around the intersection point.
			 */

			refine_space(iw, _pt1, _pt2, pt_outputs);
		}
	}

	static void refine_space_for_output(pt _pt, space::traverse st = space::traverse::root()) {
		if (!spatial_info_map.count(st.get_node())) {
			if (st.get_node()->positive)
				refine_space_for_output(_pt, st.positive());
			if (st.get_node()->negative)
				refine_space_for_output(_pt, st.negative());
			return;
		}

		point bb[2];
		st.get_view_local_bb(_pt, bb);

		if (bb[0].xy().lengthtosq(bb[1].xy()) < 2)
			return;

		if (!_pt.scaled_in_bounds(bb[0]) || !_pt.scaled_in_bounds(bb[1]))
			return;

		spatial_info_map[st.positive().get_node()];
		spatial_info_map[st.negative().get_node()];

		refine_space_for_output(_pt, st.positive());
		refine_space_for_output(_pt, st.negative());
	}

	static void refine_space_for_output(const char *d_out[], const char *v_out[],
			std::map<const char *, pt> *d3_depth_pt,
			std::map<const char *, pt> *d3_output_pt) {

		for (unsigned int f = 0; f < d2::image_rw::count(); f++) 
			if (d_out[f] || v_out[f])
				refine_space_for_output(d3::align::projective(f));

		for (std::map<const char *, pt>::iterator i = d3_depth_pt->begin(); i != d3_depth_pt->end(); i++)
			refine_space_for_output(i->second);

		for (std::map<const char *, pt>::iterator i = d3_output_pt->begin(); i != d3_output_pt->end(); i++)
			refine_space_for_output(i->second);
	}

	/*
	 * Make pt list.
	 */
	static std::vector<pt> make_pt_list(const char *d_out[], const char *v_out[],
			std::map<const char *, pt> *d3_depth_pt,
			std::map<const char *, pt> *d3_output_pt) {

		std::vector<pt> result;

		for (unsigned int n = 0; n < d2::image_rw::count(); n++) {
			if (d_out[n] || v_out[n]) {
				result.push_back(align::projective(n));
			}
		}

		for (std::map<const char *, pt>::iterator i = d3_depth_pt->begin(); i != d3_depth_pt->end(); i++) {
			result.push_back(i->second);
		}

		for (std::map<const char *, pt>::iterator i = d3_output_pt->begin(); i != d3_output_pt->end(); i++) {
			result.push_back(i->second);
		}

		return result;
	}


	/*
	 * Initialize space and identify regions of interest for the adaptive
	 * subspace model.
	 */
	static void make_space(const char *d_out[], const char *v_out[],
			std::map<const char *, pt> *d3_depth_pt,
			std::map<const char *, pt> *d3_output_pt) {

		fprintf(stderr, "Subdividing 3D space");

		std::vector<pt> pt_outputs = make_pt_list(d_out, v_out, d3_depth_pt, d3_output_pt);

		/*
		 * Initialize root space.
		 */

		space::init_root();

		/*
		 * Subdivide space to resolve intensity matches between pairs
		 * of frames.
		 */

		for (unsigned int f1 = 0; f1 < d2::image_rw::count(); f1++) {

			if (!d_out[f1] && !v_out[f1] && !d3_depth_pt->size()
			 && !d3_output_pt->size() && strcmp(pairwise_comparisons, "all"))
				continue;

			std::vector<const d2::image *> if1;

			if1.push_back(d2::image_rw::copy(f1, "3D reference image"));
			assert(if1.back());

			ale_pos decimation_index = input_decimation_exponent;
			while (decimation_index > 0
			    && if1->height() > 2
			    && if1->width() > 2) {

				if1.push_back(if1.back()->scale_by_half("3D, reduced LOD"));
				assert(if1.back());
				
				decimation_index -= 1;
			}

			if (decimation_index > 0) {
				fprintf(stderr, "Error: --di argument is too large.\n");
				exit(1);
			}


			for (unsigned int f2 = 0; f2 < d2::image_rw::count(); f2++) {

				if (f1 == f2)
					continue;

				std::vector<const d2::image *> if2;

				if2.push_back(d2::image_rw::copy(f2, "3D reference image"));
				assert(if2.back());

				ale_pos decimation_index = input_decimation_exponent;
				while (decimation_index > 0
				    && if2->height() > 2
				    && if2->width() > 2) {

					if2.push_back(if2.back()->scale_by_half("3D, reduced LOD"));
					assert(if2.back());

					decimation_index -= 1;
				}

				if (decimation_index > 0) {
					fprintf(stderr, "Error: --di argument is too large.\n");
					exit(1);
				}

				pt _pt1 = align::projective(f1);
				pt _pt2 = align::projective(f2);

				_pt1.scale(1 / _pt1.scale_2d() / pow(2, ceil(input_decimation_exponent)));
				_pt2.scale(1 / _pt2.scale_2d() / pow(2, ceil(input_decimation_exponent)));

				/*
				 * Iterate over all points in the primary frame.
				 */

				for (unsigned int i = 0; i < if1->height(); i++)
				for (unsigned int j = 0; j < if1->width();  j++) {

					total_pixels++;

					/*
					 * Generate a map from scores to 3D points for
					 * various depths in f1.
					 */

					score_map _sm = p2f_score_map(i, j, if1, if2, pt_outputs);

					/*
					 * Analyze space in a manner dependent on the score map.
					 */

					analyze_space_from_map(f1, f2, i, j, _pt1, _pt2, _sm, pt_outputs);
				}
				delete if2;
			}
			delete if1;
		}

		/*
		 * This is expensive.
		 */

		// refine_space_for_output(d_out, v_out, d3_depth_pt, d3_output_pt);

		fprintf(stderr, ".\n");
	}


	/*
	 * Update spatial information structures.
	 *
	 * XXX: the name of this function is horribly misleading.  There isn't
	 * even a 'search depth' any longer, since there is no longer any
	 * bounded DFS occurring.
	 */
	static void reduce_cost_to_search_depth(d2::exposure *exp_out, int inc_bit) {

		/*
		 * Subspace model
		 */

		for (unsigned int i = 0; i < ou_iterations; i++)
			spatial_info_update();
	}

#if 0
	/*
	 * Describe a scene to a renderer
	 */
	static void describe(render *r) {
	}
#endif
};

#endif
