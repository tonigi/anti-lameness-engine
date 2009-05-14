// Copyright 2002, 2003, 2004, 2005, 2006 David Hilvert <dhilvert@auricle.dyndns.org>, 
//                                                      <dhilvert@ugcs.caltech.edu>

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

#ifndef __input_h__
#define __input_h__

/*
 * ANSI C and POSIX include files.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
// #include <math.h>
#include <stack>
#include <map>

#include "../ale_math.h"

/*
 * Interface files
 */

#include "ui.h"
#include "accel.h"
#include "unsupported.h"
#include "implication.h"

/*
 * Configuration
 */

#if HAVE_CONFIG_H
# include <config.h>
#endif

/*
 * Types
 */

#include "../ale_pos.h"
#include "../ale_real.h"

/*
 * 2D include files
 */

#include "../d2.h"

/*
 * 3D include files
 */

#include "../d3.h"

/*
 * Thread include files
 */

#include "../thread.h"

/*
 * Device configuration files
 */

#include "../device/xvp610_320x240.h"
#include "../device/xvp610_640x480.h"
#include "../device/ov7620_raw_linear.h"
#include "../device/canon_300d_raw_linear.h"
#include "../device/canon_300d_raw_linear_85mm_1_8.h"
#include "../device/canon_300d_raw_linear_50mm_1_8.h"
#include "../device/canon_300d_raw_linear_50mm_1_4.h"
#include "../device/canon_300d_raw_linear_50mm_1_4_1_4.h"
#include "../device/nikon_d50.h"

/*
 * Help files
 */

#include "help.h"

class input {

	/*
	 * Flag for global options.
	 */

	static int global_options;

	/*
	 * Helper functions.
	 */

	/*
	 * Argument counter.
	 *
	 * Counts instances of a given option.
	 */
	static unsigned int arg_count(int argc, const char *argv[], const char *arg) {
		unsigned int count = 0;
		for (int i = 0; i < argc; i++) {
			if (!strcmp(argv[i], arg))
				count++;
			else if (!strcmp(argv[i], "--"))
				return count;
		}
		return count;
	}

	/*
	 * Argument prefix counter.
	 *
	 * Counts instances of a given option prefix.
	 */
	static unsigned int arg_prefix_count(int argc, const char *argv[], const char *pfix) {
		unsigned int count = 0;
		for (int i = 0; i < argc; i++) {
			if (!strncmp(argv[i], pfix, strlen(pfix)))
				count++;
			else if (!strcmp(argv[i], "--"))
				return count;
		}
		return count;
	}

	/*
	 * Reallocation function
	 */
	static void *local_realloc(void *ptr, size_t size) {
		void *new_ptr = realloc(ptr, size);

		if (new_ptr == NULL)
			ui::get()->memory_error_location("main()");

		return new_ptr;
	}

	/*
	 * Not enough arguments function.
	 */
	static void not_enough(const char *opt_name) {
		ui::get()->cli_not_enough(opt_name);
	}

	/*
	 * Bad argument function
	 */
	static void bad_arg(const char *opt_name) {
		ui::get()->cli_bad_arg(opt_name);
	}

	/*
	 * String comparison class.
	 */

	class compare_strings {
	public:
		int operator()(const char *A, const char *B) const {
			return strcmp(A, B) < 0;
		}
	};

	/*
	 * Environment structures.
	 *
	 * XXX: It's arguable that these should be public members of the
	 * 'input' class in order to allow passing environment values to other
	 * classes, but, since we're currently using them only to prepare state
	 * for an internal 'input' function, they can stay private for now.  A
	 * more nuanced approach will likely be required later.
	 */

	class environment {
		static std::stack<environment *> environment_stack;
		static std::set<environment *> environment_set;

		std::map<const char *, const char *, compare_strings> environment_map;

		/*
		 * Internal set operations do not protect any data.
		 */

		void internal_set(const char *name, const char *value) {
			environment_map[name] = value;
		}

		void internal_unset(const char *name) {
			environment_map.erase(name);
		}

		const char *internal_convert_pointer(const void *pointer) {
			int chars = sizeof(void *) * 2 + 3;
			char *c = (char *) malloc(sizeof(char) * chars);

			assert(c);

			if (!c)
				ui::get()->memory_error_location("environment::set_ptr");

			int count = snprintf(c, chars, "%p", pointer);

			assert (count >= 0 && count < chars);

			return c;
		}

		void internal_set_ptr(const char *name, const void *pointer) {
			internal_set(name, internal_convert_pointer(pointer));
		}

		/*
		 * Check for restricted names.
		 */

		int name_ok(const char *name) {
			if (!strcmp(name, "---chain") || !strcmp(name, "---this"))
				return 0;

			return 1;
		}

		void name_check(const char *name) {
			if (!name_ok(name)) {
				fprintf(stderr, "Bad set operation.");
				assert(0);
				exit(1);
			}
		}

	public:

		/*
		 * Get the environment map.
		 */

		std::map<const char *, const char *, compare_strings> &get_map() {
			return environment_map;
		}

		/*
		 * Public set operations restrict valid names.
		 */

		void set(const char *name, const char *value) {
			name_check(name);
			internal_set(name, value);
		}

		void unset(const char *name) {
			name_check(name);
			internal_unset(name);
		}

		void set_ptr(const char *name, const void *pointer) {
			name_check(name);
			internal_set_ptr(name, pointer);
		}
			
		const char *get(const char *name) {
			if (environment_map.count(name) == 0)
				return NULL;

			return environment_map[name];
		}

		/*
		 * Make an environment substructure.  Note that since deep
		 * structures are currently referenced rather than copied when
		 * the stack is pushed, there is no current need for any
		 * chaining mechanism.
		 */
		void make_substructure(const char *name) {
			environment *s = new environment;
			set_ptr(name, s);
			environment_set.insert(s);
		}

		static int is_env(const char *name) {
			void *ptr_value;
			sscanf(name, "%p", &ptr_value);

			/*
			 * Check for bad pointers.
			 */

			if (!environment_set.count((environment *) ptr_value)) {
				return 0;
			}

			return 1;
		}

		const char *get_option_name(const char *name) {
			if (strncmp(name, "0 ", strlen("0 ")))
				return NULL;

			name += strlen("0 ");

			if (!isdigit(name[0]))
				return NULL;

			while (isdigit(name[0]))
				name++;

			if (!isspace(name[0]))
				return NULL;

			while (isspace(name[0]))
				name++;

			if (!isalnum(name[0]))
				return NULL;

			return name;
		}

		int is_option(const char *name) {
			return (get_option_name(name) != NULL);
		}

		int is_arg(const char *name, unsigned int arg) {
			assert (is_option(name));

			int length = strlen(name) + 3 * sizeof(unsigned int);

			char *desired_string = (char *) malloc(sizeof(char) * length);

			snprintf(desired_string, length, "%u %s", arg, name + strlen("0 "));

			int result = environment_map.count(desired_string);

			free(desired_string);

			return result > 0;
		}

		void remove_arg(const char *name, unsigned int arg) {
			assert (is_option(name));
			assert (is_arg(name, arg));

			int length = strlen(name) + 3 * sizeof(unsigned int);

			char *desired_string = (char *) malloc(sizeof(char) * length);

			snprintf(desired_string, length, "%u %s", arg, name + strlen("0 "));

			environment_map.erase(desired_string);

			free(desired_string);
		}

		const char *get_string_arg(const char *name, unsigned int arg) {
			assert (is_option(name));

			int length = strlen(name) + 3 * sizeof(unsigned int);

			char *desired_string = (char *) malloc(sizeof(char) * length);

			snprintf(desired_string, length, "%u %s", arg, name + strlen("0 "));

			const char *result = environment_map[desired_string];

			assert (result);

			free(desired_string);

			return result;
		}

		long int get_long_arg(const char *name, unsigned int arg) {
			assert (is_option(name));

			const char *string = get_string_arg(name, arg);
			char *endptr;
			
			long int result = strtol(string, &endptr, 0);

			if (endptr[0] != '\0') {
				fprintf(stderr, "\n\nError: bad argument in `%s'.\n\n", get_option_name(name));
				exit(1);
			}

			return result;
		}

		int get_int_arg(const char *name, unsigned int arg) {
			return (int) get_long_arg(name, arg);
		}

		unsigned int get_unsigned_arg(const char *name, unsigned int arg) {
			long int result = get_long_arg(name, arg);

			if (result < 0) {
				fprintf(stderr, "\n\nError: bad argument in `%s'.\n\n", get_option_name(name));
				exit(1);
			}

			return (unsigned int) result;
		}

		double get_double_arg(const char *name, unsigned int arg) {
			assert (is_option(name));

			const char *string = get_string_arg(name, arg);
			char *endptr;
			
			double result = strtod(string, &endptr);

			if (endptr[0] != '\0') {
				fprintf(stderr, "\n\nError: bad argument in `%s'.\n\n", get_option_name(name));
				exit(1);
			}

			return result;
		}

		static environment *get_env(const char *name) {

			assert(name);
			
			void *ptr_value;
			sscanf(name, "%p", &ptr_value);

			/*
			 * Check for bad pointers.
			 */

			if (!environment_set.count((environment *) ptr_value)) {
				assert(0);
				fprintf(stderr, "Bad environment pointer.\n");
				exit(1);
			}

			return (environment *) ptr_value;
		}

		/*
		 * Prepend to a list.
		 */
		void prepend(const char *list, const char *element) {
			environment *d = get_env(get(list));
			make_substructure(list);
			get_env(get(list))->set("a", element);
			get_env(get(list))->set_ptr("d", d);
		}

		void prepend_ptr(const char *list, void *ptr) {
			prepend(list, internal_convert_pointer(ptr));
		}

		/*
		 * Clone the environment.
		 */
		environment *clone() {
			environment *e = new environment();

			for (std::map<const char *, const char *, compare_strings>::iterator i = environment_map.begin(); 
					i != environment_map.end(); i++) {

				if (!name_ok(i->first))
					continue;

				if (is_env(i->second)) {
					e->set_ptr(i->first, get_env(i->second)->clone());
				} else {
					e->set(i->first, i->second);
				}
			}

			return e;
		}

		static environment *top() {
			if (environment_stack.empty()) {
				environment_stack.push(new environment);
				environment_set.insert(environment_stack.top());
			}
			return environment_stack.top();
		}

		static void push() {
			environment *e = new environment;

			e->environment_map = environment_stack.top()->environment_map;

			e->internal_set_ptr("---chain", environment_stack.top());
			e->internal_set_ptr("---this", e);
			e->make_substructure("---dup");

			environment_stack.push(e);
			environment_set.insert(e);
		}

		static void dup_second() {
			environment_stack.top()->prepend_ptr("---dup", 
					environment::get_env(environment_stack.top()->get("---chain")));
		}

		static void push_and_dup_output() {
			push();
			dup_second();
		}

		static void pop() {
			assert(!environment_stack.empty());

			/*
			 * Execution environments should never be referenced by
			 * structures further up the call chain, so they can
			 * safely be deleted.  (XXX:  In particular, while
			 * lexical scoping may require copying of execution
			 * environments from lower on the call chain, there is
			 * no obvious reason that a reference should be used in
			 * this case; a shallow copy should be used instead.)
			 */

			environment_set.erase(environment_stack.top());
			delete environment_stack.top();

			environment_stack.pop();
		}

		/*
		 * Set with duplication.
		 */

		void set_with_dup(const char *name, const char *value) {
			set(name, value);

			if (!get("---dup"))
				return;

			environment *dup_item = get_env(get("---dup"));

			assert (dup_item);

			while (dup_item->get("a")) {
				get_env(dup_item->get("a"))->set_with_dup(name, value);
				assert(dup_item->get("d"));
				dup_item = get_env(dup_item->get("d"));
				assert(dup_item);
			}
		}
	};

	/*
	 * Read tokens from a stream.
	 */
	class token_reader {
	public:
		/*
		 * Get the next token
		 */
		virtual const char *get() = 0;

		/*
		 * Peek at the next token.
		 */

		virtual const char *peek() = 0;

		/*
		 * Divert the stream until the next occurrence of TOKEN.
		 */
		virtual token_reader *divert(const char *open_token, const char *close_token) = 0;

		virtual int expects_exactly_one_option(void) {
			return 0;
		}

		virtual ~token_reader() {
		}
	};

	class argument_parsing_token_reader : public token_reader {
		const char *index;
		const char *separators;
	public:
		argument_parsing_token_reader(const char *s) {
			index = s;
			separators = "=";
		}

		int expects_exactly_one_option(void) {
			return 1;
		}

		virtual const char *get() {
			int length = strcspn(index, separators);

			if (length == 0)
				return NULL;

			const char *result = strndup(index, length);
			index += length;

			if (strspn(index, separators) >= 1)
				index++;

			separators = ",";

			return result;
		}

		virtual const char *peek() {
			int length = strcspn(index, separators);

			if (length == 0)
				return NULL;

			const char *result = strndup(index, length);

			return result;
		}

		virtual token_reader *divert(const char *open_token, const char *close_token) {
			assert(0);
			return NULL;
		}
	};

	class cstring_token_reader : public token_reader {
		const char *separators;
		const char *string;
		int ephemeral;

		cstring_token_reader(const char *s, int ephemeral) {
			assert(ephemeral == 1);

			separators = "\n \t";
			string = s;
			this->ephemeral = 1;
		}

	public:
		cstring_token_reader(const char *s) {
			separators = "\n \t";
			string = s;
			ephemeral = 0;
		}

		const char *get() {

			string += strspn(string, separators);

			size_t length_to_next = strcspn(string, separators);

			if (length_to_next == 0)
				return NULL;

			const char *result = strndup(string, length_to_next);

			string += length_to_next;

			return result;
		}

		const char *peek() {
			string += strspn(string, separators);

			size_t length_to_next = strcspn(string, separators);

			if (length_to_next == 0)
				return NULL;

			return strndup(string, length_to_next);
		}

		cstring_token_reader *divert(const char *open_token, const char *close_token) {
			/*
			 * This function might be broken.
			 */

			assert(0);

			int search = 0;
			int next = strcspn(string, separators);
			int depth = 0;

			while (*(string + search) != '\0' && 
			       (depth || strcmp(close_token, (string + search)))) {
				if (!strcmp(close_token, (string + search)))
					depth--;
				if (!strcmp(open_token, (string + search)))
					depth++;
				search = next;
				next = strcspn((string + next), separators);
			}

			if (*(string + search) == '\0') {
				fprintf(stderr, "Parse error: End of scope not found.");
				exit(1);
			}

			cstring_token_reader *result = new cstring_token_reader(strndup(string, search), 1);

			string += search;

			/*
			 * Eat the closing token.
			 */

			get();

			return result;
		}

		~cstring_token_reader() {
			if (ephemeral)
				free((void *) string);
		}
	};

	class cli_token_reader : public token_reader {

		int arg_index;
		int argc;
		const char **argv;

	public:
		cli_token_reader(int c, const char *v[]) {
			argc = c;
			argv = v;
			arg_index = 0;
		}

		const char *get() {

			if (arg_index < argc)
				return argv[arg_index++];
			else
				return NULL;

		}

		const char *peek() {

			if (arg_index < argc)
				return argv[arg_index];
			else
				return NULL;

		}

		cli_token_reader *divert(const char *open_token, const char *close_token) {
			int search = 0;
			int depth = 0;

			while (arg_index + search < argc 
			    && (depth || strcmp(argv[arg_index + search], close_token))) {
				if (!strcmp(close_token, argv[arg_index + search]))
					depth--;
				if (!strcmp(open_token, argv[arg_index + search]))
					depth++;
				search++;
			}

			if (arg_index + search == argc) {
				fprintf(stderr, "Parse error: end of scope not found.\n");
				exit(1);
			}

			cli_token_reader *result = new cli_token_reader(search, argv + arg_index);

			arg_index += search;

			/*
			 * Eat the closing token.
			 */

			get();

			return result;
		}

	};

	struct simple_option {
		const char *name;
		const char *map_name;
		const char *map_value;
		int arg_count;
		int multi;
	};

	static const char *supported_nonglobal_option_table[];
	static const char *focus_prefixes[];
	static simple_option simple_option_table[];

	static int option_name_match(const char *unadorned, const char *token, int require_ornamentation = 1) {
		int strip_max = 2;

		if (!strcmp(unadorned, token) && !require_ornamentation)
			return 1;

		while (token[0] == '-' && strip_max) {
			token++;
			strip_max--;
			if (!strcmp(unadorned, token))
				return 1;
		}

		return 0;
	}

	static int is_scope_operator(const char *string) {
		if (!strcmp("{", string)
		 || !strcmp("}", string)
		 || !strcmp("[", string)
		 || !strcmp("]", string)
		 || !strcmp("<", string)
		 || !strcmp(">", string))
			return 1;

		return 0;
	}

	static const char *option_name_gen(const char *unadorned, const char *map_name, int arg_num, int multi) {
		static unsigned int multi_counter = 0;

		if (map_name) {
			unadorned = map_name;
		}

		int length = (strlen(unadorned) + sizeof(unsigned int) * 3 + sizeof(int) * 3 + 2) + 1;

		char *result = (char *) malloc(sizeof(char) * length);

		assert (result);

		if (!multi) {
			snprintf(result, length, "%u 0 %s", arg_num, unadorned);
		} else {

			/*
			 * XXX: This assumes that generating calls for
			 * options other than 0 exist in the same
			 * multiplicity group as the most recently 
			 * generated 0-option multiplicity.
			 */

			if (arg_num == 0)
				multi_counter++;

			snprintf(result, length, "%u %u %s", arg_num, multi_counter, unadorned);
		}

		return result;
	}

	static environment *genv;

	static const char *get_next(token_reader *tr, const char *option_name) {
		const char *argument = tr->get();

		if (argument == NULL) {
			fprintf(stderr, "\n\nError: not enough arguments for `%s'.\n\n", option_name);
			exit(1);
		}

		return argument;
	}

	static int table_contains(const char **haystack, const char *needle, int prefix_length = 0) {

		if (needle == NULL)
			return 0;

		while (*haystack != NULL) {
			if (prefix_length == 0 && !strcmp(*haystack, needle))
				return 1;
			if (prefix_length > 0 && !strncmp(*haystack, needle, prefix_length))
				return 1;
			haystack++;
		}

		return 0;
	}

	static int option_is_identical(environment *a, environment *b, const char *option_name) {
		if (!a->is_option(option_name) || !b->is_option(option_name))
			return 0;
		
		int option_number = 0;

		while (a->is_arg(option_name, option_number) || b->is_arg(option_name, option_number)) {
			if (!a->is_arg(option_name, option_number)
			 || !b->is_arg(option_name, option_number))
				return 0;

			const char *a_str = a->get_string_arg(option_name, option_number);
			const char *b_str = b->get_string_arg(option_name, option_number);

			if (strcmp(a_str, b_str))
				return 0;

			option_number++;
		}

		return 1;
	}

	static void remove_option(environment *a, const char *option_name) {
		assert(a->is_option(option_name));

		int option_number = 0;

		while (a->is_arg(option_name, option_number)) {
			a->remove_arg(option_name, option_number);
			option_number++;
		}
	}

	static void remove_nonglobals(environment *a) {
		assert(a);

		std::stack<const char *> removal_stack;

		for (std::map<const char *, const char *, compare_strings>::iterator i = a->get_map().begin();
				i != a->get_map().end(); i++) {

			if (!a->is_option(i->first))
				continue;

			if (!table_contains(supported_nonglobal_option_table, a->get_option_name(i->first)))
				continue;

			removal_stack.push(i->first);
		}

		while (!removal_stack.empty()) {
			remove_option(a, removal_stack.top());
			removal_stack.pop();
		}
	}

	static void option_intersect(environment *a, environment *b) {
		assert(a);
		assert(b);

		std::stack<const char *> removal_stack;

		for (std::map<const char *, const char *, compare_strings>::iterator i = a->get_map().begin();
				i != a->get_map().end(); i++) {

			if (!a->is_option(i->first))
				continue;

			if (option_is_identical(a, b, i->first))
				continue;

			removal_stack.push(i->first);
		}

		while (!removal_stack.empty()) {
			remove_option(a, removal_stack.top());
			removal_stack.pop();
		}
	}

	static void option_difference(environment *a, environment *b) {
		assert(a);
		assert(b);

		std::stack<const char *> removal_stack;

		for (std::map<const char *, const char *, compare_strings>::iterator i = a->get_map().begin();
				i != a->get_map().end(); i++) {

			if (!a->is_option(i->first))
				continue;

			if (!option_is_identical(a, b, i->first))
				continue;

			removal_stack.push(i->first);
		}

		while (!removal_stack.empty()) {
			remove_option(a, removal_stack.top());
			removal_stack.pop();
		}
	}

	static void evaluate_stream(token_reader *tr, 
			std::vector<std::pair<const char *, environment *> > *files) {
		const char *token;
		int end_of_options = 0;

		while ((token = tr->get())) {

			/*
			 * Check for nesting
			 */

			if (!strcmp(token, "{") && !end_of_options) {
				environment::push_and_dup_output();
				token_reader *tr_nest = tr->divert("{", "}");
				evaluate_stream(tr_nest, files);
				delete tr_nest;
				environment::pop();
			} else if (!strcmp(token, "[") && !end_of_options) {
				global_options = 0;
				environment::push();
				token_reader *tr_nest = tr->divert("[", "]");
				evaluate_stream(tr_nest, files);
				delete tr_nest;
				environment::pop();
			} else if (!strcmp(token, "<") && !end_of_options) {
				environment *dup_list = environment::get_env(environment::top()->get("---dup"));
				assert (dup_list != NULL);
				dup_list = dup_list->clone();

				environment::dup_second();
				token_reader *tr_nest = tr->divert("<", ">");
				evaluate_stream(tr_nest, files);
				delete tr_nest;

				environment::top()->set_ptr("---dup", dup_list);
			}

			/*
			 * Check for non-whitespace argument separators
			 */

			else if (!end_of_options && token && token[0] == '-' && strchr(token, '=')) {
				environment::push_and_dup_output();
				token_reader *tr_nest = new argument_parsing_token_reader(token);
				evaluate_stream(tr_nest, files);
				delete tr_nest;
				environment::pop();
			}

			/*
			 * Trap the end-of-option indicator.
			 */
			
			else if (!strcmp(token, "--")) {
				global_options = 0;
				end_of_options = 1;
			}

			/*
			 * Check for options and filenames
			 */
			
			else {
				/*
				 * Handle filenames.
				 */

				if (strncmp("-", token, strlen("-")) || end_of_options) {

					assert(files);

					global_options = 0;
					files->push_back(std::pair<const char *, environment *>(strdup(token), 
								environment::top()->clone()));

					if (tr->expects_exactly_one_option() && tr->get()) {
						fprintf(stderr, "\n\nError: Too many arguments for `%s'.\n\n", token);
						exit(1);
					}

					continue;
				}

				/*
				 * Handle focus option.
				 */

				if (option_name_match("focus", token)) {

					environment *target;
					target = environment::top();

					target->set_with_dup(option_name_gen("focus", NULL, 0, 0), "1");

					const char *option = get_next(tr, "focus");

					target->set_with_dup(option_name_gen("focus", NULL, 1, 0), option);

					if (!strcmp(option, "d")) {
						target->set_with_dup(option_name_gen("focus", NULL, 2, 0), 
								get_next(tr, "focus"));
					} else if (!strcmp(option, "p")) {
						target->set_with_dup(option_name_gen("focus", NULL, 2, 0), 
								get_next(tr, "focus"));
						target->set_with_dup(option_name_gen("focus", NULL, 3, 0), 
								get_next(tr, "focus"));
					} else 
						bad_arg("focus");

					int arg = 0;

					while (table_contains(focus_prefixes, tr->peek(), 3)) {
						target->set_with_dup(option_name_gen("focus", NULL, 4 + arg, 0), 
								get_next(tr, "focus"));
						arg++;
					}

					continue;
				}

				/*
				 * Handle simple options.
				 */

				int found_option = 0;
				for (int i = 0; simple_option_table[i].name; i++) {
					if (!option_name_match(simple_option_table[i].name, token))
						continue;
					
					/*
					 * Handle the match case.
					 */

					found_option = 1;

					/*
					 * Determine which environment should be modified
					 */

					environment *target;
					target = environment::top();
					
					/*
					 * Store information required for
					 * handling the local case later.
					 */

					const char *map_value = "1";

					if (simple_option_table[i].map_value) {
						map_value = simple_option_table[i].map_value;
					} else if (simple_option_table[i].map_name) {
						map_value = simple_option_table[i].name;
					}

					target->set_with_dup(option_name_gen(simple_option_table[i].name,
								simple_option_table[i].map_name,
								0,
								simple_option_table[i].multi),
							map_value);

					for (int j = 0; j < simple_option_table[i].arg_count; j++) {
						const char *option = tr->get();

						if (option == NULL) {
							fprintf(stderr, "\n\nError: not enough options for `%s'.\n\n", token);
							exit(1);
						}

						/*
						 * Reject scope operators as options,
						 * at least for now.
						 */

						if (is_scope_operator(option)) {
							fprintf(stderr, "\n\nError: illegal argument to `%s'.\n\n", token);
							exit(1);
						}

						target->set_with_dup(option_name_gen(simple_option_table[i].name,
									simple_option_table[i].map_name,
									j + 1,
									simple_option_table[i].multi),
								option);
					}
				}

				/*
				 * Trap illegal options.
				 */

				if (!found_option) 
					ui::get()->illegal_option(token);
			}

			if (tr->expects_exactly_one_option() && tr->get()) {
				fprintf(stderr, "\n\nError: Too many arguments for `%s'.\n\n", token);
				exit(1);
			}
		}
	}

public:

	/*
	 * Libale Sequence structure and functions.
	 */

	struct seq_struct {
		ale_exclusion_list ex;
		ale_render *ochain;
		int oc_count;
	};

	ale_image seq_file(int n, ale_sequence s) {
		return d2::image_rw::open_simple(n);
	}

	ale_trans seq_trans(int n, ale_sequence s) {
		return align::of(n);
	}

	ale_exclusion_list seq_ex(int n, ale_sequence s) {
		return ((seq_struct *) ale_sequence_data(s))->ex;
	}

	void seq_step(int n, ale_sequence s) {
		seq_struct *seq_data = (seq_struct *) ale_sequence_data(s);

	}
	
	/*
	 * Input handler.
	 *
	 * Does one of two things:
	 *
	 * (1) Output version information if called with '--version'
	 *
	 * (2) Read options and file arguments, and if the arguments are correct, 
	 * write output.  If an error is detected, print the usage statement.
	 *
	 */

	static void handle(int argc, const char *argv[], const char *package, const char *short_version, const char *version) {

		/*
		 * Initialize help object
		 */
		
		help hi(package, argv[0], short_version);

		/*
		 * Output version information if --version appears
		 * on the command line.
		 */

		if (arg_count(argc, argv, "--version")) {
			/*
			 * Output the version
			 */

			fprintf(stdout, "%s", version);

			/*
			 * Output relevant environment variables
			 */

			fprintf(stdout, "Environment:\n");

			const char *env_names[] = {
				"ALE_BIN",
				"DCRAW",
				"EXIF_UTILITY",
				"ALE_COUNT_THREADS",
				"PAGER",
				NULL
			};

			for (int i = 0; env_names[i]; i++) {
				char *value = getenv(env_names[i]);

				fprintf(stdout, "   %s=%s\n",
					env_names[i], value ? value : "");
			}

			return;
		}

		/*
		 * Handle help options
		 */

		if  (arg_prefix_count(argc, argv, "--h"))
		for (int i = 1; i < argc; i++) {
			int all = !strcmp(argv[i], "--hA");
			int is_help_option = !strncmp(argv[i], "--h", strlen("--h"));
			int found_help = 0;

			if (!strcmp(argv[i], "--hu") || all)
				hi.usage(), found_help = 1;
			if (!strcmp(argv[i], "--hq") || all)
				hi.defaults(), found_help = 1;
			if (!strcmp(argv[i], "--hf") || all)
				hi.file(), found_help = 1;
			if (!strcmp(argv[i], "--he") || all)
				hi.exclusion(), found_help = 1;
			if (!strcmp(argv[i], "--ha") || all)
				hi.alignment(), found_help = 1;
			if (!strcmp(argv[i], "--hr") || all)
				hi.rendering(), found_help = 1;
			if (!strcmp(argv[i], "--hx") || all)
				hi.exposure(), found_help = 1;
			if (!strcmp(argv[i], "--ht") || all)
				hi.tdf(), found_help = 1;
			if (!strcmp(argv[i], "--hl") || all)
				hi.filtering(), found_help = 1;
			if (!strcmp(argv[i], "--hd") || all)
				hi.device(), found_help = 1;
			if (!strcmp(argv[i], "--hi") || all)
				hi.interface(), found_help = 1;
			if (!strcmp(argv[i], "--hv") || all)
				hi.visp(), found_help = 1;
			if (!strcmp(argv[i], "--hc") || all)
				hi.cp(), found_help = 1;
			if (!strcmp(argv[i], "--h3") || all)
				hi.d3(), found_help = 1;
			if (!strcmp(argv[i], "--hs") || all)
				hi.scope(), found_help = 1;
			if (!strcmp(argv[i], "--hp") || all)
				hi.process(), found_help = 1;
			if (!strcmp(argv[i], "--hz") || all)
				hi.undocumented(), found_help = 1;

			if (is_help_option && !found_help)
				hi.usage();

			/*
			 * Check for the end-of-options marker, a non-option argument,
			 * or the end of arguments.  In all of these cases, we exit.
			 */

			if (!strcmp(argv[i], "--")
			 || strncmp(argv[i], "--", strlen("--"))
			 || i == argc - 1)
				return;
		}

		/*
		 * Undocumented projective transformation utility
		 */

		if (arg_count(argc, argv, "--ptcalc") > 0) {
			fprintf(stderr, "\n\n*** Warning: this feature is not documented ***\n\n");
			printf("Enter: w h tlx tly blx bly brx bry trx try x y\n\n");
		
			double w, h, tlx, tly, blx, bly, brx, bry, trx, tr_y, x, y;

			printf("> ");

			if (scanf("%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
				&w, &h, &tlx, &tly, &blx, &bly, &brx, &bry, &trx, &tr_y, &x, &y) != 12) {

				fprintf(stderr, "Error reading input.\n");
				exit(1);
			}

			d2::image *i = d2::new_image_ale_real((int)h, (int)w, 3);
			d2::transformation t = d2::transformation::gpt_identity(i, 1);
			d2::point q[4] = {
				d2::point(tly, tlx),
				d2::point(bly, blx),
				d2::point(bry, brx),
				d2::point(tr_y, trx)
			};
			t.gpt_set(q);

			d2::point a(y, x), b;

			b = t.transform_scaled(a);

			printf("TRANSFORM t(a): (%f, %f)\n", (double) b[1], (double) b[0]);

			b = t.scaled_inverse_transform(a);

			printf("INVERSE t^-1(a): (%f, %f)\n", (double) b[1], (double) b[0]);

			exit(0);
		}

		/*
		 * Thread initialization.
		 */

		thread::init();

		/*
		 * Flags and variables
		 */

		double scale_factor = 1;
		double vise_scale_factor = 1;
#if 0
		double usm_multiplier = 0.0;
#endif
		int extend = 0;
		struct d2::tload_t *tload = NULL;
		struct d2::tsave_t *tsave = NULL;
		struct d3::tload_t *d3_tload = NULL;
		struct d3::tsave_t *d3_tsave = NULL;
		int ip_iterations = 0;
		int ip_use_median = 0;
		double ipwl = 0;
		enum { psf_linear, psf_nonlinear, psf_N };
		const char *psf[psf_N] = {NULL, NULL};
		const char *device = NULL;
		int psf_match = 0;
		double psf_match_args[6];
		int inc = 0;
		int exposure_register = 1;
		const char *wm_filename = NULL;
		int wm_offsetx = 0, wm_offsety = 0;
		double cx_parameter = 1;
		double *d3px_parameters = NULL;
		int d3px_count = 0;
		d2::exclusion *ex_parameters = NULL;
		int ex_count = 0;
		int ex_show = 0;
		d2::render *achain;
		const char *achain_type = "triangle:2";
		const char *afilter_type = "internal";
		d2::render **ochain = NULL;
		const char **ochain_names = NULL;
		const char **ochain_types = NULL;
		const char *d3chain_type = NULL;
		int oc_count = 0;
		const char **visp = NULL;
		int vise_count = 0;
		const char **d3_output = NULL;
		const char **d3_depth = NULL;
		unsigned int d3_count = 0;
		double user_view_angle = 0;
		int user_bayer = IMAGE_BAYER_DEFAULT;
		d2::pixel exp_mult = d2::pixel(1, 1, 1);
		std::map<const char *, d3::pt> d3_output_pt;
		std::map<const char *, d3::pt> d3_depth_pt;
		double cache = 256; /* MB */

		/*
		 * dchain is ochain[0].
		 */

		ochain = (d2::render **) local_realloc(ochain, 
					(oc_count + 1) * sizeof(d2::render *));
		ochain_names = (const char **) local_realloc((void *)ochain_names, 
					(oc_count + 1) * sizeof(const char *));
		ochain_types = (const char **) local_realloc((void *)ochain_types, 
					(oc_count + 1) * sizeof(const char *));

		ochain_types[0] = "sinc*lanc:8";

		oc_count = 1;

		/*
		 * Handle default settings
		 */

		if (arg_prefix_count(argc, argv, "--q") > 0)
			ui::get()->error("Default settings --q* are no longer recognized.");

#define FIXED16 4
#if ALE_COLORS == FIXED16
		const char *defaults = 
   			"--dchain auto:triangle:2,fine:box:1,triangle:2 "
			"--achain triangle:2 "
			"--ips 0 "
			"--3d-chain fine:triangle:2,fine:gauss:0.75,triangle:2 ";
#else
		const char *defaults = 
   			"--dchain auto:triangle:2,fine:box:1,triangle:2 "
			"--achain triangle:2 "
			"--ips 1 "
			"--3d-chain fine:triangle:2,fine:gauss:0.75,triangle:2 ";
#endif
#undef FIXED16

		token_reader *default_reader = new cstring_token_reader(defaults);

		evaluate_stream(default_reader, NULL);

		/*
		 * Set basic program information in the environment.
		 */

		environment::top()->set_with_dup("---package", package);
		environment::top()->set_with_dup("---short-version", short_version);
		environment::top()->set_with_dup("---version", version);
		environment::top()->set_with_dup("---invocation", argv[0]);

		/*
		 * Initialize the top-level token-reader and generate
		 * an environment variable for it.
		 */

		token_reader *tr = new cli_token_reader(argc - 1, argv + 1);
		environment::top()->set_ptr("---token-reader", tr);

		/*
		 * Evaluate the command-line arguments to generate environment
		 * structures.
		 */

		std::vector<std::pair<const char *, environment *> > files;

		evaluate_stream(tr, &files);

		/*
		 * If there are fewer than two files, then output usage information.
		 */

		if (files.size() < 2) {
			hi.usage();
			exit(1);
		}

		/*
		 * Extract the global environment and check non-globals
		 * against a list of supported non-global options.
		 */

		genv = files[0].second->clone();

		remove_nonglobals(genv);
		
		for (unsigned int i = 0; i < files.size(); i++) {
			option_intersect(genv, files[i].second);
		}

		for (unsigned int i = 0; i < files.size(); i++) {
			option_difference(files[i].second, genv);

			for (std::map<const char *, const char *>::iterator j = files[i].second->get_map().begin();
					j != files[i].second->get_map().end(); j++) {

				environment *env = files[i].second;

				if (!env->is_option(j->first))
					continue;

				const char *option_name = env->get_option_name(j->first);

				if (!table_contains(supported_nonglobal_option_table, option_name)) {
					fprintf(stderr, "\n\nError: option `%s' must be applied globally.", option_name);
					fprintf(stderr, "\n\nHint:  Move option `%s' prior to file and scope operators.\n\n", 
							option_name);
					exit(1);
				}
			}
		}

		/*
		 * Iterate through the global environment,
		 * looking for options.
		 */

		for (std::map<const char *, const char *>::iterator i = genv->get_map().begin();
				i != genv->get_map().end(); i++) {

			environment *env = genv;

			if (!env->is_option(i->first))
				continue;

			const char *option_name = env->get_option_name(i->first);

			if (!strcmp(option_name, "default")) {
				/*
				 * Do nothing.  Defaults have already been set.
				 */
			} else if (!strcmp(option_name, "bpc")) {
				if (!strcmp(env->get_string_arg(i->first, 0), "8bpc"))
					d2::image_rw::depth8();
				else if (!strcmp(env->get_string_arg(i->first, 0), "16bpc"))
					d2::image_rw::depth16();
				else
					assert(0);
			} else if (!strcmp(option_name, "format")) {
				if (!strcmp(env->get_string_arg(i->first, 0), "plain"))
					d2::image_rw::ppm_plain();
				else if (!strcmp(env->get_string_arg(i->first, 0), "raw"))
					d2::image_rw::ppm_raw();
				else if (!strcmp(env->get_string_arg(i->first, 0), "auto"))
					d2::image_rw::ppm_auto();
				else
					assert(0);
			} else if (!strcmp(option_name, "align")) {
				if (!strcmp(env->get_string_arg(i->first, 0), "align-all"))
					d2::align::all();
				else if (!strcmp(env->get_string_arg(i->first, 0), "align-green"))
					d2::align::green();
				else if (!strcmp(env->get_string_arg(i->first, 0), "align-sum"))
					d2::align::sum();
				else
					assert(0);
			} else if (!strcmp(option_name, "transformation")) {
				if (!strcmp(env->get_string_arg(i->first, 0), "translation"))
					d2::align::class_translation();
				else if (!strcmp(env->get_string_arg(i->first, 0), "euclidean"))
					d2::align::class_euclidean();
				else if (!strcmp(env->get_string_arg(i->first, 0), "projective"))
					d2::align::class_projective();
				else
					assert(0);
			} else if (!strcmp(option_name, "transformation-default")) {
				if (!strcmp(env->get_string_arg(i->first, 0), "identity"))
					d2::align::initial_default_identity();
				else if (!strcmp(env->get_string_arg(i->first, 0), "follow"))
					d2::align::initial_default_follow();
				else
					assert(0);
			} else if (!strcmp(option_name, "perturb")) {
				if (!strcmp(env->get_string_arg(i->first, 0), "perturb-output"))
					d2::align::perturb_output();
				else if (!strcmp(env->get_string_arg(i->first, 0), "perturb-source"))
					d2::align::perturb_source();
				else
					assert(0);
			} else if (!strcmp(option_name, "fail")) {
				if (!strcmp(env->get_string_arg(i->first, 0), "fail-optimal"))
					d2::align::fail_optimal();
				else if (!strcmp(env->get_string_arg(i->first, 0), "fail-default"))
					d2::align::fail_default();
				else
					assert(0);
			} else if (!strcmp(option_name, "profile")) {
				ui::set_profile();
			} else if (!strcmp(option_name, "extend")) {
				if (env->get_int_arg(i->first, 0))
					extend = 1;
				else
					extend = 0;
			} else if (!strcmp(option_name, "oc")) {
				if (env->get_int_arg(i->first, 0))
					d3::scene::oc();
				else 
					d3::scene::no_oc();
			} else if (!strcmp(option_name, "focus")) {

				double one = +1;
				double zero = +0;
				double inf = one / zero;

				assert (isinf(inf) && inf > 0);

				/*
				 * Focus type
				 */

				unsigned int type = 0;
				double distance = 0;
				double px = 0, py = 0;

				if (!strcmp(env->get_string_arg(i->first, 1), "d")) {

					type = 0;

					distance = env->get_double_arg(i->first, 2);

				} else if (!strcmp(env->get_string_arg(i->first, 1), "p")) {

					type = 1;

					px = env->get_double_arg(i->first, 2);
					py = env->get_double_arg(i->first, 3);

				} else {
					bad_arg(option_name);
				}

				/*
				 * Options
				 */

				unsigned int ci = 0;
				double fr = 0;
				double ht = 0;
				double vt = 0;
				double sd = 0;
				double ed = inf;
				double sx = -inf;
				double ex = inf;
				double sy = -inf;
				double ey = inf;
				double ap = 3;
				unsigned int sc = 3;
				unsigned int fs = 0;
				unsigned int sr = 0;

				for (int arg_num = 4; env->is_arg(i->first, arg_num); arg_num++) {
					const char *option = env->get_string_arg(i->first, arg_num);
					if (!strncmp(option, "ci=", 3)) {
						if(sscanf(option + 3, "%u", &ci) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "fr=", 3)) {
						if(sscanf(option + 3, "%lf", &fr) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "ht=", 3)) {
						if(sscanf(option + 3, "%lf", &ht) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "vt=", 3)) {
						if(sscanf(option + 3, "%lf", &vt) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "sy=", 3)) {
						if(sscanf(option + 3, "%lf", &sy) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "ey=", 3)) {
						if(sscanf(option + 3, "%lf", &ey) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "sx=", 3)) {
						if(sscanf(option + 3, "%lf", &sx) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "ex=", 3)) {
						if(sscanf(option + 3, "%lf", &ex) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "sd=", 3)) {
						if(sscanf(option + 3, "%lf", &sd) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "ed=", 3)) {
						if(sscanf(option + 3, "%lf", &ed) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "ap=", 3)) {
						if(sscanf(option + 3, "%lf", &ap) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "sc=", 3)) {
						if(sscanf(option + 3, "%u", &sc) != 1)
							bad_arg("--focus");
					} else if (!strncmp(option, "sr=", 3)) {
						if (!strcmp(option, "sr=aperture")) {
							sr = 0;
						} else if (!strcmp(option, "sr=pixel")) {
							sr = 1;
						} else
							bad_arg("--focus");

					} else if (!strncmp(option, "fs=", 3)) {
						if (!strcmp(option, "fs=mean")) {
							fs = 0;
						} else if (!strcmp(option, "fs=median")) {
							fs = 1;
						} else 
							bad_arg("--focus");
					} else
						bad_arg("--focus");
				}

				d3::focus::add_region(type, distance, px, py, ci, fr, ht, vt, sd, ed, sx, ex, sy, ey, ap, sc, fs, sr);

			} else if (!strcmp(option_name, "3ddp") || !strcmp(option_name, "3dvp")) {
				d2::align::keep();

				/*
				 * Unsupported configurations
				 */

				if (ip_iterations)
					unsupported::fornow("3D modeling with Irani-Peleg rendering");

#if 0
				if (usm_multiplier)
					unsupported::fornow("3D modeling with unsharp mask");
#endif

				/*
				 * Initialize if necessary
				 *
				 * Note: because their existence is checked as an
				 * indicator of the presence of 3D arguments, we
				 * initialize these structures here.
				 */

				if (d3_output == NULL) {
					d3_count  = argc;
					d3_output = (const char **) calloc(d3_count, sizeof(char *));
					d3_depth = (const char **) calloc(d3_count, sizeof(char *));
				}

				unsigned int width, height;
				double view_angle;
				double x, y, z;
				double P, Y, R;

				width = env->get_unsigned_arg(i->first, 1);
				height = env->get_unsigned_arg(i->first, 2);
				view_angle = env->get_double_arg(i->first, 3);
				x = env->get_double_arg(i->first, 4);
				y = env->get_double_arg(i->first, 5);
				z = env->get_double_arg(i->first, 6);
				P = env->get_double_arg(i->first, 7);
				Y = env->get_double_arg(i->first, 8);
				R = env->get_double_arg(i->first, 9);

				view_angle *= M_PI / 180;
				P *= M_PI / 180;
				Y *= M_PI / 180;
				R *= M_PI / 180;

				d2::transformation t = 
					d2::transformation::eu_identity();
				t.set_domain(height, width);
				d3::pt _pt(t, d3::et(y, x, z, Y, P, R), view_angle);
				
				if (!strcmp(option_name, "3dvp")) {
					d3_output_pt[env->get_string_arg(i->first, 10)] = _pt;
				} else if (!strcmp(option_name, "3ddp")) {
					d3_depth_pt[env->get_string_arg(i->first, 10)] = _pt;
				} else {
					assert(0);
				}
			} else if (!strcmp(option_name, "3dv")) {
				d2::align::keep();

				unsigned int frame_no;

				/*
				 * Unsupported configurations
				 */

				if (ip_iterations)
					unsupported::fornow("3D modeling with Irani-Peleg rendering");

#if 0
				if (usm_multiplier)
					unsupported::fornow("3D modeling with unsharp mask");
#endif

				/*
				 * Initialize if necessary
				 */

				if (d3_output == NULL) {
					d3_count  = argc;
					d3_output = (const char **) calloc(d3_count, sizeof(char *));
					d3_depth = (const char **) calloc(d3_count, sizeof(char *));
				}

				frame_no = env->get_int_arg(i->first, 1);

				if (frame_no >= d3_count)
					ui::get()->error("--3dv argument 0 is too large");

				if (d3_output[frame_no] != NULL) {
					unsupported::fornow ("Writing a single 3D view to more than one output file");
				}

				d3_output[frame_no] = env->get_string_arg(i->first, 2);

			} else if (!strcmp(option_name, "3dd")) {
				d2::align::keep();

				unsigned int frame_no;

				/*
				 * Unsupported configurations
				 */

				if (ip_iterations)
					unsupported::fornow("3D modeling with Irani-Peleg rendering");

#if 0
				if (usm_multiplier)
					unsupported::fornow("3D modeling with unsharp mask");
#endif

				/*
				 * Initialize if necessary
				 */

				if (d3_output == NULL) {
					d3_count  = argc;
					d3_output = (const char **) calloc(d3_count, sizeof(char *));
					d3_depth = (const char **) calloc(d3_count, sizeof(char *));
				}

				frame_no = env->get_int_arg(i->first, 1);

				if (frame_no >= d3_count)
					ui::get()->error("--3dd argument 0 is too large");

				if (d3_depth[frame_no] != NULL) {
					unsupported::fornow ("Writing a single frame's depth info to more than one output file");
				}

				d3_depth[frame_no] = env->get_string_arg(i->first, 2);

			} else if (!strcmp(option_name, "view-angle")) {
				user_view_angle = env->get_double_arg(i->first, 1) * M_PI / 180;
			} else if (!strcmp(option_name, "cpf-load")) {
				d3::cpf::init_loadfile(env->get_string_arg(i->first, 1));
			} else if (!strcmp(option_name, "accel")) {
				if (!strcmp(env->get_string_arg(i->first, 1), "gpu"))
					accel::set_gpu();
				else if (!strcmp(env->get_string_arg(i->first, 1), "cpu"))
					accel::set_cpu();
				else if (!strcmp(env->get_string_arg(i->first, 1), "accel"))
					accel::set_accel();
				else if (!strcmp(env->get_string_arg(i->first, 1), "auto"))
					accel::set_auto();
				else {
					fprintf(stderr, "Error: Unknown acceleration type '%s'\n", 
							env->get_string_arg(i->first, 1));
					exit(1);
				}
			} else if (!strcmp(option_name, "ui")) {
				if (!strcmp(env->get_string_arg(i->first, 1), "stream"))
					ui::set_stream();
				else if (!strcmp(env->get_string_arg(i->first, 1), "tty"))
					ui::set_tty();
				else if (!strcmp(env->get_string_arg(i->first, 1), "log"))
					ui::set_log();
				else if (!strcmp(env->get_string_arg(i->first, 1), "quiet"))
					ui::set_quiet();
				else {
					fprintf(stderr, "Error: Unknown user interface type '%s'\n", 
							env->get_string_arg(i->first, 1));
					exit(1);
				}
			} else if (!strcmp(option_name, "3d-fmr")) {
				d3::scene::fmr(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "3d-dmr")) {
				d3::scene::dmr(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "et")) {
				d3::scene::et(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "st")) {
				d3::cpf::st(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "di-lower")) {
				d3::scene::di_lower(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "rc")) {
				d3::scene::rc(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "do-try")) {
				d3::scene::do_try(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "di-upper")) {
				d3::scene::di_upper(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "fc")) {
				d3::scene::fc(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "ecm")) {
				unsupported::discontinued("--ecm <x>");
			} else if (!strcmp(option_name, "acm")) {
				unsupported::discontinued("--acm <x>");
			} else if (!strcmp(option_name, "def-nn")) {
				d2::image_rw::def_nn(env->get_double_arg(i->first, 1));

				if (env->get_double_arg(i->first, 1) > 2) {
					fprintf(stderr, "\n\n*** Warning: --def-nn implementation is currently "
							     "inefficient for large radii. ***\n\n");
				}

			} else if (!strcmp(option_name, "fx")) {
				d3::scene::fx(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "tcem")) {
				d3::scene::tcem(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "oui")) {
				d3::scene::oui(env->get_unsigned_arg(i->first, 1));
			} else if (!strcmp(option_name, "pa")) {
				d3::scene::pa(env->get_unsigned_arg(i->first, 1));
			} else if (!strcmp(option_name, "pc")) {
				d3::scene::pc(env->get_string_arg(i->first, 1));
			} else if (!strcmp(option_name, "cw")) {
				d2::align::certainty_weighted(env->get_unsigned_arg(i->first, 0));
			} else if (!strcmp(option_name, "wm")) {
				if (wm_filename != NULL)
					ui::get()->error("only one weight map can be specified");

				wm_filename = env->get_string_arg(i->first, 1);
				wm_offsetx = env->get_int_arg(i->first, 2);
				wm_offsety = env->get_int_arg(i->first, 3);
				
			} else if (!strcmp(option_name, "fl")) {
#ifdef USE_FFTW
				d2::align::set_frequency_cut(env->get_double_arg(i->first, 1),
						             env->get_double_arg(i->first, 2),
						             env->get_double_arg(i->first, 3));

#else
				ui::get()->error_hint("--fl is not supported", "rebuild ALE with FFTW support");
#endif
			} else if (!strcmp(option_name, "wmx")) {
#ifdef USE_UNIX
				d2::align::set_wmx(env->get_string_arg(i->first, 1),
						   env->get_string_arg(i->first, 2),
						   env->get_string_arg(i->first, 3));
#else
				ui::get()->error_hint("--wmx is not supported", "rebuild ALE with support for --wmx");
#endif
			} else if (!strcmp(option_name, "flshow")) {
				d2::align::set_fl_show(env->get_string_arg(i->first, 1));
			} else if (!strcmp(option_name, "3dpx")) {

				d3px_parameters = (double *) local_realloc(d3px_parameters, (d3px_count + 1) * 6 * sizeof(double));

				for (int param = 0; param < 6; param++)
					d3px_parameters[6 * d3px_count + param] = env->get_double_arg(i->first, param + 1);

				/*
				 * Swap x and y, since their internal meanings differ from their external meanings.
				 */

				for (int param = 0; param < 2; param++) {
					double temp = d3px_parameters[6 * d3px_count + 2 + param];
					d3px_parameters[6 * d3px_count + 2 + param] = d3px_parameters[6 * d3px_count + 0 + param];
					d3px_parameters[6 * d3px_count + 0 + param] = temp;
				}


				/*
				 * Increment counters
				 */

				d3px_count++;
				
			} else if (!strcmp(option_name, "ex") || !strcmp(option_name, "fex")) {

				ex_parameters = (d2::exclusion *) local_realloc(ex_parameters, 
						(ex_count + 1) * sizeof(d2::exclusion));

				ex_parameters[ex_count].type = (!strcmp(option_name, "ex"))
						             ? d2::exclusion::RENDER
							     : d2::exclusion::FRAME;

				/*
				 * Get parameters, swapping x and y coordinates
				 */

				ex_parameters[ex_count].x[0] = env->get_int_arg(i->first, 1 + 2);
				ex_parameters[ex_count].x[1] = env->get_int_arg(i->first, 1 + 3);
				ex_parameters[ex_count].x[2] = env->get_int_arg(i->first, 1 + 0);
				ex_parameters[ex_count].x[3] = env->get_int_arg(i->first, 1 + 1);
				ex_parameters[ex_count].x[4] = env->get_int_arg(i->first, 1 + 4);
				ex_parameters[ex_count].x[5] = env->get_int_arg(i->first, 1 + 5);

				/*
				 * Increment counters
				 */

				ex_count++;

			} else if (!strcmp(option_name, "crop") || !strcmp(option_name, "fcrop")) {

				ex_parameters = (d2::exclusion *) local_realloc(ex_parameters, 
						(ex_count + 4) * sizeof(d2::exclusion));

				for (int r = 0; r < 4; r++) 
					ex_parameters[ex_count + r].type = (!strcmp(option_name, "crop"))
								         ? d2::exclusion::RENDER
								         : d2::exclusion::FRAME;


				int crop_args[6];

				for (int param = 0; param < 6; param++)
					crop_args[param] = env->get_int_arg(i->first, param + 1);

				/*
				 * Construct exclusion regions from the crop area,
				 * swapping x and y, since their internal meanings
				 * differ from their external meanings.
				 */

				/*
				 * Exclusion region 1: low x
				 */

				ex_parameters[ex_count + 0].x[0] = INT_MIN;
				ex_parameters[ex_count + 0].x[1] = crop_args[2] - 1;
				ex_parameters[ex_count + 0].x[2] = INT_MIN;
				ex_parameters[ex_count + 0].x[3] = INT_MAX;
				ex_parameters[ex_count + 0].x[4] = crop_args[4];
				ex_parameters[ex_count + 0].x[5] = crop_args[5];

				/*
				 * Exclusion region 2: low y
				 */

				ex_parameters[ex_count + 1].x[0] = INT_MIN;
				ex_parameters[ex_count + 1].x[1] = INT_MAX;
				ex_parameters[ex_count + 1].x[2] = INT_MIN;
				ex_parameters[ex_count + 1].x[3] = crop_args[0] - 1;
				ex_parameters[ex_count + 1].x[4] = crop_args[4];
				ex_parameters[ex_count + 1].x[5] = crop_args[5];

				/*
				 * Exclusion region 3: high y
				 */

				ex_parameters[ex_count + 2].x[0] = INT_MIN;
				ex_parameters[ex_count + 2].x[1] = INT_MAX;
				ex_parameters[ex_count + 2].x[2] = crop_args[1] + 1;
				ex_parameters[ex_count + 2].x[3] = INT_MAX;
				ex_parameters[ex_count + 2].x[4] = crop_args[4];
				ex_parameters[ex_count + 2].x[5] = crop_args[5];

				/*
				 * Exclusion region 4: high x
				 */

				ex_parameters[ex_count + 3].x[0] = crop_args[3] + 1;
				ex_parameters[ex_count + 3].x[1] = INT_MAX;
				ex_parameters[ex_count + 3].x[2] = INT_MIN;
				ex_parameters[ex_count + 3].x[3] = INT_MAX;
				ex_parameters[ex_count + 3].x[4] = crop_args[4];
				ex_parameters[ex_count + 3].x[5] = crop_args[5];

				/*
				 * Increment counters
				 */

				ex_count += 4;
				
			} else if (!strcmp(option_name, "exshow")) {
				ex_show = 1;
			} else if (!strcmp(option_name, "wt")) {
				d2::render::set_wt(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "3d-chain")) {
				d3chain_type = env->get_string_arg(i->first, 1);
			} else if (!strcmp(option_name, "dchain")) {
				ochain_types[0] = env->get_string_arg(i->first, 1);
			} else if (!strcmp(option_name, "achain")) {
				achain_type = env->get_string_arg(i->first, 1);
			} else if (!strcmp(option_name, "afilter")) {
				afilter_type = env->get_string_arg(i->first, 1);
			} else if (!strcmp(option_name, "ochain")) {

				ochain = (d2::render **) local_realloc(ochain, 
							(oc_count + 1) * sizeof(d2::render *));
				ochain_names = (const char **) local_realloc((void *)ochain_names, 
							(oc_count + 1) * sizeof(const char *));
				ochain_types = (const char **) local_realloc((void *)ochain_types, 
							(oc_count + 1) * sizeof(const char *));

				ochain_types[oc_count] = env->get_string_arg(i->first, 1);
				ochain_names[oc_count] = env->get_string_arg(i->first, 2);

				oc_count++;

			} else if (!strcmp(option_name, "visp")) {

				visp = (const char **) local_realloc((void *)visp, 4 *
							(vise_count + 1) * sizeof(const char *));

				for (int param = 0; param < 4; param++)
					visp[vise_count * 4 + param] = env->get_string_arg(i->first, param + 1);

				vise_count++;

			} else if (!strcmp(option_name, "cx")) {
				cx_parameter = env->get_int_arg(i->first, 0) ? env->get_double_arg(i->first, 1) : 0;
			} else if (!strcmp(option_name, "ip")) {
				unsupported::discontinued("--ip <r> <i>", "--lpsf box=<r> --ips <i>");
			} else if (!strcmp(option_name, "cache")) {
				cache = env->get_double_arg(i->first, 1);
			} else if (!strcmp(option_name, "resident")) {
				double resident = env->get_double_arg(i->first, 1);

				d2::image::set_resident(resident);

			} else if (!strcmp(option_name, "bayer")) {

				/*
				 * External order is clockwise from top-left.  Internal
				 * order is counter-clockwise from top-left.
				 */

				const char *option = env->get_string_arg(i->first, 1);

				if (!strcmp(option, "rgbg")) {
					user_bayer = IMAGE_BAYER_RGBG;
				} else if (!strcmp(option, "bgrg")) {
					user_bayer = IMAGE_BAYER_BGRG;
				} else if (!strcmp(option, "gbgr")) {
					user_bayer = IMAGE_BAYER_GRGB;
				} else if (!strcmp(option, "grgb")) {
					user_bayer = IMAGE_BAYER_GBGR;
				} else if (!strcmp(option, "none")) {
					user_bayer = IMAGE_BAYER_NONE;
				} else {
					bad_arg("--bayer");
				}
				
			} else if (!strcmp(option_name, "lpsf")) {
				psf[psf_linear] = env->get_string_arg(i->first, 1);
			} else if (!strcmp(option_name, "nlpsf")) {
				psf[psf_nonlinear] = env->get_string_arg(i->first, 1);
			} else if (!strcmp(option_name, "psf-match")) {

				psf_match = 1;

				for (int index = 0; index < 6; index++) {
					psf_match_args[index] = env->get_double_arg(i->first, index + 1);
				}

			} else if (!strcmp(option_name, "device")) {
				device = env->get_string_arg(i->first, 1);
#if 0
			} else if (!strcmp(option_name, "usm")) {

				if (d3_output != NULL)
					unsupported::fornow("3D modeling with unsharp mask");

				usm_multiplier = env->get_double_arg(i->first, 1);
#endif

			} else if (!strcmp(option_name, "ipr")) {

				ip_iterations = env->get_int_arg(i->first, 1);

				ui::get()->warn("--ipr is deprecated.  Use --ips instead");

			} else if (!strcmp(option_name, "cpp-err")) {
				if (!strcmp(env->get_string_arg(i->first, 0), "median"))
					d3::cpf::err_median();
				else if (!strcmp(env->get_string_arg(i->first, 0), "mean"))
					d3::cpf::err_mean();
			} else if (!strcmp(option_name, "vp-adjust")) {
				if (env->get_int_arg(i->first, 0))
					d3::align::vp_adjust();
				else
					d3::align::vp_noadjust();
			} else if (!strcmp(option_name, "vo-adjust")) {
				if (env->get_int_arg(i->first, 0))
					d3::align::vo_adjust();
				else
					d3::align::vo_noadjust();
			} else if (!strcmp(option_name, "ip-statistic")) {
				if (!strcmp(env->get_string_arg(i->first, 0), "mean"))
					ip_use_median = 0;
				else if (!strcmp(env->get_string_arg(i->first, 0), "median"))
					ip_use_median = 1;
			} else if (!strcmp(option_name, "ips")) {
				ip_iterations = env->get_int_arg(i->first, 1);
			} else if (!strcmp(option_name, "ip-wl")) {
				int limited = env->get_int_arg(i->first, 0);
				if (limited) {
					ipwl = env->get_double_arg(i->first, 1);
				} else {
					ipwl = 0;
				}
			} else if (!strcmp(option_name, "ipc")) {
				unsupported::discontinued("--ipc <c> <i>", "--ips <i> --lpsf <c>", "--ips <i> --device <c>");
			} else if (!strcmp(option_name, "exp-extend")) {
				if (env->get_int_arg(i->first, 0))
					d2::image_rw::exp_scale();
				else
					d2::image_rw::exp_noscale();
			} else if (!strcmp(option_name, "exp-register")) {
				if (env->get_int_arg(i->first, 0) == 1) {
					exposure_register = 1;
					d2::align::exp_register();
				} else if (env->get_int_arg(i->first, 0) == 0) {
					exposure_register = 0;
					d2::align::exp_noregister();
				} else if (env->get_int_arg(i->first, 0) == 2) {
					exposure_register = 2;
					d2::align::exp_meta_only();
				}
			} else if (!strcmp(option_name, "drizzle-only")) {
				unsupported::discontinued("--drizzle-only", "--dchain box:1");
			} else if (!strcmp(option_name, "subspace-traverse")) {
				unsupported::undocumented("--subspace-traverse");
				d3::scene::set_subspace_traverse();
			} else if (!strcmp(option_name, "3d-filter")) {
				if (env->get_int_arg(i->first, 0))
					d3::scene::filter();
				else
					d3::scene::nofilter();
			} else if (!strcmp(option_name, "occ-norm")) {
				if (env->get_int_arg(i->first, 0))
					d3::scene::nw();
				else
					d3::scene::no_nw();
			} else if (!strcmp(option_name, "inc")) {
				inc = env->get_int_arg(i->first, 0);
			} else if (!strcmp(option_name, "exp-mult")) {
				double exp_c, exp_r, exp_b;

				exp_c = env->get_double_arg(i->first, 1);
				exp_r = env->get_double_arg(i->first, 2);
				exp_b = env->get_double_arg(i->first, 3);

				exp_mult = d2::pixel(1/(exp_r * exp_c), 1/exp_c, 1/(exp_b * exp_c));

			} else if (!strcmp(option_name, "visp-scale")) {

				vise_scale_factor = env->get_double_arg(i->first, 1);

				if (vise_scale_factor <= 0.0)
					ui::get()->error("VISP scale must be greater than zero");

				if (!finite(vise_scale_factor))
					ui::get()->error("VISP scale must be finite");

			} else if (!strcmp(option_name, "scale")) {

				scale_factor = env->get_double_arg(i->first, 1);

				if (scale_factor <= 0)
					ui::get()->error("Scale factor must be greater than zero");

				if (!finite(scale_factor))
					ui::get()->error("Scale factor must be finite");

			} else if (!strcmp(option_name, "metric")) {
				d2::align::set_metric_exponent(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "threshold")) {
				d2::align::set_match_threshold(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "drizzle-diam")) {
				unsupported::discontinued("--drizzle-diam=<x>", "--dchain box:1");
			} else if (!strcmp(option_name, "perturb-lower")) {
				const char *option = env->get_string_arg(i->first, 1);
				double perturb_lower;
				int characters;
				sscanf(option, "%lf%n", &perturb_lower, &characters);
				if (perturb_lower <= 0)
					ui::get()->error("--perturb-lower= value is non-positive");

				if (*(option + characters) == '%')
					d2::align::set_perturb_lower(perturb_lower, 1);
				else
					d2::align::set_perturb_lower(perturb_lower, 0);
			} else if (!strcmp(option_name, "stepsize")) {
				ui::get()->warn("--stepsize is deprecated.  Use --perturb-lower instead");
				d2::align::set_perturb_lower(env->get_double_arg(i->first, 1), 0);
			} else if (!strcmp(option_name, "va-upper")) {
				const char *option = env->get_string_arg(i->first, 1);
				double va_upper;
				int characters;
				sscanf(option, "%lf%n", &va_upper, &characters);
				if (*(option + characters) == '%')
					ui::get()->error("--va-upper= does not accept '%' arguments\n");
				else
					d3::cpf::set_va_upper(va_upper);
			} else if (!strcmp(option_name, "cpp-upper")) {
				const char *option = env->get_string_arg(i->first, 1);
				double perturb_upper;
				int characters;
				sscanf(option, "%lf%n", &perturb_upper, &characters);
				if (*(option + characters) == '%')
					ui::get()->error("--cpp-upper= does not currently accept '%' arguments\n");
				else
					d3::cpf::set_cpp_upper(perturb_upper);
			} else if (!strcmp(option_name, "cpp-lower")) {
				const char *option = env->get_string_arg(i->first, 1);
				double perturb_lower;
				int characters;
				sscanf(option, "%lf%n", &perturb_lower, &characters);
				if (*(option + characters) == '%')
					ui::get()->error("--cpp-lower= does not currently accept '%' arguments\n");
				else
					d3::cpf::set_cpp_lower(perturb_lower);
			} else if (!strcmp(option_name, "hf-enhance")) {
				unsupported::discontinued("--hf-enhance=<x>");
			} else if (!strcmp(option_name, "multi")) {
				d2::trans_multi::set_multi(env->get_string_arg(i->first, 1));
			} else if (!strcmp(option_name, "track")) {
				if (!strcmp(env->get_string_arg(i->first, 0), "none")) {
					d2::trans_multi::track_none();
				} else if (!strcmp(env->get_string_arg(i->first, 0), "primary")) {
					d2::trans_multi::track_primary();
				} else if (!strcmp(env->get_string_arg(i->first, 0), "point")) {
					d2::trans_multi::track_point(env->get_double_arg(i->first, 1),
					                             env->get_double_arg(i->first, 2));
				} else {
					assert(0);
				}
			} else if (!strcmp(option_name, "gs")) {
				d2::align::gs(env->get_string_arg(i->first, 1));
			} else if (!strcmp(option_name, "rot-upper")) {
				d2::align::set_rot_max((int) floor(env->get_double_arg(i->first, 1)));
			} else if (!strcmp(option_name, "bda-mult")) {
				d2::align::set_bda_mult(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "bda-rate")) {
				d2::align::set_bda_rate(env->get_double_arg(i->first, 1));
			} else if (!strcmp(option_name, "lod-preferred")) {
				d2::align::set_lod_preferred((int) floor(env->get_double_arg(i->first, 1)));
			} else if (!strcmp(option_name, "min-dimension")) {
				d2::align::set_min_dimension((int) ceil(env->get_double_arg(i->first, 1)));
			} else if (!strcmp(option_name, "cpf-load")) {
				d3::cpf::init_loadfile(env->get_string_arg(i->first, 1));
#if 0
			} else if (!strcmp(option_name, "model-load")) {
				d3::scene::load_model(env->get_string_arg(i->first, 1));
			} else if (!strcmp(option_name, "model-save")) {
				d3::scene::save_model(env->get_string_arg(i->first, 1));
#endif
			} else if (!strcmp(option_name, "trans-load")) {
				d2::tload_delete(tload);
				tload = d2::tload_new(env->get_string_arg(i->first, 1));
				d2::align::set_tload(tload);
			} else if (!strcmp(option_name, "trans-save")) {
				tsave_delete(tsave);
				tsave = d2::tsave_new(env->get_string_arg(i->first, 1));
				d2::align::set_tsave(tsave);
			} else if (!strcmp(option_name, "3d-trans-load")) {
				d3::tload_delete(d3_tload);
				d3_tload = d3::tload_new(env->get_string_arg(i->first, 1));
				d3::align::set_tload(d3_tload);
			} else if (!strcmp(option_name, "3d-trans-save")) {
				d3::tsave_delete(d3_tsave);
				d3_tsave = d3::tsave_new(env->get_string_arg(i->first, 1));
				d3::align::set_tsave(d3_tsave);
			} else {
				assert(0);
			}
		}

		/*
		 * Initialize the interface.
		 */

		ui::get();

		/*
		 * Apply implication logic.
		 */

		if (extend == 0 && vise_count != 0) {
			implication::changed("VISP requires increased image extents.",
					     "Image extension is now enabled.",
					     "--extend");
			extend = 1;
		}

		if (psf_match && ex_count)
			unsupported::fornow("PSF calibration with exclusion regions.");

		
		if (d3_output != NULL && ip_iterations != 0) 
			unsupported::fornow("3D modeling with Irani-Peleg rendering");

#if 0
		if (extend == 0 && d3_output != NULL) {
			implication::changed("3D modeling requires increased image extents.",
					     "Image extension is now enabled.",
					     "--extend");
			extend = 1;
		}
#endif

#if 0
		if (cx_parameter != 0 && !exposure_register) {
			implication::changed("Certainty-based rendering requires exposure registration.",
					     "Exposure registration is now enabled.",
					     "--exp-register");
			d2::align::exp_register();
			exposure_register = 1;
		}
#endif

		/*
		 * Set alignment class exclusion region static variables
		 */

		d2::align::set_exclusion(ex_parameters, ex_count);

		/*
		 * Initialize renderer class statics.
		 */

		d2::render::render_init(ex_count, ex_parameters, ex_show, extend, scale_factor);

		/*
		 * Set confidence
		 */

		d2::exposure::set_confidence(cx_parameter);

		/*
		 * Keep transformations for Irani-Peleg, psf-match, and
		 * VISE
		 */

		if (ip_iterations > 0 || psf_match || vise_count > 0) {
			d2::align::keep();
		}

		/*
		 * Initialize device-specific variables
		 */

		// int input_file_count = argc - i - 1;
		int input_file_count = files.size() - 1;

		d2::psf *device_response[psf_N] = { NULL, NULL };
		d2::exposure **input_exposure = NULL;
		ale_pos view_angle = 43.7 * M_PI / 180;  
		// ale_pos view_angle = 90 * M_PI / 180;  
		input_exposure = (d2::exposure **)
			// malloc((argc - i - 1) * sizeof(d2::exposure *));
			malloc(input_file_count * sizeof(d2::exposure *));

		if (device != NULL) {
			if (!strcmp(device, "xvp610_640x480")) {
				device_response[psf_linear] = new xvp610_640x480::lpsf();
				device_response[psf_nonlinear] = new xvp610_640x480::nlpsf();
				for (int ii = 0; ii < input_file_count; ii++)
					input_exposure[ii] = new xvp610_640x480::exposure();
				view_angle = xvp610_640x480::view_angle();
			} else if (!strcmp(device, "xvp610_320x240")) {
				device_response[psf_linear] = new xvp610_320x240::lpsf();
				device_response[psf_nonlinear] = new xvp610_320x240::nlpsf();
				for (int ii = 0; ii < input_file_count; ii++)
					input_exposure[ii] = new xvp610_320x240::exposure();
				view_angle = xvp610_320x240::view_angle();
			} else if (!strcmp(device, "ov7620")) {
				device_response[psf_linear] = new ov7620_raw_linear::lpsf();
				device_response[psf_nonlinear] = NULL;
				for (int ii = 0; ii < input_file_count; ii++)
					input_exposure[ii] = new ov7620_raw_linear::exposure();
				d2::image_rw::set_default_bayer(IMAGE_BAYER_BGRG);
			} else if (!strcmp(device, "canon_300d")) {
				device_response[psf_linear] = new canon_300d_raw_linear::lpsf();
				device_response[psf_nonlinear] = NULL;
				for (int ii = 0; ii < input_file_count; ii++)
					input_exposure[ii] = new canon_300d_raw_linear::exposure();
				d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
			} else if (!strcmp(device, "nikon_d50")) {
				device_response[psf_linear] = nikon_d50::lpsf();
				device_response[psf_nonlinear] = nikon_d50::nlpsf();
				for (int ii = 0; ii < input_file_count; ii++)
					input_exposure[ii] = new nikon_d50::exposure();
				d2::image_rw::set_default_bayer( nikon_d50::bayer() );
			} else if (!strcmp(device, "canon_300d+85mm_1.8")) {
				device_response[psf_linear] = new canon_300d_raw_linear_85mm_1_8::lpsf();
				device_response[psf_nonlinear] = NULL;
				for (int ii = 0; ii < input_file_count; ii++)
					input_exposure[ii] = new canon_300d_raw_linear_85mm_1_8::exposure();
				d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
				view_angle = canon_300d_raw_linear_85mm_1_8::view_angle();
			} else if (!strcmp(device, "canon_300d+50mm_1.8")) {
				device_response[psf_linear] = new canon_300d_raw_linear_50mm_1_8::lpsf();
				device_response[psf_nonlinear] = NULL;
				for (int ii = 0; ii < input_file_count; ii++)
					input_exposure[ii] = new canon_300d_raw_linear_50mm_1_8::exposure();
				d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
				view_angle = canon_300d_raw_linear_50mm_1_8::view_angle();
			} else if (!strcmp(device, "canon_300d+50mm_1.4")) {
				device_response[psf_linear] = new canon_300d_raw_linear_50mm_1_4::lpsf();
				device_response[psf_nonlinear] = NULL;
				for (int ii = 0; ii < input_file_count; ii++)
					input_exposure[ii] = new canon_300d_raw_linear_50mm_1_4::exposure();
				d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
				view_angle = canon_300d_raw_linear_50mm_1_4::view_angle();
			} else if (!strcmp(device, "canon_300d+50mm_1.4@1.4")) {
				device_response[psf_linear] = new canon_300d_raw_linear_50mm_1_4_1_4::lpsf();
				device_response[psf_nonlinear] = NULL;
				for (int ii = 0; ii < input_file_count; ii++)
					input_exposure[ii] = new canon_300d_raw_linear_50mm_1_4_1_4::exposure();
				d2::image_rw::set_default_bayer(IMAGE_BAYER_RGBG);
				view_angle = canon_300d_raw_linear_50mm_1_4_1_4::view_angle();
			} else {
				ui::get()->unknown_device(device);
			}
		} else {
			for (int ii = 0; ii < input_file_count; ii++)
				input_exposure[ii] = new d2::exposure_default();
		}

		/*
		 * User-specified variables.
		 */

		if (user_view_angle != 0) {
			view_angle = user_view_angle;
		}

		if (user_bayer != IMAGE_BAYER_DEFAULT) {
			d2::image_rw::set_default_bayer(user_bayer);
		}

		/*
		 * PSF-match exposure.
		 */
		if (psf_match) {
			delete input_exposure[input_file_count - 1];
			input_exposure[input_file_count - 1] = new d2::exposure_default();
		}

		/*
		 * Initialize output exposure
		 */

		d2::exposure *output_exposure = new d2::exposure_default();
		output_exposure->set_multiplier(exp_mult);

		/*
		 * Configure the response function.
		 */

		d2::psf *response[2] = {NULL, NULL};

		for (int n = 0; n < psf_N; n++ ) {
			if (psf[n] != NULL) {

				response[n] = d2::psf_parse::get((n == psf_linear), psf[n]);

			} else if (device_response[n] != NULL) {

				/*
				 * Device-specific response
				 */

				response[n] = device_response[n];

			} else {

				/*
				 * Default point-spread function.
				 */

				if (n == psf_linear) {

					/*
					 * Default lpsf is a box filter
					 * of diameter 1.0 (radius
					 * 0.5).
					 */

					response[n] = new d2::box(0.5);

				} else if (n == psf_nonlinear) {

					/*
					 * nlpsf is disabled by default.
					 */

					 response[n] = NULL;
				}
			}
		}

		/* 
		 * First file argument.  Print general file information as well
		 * as information specific to this argument.  Initialize image
		 * file handler.
		 */

		// d2::image_rw::init(argc - i - 1, argv + i, argv[argc - 1], input_exposure, output_exposure);
		// ochain_names[0] = argv[argc - 1];
		
		const char **input_files = (const char **) malloc(sizeof(const char *) * input_file_count);
		for (int i = 0; i < input_file_count; i++) {
			input_files[i] = files[i].first;
		}

		d2::image_rw::init(input_file_count, input_files, files[files.size() - 1].first,
				input_exposure, output_exposure);

		seq_struct seq;

		ale_sequence sequence = ale_new_sequence(seq_file, seq_trans, seq_ex, seq_step, seq_ui, (void *) &seq, cache);

		ochain_names[0] = files[files.size() - 1].first;

		/*
		 * Handle control point data for alignment
		 */
		d2::align::set_cp_count(d3::cpf::count());
		for (unsigned int ii = 0; ii < d3::cpf::count(); ii++)
			d2::align::set_cp(ii, d3::cpf::get_2d(ii));

		/*
		 * PSF-match bayer patterns.
		 */

		if (psf_match) {
			// d2::image_rw::set_specific_bayer(argc - i - 2, IMAGE_BAYER_NONE);
			d2::image_rw::set_specific_bayer(input_file_count - 1, IMAGE_BAYER_NONE);
		}

		/*
		 * Handle alignment weight map, if necessary
		 */

		if (wm_filename != NULL) {
			d2::image *weight_map;
			weight_map = d2::image_rw::read_image(wm_filename, new d2::exposure_linear());
			weight_map->set_offset(wm_offsety, wm_offsetx);
			d2::align::set_weight_map(weight_map);
		}

		/*
		 * Initialize alignment interpolant.
		 */

		if (strcmp(afilter_type, "internal"))
			d2::align::set_interpolant(d2::render_parse::get_SSF(afilter_type));

		/*
		 * Initialize achain and ochain.
		 */

		achain = d2::render_parse::get(achain_type);
		
		for (int chain = 0; chain < oc_count; chain++)
			ochain[chain] = d2::render_parse::get(ochain_types[chain]);

		/*
		 * Use merged renderings as reference images in
		 * alignment.
		 */

		d2::align::set_reference(achain);

		/*
		 * Tell the alignment class about the scale factor.
		 */

		d2::align::set_scale(scale_factor);

		/*
		 * Initialize visp.
		 */

		d2::vise_core::set_scale(vise_scale_factor);

		for (int opt = 0; opt < vise_count; opt++) {
			d2::vise_core::add(d2::render_parse::get(visp[opt * 4 + 0]),
					   visp[opt * 4 + 1],
					   visp[opt * 4 + 2],
					   visp[opt * 4 + 3]);
		}

		/*
		 * Initialize non-incremental renderers
		 */

#if 0
		if (usm_multiplier != 0) {

			/*
			 * Unsharp Mask renderer
			 */

			ochain[0] = new d2::usm(ochain[0], scale_factor,
					usm_multiplier, inc, response[psf_linear],
					response[psf_nonlinear], &input_exposure[0]);
		}
#endif

		if (psf_match) {
			
			/*
			 * Point-spread function calibration renderer.
			 * This renderer does not produce image output.
			 * It is reserved for use with the point-spread
			 * function calibration script
			 * ale-psf-calibrate.
			 */

			ochain[0] = new d2::psf_calibrate(ochain[0],
					1, inc, response[psf_linear],
					response[psf_nonlinear],
					psf_match_args);

		} else if (ip_iterations != 0) {

			/*
			 * Irani-Peleg renderer
			 */

			ochain[0] = new d2::ipc( ochain[0], ip_iterations,
					inc, response[psf_linear],
					response[psf_nonlinear],
					(exposure_register == 1), ip_use_median, ipwl);
		}

		/*
		 * Iterate through all files.
		 */

		ale_sequence_run(sequence);

		for (unsigned int j = 0; j < d2::image_rw::count(); j++) {

			/*
			 * Iterate through non-global options
			 */

			environment *env = files[j].second;

			for (std::map<const char *, const char *>::iterator i = env->get_map().begin();
					i != env->get_map().end(); i++) {

				if (!env->is_option(i->first))
					continue;

				const char *option_name = env->get_option_name(i->first);

				if (!strcmp(option_name, "mc")) {
					d2::align::mc(env->get_double_arg(i->first, 1));
				} else if (!strcmp(option_name, "md")) {
					d2::trans_multi::set_md(env->get_double_arg(i->first, 1));
				} else if (!strcmp(option_name, "ma-cert")) {
					d2::align::set_ma_cert(env->get_double_arg(i->first, 1));
				} else if (!strcmp(option_name, "mi")) {
					d2::trans_multi::set_mi(env->get_double_arg(i->first, 1));
				} else if (!strcmp(option_name, "gs-mo")) {
					const char *option = env->get_string_arg(i->first, 1);
					double gs_mo;
					int characters;
					sscanf(option, "%lf%n", &gs_mo, &characters);
					if (*(option + characters) == '%')
						d2::align::gs_mo(gs_mo, 1);
					else
						d2::align::gs_mo(gs_mo, 0);
				} else if (!strcmp(option_name, "ev")) {
					double ev = env->get_double_arg(i->first, 1);
					double gain_value = pow(2, -ev);

					if (j == 0)
						d2::exposure::set_gain_reference(gain_value);
					else 
						input_exposure[j]->set_gain_multiplier(
								(double) d2::exposure::get_gain_reference()
							      / gain_value);

				} else if (!strcmp(option_name, "black")) {
					double black = env->get_double_arg(i->first, 1);
					input_exposure[j]->set_black_level(black);
				} else if (!strcmp(option_name, "perturb-upper")) {
					const char *option = env->get_string_arg(i->first, 1);
					double perturb_upper;
					int characters;
					sscanf(option, "%lf%n", &perturb_upper, &characters);
					if (*(option + characters) == '%')
						d2::align::set_perturb_upper(perturb_upper, 1);
					else
						d2::align::set_perturb_upper(perturb_upper, 0);
				} else if (!strcmp(option_name, "threads")) {
					thread::set_count((unsigned int) env->get_int_arg(i->first, 1));
				} else if (!strcmp(option_name, "per-cpu")) {
					thread::set_per_cpu((unsigned int) env->get_int_arg(i->first, 1));
				} else {
					/*
					 * This error should be encountered earlier.
					 */

					assert(0);

					fprintf(stderr, "\n\nError: option `%s' must be applied globally.", option_name);
					fprintf(stderr, "\n\nHint:  Move option `%s' prior to file and scope operators.\n\n", 
							option_name);
					exit(1);
				}
			}



			if (j == 0) {

				/*
				 * Write comment information about original frame and
				 * target image to the transformation save file, if we
				 * have one.
				 */

				tsave_orig(tsave, files[0].first);
				tsave_target(tsave, files[files.size() - 1].first);

				/*
				 * Handle the original frame.
				 */

				// ui::get()->original_frame_start(argv[i]);
				ui::get()->original_frame_start(files[0].first);

				for (int opt = 0; opt < oc_count; opt++) {
					ui::get()->set_orender_current(opt);
					ochain[opt]->sync(0);
					if  (inc) {
						ui::get()->writing_output(opt);
						d2::image_rw::write_image(ochain_names[opt], 
							ochain[opt]->get_image(0));
					}
				}

				d2::vise_core::frame_queue_add(0);

				ui::get()->original_frame_done();

				continue;
			}

			/*
			 * Handle supplemental frames.
			 */

			const char *name = d2::image_rw::name(j);

			ui::get()->supplemental_frame_start(name);

			/*
			 * Write comment information about the
			 * supplemental frame to the transformation
			 * save file, if we have one.
			 */

			tsave_info (tsave, name);

			for (int opt = 0; opt < oc_count; opt++) {
				ui::get()->set_orender_current(opt);
				ochain[opt]->sync(j);
				if (inc) {
					ui::get()->writing_output(opt);
					d2::image_rw::write_image(ochain_names[opt], 
						ochain[opt]->get_image(j));
				}
			}

			d2::vise_core::frame_queue_add(j);

			ui::get()->supplemental_frame_done();
		}

		/*
		 * Do any post-processing and output final image
		 *
		 * XXX: note that the Irani-Peleg renderer currently
		 * returns zero for ochain[0]->sync(), since it writes
		 * output internally when inc != 0.
		 *
		 * XXX: Leave ochain[0] last, since this may allow disposal of
		 * rendering structures not required by the Irani-Peleg
		 * renderer.
		 */

		for (int opt = 1; opt < oc_count; opt++) 
		if  ((ochain[opt]->sync() || !inc) && !psf_match) {
			ui::get()->writing_output(opt);
			d2::image_rw::write_image(ochain_names[opt], ochain[opt]->get_image());
		}

		if (oc_count > 0)
		if ((ochain[0]->sync() || !inc) && !psf_match) {
			ui::get()->writing_output(0);
			d2::image_rw::write_image(ochain_names[0], ochain[0]->get_image());
		}

		/*
		 * Output a summary match statistic.
		 */

		ui::get()->ale_2d_done((double) d2::align::match_summary());

		/*
		 * Perform any 3D tasks
		 */

		optimizations::begin_3d_work();

		if (d3_count > 0) {

			ui::get()->d3_start();

			d3::align::init_angle(view_angle);

			ui::get()->d3_init_view_angle((double) view_angle / M_PI * 180);

			d3::align::init_from_d2();

			if (d3::cpf::count() > 0) {
				ui::get()->d3_control_point_solve();
				d3::cpf::solve_3d();
				ui::get()->d3_control_point_solve_done();
			}

			ui::get()->d3_final_view_angle(d3::align::angle_of(0) / M_PI * 180);

			d3::align::write_alignments();

			d3::scene::set_filter_type(d3chain_type);

			d3::scene::init_from_d2();

			ui::get()->d3_subdividing_space();
			d3::scene::make_space(d3_depth, d3_output, &d3_depth_pt, &d3_output_pt);
			ui::get()->d3_subdividing_space_done();

			ui::get()->d3_updating_occupancy();
			d3::scene::reduce_cost_to_search_depth(output_exposure, inc);
			ui::get()->d3_updating_occupancy_done();

			d3::scene::d3px(d3px_count, d3px_parameters);
			int view_count = 0;
			for (unsigned int i = 0; i < d2::image_rw::count(); i++) {
				assert (i < d3_count);

				if (d3_depth[i] != NULL) {
					ui::get()->d3_writing_output(d3_depth[i]);
					ui::get()->d3_render_status(0, 0, -1, -1, -1, -1, 0);
					const d2::image *im = d3::scene::depth(i);
					d2::image_rw::write_image(d3_depth[i], im, output_exposure, 1, 1);
					delete im;
					ui::get()->d3_writing_output_done();
				}

				if (d3_output[i] != NULL) {
					ui::get()->d3_writing_output(d3_output[i]);
					const d2::image *im = d3::scene::view(i);
					d2::image_rw::write_image(d3_output[i], im, output_exposure);
					delete im;
					d3::focus::set_camera(view_count++);
					ui::get()->d3_writing_output_done();
				}

				for (std::map<const char *, d3::pt>::iterator i = d3_output_pt.begin();
						i != d3_output_pt.end(); i++) {

					ui::get()->d3_writing_output(i->first);
					const d2::image *im = d3::scene::view(i->second);
					d2::image_rw::write_image(i->first, im, output_exposure);
					delete im;
					d3::focus::set_camera(view_count++);
					ui::get()->d3_writing_output_done();
				}

				for (std::map<const char *, d3::pt>::iterator i = d3_depth_pt.begin();
						i != d3_depth_pt.end(); i++) {

					ui::get()->d3_writing_output(i->first);
					ui::get()->d3_render_status(0, 0, -1, -1, -1, -1, 0);
					const d2::image *im = d3::scene::depth(i->second);
					d2::image_rw::write_image(i->first, im, output_exposure, 1, 1);
					delete im;
					ui::get()->d3_writing_output_done();
				}
			}

			for (unsigned int i = d2::image_rw::count(); i < d3_count; i++) {
				if (d3_depth[i] != NULL) {
					fprintf(stderr, "\n\n*** Frame number for --3dd too high. ***\n\n");
				}
				if (d3_output[i] != NULL) {
					fprintf(stderr, "\n\n*** Frame number for --3dv too high. ***\n\n");
				}
			}
		}

		/*
		 * Destroy the image file handler
		 */

		d2::image_rw::destroy();

		/*
		 * Delete the transformation file structures, if any
		 * exist.
		 */

		tsave_delete(tsave);
		tload_delete(tload);

		/*
		 * We're done.
		 */

		exit(0);
	}
};

#endif
