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

#ifndef __ui_tty_h__
#define __ui_tty_h__

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "ui.h"
#include "util.h"
#include "../d2.h"

/*
 * TTY user interface
 */

class ui_tty : public ui {
private:
	int terminal_width; 
	char *buffer;
	int buffer_index;
	int status_index;

	void clear_buffer() {
		buffer[0] = '\0';
		buffer_index = 0;
	}

	void write_buffer() {
		if (buffer_index < 0)
			return;
		fputc('\r', ui_stream);
		fputs(buffer, ui_stream);
		status_index = buffer_index;
	}

	void status_clear() {
		while (status_index < terminal_width) {
			status_index++;
			fprintf(ui_stream, " ");
		}
	}


	void status_printf(int count, ...) {

		if (buffer_index < 0)
			return;

		int n;

		char *local_buffer = (char *) calloc(terminal_width - status_index, sizeof(char));

		assert(local_buffer);

		if (!local_buffer)
			return;

		for (int i = 0; i < count; i++) {
			char *format = NULL;
			va_list ap;

			va_start(ap, count);
			
			for (int arg = 0; arg < i + 1; arg++)
				format = va_arg(ap, char *);
			for (int arg = i + 1; arg < count; arg++)
				va_arg(ap, char *);

			assert (format);

			n = vsnprintf(local_buffer, terminal_width - status_index, format, ap);
			va_end(ap);

			if (n < 0 || n > terminal_width - status_index - 1)
				continue;

			fputs(local_buffer, ui_stream);

			status_index += n;

			free(local_buffer);

			return;
		}

		/*
		 * If we reach this point, then there was no valid string produced.
		 */

		status_clear();

		free(local_buffer);
	}

	void pad_align_status() {
		for (int i = 0; i < status.perturb_steps - status.perturb_steps_completed; i++) {
			status_printf(1, " ");
		}
	}
	void pad_match_status() {
		status_printf(1, "               ");
	}

	void write_match_status() {
		status_printf(1, format_string_working(), status.match_value);
	}

	void write_status() {

		if (status.code == status.UNDEFINED
		 || status.code == status.FRAME_DONE
		 || status.code == status.SET_DONE
		 || status.code == status.IP_STEP_DONE) {
			status_clear();
			return;
		}

		if (status.code == status.ALIGN
		 || status.code == status.POSTMATCH
		 || status.code == status.EXPOSURE_PASS_2) {
			pad_align_status();
			write_match_status();
		}

		status_printf(1, " | ");

		switch (status.code) {
			case status.LOAD_FILE:
				status_printf(3, "Loading Image", "Loading", "load");
				break;
			case status.EXPOSURE_PASS_1:
				status_printf(3, "Registering exposure (first pass)", "Registering exposure", "regexp1");
				break;
			case status.LODCLUSTER_CREATE:
				status_printf(3, "Creating LOD cluster, scale %f", "Creating LOD clusters", "lodcluster", 
						pow(2, -status.align_lod));
				break;
			case status.PREMATCH:
				status_printf(3, "Calculating pre-alignment match", "Calculating match", "prematch");
				break;
			case status.ALIGN:
				status_printf(4, "Aligning [perturb=%f] [lod=%f] [exp_mult=%f %f %f]", 
						 "Aligning [perturb=%f] [lod=%f]", "Aligning...", "align",
						 status.perturb_size,
						 pow(2, -status.align_lod),
						 status.exp_multiplier[0],
						 status.exp_multiplier[1],
						 status.exp_multiplier[2]);
				break;
			case status.POSTMATCH:
				status_printf(3, "Calculating post-alignment match", "Calculating match", "postmatch");
				break;
			case status.EXPOSURE_PASS_2:
				status_printf(3, "Registering exposure (second pass)", "Registering exposure", "regexp2");
				break;
			case status.RENDERA:
				status_printf(3, "Rendering alignment reference image", "Rendering", "render-a");
				break;
			case status.RENDERD:
				status_printf(3, "Rendering default chain", "Rendering", "render-d");
				break;
			case status.RENDERO:
				status_printf(3, "Rendering chain %d", "Rendering", "render-o%d", status.onum);
				break;
			case status.WRITED:
				status_printf(4, "Writing default chain to '%s'", 
						 "Writing '%s'", "Writing", "write-d", d2::image_rw::output_name());
				break;
			case status.WRITEO:
				status_printf(3, "Writing image for chain %d", "Writing", "write-o%d", status.onum);
				break;
			case status.IP_RENDER:
				status_printf(2, "Processing frame '%s'", "Processing", d2::image_rw::name(status.ip_frame_num));
				break;
			case status.IP_UPDATE:
				status_printf(3, "Updating approximation", "Updating", "update");
				break;
			case status.IP_WRITE:
				status_printf(3, "Writing '%s'", "Writing", "write", d2::image_rw::output_name());
				break;
			default:
				break;
		}

		status_clear();
	}

	void write_all() {
		write_buffer();
		write_status();
	}


	void printf(char *format, ...) {
		va_list ap;
		int n = -1;

		if (buffer_index >= 0 && buffer_index < terminal_width && format[strlen(format) - 1] != '\n') {
			va_start(ap, format);
			n = vsnprintf(buffer + buffer_index, terminal_width - buffer_index, format, ap);
			va_end(ap);
		}

		if (n >= 0 && n < terminal_width - buffer_index) {
			/*
			 * The message fits in the buffer, so update the index
			 * and write buffer and status information.
			 */
			buffer_index += n;
			write_all();
		} else {
			/*
			 * The message does not fit in the buffer, so write any
			 * existing buffer and append the new text to the stream.
			 */
			if (buffer_index >= 0) {
				assert(buffer_index < terminal_width);
				buffer[buffer_index] = '\0';
				write_buffer();
			}
			buffer_index = -1;
			va_start(ap, format);
			vfprintf(ui_stream, format, ap);
			va_end(ap);
		}

		/*
		 * This is not the only case that produces a newline,
		 * but ignoring other cases should be safe.
		 */
		if (format[strlen(format) - 1] == '\n') {
			buffer_index = 0;
			buffer[0] = '\0';
		}
	}

	void update() {
		if (status.code == status.FRAME_DONE
		 || status.code == status.SET_DONE) {
			printf(".");
			fputc('\n', ui_stream);
			buffer_index = 0;
			buffer[0] = '\0';
		} else {
			write_all();
		}
	}

public:
	/*
	 * Constructor may throw an exception to signal that using ui_wo would
	 * be more appropriate.
	 */
	ui_tty() {
		int exception_value = 1;

		if (!isatty(fileno(ui_stream)))
			throw exception_value;

		terminal_width = get_terminal_width(ui_stream);

		if (terminal_width < 0)
			throw exception_value;

		buffer = (char *) calloc(terminal_width + 1, sizeof(char));

		assert (buffer);

		if (!buffer)
			throw exception_value;

		buffer[0] = '\0';
		buffer_index = 0;
	}

	~ui_tty() {
		free(buffer);
	}
};

#endif
