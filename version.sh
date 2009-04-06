#!/bin/sh

# ALE package version processor
#
# Copyright (C) 2009 David Hilvert
#
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

if [ "$1" = "-n" ]; then
	ECHO_OPTS="-n"
	shift
elif [ "$1" = "-x" ]; then
	XSLT_PROCESSOR="1"
	shift
elif [ "$1" = "-a" ]; then
	ASCIIDOC_PROCESSOR="1"
	shift
fi

#
# If no tag is specified for display, display all
#

if [ -z "$1" ]; then

	# XXX: surely there's a better way to do this than sed and un-sed.
	# Setting IFS seems to cause various disaster.

	VARS=`cat VERSION | grep "^:.*:" | sed -e 's/ /%20/' | cut -d: -s -f 2`
	for i in $VARS; do
		i=`echo -n $i | sed -e 's/%20/ /'`
		j=`$0 "$i"`
		if [ -n "$XSLT_PROCESSOR" ]; then
			
			# XXX: xsltproc reports a problem with too many
			# parameters (!!!).  Get around this by excluding
			# certain parameters.  We choose man page parameters,
			# since we probably don't need these.

			if ! echo $i | grep -q '^man '; then
				echo -n "--stringparam $i \"$j\" "
			fi

		elif [ -n "$ASCIIDOC_PROCESSOR" ]; then

			# XXX: Using include:: is better, but this option is
			# here for completeness.  Unfortunately, I can't get
			# asciidoc -a values having spaces to work with either
			# double or single quotes; do this instead, as it
			# allows values without spaces to continue to work.
			#
			# - dhilvert

			j=`echo -n $j | sed -e 's/ /\\\\ /g'`
			i=`echo -n $i | sed -e 's/ /\\\\ /g'`
			echo -n "-a $i=$j "
		else
			echo ":$i:	$j"
		fi
	done
	exit
fi

#
# Display a particular tag.
#

RESULT=`grep "^:$1:" VERSION | sed -e 's/:[^:]*:\s*//' | sed -e 's/\n//'`

while TAG=`echo -n $RESULT | grep "{.*}"` && TAG=`echo -n $TAG | sed -e "s/.*{\(.*\)}.*/\1/"`; do
	TAGVALUE=`$0 -n $TAG`
	RESULT=`echo $RESULT | sed -e "s*{$TAG}*$TAGVALUE*g"`
done

echo $ECHO_OPTS $RESULT
