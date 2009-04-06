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
	VARS=`cat VERSION | grep "^:.*:" | cut -d: -s -f 2`
	for i in $VARS; do
		j=`$0 $i`;
		if [ -n "$XSLT_PROCESSOR" ]; then
			echo -n "--stringparam $i \"$j\" "
		elif [ -n "$ASCIIDOC_PROCESSOR" ]; then

			# XXX: Using include:: is better, but this option is
			# here for completeness.  Unfortunately, I can't get
			# asciidoc -a values having spaces to work with either
			# double or single quotes; do this instead, as it
			# allows values without spaces to continue to work.
			#
			# - dhilvert

			j=`echo -n $j | sed -e 's/ /\\\\ /g'`
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
