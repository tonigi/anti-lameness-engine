TARGETS=$(shell find -name index.txt | sed -e 's/.txt/.html/')

all: $(TARGETS)

VERSION.txt:
	ln -s ../../VERSION ./VERSION.txt

%.html: %.txt VERSION.txt footer-html.txt
	# XXX: asciidoc considers out-of-directory includes unsafe.
	asciidoc --unsafe -f footer-html.txt -b xhtml11 $<

clean:
	rm $(TARGETS)
