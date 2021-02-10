PYPATH=${PWD}/..
MDFILES:=$(shell find . -type f -name '*.md')
INDEX:=_build/html/index.html

all: $(INDEX)
	PYTHONPATH=$(PYPATH) bash build.sh
clean:
	rm -rf _build

$(INDEX): $(MDFILE) _config.yml _toc.yml build.sh references.bib
