# Homebrew's Boost (/opt/homebrew on Apple Silicon) is not on the
# compiler's default include path, so add it when brew is installed.
BREW_PREFIX := $(shell brew --prefix 2>/dev/null)
ifneq ($(BREW_PREFIX),)
export CPLUS_INCLUDE_PATH := $(if $(CPLUS_INCLUDE_PATH),$(CPLUS_INCLUDE_PATH):)$(BREW_PREFIX)/include
endif

JARS:
	make -C Java jarfiles

CBINS: 
	make -C src install
	# also call the suspect perl script
	# Writing log to install_perl.log
	./setup.pl -i &> install_perl.log

# BGP dump parser (Rust), used by scripts/parsebview.sh
BVIEWPAIRS:
	cargo build --release --locked --manifest-path tools/bviewpairs/Cargo.toml
	mkdir -p bin
	cp tools/bviewpairs/target/release/bviewpairs bin/

all: JARS BVIEWPAIRS
	make -C src all

install: all CBINS

.PHONY: all install JARS CLASS CBINS BVIEWPAIRS