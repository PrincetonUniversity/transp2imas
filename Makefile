SHELL := /bin/bash

PYTHON ?= python3
VENV ?= .venv
VENV_PYTHON := $(VENV)/bin/python
DATA_ROOT ?= $(CURDIR)/examples/data/nstx
CASE ?= 141711Y01
CASES ?= 141711Y01 141711Y02
CASE_DIR ?= $(DATA_ROOT)/$(CASE)
CDF ?= $(CASE_DIR)/$(CASE).CDF
TRDAT ?= $(CASE_DIR)/$(CASE)TR.DAT
BACKEND ?= HDF5
DD_VERSION ?= auto
OUTPUT_DIR ?= $(CURDIR)/build/imas
CONVERT_FLAGS ?= -new
PYTHON_FILES := $(shell find . -type f -name '*.py' -not -path './.venv/*' -not -path './build/*' | sort)

.DEFAULT_GOAL := help

.PHONY: help venv install convert check check-shebang check-literals test verify verify-all verify-dd audit inspect plot-ids fortran clean

help:
	@echo "TRANSP-to-IMAS developer commands"
	@echo
	@echo "  make venv         Create a local Python virtual environment"
	@echo "  make install      Install requirements into the local environment"
	@echo "  make convert      Convert CDF into OUTPUT_DIR (fresh HDF5 entry by default)"
	@echo "  make test         Verify and convert both NSTX regression cases"
	@echo "  make check        Enforce Python headers, compile-check Python, and run fast regressions"
	@echo "  make verify       Run source diagnostics for CASE"
	@echo "  make verify-all   Run source diagnostics for both NSTX cases"
	@echo "  make verify-dd    Check supported IMAS Data Dictionary layouts"
	@echo "  make audit        Write CDF mapping-coverage reports for CASE"
	@echo "  make inspect      Plot the final six q profiles from CDF"
	@echo "  make plot-ids     Print usage for the converted-IDS plotting example"
	@echo "  make fortran      Build the optional Fortran Access Layer verifier"
	@echo "  make clean        Remove generated local build/example artifacts"
	@echo
	@echo "Overrides: CASE, CASES, DATA_ROOT, CDF, TRDAT, OUTPUT_DIR, BACKEND, DD_VERSION, PYTHON"

venv:
	$(PYTHON) -m venv "$(VENV)"
	$(VENV_PYTHON) -m pip install --upgrade pip setuptools wheel

install: venv
	$(VENV_PYTHON) -m pip install -r requirements.txt

convert:
	@test -f "$(CDF)" || { echo "CDF not found: $(CDF)" >&2; exit 1; }
	@mkdir -p "$(OUTPUT_DIR)"
	$(PYTHON) convert.py -i "$(CDF)" -b "$(BACKEND)" -u "$(OUTPUT_DIR)" \
		--dd-version "$(DD_VERSION)" $(CONVERT_FLAGS)

check-shebang:
	@bad=0; \
	for file in $(PYTHON_FILES); do \
		first=$$(head -n 1 "$$file"); \
		if [ "$$first" != "#!/usr/bin/python" ]; then \
			echo "Invalid Python header: $$file: $$first" >&2; \
			bad=1; \
		fi; \
	done; \
	if [ $$bad -ne 0 ]; then exit 1; fi
	@echo "Python header policy OK: $$(printf '%s\n' $(PYTHON_FILES) | wc -l) files"

check-literals:
	$(PYTHON) examples/verify_python_literals.py

check: check-shebang check-literals
	$(PYTHON) -m compileall -q -x '(^|/)(\.venv|build)/' .
	$(MAKE) -C examples quick PYTHON="$(PYTHON)"

test: check verify-all
	@set -e; for case in $(CASES); do \
		echo "Converting NSTX case $$case"; \
		$(MAKE) --no-print-directory convert CASE="$$case"; \
	done

verify:
	$(MAKE) -C examples verify PYTHON="$(PYTHON)" CDF="$(CDF)" TRDAT="$(TRDAT)"

verify-all:
	$(MAKE) -C examples verify-all PYTHON="$(PYTHON)" CASES="$(CASES)"

verify-dd:
	$(MAKE) -C examples verify-dd PYTHON="$(PYTHON)"

audit:
	$(MAKE) -C examples audit PYTHON="$(PYTHON)" CDF="$(CDF)" CASE="$(CASE)"

inspect:
	$(MAKE) -C examples inspect PYTHON="$(PYTHON)" CDF="$(CDF)"

plot-ids:
	$(MAKE) -C examples plot-ids PYTHON="$(PYTHON)"

fortran:
	$(MAKE) -C examples fortran

clean:
	$(RM) -r __pycache__ build
	$(MAKE) -C examples clean
