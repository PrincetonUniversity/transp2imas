SHELL := /bin/bash

ROOT := $(abspath ..)
PYTHON ?= python3
DATA_ROOT ?= $(CURDIR)/data/nstx
CASE ?= 141711Y01
CASES ?= 141711Y01 141711Y02
CASE_DIR ?= $(DATA_ROOT)/$(CASE)
CDF ?= $(CASE_DIR)/$(CASE).CDF
TRDAT ?= $(CASE_DIR)/$(CASE)TR.DAT
PLOT ?= $(CURDIR)/q-$(CASE).png
AUDIT_DIR ?= $(ROOT)/build/audit
AUDIT_MD ?= $(AUDIT_DIR)/$(CASE)-mapping-coverage.md
AUDIT_CSV ?= $(AUDIT_DIR)/$(CASE)-mapping-coverage.csv
AUDIT_JSON ?= $(AUDIT_DIR)/$(CASE)-mapping-coverage.json
WAVES_CDFS = $(shell find data -type f -name "*.CDF" | sort)

F90 ?= gfortran
FFLAGS ?= -cpp -g -O0 -fPIE
CPPFLAGS += -I$(IMAS_PREFIX)/include/gfortran
LDLIBS += -L$(IMAS_PREFIX)/lib -limas-gfortran -limas

.DEFAULT_GOAL := help

.PHONY: help quick verify verify-all verify-deep verify-waves verify-neoclassical verify-neutrals verify-theory verify-mmm verify-pt-solver verify-sawteeth verify-pf-passive verify-pf-active verify-isolver verify-core-profiles verify-equilibrium verify-dd verify-end-to-end audit inspect plot-ids fortran clean

help:
	@echo "Diagnostic and example commands"
	@echo
	@echo "  make quick       Run checks that do not scan the full CDF"
	@echo "  make verify      Run source-side checks for CASE"
	@echo "  make verify-all  Discover, validate, and report all bundled CDF examples"
	@echo "  make verify-deep Run optional slow per-CDF source diagnostics on ALL CDFs"
	@echo "  make verify-waves Validate every bundled CDF for EC/IC wave deposition"
	@echo "  make verify-core-profiles  Verify Te/Ti peaking on TIME/TIME3 for every CDF"
	@echo "  make verify-equilibrium    Verify equilibrium TIME3 slices and interpolated TIME boundary moments"
	@echo "  make verify-neoclassical  Check neoclassical fluxes separately on every CDF"
	@echo "  make verify-neutrals  Check neutral profiles only, separately on every CDF"
	@echo "  make verify-theory  Check MMM07, GLF23 and TGLF family status for every CDF"
	@echo "  make verify-mmm  Alias for theory transport verifier"
	@echo "  make verify-pt-solver  Check PT_SOLVER regional model choices, equation limits and updates"
	@echo "  make verify-sawteeth  Verify TLASTSAW crash histories and q=1 native TIME diagnostics"
	@echo "  make verify-inputs Validate embedded TR.DAT and model namelists"
	@echo "  make verify-pf-passive Validate synthetic PS_n currents and DD3/DD4 passive geometries"
	@echo "  make verify-pf-active  Validate PF coil/supply/circuit source mappings"
	@echo "  make verify-isolver    Require BOTH real ITER/NSTX ISolver sources; write PF reports"
	@echo "  make verify-dd   Inspect DD 3.39.0 and DD 4.1.1 schemas"
	@echo "  make audit       Write full CDF mapping-coverage reports for CASE"
	@echo "  make inspect     Write a q-profile plot to PLOT"
	@echo "  make plot-ids    Show the converted-IDS plotting example usage"
	@echo "  make fortran     Build build/verify using IMAS_PREFIX"

quick:
	$(PYTHON) -m py_compile $(wildcard *.py)
	$(PYTHON) verify_writer_api.py
	$(PYTHON) verify_pf_passive_mapping.py
	$(PYTHON) verify_pf_active_mapping.py
	$(PYTHON) verify_isolver_fixture_policy.py
	$(PYTHON) verify_input_provenance.py
	$(PYTHON) verify_plot_ids.py
	$(PYTHON) verify_core_profiles_timegrids.py --quick
	$(PYTHON) verify_equilibrium_timegrids.py --quick
	$(PYTHON) verify_neoclassical_transport.py --quick
	$(PYTHON) verify_neutral_profiles.py --quick
	$(PYTHON) verify_theory_transport.py --quick
	$(PYTHON) verify_pt_solver_regions.py --quick
	$(PYTHON) verify_release_config.py
	$(PYTHON) verify_sawteeth_mapping.py --quick
	$(PYTHON) verify_end_to_end_contract.py

verify: quick
	@test -f "$(CDF)" || { echo "CDF not found: $(CDF)" >&2; exit 1; }
	@test -f "$(TRDAT)" || { echo "TR.DAT not found: $(TRDAT)" >&2; exit 1; }
	$(PYTHON) verify_distribution_sources_namelist.py "$(TRDAT)"
	$(PYTHON) verify_cdf_extension.py "$(CDF)" "$(TRDAT)"
	$(PYTHON) verify_cocos17.py "$(CDF)" --dat "$(TRDAT)"
	$(PYTHON) verify_nbi_distribution_sources.py "$(CDF)"
	$(PYTHON) verify_sources_transport.py "$(CDF)"
	$(PYTHON) verify_next_ids_inventory.py "$(CDF)" "$(TRDAT)"
	$(PYTHON) verify_metadata_provenance_inputs.py "$(CDF)" "$(TRDAT)"

verify-all: quick
	$(PYTHON) verify_all_examples.py --quick --report-dir "$(ROOT)/build/verification" --python "$(PYTHON)"
	$(PYTHON) verify_pt_solver_regions.py --data-root data --report-dir "$(ROOT)/build/pt-solver-regions"
	$(PYTHON) verify_isolver_mapping.py --allow-missing --report-dir "$(ROOT)/build/pf-isolver"

verify-deep: quick
	$(PYTHON) verify_all_examples.py --report-dir "$(ROOT)/build/verification-deep" --python "$(PYTHON)"

verify-waves:
	@test -n "$(WAVES_CDFS)" || { echo "No bundled CDF files available" >&2; exit 1; }
	$(PYTHON) verify_waves_mapping.py $(WAVES_CDFS)

verify-inputs:
	$(PYTHON) verify_input_provenance.py

verify-pf-passive:
	$(PYTHON) verify_pf_passive_mapping.py

verify-pf-active:
	$(PYTHON) verify_pf_active_mapping.py

verify-core-profiles:
	$(PYTHON) verify_core_profiles_timegrids.py --data-root data

verify-equilibrium:
	$(PYTHON) verify_equilibrium_timegrids.py --data-root data

verify-neoclassical:
	$(PYTHON) verify_neoclassical_transport.py --data-root data --report-dir "$(ROOT)/build/neoclassical-transport"

verify-neutrals:
	$(PYTHON) verify_neutral_profiles.py --data-root data --report-dir "$(ROOT)/build/neutral-profiles"

verify-theory verify-mmm:
	$(PYTHON) verify_theory_transport.py --data-root data --report-dir "$(ROOT)/build/theory-transport"

verify-pt-solver:
	$(PYTHON) verify_pt_solver_regions.py --data-root data --report-dir "$(ROOT)/build/pt-solver-regions"

verify-sawteeth:
	$(PYTHON) verify_sawteeth_mapping.py --data-root data --report-dir "$(ROOT)/build/sawteeth"

verify-isolver:
	$(PYTHON) verify_isolver_mapping.py --report-dir "$(ROOT)/build/pf-isolver"

verify-dd:
	$(PYTHON) verify_dd_compat.py 3.39.0 4.1.1

verify-end-to-end:
	$(PYTHON) verify_end_to_end.py --python "$(PYTHON)" --dd-version "$(if $(filter auto,$(DD_VERSION)),4.1.1,$(DD_VERSION))" \
		--report-dir "$(ROOT)/build/end-to-end"

audit:
	@test -f "$(CDF)" || { echo "CDF not found: $(CDF)" >&2; exit 1; }
	@mkdir -p "$(AUDIT_DIR)"
	$(PYTHON) audit_mapping_coverage.py "$(CDF)" --repo-root "$(ROOT)" \
		--markdown "$(AUDIT_MD)" --csv "$(AUDIT_CSV)" --json "$(AUDIT_JSON)"

inspect:
	$(PYTHON) inspect_cdf.py "$(CDF)" --output "$(PLOT)"

plot-ids:
	@echo "Example: $(PYTHON) examples/plot_ids.py 38530P94 ITER HDF5 psirz --time-index 0"
	@echo "Save:    $(PYTHON) examples/plot_ids.py 38530P94 ITER HDF5 q --time-index 0 --save build/iter_q.png"
	@echo "Help:    $(PYTHON) examples/plot_ids.py --help"

build:
	@mkdir -p $@

build/verify: verify.f90 | build
	$(F90) $(FFLAGS) $(CPPFLAGS) $< $(LDLIBS) -o $@

fortran: build/verify

clean:
	$(RM) -r __pycache__ build q-*.png
