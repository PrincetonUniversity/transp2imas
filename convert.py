#!/usr/bin/python
"""Command-line TRANSP CDF -> IMAS converter with selectable IDS writers."""

import sys
from pathlib import Path

import mapping_status
from conversion_plan import ConversionJournal, ConversionRunner, select_ids
from writer_base import run_writer
from nc2imas import nc2imas, get_args, TRANSP2IMAS_VERSION
from utilities import find_all_files
from imas_summary import imas_summary
from imas_equilibrium import imas_equilibrium
from imas_wall import imas_wall
from imas_core_profiles import imas_core_profiles
from imas_ec_launchers import imas_ec_launchers
from imas_lh_antennas import imas_lh_antennas
from imas_ic_antennas import imas_ic_antennas
from imas_nbi import imas_nbi
from imas_tf import imas_tf
from imas_pf_active import imas_pf_active
from imas_pf_passive import imas_pf_passive
from imas_edge_sources import imas_edge_sources
from imas_radiation import imas_radiation
from imas_core_sources import imas_core_sources
from imas_waves import imas_waves
from imas_workflow import imas_workflow
from imas_dataset_fair import imas_dataset_fair
from imas_distribution_sources import imas_distribution_sources
from imas_distributions import imas_distributions
from imas_core_transport import imas_core_transport
from imas_transport_solver_numerics import imas_transport_solver_numerics
from imas_sawteeth import imas_sawteeth

TRANSP2IMAS_DD4_PATCHLEVEL = TRANSP2IMAS_VERSION


def make_writers(nci, moments, nzones, ec_nml_launch, pf_passive_geometry,
                 pf_active_geometry):
    """Map canonical IDS names to writers; preserve existing execution order."""
    return {
        'workflow': lambda: run_writer(imas_workflow(nci)),
        'wall': lambda: run_writer(imas_wall(nci)),
        'tf': lambda: run_writer(imas_tf(nci)),
        'lh_antennas': lambda: run_writer(imas_lh_antennas(nci)),
        'ic_antennas': lambda: run_writer(imas_ic_antennas(nci)),
        'ec_launchers': lambda: run_writer(imas_ec_launchers(nci, ec_nml_launch=ec_nml_launch)),
        'nbi': lambda: run_writer(imas_nbi(nci)),
        'distribution_sources': lambda: run_writer(imas_distribution_sources(nci)),
        'distributions': lambda: run_writer(imas_distributions(nci)),
        'pf_active': lambda: run_writer(imas_pf_active(nci, geometry_path=pf_active_geometry)),
        'core_profiles': lambda: run_writer(imas_core_profiles(nci)),
        'pf_passive': lambda: run_writer(imas_pf_passive(nci, geometry_path=pf_passive_geometry)),
        'radiation': lambda: run_writer(imas_radiation(nci)),
        'equilibrium': lambda: run_writer(imas_equilibrium(nci, conv_moments=moments, nzones_in=nzones)),
        'summary': lambda: run_writer(imas_summary(nci)),
        'core_sources': lambda: run_writer(imas_core_sources(nci)),
        'waves': lambda: run_writer(imas_waves(nci)),
        'core_transport': lambda: run_writer(imas_core_transport(nci)),
        'transport_solver_numerics': lambda: run_writer(imas_transport_solver_numerics(nci)),
        'sawteeth': lambda: run_writer(imas_sawteeth(nci)),
        'edge_sources': lambda: run_writer(imas_edge_sources(nci)),
        'dataset_fair': lambda: run_writer(imas_dataset_fair(nci)),
    }


def do_convert(argvs):
    print(f'transp2imas patchlevel {TRANSP2IMAS_DD4_PATCHLEVEL}')
    parsed_args = get_args(argvs)
    if len(parsed_args) == 20:
        (backend, run, occur, user, machine, newdb, ncfile, nml, log,
         moments, nzones, deb, ec_nml_launch, ip_direction, bt_direction,
         dd_version, pf_passive_geometry, pf_active_geometry, requested,
         report_json) = parsed_args
    elif len(parsed_args) == 18:
        (backend, run, occur, user, machine, newdb, ncfile, nml, log,
         moments, nzones, deb, ec_nml_launch, ip_direction, bt_direction,
         dd_version, pf_passive_geometry, pf_active_geometry) = parsed_args
        requested, report_json = None, None
    else:
        raise ValueError(f'Unexpected get_args return length: {len(parsed_args)}')

    selected = select_ids(requested)
    journal = ConversionJournal(ncfile, selected)
    try:
        with nc2imas(backend=backend, run=run, occur=occur, user=user,
                     machine=machine, new_db=newdb, ncfile=ncfile,
                     nml_file=nml, log=log, debug=deb,
                     ip_direction=ip_direction, bt_direction=bt_direction,
                     dd_version=dd_version) as nci:
            journal.attach_converter(nci)
            journal.data['setup_reads'] = nci.mapping.consumed(None)
            journal.data['mapping_policy'] = 'strict' if mapping_status.is_strict() else 'lenient'
            print(f'IMAS target DD: {nci.dd_version} (schema family DD{nci.dd_major_version})')
            writers = make_writers(nci, moments, nzones, ec_nml_launch,
                                   pf_passive_geometry, pf_active_geometry)
            ConversionRunner(nci, writers, journal).run(selected)
        journal.finish()
    except Exception as exc:
        journal.fail(exc)
        raise
    finally:
        if report_json:
            journal.save(report_json)
            print(f'conversion report: {report_json}')
    return journal.data


if __name__ == '__main__':
    argv = sys.argv
    indices = [i for i, value in enumerate(argv)
               if value == 'all' and i > 0 and argv[i - 1] in ('-i', '--input-cdf')]
    if not indices:
        do_convert(argv)
    else:
        i = indices[0]
        paths = [j + 1 for j, value in enumerate(argv) if value in ('-d', '--input-dir')]
        if not paths or paths[0] >= len(argv):
            raise ValueError('Use -i all -d <directory> to convert all CDFs')
        directory = argv[paths[0]]
        files = find_all_files('.CDF', directory)
        for file in files:
            args = list(argv)
            args[i] = file
            # A batch should not overwrite a report intended for a single run.
            for j, value in enumerate(args[:-1]):
                if value == '--report-json':
                    dest = Path(args[j + 1])
                    args[j + 1] = str(dest.with_name(dest.stem + '-' + Path(file).stem + dest.suffix))
                    break
            do_convert(args)
