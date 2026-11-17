#!/usr/bin/python
"""Real-browser checks of the workbench (opt-in: needs Playwright and a Chromium build).

Phase 1 serves a report whose analyst-written text is hostile (script and image tags, an ``onerror`` handler) and asserts
that nothing executes, the text is shown literally, no request leaves the loopback address, no console error or
Content-Security-Policy violation occurs, and the canvas and SVG panels draw.

Phase 2 drives the whole conversion workflow in the page against the bundled NSTX case: choose the CDF and IDSs, see the exact
command, get a validation error for a bad value, convert for real, follow the job, open the report, meet the overwrite guard
for the now-existing entry, then explore the converted entry: draw profiles, scrub and play time, overlay times, log axis,
a flux map, and export PNG and CSV.

    python examples/verify_gui_browser.py [--server-python /path/to/python]   # run with the Playwright interpreter
"""

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
HOSTILE = '<img src=x onerror="window.__xss=1"><script>window.__xss=2</script>'


def make_report(directory, server_python):
    """A real curation report for the synthetic fixture, with hostile text injected afterwards."""
    code = ("import sys; sys.path.insert(0, %r); sys.path.insert(0, %r)\n"
            "from verify_curation import run_screen\nfrom curation.reporting import write_report\n"
            "from curation import MANIFEST_SCHEMA\n"
            "manifest = {'schema_version': MANIFEST_SCHEMA, 'notes': 'n', 'events': [{'id': 'p', 'type': 'pellet', 'start': 0.4, "
            "'end': 0.6, 'rules': ['PROFILE.TE.FINITE']}], 'evidence_sources': [{'id': 'raw', 'description': 'd'}, {'id': 'raw2', "
            "'description': 'd', 'parents': ['raw']}]}\n"
            "report, _ = run_screen(manifest=manifest)\nwrite_report(report, %r)\n"
            % (str(ROOT), str(ROOT / 'examples'), str(directory)))
    subprocess.run([server_python, '-c', code], check=True, cwd=ROOT)
    path = Path(directory) / 'curation.json'
    document = json.loads(path.read_text())
    document['findings'][0]['message'] = HOSTILE
    document['findings'][0]['recommendation'] = HOSTILE
    document['findings'][0]['evidence'] = {'note': HOSTILE}
    document['source']['path'] = HOSTILE
    path.write_text(json.dumps(document))
    return path


def read_only_page(playwright, python, check):
    """A server without --enable-run: Convert still plans, says plainly why Run is off, and the option wording is plain."""
    with tempfile.TemporaryDirectory() as tmp:
        server = subprocess.Popen([python, '-m', 'gui', '--root', str(ROOT / 'examples' / 'data'), '--output-root', tmp], cwd=ROOT,
                                  stdout=subprocess.PIPE, text=True)
        try:
            url = next(m.group(1) for m in (re.search(r'Open: (\S+)', line) for line in server.stdout) if m)
            browser = playwright.chromium.launch()
            page = browser.new_page(viewport={'width': 1400, 'height': 1000})
            page.goto(url)
            page.wait_for_selector('select[aria-label="CDF file"] option:nth-child(2)', state='attached')
            options = page.eval_on_selector_all('select[aria-label="CDF file"] option', 'els => els.map(e => e.value)')
            page.select_option('select[aria-label="CDF file"]', value=[o for o in options if o.endswith('.CDF')][0])
            page.wait_for_selector('pre.command')
            check(page.inner_text('#mode-badge') == 'conversion off', 'the badge does not say conversion is off')
            check(page.is_disabled('#run-button'), 'Run is enabled on a server without --enable-run')
            plan = page.inner_text('.plancol')
            check('--enable-run' in plan and 'switched off' in plan, f'the page does not say why Run is off: {plan[:300]}')
            body = page.inner_text('.convert')
            check('Skip variables that cannot be mapped' in body, 'the renamed option is not shown')
            labels = page.eval_on_selector_all('.convert .check span:first-of-type, .convert .flabel', 'els => els.map(e => e.textContent)')
            check(not [label for label in labels if 'enient' in label], f'the old wording is still a label: {labels}')
            check('--lenient-mapping' in body, 'the converter flag is no longer mentioned in the help')
            # --- curation runs on this server although conversion is off: plan, run, outcome, open the report
            page.click('.view-button[data-view=curate]')
            page.wait_for_selector('select[aria-label="CDF file"] option:nth-child(2)', state='attached')
            iter_cdf = [o for o in page.eval_on_selector_all('select[aria-label="CDF file"] option', 'els => els.map(e => e.value)') if o.endswith('38530P94.CDF')]
            check(len(iter_cdf) == 1, 'the ITER CDF is not offered by the Curate view')
            check(not [o for o in page.eval_on_selector_all('select[aria-label="CDF file"] option', 'els => els.map(e => e.value)') if o.upper().endswith('PH.CDF')],
                  'a PH.CDF companion file is offered')
            page.select_option('select[aria-label="CDF file"]', value=iter_cdf[0])
            page.wait_for_selector('pre.command')
            check('curate.py' in page.inner_text('pre.command') and '--output-dir' in page.inner_text('pre.command'), 'the curation command is not shown')
            check('Plugins (--plugin) are not offered' in page.inner_text('.formcol'), 'the plugin note is missing')
            check(page.is_enabled('#run-button'), 'Run curation is disabled on a server without --enable-run')
            page.click('#run-button')
            for _ in range(120):
                if 'completed' in page.inner_text('#job-state') or 'failed' in page.inner_text('#job-state'):
                    break
                page.wait_for_timeout(250)
            check('completed' in page.inner_text('#job-state'), f'the curation did not complete: {page.inner_text("#job-state")}')
            check('NOT performed' in page.inner_text('#job-outcome'), 'the outcome does not say certification was not performed')
            check('Curation screening:' in page.inner_text('#job-log'), 'the curation log is not shown')
            page.click('#open-report')
            page.wait_for_selector('.tab')
            check(page.get_attribute('.view-button[data-view=reports]', 'aria-selected') == 'true', 'Open curation report did not open Reports')
            check('conversion report' not in page.inner_text('#view').lower(), 'the curation report opened as a conversion report')
            page.click('.view-button[data-view=convert]')
            page.wait_for_selector('#run-button')
            check(page.is_disabled('#run-button'), 'running curation enabled conversion')
            page.goto(url + '&view=reports')
            page.wait_for_selector('.tab, #view .empty')
            check(page.get_attribute('.view-button[data-view=reports]', 'aria-selected') == 'true', '?view=reports does not open Reports')
            browser.close()
        finally:
            server.terminate()


def workflow(playwright, python, check, fast_ion_writer=None):
    """Phase 2: convert and plot through the page."""
    import time
    nstx = ROOT / 'examples' / 'data' / 'nstx' / '141711Y05' / '141711Y05.CDF'
    if not nstx.is_file():
        check(False, f'bundled fixture missing: {nstx}')
        return
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / 'out'
        server = subprocess.Popen([python, '-m', 'gui', '--root', str(ROOT / 'examples' / 'data'), '--enable-run', '--imas-python', python,
                                   '--output-root', str(out)], cwd=ROOT, stdout=subprocess.PIPE, text=True)
        try:
            url = next(m.group(1) for m in (re.search(r'Open: (\S+)', line) for line in server.stdout) if m)
            base = url.split('/?')[0]
            problems, requests = [], []
            browser = playwright.chromium.launch()
            page = browser.new_page(viewport={'width': 1400, 'height': 1000}, accept_downloads=True)
            page.on('console', lambda m: problems.append((m.type, m.text)) if m.type in ('error', 'warning') else None)
            page.on('pageerror', lambda e: problems.append(('pageerror', str(e))))
            page.on('request', lambda r: requests.append(r.url))
            page.goto(url)
            page.wait_for_selector('.view-button')
            check(page.inner_text('#mode-badge') == 'conversion enabled', 'a server with --enable-run does not say so')
            page.click('.view-button[data-view=convert]')
            page.wait_for_selector('select[aria-label="CDF file"] option:nth-child(2)', state='attached')
            options = page.eval_on_selector_all('select[aria-label="CDF file"] option', 'els => els.map(e => e.value)')
            page.select_option('select[aria-label="CDF file"]', value=[o for o in options if o.endswith('141711Y05.CDF')][0])
            page.wait_for_selector('pre.command')
            page.click('button:text-is("None")')
            for name in ('wall', 'tf', 'core_profiles', 'equilibrium', 'summary'):
                page.check(f'input[aria-label="{name}"]')
            page.wait_for_timeout(800)
            command = page.inner_text('pre.command')
            check('--ids wall,tf,core_profiles,equilibrium,summary' in command, f'the command does not list the chosen IDSs: {command}')
            check(f'-u {out}/imasdb' in command and '--dd-version 4.1.1' in command and command.rstrip().endswith('-new'), f'unexpected command: {command}')
            check('entry exists' in page.inner_text('.plancol') and page.is_enabled('#run-button'), 'the plan card is incomplete or Run is disabled')
            # typing digit by digit, with a plan refresh between keystrokes, must keep the order (a rebuilt form once put each
            # new digit on the left) and the caret where the user put it
            page.click('input[aria-label="Run number"]')
            for digit in '4321':
                page.keyboard.type(digit)
                page.wait_for_timeout(700)
            run_number = page.input_value('input[aria-label="Run number"]')
            check(run_number == '4321', f'digits were not kept in order: {run_number}')
            check('4321' in page.inner_text('.plancol'), 'the plan does not show the typed run number')
            page.fill('input[aria-label="Run number"]', '')
            page.click('input[aria-label="Machine / database name"]')
            page.keyboard.type('ac')
            page.wait_for_timeout(700)
            page.keyboard.press('ArrowLeft')
            page.keyboard.type('b')
            page.wait_for_timeout(700)
            check(page.input_value('input[aria-label="Machine / database name"]') == 'abc', 'the caret was moved while typing in a text field')
            page.fill('input[aria-label="Machine / database name"]', '')
            page.wait_for_timeout(900)
            page.fill('input[aria-label="Run number"]', '-5')
            page.wait_for_timeout(900)
            check(page.query_selector('#run-button') is None and 'cannot run yet' in page.inner_text('.plancol'), 'a bad run number did not block the plan')
            page.fill('input[aria-label="Run number"]', '')
            page.wait_for_timeout(900)
            check(page.query_selector('#run-button') is not None, 'clearing the bad value did not restore the plan')
            page.click('#run-button')
            deadline = time.time() + 240
            state = ''
            while time.time() < deadline:
                page.wait_for_selector('#job-state')
                state = page.inner_text('#job-state')
                if any(word in state for word in ('succeeded', 'failed', 'cancelled')):
                    break
                page.wait_for_timeout(1000)
            check('succeeded' in state, f'the conversion did not succeed: {state}')
            check('ids stored in' in page.inner_text('#job-log'), 'the live log does not show the converter output')
            page.click('button:text-is("Open conversion report")')
            page.wait_for_selector('.tab')
            check('PASS' in page.inner_text('#view') or 'pass' in page.inner_text('#view'), 'the report view does not show a passing conversion')
            check('complete' in page.inner_text('#view'), 'the publication status is not shown')
            page.click('.view-button[data-view=convert]')
            page.wait_for_selector('.confirm')
            check(page.is_disabled('#run-button'), 'Run is enabled for an existing entry before confirmation')
            phrase = page.get_attribute('input[aria-label="confirm entry path"]', 'placeholder')
            page.fill('input[aria-label="confirm entry path"]', phrase + 'x')
            check(page.is_disabled('#run-button'), 'a wrong confirmation enabled Run')
            page.fill('input[aria-label="confirm entry path"]', phrase)
            check(page.is_enabled('#run-button'), 'the exact confirmation did not enable Run')
            page.click('.view-button[data-view=explore]')
            page.wait_for_selector('.recipes button')
            page.click('.recipes button:text-is("Electron and ion temperature")')
            page.wait_for_selector('svg.plot path.line')
            check(len(page.query_selector_all('svg.plot path.line')) >= 2, 'the temperature plot has fewer than two lines')
            first = page.inner_text('.timelabel')
            page.locator('input[type=range][aria-label="time"]').fill('4')
            check(page.inner_text('.timelabel') != first, 'moving the time slider did not change the time')
            after_slider = page.inner_text('.timelabel')
            page.click('button:text-is("Play")')
            page.wait_for_timeout(700)
            page.click('button:text-is("Pause")')
            check(page.inner_text('.timelabel') != after_slider, 'Play did not advance time')
            page.check('input[aria-label="overlay several times"]')
            check(len(page.query_selector_all('svg.plot path.line')) > 2, 'overlaying times did not add lines')
            page.uncheck('input[aria-label="overlay several times"]')
            page.check('input[aria-label="logarithmic y axis"]')
            check('(log)' in page.text_content('svg.plot'), 'the log axis is not labelled')
            with page.expect_download() as download:
                page.click('button:text-is("CSV")')
            csv = Path(download.value.path()).read_text()
            check(csv.startswith('series,x,y') and len(csv.splitlines()) > 20, 'the CSV export is not the plotted data')
            with page.expect_download() as download:
                page.click('button:text-is("PNG")')
            check(Path(download.value.path()).read_bytes()[:8] == b'\x89PNG\r\n\x1a\n', 'the PNG export is not a PNG')
            page.click('.recipes button:text-is("Poloidal flux map with flux surfaces, boundary, wall and magnetic axis")')
            page.wait_for_selector('canvas.heat')
            check(page.eval_on_selector('canvas.heat', 'c => c.width > 100 && c.height > 100'), 'the flux map canvas is empty')
            # contours: ten flux surfaces are drawn by default; the count, colouring and wall are controls that redraw the picture
            check(page.input_value('input[aria-label="number of contour lines"]') == '10', 'the flux map does not start with flux surfaces')
            shot = lambda: page.locator('canvas.heat').screenshot()
            with_ten = shot()
            page.fill('input[aria-label="number of contour lines"]', '0')
            without = shot()
            check(without != with_ten, 'removing the flux surfaces did not change the map')
            page.fill('input[aria-label="number of contour lines"]', '25')
            check(shot() not in (with_ten, without), '25 flux surfaces look like 10 or none')
            page.uncheck('input[aria-label="show wall"]')
            without_wall = shot()
            page.check('input[aria-label="show wall"]')
            check(shot() != without_wall, 'hiding the wall did not change the map')
            page.uncheck('input[aria-label="colour the map"]')
            check(shot() != without_wall, 'switching the colouring off did not change the map')
            page.check('input[aria-label="colour the map"]')
            page.click('.recipes button:text-is("Wall: limiter outline in R,Z")')
            page.wait_for_selector('canvas.shapes')
            check(page.eval_on_selector('canvas.shapes', 'c => c.width > 100 && c.height > 100'), 'the wall outline canvas is empty')
            check(page.query_selector('canvas.colorbar:visible') is None, 'a colour bar is shown for an outline without values')
            page.click('.recipes button:text-is("Density of each ion species")')
            page.wait_for_selector('svg.plot path.line')
            check(len(page.query_selector_all('svg.plot path.line')) >= 2, 'the ion density plot has fewer than two lines')
            page.click('.recipes button:text-is("Toroidal angular momentum density")')
            page.wait_for_selector('svg.plot path.line')
            # fast-ion distribution function on a small synthetic entry placed in the output database
            if fast_ion_writer is not None:
                fast_ion_writer(out / 'imasdb' / 'TESTDF' / '4' / '1' / '1')
                keep = page.input_value('select[aria-label="entry"]')
                page.click('.explorecol button:text-is("Rescan")')
                for _ in range(40):
                    if page.query_selector(f'select[aria-label="entry"] option[value$="TESTDF/4/1/1"]'):
                        break
                    page.wait_for_timeout(250)
                page.select_option('select[aria-label="entry"]', value=str(out / 'imasdb' / 'TESTDF' / '4' / '1' / '1'))
                page.wait_for_selector('#df-ep')
                page.click('#df-ep')
                page.wait_for_selector('canvas.heat')
                check(page.query_selector('input[aria-label="number of contour lines"]') is not None, 'the distribution map has no contour control')
                page.click('#df-energy')
                page.wait_for_selector('svg.plot path.line')
                page.click('#df-zones')
                page.wait_for_selector('canvas.shapes')
                canvas = page.locator('canvas.shapes')
                box = canvas.bounding_box()
                canvas.click(position={'x': box['width'] * 0.35, 'y': box['height'] * 0.5})
                for _ in range(40):
                    if page.input_value('input[aria-label="Zone"]') != '':
                        break
                    page.wait_for_timeout(250)
                check(page.input_value('input[aria-label="Zone"]') in ('0', '1'), 'clicking a zone did not select it')
                page.click('#df-pitch')
                page.wait_for_selector('svg.plot path.line')
                page.select_option('select[aria-label="entry"]', value=keep)
                page.wait_for_selector('#df-ep', state='detached')
            page.click('button:text-is("Browse nodes") >> nth=0')
            page.wait_for_selector('.nodelist .node')
            page.fill('input[aria-label="filter nodes"]', 'electrons/temperature')
            page.check('.nodelist input[type=checkbox] >> nth=0')
            page.click('button:text-is("Plot selected")')
            page.wait_for_selector('svg.plot path.line')
            # --- the existing entry: pick a free run number instead, or overwrite deliberately and see the new contents
            page.click('.view-button[data-view=convert]')
            page.wait_for_selector('.confirm')
            suggestion = page.locator('button', has_text='instead (nothing exists there)')
            check(suggestion.count() == 1, 'no free run number is offered for an existing entry')
            suggestion.click()
            page.wait_for_timeout(1200)
            check(page.input_value('input[aria-label="Run number"]') == '2506', 'the suggested run number was not filled in')
            check(page.query_selector('.confirm') is None and page.is_enabled('#run-button'), 'a free run number still asks for confirmation')
            page.fill('input[aria-label="Run number"]', '')
            page.wait_for_selector('.confirm')
            page.click('button:text-is("None")')
            for name in ('wall', 'tf'):
                page.check(f'input[aria-label="{name}"]')
            page.wait_for_selector('.confirm')
            check(page.is_disabled('#run-button'), 'Run is enabled before the overwrite is confirmed')
            page.click('button[aria-label="fill in the entry path"]')
            check(page.is_enabled('#run-button'), 'Fill in the path did not allow the overwrite')
            first_job = page.inner_text('#job-id')
            page.click('#run-button')
            for _ in range(60):                                    # wait for the NEW job to replace the previous one on the page
                if page.inner_text('#job-id') != first_job:
                    break
                page.wait_for_timeout(250)
            check(page.inner_text('#job-id') != first_job, 'starting the overwrite did not start a new job')
            deadline = time.time() + 240
            state = ''
            while time.time() < deadline:
                state = page.inner_text('#job-state')
                if any(word in state for word in ('succeeded', 'failed', 'cancelled')):
                    break
                page.wait_for_timeout(1000)
            check('succeeded' in state, f'the overwrite did not succeed: {state}')
            page.click('.view-button[data-view=explore]')
            page.wait_for_selector('.idsblock')
            page.wait_for_selector('#refreshing', state='attached', timeout=5000)           # it says it is re-reading the entry
            page.wait_for_selector('#refreshing', state='detached', timeout=60000)
            names = sorted(page.eval_on_selector_all('.idshead strong', 'els => els.map(e => e.textContent)'))
            check(names == ['tf', 'wall'], f'the overwritten entry still shows old contents: {names}')
            check(all(u.startswith(base) for u in requests), f'requests left the viewer: {[u for u in requests if not u.startswith(base)]}')
            # the only console messages allowed are the browser's note about the 400 answers to the deliberately invalid plan
            expected = [p for p in problems if p[0] == 'error' and 'status of 400 (Bad Request)' in p[1]]
            others = [p for p in problems if p not in expected]
            check(1 <= len(expected) <= 3 and not others, f'unexpected console messages: {others} (400 answers: {len(expected)})')
            browser.close()
        finally:
            server.terminate()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--server-python', default=sys.executable, help='interpreter that runs the viewer (default: this one)')
    args = parser.parse_args(argv)
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print('Playwright is not installed in this interpreter; the browser check was NOT run.', file=sys.stderr)
        return 2
    failures = []

    def check(condition, message):
        if not condition:
            failures.append(message)

    with tempfile.TemporaryDirectory() as tmp:
        make_report(tmp, args.server_python)
        server = subprocess.Popen([args.server_python, '-m', 'gui', '--root', tmp, '--output-root', str(Path(tmp) / 'out')], cwd=ROOT,
                                  stdout=subprocess.PIPE, text=True)      # its own output root: nothing else may show up as a report
        try:
            url = next(m.group(1) for m in (re.search(r'Open: (\S+)', line) for line in server.stdout) if m)
            problems, requests, dialogs = [], [], []
            with sync_playwright() as playwright:
                browser = playwright.chromium.launch()
                page = browser.new_page(viewport={'width': 1280, 'height': 900})
                page.on('console', lambda m: problems.append((m.type, m.text)) if m.type in ('error', 'warning') else None)
                page.on('pageerror', lambda e: problems.append(('pageerror', str(e))))
                page.on('request', lambda r: requests.append(r.url))
                page.on('dialog', lambda d: (dialogs.append(d.message), d.dismiss()))
                page.goto(url)
                page.wait_for_selector('.view-button')
                labels = page.eval_on_selector_all('.view-button', 'els => els.map(e => e.textContent)')
                check(labels == ['Convert', 'Curate', 'Explore', 'Reports'], f'view order: {labels}')
                check(page.get_attribute('.view-button[data-view=convert]', 'aria-selected') == 'true', 'the page does not land on Convert')
                page.click('.view-button[data-view=reports]')
                page.wait_for_selector('.tab')                       # one report under the root: opened automatically
                tabs = page.eval_on_selector_all('.tab', 'els => els.map(e => e.textContent)')
                check(tabs == ['Overview', 'Findings', 'Rules', 'Evidence', 'Events', 'Masks', 'Inputs'], f'tabs: {tabs}')
                check(page.inner_text('#mode-badge') == 'conversion off', 'a default server does not say conversion is off')
                for name in tabs:
                    page.click(f'.tab:text-is("{name}")')
                    page.wait_for_timeout(150)
                    if name == 'Findings':
                        page.fill('input[aria-label="search findings"]', '')
                        for box in page.query_selector_all('.toolbar input[type=checkbox]'):
                            if not box.is_checked():
                                box.check()
                        page.click('details.finding >> nth=0 >> summary')
                        body = page.inner_text('#view')
                        check(HOSTILE in body, 'hostile text was not shown literally in Findings')
                    if name == 'Overview':
                        check(HOSTILE in page.inner_text('#view'), 'hostile source path was not shown literally in Overview')
                    if name == 'Evidence':
                        check(page.query_selector('svg g.node') is not None, 'no evidence graph node drawn')
                    if name == 'Events':
                        check(page.query_selector('svg rect.ev-declared') is not None, 'no declared event drawn')
                    if name == 'Masks':
                        page.wait_for_selector('canvas')
                        check(page.eval_on_selector('canvas', 'c => c.width > 0 && c.height > 0'), 'mask canvas has no size')
                check(page.evaluate('window.__xss') is None, 'hostile markup executed script')
                check(page.query_selector('img[src="x"]') is None, 'hostile text was parsed into an element')
                check(not dialogs, f'dialogs opened: {dialogs}')
                check(all(u.startswith(url.split('/?')[0]) for u in requests), f'requests left the viewer: {requests}')
                check(not problems, f'console problems: {problems}')
                browser.close()
        finally:
            server.terminate()
    if not failures:
        with sync_playwright() as playwright:
            read_only_page(playwright, args.server_python, check)
            def write_fast_ion_entry(directory):
                code = ('import sys; sys.path[:0] = [%r, %r]; from verify_gui_plots import FastIonSyntheticEntryTests as T; T.write_entry(%r)'
                        % (str(ROOT), str(ROOT / 'examples'), str(directory)))
                subprocess.run([args.server_python, '-c', code], check=True, cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            workflow(playwright, args.server_python, check, write_fast_ion_entry)
    if failures:
        for failure in failures:
            print('FAIL:', failure)
        return 1
    print('Browser check: OK (reports render, hostile text inert, convert-and-plot workflow works, no external requests, no unexpected console messages)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
