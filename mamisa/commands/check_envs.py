#!/usr/bin/env python3
"""
MaMISA - check-envs command

Sanity-check the conda environments a user already has: for each env it reports
the version of every relevant tool, whether MaMISA is importable there, and
whether the versions are compatible with MaMISA's wrappers. On any problem it
prints exactly what to install.

Roles and what each env must provide:
  mamisa   light CLI tools (samtools, blastn, bedtools, kraken2, meryl) — optional
  checkm2  checkm2 (+ diamond)
  gtdbtk   gtdbtk
  gunc     gunc (runnable) + diamond + prodigal, with pandas 2.x / numpy < 2

Checks run via `conda run -n <env> ...`; no env is modified.
"""

import os
import re
import sys
import shutil
import argparse
import subprocess
from pathlib import Path

from .. import PACKAGE_ROOT
from ..utils.logging import log_info, log_warning, log_error, print_header, print_section


PY_MIN = (3, 9)

# command to print each tool's version
TOOL_VERSION_CMD = {
    'samtools': ['samtools', '--version'],
    'blastn':   ['blastn', '-version'],
    'bedtools': ['bedtools', '--version'],
    'kraken2':  ['kraken2', '--version'],
    'meryl':    ['meryl', '--version'],
    'checkm2':  ['checkm2', '--version'],
    'gtdbtk':   ['gtdbtk', '--version'],
    'gunc':     ['gunc', '--version'],
    'diamond':  ['diamond', 'version'],
    'prodigal': ['prodigal', '-v'],
    'barrnap':  ['barrnap', '--version'],
    'tRNAscan-SE': ['tRNAscan-SE', '-h'],
}

# per role: required tools (FAIL if missing) and optional (WARN if missing)
ROLE_TOOLS = {
    'mamisa':  {'required': [], 'optional': ['samtools', 'blastn', 'bedtools', 'kraken2', 'meryl']},
    'checkm2': {'required': ['checkm2'], 'optional': ['diamond']},
    'gtdbtk':  {'required': ['gtdbtk'], 'optional': []},
    'gunc':    {'required': ['gunc', 'diamond', 'prodigal'], 'optional': []},
    'mimag':   {'required': ['barrnap', 'tRNAscan-SE'], 'optional': []},
}

INSTALL_HINT = {
    'samtools': 'conda install -n {env} -c bioconda samtools',
    'blastn':   'conda install -n {env} -c bioconda blast',
    'bedtools': 'conda install -n {env} -c bioconda bedtools',
    'kraken2':  'conda install -n {env} -c bioconda kraken2',
    'meryl':    'conda install -n {env} -c bioconda meryl',
    'checkm2':  'conda install -n {env} -c bioconda checkm2',
    'gtdbtk':   'conda install -n {env} -c bioconda gtdbtk',
    'diamond':  'conda install -n {env} -c bioconda diamond',
    'prodigal': 'conda install -n {env} -c bioconda prodigal',
    'gunc':     ('recreate the env pinned: mamba create -n {env} -c conda-forge -c bioconda '
                 '"gunc=1.1.1" "python=3.10" "pandas>=2,<3" "numpy<2" diamond prodigal'),
    'barrnap':  'conda install -n {env} -c bioconda barrnap',
    'tRNAscan-SE': 'conda install -n {env} -c bioconda trnascan-se',
}

_VER_RE = re.compile(r'(\d+)\.(\d+)(?:\.(\d+))?')
_BROKEN = ('Traceback', 'ImportError', 'ModuleNotFoundError', 'AttributeError',
           'cannot import name')


def _conda_run(env, cmd, timeout=120):
    """Run `conda run -n env cmd`; return (rc, combined_output)."""
    full = ['conda', 'run', '-n', env] + cmd
    try:
        r = subprocess.run(full, capture_output=True, text=True, timeout=timeout)
        return r.returncode, (r.stdout or '') + (r.stderr or '')
    except subprocess.TimeoutExpired:
        return 124, 'timeout'
    except FileNotFoundError:
        return 127, 'conda not found'


def _parse_ver(text):
    m = _VER_RE.search(text or '')
    if not m:
        return None
    return tuple(int(x) for x in m.groups(default='0'))


def _fmt_ver(t):
    return '.'.join(str(x) for x in t) if t else '?'


def check_env(env, role, repo):
    """
    Return {'env', 'role', 'rows': [(check, status, detail, fix)], 'ok': bool}.
    status in {'OK','WARN','FAIL'}; ok is False if any FAIL.
    """
    rows = []

    def add(check, status, detail='', fix=''):
        rows.append((check, status, detail, fix))

    # 0. env reachable
    rc, out = _conda_run(env, ['python', '-c', 'import sys;print(sys.version.split()[0])'])
    if rc == 127:
        add('conda', 'FAIL', 'conda not found in PATH', 'install conda/mamba')
        return {'env': env, 'role': role, 'rows': rows, 'ok': False}
    if rc != 0:
        add('env', 'FAIL', f"env '{env}' not reachable", f'conda create -n {env} ...')
        return {'env': env, 'role': role, 'rows': rows, 'ok': False}

    pyver = _parse_ver(out)
    if pyver and pyver >= PY_MIN:
        add('python', 'OK', _fmt_ver(pyver))
    else:
        add('python', 'FAIL', _fmt_ver(pyver),
            f'needs >= {_fmt_ver(PY_MIN)} (MaMISA requires it)')

    # 1. MaMISA importable in this env
    rc, out = _conda_run(env, ['python', '-c', 'import mamisa;print(mamisa.__version__)'])
    if rc == 0:
        add('mamisa', 'OK', _parse_ver(out) and _fmt_ver(_parse_ver(out)) or 'installed')
    else:
        add('mamisa', 'FAIL', 'not importable',
            f'conda run -n {env} pip install -e {repo}')

    # 2. tools
    spec = ROLE_TOOLS.get(role, {'required': [], 'optional': []})
    for tool in spec['required'] + spec['optional']:
        required = tool in spec['required']
        rc, out = _conda_run(env, TOOL_VERSION_CMD[tool])
        if rc in (0,) and _parse_ver(out):
            add(tool, 'OK', _fmt_ver(_parse_ver(out)))
        elif rc == 0:
            add(tool, 'OK', 'present')
        else:
            add(tool, 'FAIL' if required else 'WARN',
                'missing', INSTALL_HINT.get(tool, '').format(env=env))

    # 3. role-specific library constraints
    if role == 'gunc':
        # gunc must actually run (old code vs modern pandas)
        rc, out = _conda_run(env, ['gunc', '--version'])
        if rc != 0 or any(m in out for m in _BROKEN):
            add('gunc runnable', 'FAIL', 'gunc errors on startup',
                INSTALL_HINT['gunc'].format(env=env))
        else:
            add('gunc runnable', 'OK')

        rc, out = _conda_run(env, [
            'python', '-c',
            'import pandas,numpy;print(pandas.__version__);print(numpy.__version__)'])
        lines = [l for l in out.splitlines() if _VER_RE.match(l.strip())]
        pv = _parse_ver(lines[0]) if len(lines) >= 1 else None
        nv = _parse_ver(lines[1]) if len(lines) >= 2 else None
        if pv and (2, 0) <= pv < (3, 0):
            add('pandas', 'OK', _fmt_ver(pv))
        else:
            add('pandas', 'FAIL', _fmt_ver(pv),
                f'GUNC needs pandas 2.x: conda install -n {env} "pandas>=2,<3"')
        if nv and nv < (2, 0):
            add('numpy', 'OK', _fmt_ver(nv))
        elif nv:
            add('numpy', 'WARN', _fmt_ver(nv),
                f'verified combo uses numpy<2: conda install -n {env} "numpy<2"')

    ok = not any(s == 'FAIL' for _, s, _, _ in rows)
    return {'env': env, 'role': role, 'rows': rows, 'ok': ok}


_MARK = {'OK': '✓', 'WARN': '!', 'FAIL': '✗'}


def render(result):
    """Print one env's report. Returns result['ok']."""
    env, role = result['env'], result['role']
    status = 'PASS' if result['ok'] else 'FAIL'
    print(f"\n  [{role}]  env='{env}'  →  {status}")
    for check, st, detail, fix in result['rows']:
        line = f"    {_MARK.get(st,'?')} {check:<16} {st:<4} {detail}"
        print(line)
        if st in ('FAIL', 'WARN') and fix:
            print(f"        fix: {fix}")
    return result['ok']


def run_checks(envs: dict, repo: Path):
    """
    envs: {role: env_name} for roles to check. Returns overall ok bool.
    """
    overall = True
    for role in ('mamisa', 'checkm2', 'gtdbtk', 'gunc', 'mimag'):
        env = envs.get(role)
        if not env:
            continue
        result = check_env(env, role, repo)
        overall = render(result) and overall
    return overall


def register_parser(subparsers):
    parser = subparsers.add_parser(
        'check-envs',
        help='Sanity-check existing conda envs: tool versions + MaMISA compatibility',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  mamisa check-envs --mamisa-env mamisa --checkm2-env checkm2 \\
                    --gtdbtk-env gtdbtk --gunc-env gunc

  # just one env
  mamisa check-envs --gunc-env my_gunc
"""
    )
    parser.add_argument('--repo', type=Path, default=None,
                        help='MaMISA checkout used in fix hints (default: auto-detected)')
    parser.add_argument('--mamisa-env', help='env with the light CLI tools + MaMISA')
    parser.add_argument('--checkm2-env', help='env with CheckM2')
    parser.add_argument('--gtdbtk-env', help='env with GTDB-Tk')
    parser.add_argument('--gunc-env', help='env with GUNC')
    parser.add_argument('--mimag-env', help='env with barrnap + tRNAscan-SE')
    parser.set_defaults(func=run)
    return parser


def run(args):
    print_header("MaMISA - Check Environments")

    if shutil.which('conda') is None:
        log_error("conda not found in PATH — required to inspect envs")
        sys.exit(1)

    repo = (args.repo or PACKAGE_ROOT.parent).resolve()
    envs = {
        'mamisa': args.mamisa_env,
        'checkm2': args.checkm2_env,
        'gtdbtk': args.gtdbtk_env,
        'gunc': args.gunc_env,
        'mimag': args.mimag_env,
    }
    envs = {r: e for r, e in envs.items() if e}
    if not envs:
        log_error("No envs given. Pass at least one of "
                  "--mamisa-env/--checkm2-env/--gtdbtk-env/--gunc-env/--mimag-env")
        sys.exit(1)

    print_section("Sanity check")
    ok = run_checks(envs, repo)

    print_section("SUMMARY")
    if ok:
        log_info("✓ All checked envs are compatible with MaMISA")
    else:
        log_error("Some checks FAILED — see the fix hints above")
        sys.exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    register_parser(parser.add_subparsers(dest='command'))
    args = parser.parse_args()
    run(args)
