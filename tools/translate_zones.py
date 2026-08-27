#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Translate the MATLAB zone solvers to Python, mechanically.

The zone expressions are Maple output, thousands of characters each. They are not
retyped by hand anywhere. This script does operator substitution only, so the
algebra that reaches Python is character for character the algebra in the .m file.
Correctness is then established numerically by tests/test_zones.py, which compares
against MATLAB over a random parameter sweep.

Usage:  python tools/translate_zones.py <matlab_dir> <out_py>
"""
import os, re, sys

ZONES = ["111", "211", "212", "221", "222", "311", "312", "321", "322",
         "411", "412", "421", "422", "4222"]

HEADER = '''"""
Fourteen stage closed-form solvers, neutral axis and moment.

Generated from the MATLAB sources by tools/translate_zones.py. The expressions are
Maple output and are reproduced by operator substitution only, never retyped. Every
zone is checked against the MATLAB original in tests/test_zones.py.

sqrt is emitted as np.emath.sqrt, because MATLAB returns a complex root for a
negative argument and these expressions rely on it, taking real() at the end.

Zone naming is zone<T><C><R>, the tension, compression and reinforcement states.
zone4222 carries the extra compression-steel state.

The section carries a quad-linear tension law, a bilinear compression law, bilinear
or linear-elastic bars, and a bonded FRP skin that activates at iota. The skin enters
through psi, rho_x and iota, and creates no new zone.
"""

import numpy as np

__all__ = [%s]


def _as_array(x):
    a = np.atleast_1d(np.asarray(x, dtype=float))
    return a
'''


def strip_continuations(src):
    """Join MATLAB ... line continuations."""
    return re.sub(r'\.\.\.\s*\n\s*', ' ', src)


def translate_expr(e):
    """MATLAB expression -> Python. Operator substitution only."""
    e = e.replace('.^', '**').replace('^', '**')
    e = e.replace('.*', '*').replace('./', '/')
    # function names that differ
    # MATLAB returns a complex root for a negative argument and these expressions
    # rely on it, taking real() at the end. np.sqrt would give nan, np.emath.sqrt
    # matches MATLAB.
    e = re.sub(r'\bsqrt\s*\(', 'np.emath.sqrt(', e)
    e = re.sub(r'\babs\s*\(', 'np.abs(', e)
    e = re.sub(r'\breal\s*\(', 'np.real(', e)
    e = re.sub(r'\bexp\s*\(', 'np.exp(', e)
    e = re.sub(r'\blog\s*\(', 'np.log(', e)
    return e


def parse_zone(path):
    """Pull the signature and the k and M expressions out of one zone file."""
    src = strip_continuations(open(path, encoding='utf-8', errors='replace').read())

    m = re.search(r'function\s*\[([^\]]*)\]\s*=\s*(\w+)\s*\(([^)]*)\)', src)
    if not m:
        raise ValueError(f'no function signature in {path}')
    name = m.group(2)
    args = [a.strip() for a in m.group(3).split(',')]

    k_expr = re.search(r'^\s*k\(i,1\)\s*=\s*(.+?);\s*$', src, re.M)
    M_expr = re.search(r'^\s*M\(i,1\)\s*=\s*(.+?);\s*$', src, re.M)
    if not (k_expr and M_expr):
        raise ValueError(f'could not find k(i,1) or M(i,1) in {path}')

    # how the file de-normalises the moment at the end
    per_step = bool(re.search(r'M_final\s*=\s*real\(M\)\s*\.\*', src))
    return name, args, translate_expr(k_expr.group(1)), translate_expr(M_expr.group(1)), per_step


def emit(name, args, k_src, M_src, per_step):
    sig = ', '.join(args)
    driver = 'beta_array'
    body = [f'def {name}({sig}):',
            f'    """Zone {name[4:]}. Returns (k, M) as arrays over beta_array."""',
            f'    {driver} = _as_array({driver})',
            f'    k = np.zeros({driver}.shape, dtype=complex)',
            f'    M = np.zeros({driver}.shape, dtype=complex)',
            f'    for i in range({driver}.size):',
            f'        beta = {driver}[i]',
            f'        k[i] = {k_src}',
            f'        M[i] = {M_src}',
            f'    k_final = np.real(k)']
    if per_step:
        body.append('    M_final = np.real(M) * (epcr * E * b * h**2 / (12 * (1 - k_final)))')
    else:
        body.append('    M_cr_local = epcr * E * b * h**2 / (12 * (1 - k_final[-1]))')
        body.append('    M_final = np.real(M) * M_cr_local')
    body.append('    return k_final, M_final')
    return '\n'.join(body)


def main():
    mdir, out = sys.argv[1], sys.argv[2]
    blocks, names = [], []
    for z in ZONES:
        p = os.path.join(mdir, f'zone{z}.m')
        name, args, k_src, M_src, per_step = parse_zone(p)
        names.append(name)
        blocks.append(emit(name, args, k_src, M_src, per_step))
        print(f'  {name:<10} args={len(args):<3} per-step de-norm={per_step}  '
              f'k {len(k_src)} chars, M {len(M_src)} chars')

    exports = ', '.join(f'"{n}"' for n in names)
    with open(out, 'w', encoding='utf-8') as fh:
        fh.write(HEADER % exports)
        fh.write('\n\n')
        fh.write('\n\n\n'.join(blocks))
        fh.write('\n')
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
