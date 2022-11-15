# Copyright 2020 Tom Eulenfeld, MIT license

import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter
from qopen.core import collect_results
from qopen.imaging import calc_dependent
from qopen.util import gerr


fin = '../qopen/01_go/results.json'


def Q(q, obs='sc', v0=None, error=False):
    if v0 is None:
        v0 = q['config']['v0']
    freq = np.array(q['freq'])
    if obs == 'sc':
        mean = np.array(q['g0'], dtype=float) * v0 / (2 * np.pi * freq)
    else:
        mean = np.array(q['b'], dtype=float) / (2 * np.pi * freq)
    if not error:
        return freq, mean
    q = collect_results(q)
    if obs == 'sc':
        vals = np.array(q['g0']) * v0 / (2 * np.pi * freq)
    else:
        vals = np.array(q['b']) / (2 * np.pi * freq)
    mean2, err1, err2 = gerr(vals, axis=0, robust=True)
    np.testing.assert_allclose(mean, mean2)
    return freq, mean, (err1, err2)


def printQ_json(freq, Qsc, Qi, label=None):
    freqx = np.round(freq, 3).tolist()
    Qscx = np.round(Qsc*1e3, 3).tolist()
    Qix = np.round(Qi*1e3, 3).tolist()
    print(f""""{label}": {{
              "f": {freqx},
              "Qsc": {Qscx},
              "Qi": {Qix}}}""")


def plotQ():
    with open(fin) as f:
        q1 = json.load(f)
    fig = plt.figure(figsize=(4, 4))
    ax1 = fig.add_subplot(111)

    freq, Qsc = Q(q1, 'sc')
    _, Qi = Q(q1, 'i')

    printQ_json(freq, Qsc, Qi, label='Qopen_Ridgecrest')
    ax1.errorbar(*Q(q1, 'sc', error=True), color='C9', marker='.', ms=6, zorder=5, label='scattering')
    ax1.errorbar(*Q(q1, 'i', error=True), color='C6', marker='.', ms=6, zorder=5, label='intrinsic')
    ax1.plot(freq, Qsc+Qi, color='0.5', marker='.', ms=6, zorder=5, label='total')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('frequency (Hz)')
    ax1.set_ylabel(r'$Q^{-1}$')
    ax1.legend()
    ax1.set_xticks((1, 10, 100))
    ax1.set_xticklabels(('1', '10', '100'))
    fig.savefig('../figs/Q.pdf')


def plotl(ax=None, marker=None):
    with open(fin) as f:
        result = json.load(f)
    fig = None
    if ax is None:
        fig = plt.figure(figsize=(5, 3))
        ax = fig.add_subplot(111)
    kw = {'marker': marker}
    if marker == 'x':
        kw['mew'] = 1
        kw['ms'] = 4
    freq = np.array(result['freq'])
    v0 = result['config']['v0']
    col = collect_results(result, only=('g0', 'b', 'error'))
    lsc_all = calc_dependent('lsc', col['g0'], freq=freq, v0=v0)
    li_all = calc_dependent('li', col['b'], freq=freq, v0=v0)
    lsc, err1, err2 = gerr(lsc_all, axis=0, robust=True)
    lscerrs = (err1, err2)
    li, err1, err2 = gerr(li_all, axis=0, robust=True)
    lierrs = (err1, err2)
    ax.errorbar(freq*1.02, lsc, color='C9', ls='-', yerr=lscerrs, label='transport\nmean free path', **kw)
    ax.errorbar(freq/1.02, li, color='C6', ls='-', yerr=lierrs, label='intrinsic\nabsorption length', **kw)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_formatter(ScalarFormatter())

    ax.set_xlabel('frequency (Hz)')
    ax.set_ylabel('path length (km)')
    ax.legend()
    if fig is not None:
        fig.savefig('../figs/l.pdf', bbox_inches='tight', dpi=600)


if __name__ == '__main__':
    # plotQ()
    plotl()
