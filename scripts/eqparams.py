# Copyright 2021 Tom Eulenfeld, MIT license

import json
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
from obspy import read_events
from scipy.optimize import minimize
from scipy.stats import linregress, median_abs_deviation as mad
from qopen.imaging import _secondary_yaxis_seismic_moment, plot_all_sds
from qopen.source import moment_magnitude
from qopen.util import linear_fit

Mls = r'$M_{\rm{L}}$'
Mws = r'$M_{\rm w}$'
fcs = r'$f_{\rm c}$'


QOPENEVENTRESULTS = '../qopen/03_source/results.json'
QOPENEVENTRESULTS2 = '../qopen/04_source_nconst/results.json'
QOPENSITERESULTS = '../qopen/02_sites/results.json'

EQ_PARAMS = '../data/eq_params{}.csv'
EPHEADER = 'evid,Ml,Mw,fc,n,stressdrop'
EPDTYPE = '<U20,f4,f4,f4,f4,f4'
EPFMT = '%s,%.3f,%.3f,%.3f,%.3f,%.3f'


EVIDS8 = (  # ev ids of 8 events for detailed study
    '38445975',
    '38451079',
    '38471103',
    '38483215',
    '38450263',
    '38538991',
    '38489543',
    '38496551'
    )


VSMODEL = """0   4.96     2.93
1    5.14     3.01
2    5.45     3.14
4    6.07     3.52
8    6.12     3.62
16    6.24     3.72
32    7.12     4.01"""
VSMODEL = [list(map(float, l.split())) for l in VSMODEL.splitlines()]

def _vs(depth):
    llay = None
    for lay in VSMODEL:
        if depth > lay[0]:
            llay = lay
        else:
            return llay[2]


def print_submission_55events():
    with open(QOPENEVENTRESULTS2) as f:
            results = json.load(f)
    events = read_events('../data/cat55.csz', check_compression=False)
    for event in events:
        ori = event.origins[0]
        t = str(ori.time).replace('-', ',').replace('T', ',').replace(':', ',').strip('Z')
        evid = str(event.resource_id).split('/')[-1]
        evres = results['events'][evid]
        M0orig = evres['M0']
        dep = ori.depth/1000
        v = _vs(dep)
        M0 = M0orig * (v / 3.2) ** 2.5  # velocity correction for M0
        fc = evres['fc']
        Ml = evres['Mcat']
        sd = fc2stress_drop(fc, M0, beta=v*1000) / 1e6
        print(f'Eulenfeld,4,2700,{v*1000},,,,'
              f'{evid},{ori.latitude},{ori.longitude},{dep},{t},{Ml},'
              f'4,{M0:.0f},{sd:.3f},,{fc:.2f},,,')


def correct_results_velocity(results):
    events = read_events('../data/cat55+1.csz', check_compression=False)
    for event in events:
        ori = event.origins[0]
        evid = str(event.resource_id).split('/')[-1]
        if evid not in results['events']:
            import warnings
            warnings.warn(f'event {evid} missing in results')
            continue
        evres = results['events'][evid]
        M0orig = evres['M0']
        dep = ori.depth/1000
        v = _vs(dep)
        evres['M0'] = M0 = M0orig * (v / 3.2) ** 2.5  # velocity correction for M0
        evres['Mw'] = moment_magnitude(M0)
        evres['v1'] = v * 1000


def _linear_fit_L1(y, x, m0, b0):
    def cost_function(params, x, y):
        m, b = params
        return np.sum(np.abs(y - m * x - b))
    out = minimize(cost_function, (m0, b0), args=(x, y))
    return out.x


def _load_json(fname):
    with open(fname) as f:
        return json.load(f)


def plot_sds_in_one():
    fixn = True
    results = _load_json(QOPENEVENTRESULTS2 if fixn else QOPENEVENTRESULTS)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    freq = np.array(results['freq'])
    sdsa = [r['sds'] for r in results['events'].values()]
    ax.loglog(freq, np.transpose(sdsa), color='k', alpha=0.2)
    ax.set_xlabel('frequency (Hz)')
    ax.set_xticks([10, 20, 30, 40, 50])
    ax.set_xticklabels(['10', '20', '30', '40', ''])
    ax.set_ylabel(r'source displacement spectrum $\omega M$ (Nm)')
    fig.savefig('../figs/all_sds_in_one.pdf', bbox_inches='tight')


def calc_stressdrops(fixn=True, correct_velocity=False, only8=False):
    results = _load_json(QOPENEVENTRESULTS2 if fixn else QOPENEVENTRESULTS)
    if only8:
        results['events'] = {evid: evres for evid, evres in results['events'].items() if evid in EVIDS8}
    if correct_velocity:
        correct_results_velocity(results)
    vals = [(evid, evres['Mcat'], evres['Mw'], evres['M0'], evres['fc'], evres.get('v1', 3200)) for evid, evres in results['events'].items() if 'Mw' in evres]
    evid, Ml, Mw, M0, fc, v1 = map(np.array, zip(*vals))
    # load n values
    results = _load_json(QOPENEVENTRESULTS)  # 2018
    nd = {evid: evres['n'] for evid, evres in results['events'].items() if 'Mw' in evres}
    n = np.array([nd.get(id_) for id_ in evid], dtype=float)
    print(f'high frequency fall-of mean: {np.nanmean(n):.3f}  '
          f'median: {np.nanmedian(n):.3f}  std: {np.nanstd(n):.3f}')
    sd = fc2stress_drop(fc, M0, beta=v1) / 1e6  # MPa

    arr = np.asarray(list(zip(*[evid, Ml, Mw, fc, n, sd])), dtype=EPDTYPE)
    fname = EQ_PARAMS.format('8' * only8 + '_freen' * (not fixn) + '_cv' * correct_velocity)
    np.savetxt(fname, arr, fmt=EPFMT, header=EPHEADER)


def compare_mags():
    t1 = np.genfromtxt(EQ_PARAMS.format(''), dtype=EPDTYPE, names=True, delimiter=',')
    Ml = t1['Ml']
    Mw = t1['Mw']
    ind8 = np.isin(t1['evid'], EVIDS8)
    mmin, mmax = np.min(Ml), np.max(Ml)
    m = np.linspace(mmin, mmax, 100)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, aspect='equal')

    a, b = linear_fit(Mw, Ml)
    label = f'{Mws}={a:.2f}{Mls} {b:+.2f}'
    print(label)
    ax.plot(m, a * m + b, '-C0', label=label, alpha=0.6)

    ax.plot(Ml[~ind8], Mw[~ind8], 'x', ms=5)
    ax.plot(Ml[ind8], Mw[ind8], 'xC6', ms=5)
    ax.legend(loc='upper left', frameon=False)
    _secondary_yaxis_seismic_moment(ax)
    ax.set_ylabel(r'moment magnitude $M_{\rm w}$')
    ax.set_xlabel(r'local magnitude  $M_{\rm l}$')
    fig.savefig('../figs/mags2.pdf', bbox_inches='tight')


def fc_vs_Mw_vs_n():
    t1 = np.genfromtxt(EQ_PARAMS.format(''), dtype=EPDTYPE, names=True, delimiter=',')
    Mw, fc, n = t1['Mw'], t1['fc'], t1['n']
    fig = plt.figure(figsize=(12,5))
    ax1 = fig.add_subplot(131)
    ax1.plot(Mw, fc, 'x')
    ax1.set_xlabel(f'moment magnitude {Mws}')
    ax1.set_ylabel(f'corner frequency {fcs} (Hz)')
    ax1.set_yscale('log')
    ax2 = fig.add_subplot(132)
    ax2.hist(n, rwidth=0.9)
    ax2.set_xlabel('high frequency fall-off $n$')
    ax3 = fig.add_subplot(133)
    ax3.plot(Mw, n, 'x')
    ax3.set_xlabel(f'moment magnitude {Mws}')
    ax3.set_ylabel('high frequency fall-off $n$')
    ax1.yaxis.set_minor_formatter(ScalarFormatter())
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    plt.tight_layout()
    fig.savefig('../figs/eqparams2.pdf', bbox_inches='tight', pad_inches=0.1)


def fc2Mw(stress_drop, fc, beta=3200):
    r = beta * 0.21 / np.array(fc)  # Madariaga (1976) for S waves
    M0 = 16 * r ** 3 * stress_drop / 7
    return moment_magnitude(M0)


def fc2stress_drop(fc, M0, beta=3200):
    r = beta * 0.21 / np.array(fc)
    stress_drop = 7 * M0 / (16 * r ** 3)
    return stress_drop


def geomad(vals):
    lmad = mad(np.log(vals))
    med = np.median(vals)
    return med-med*np.exp(-lmad), med*np.exp(lmad)-med


def _logregr(Mw, fc):
    m2, b2, _, _, m2_stderr = linregress(Mw, np.log10(fc))
    m, b = 1 / m2, -b2 / m2
    m_stderr = m2_stderr / m2 ** 2
    print(f'L2 fit, Mw independent variable: M0 ∝ fc^{1.5*m:.2f}+-{1.5*m_stderr:.2f}')
    label = f'$M_0$ ∝ {fcs}$^{{{1.5*m:.2f} \pm {1.5*m_stderr:.2f}}}$'
    return m, b, m_stderr, label



def fc_vs_Mw(fixn=True, correct_velocity=False):
    dxxx = 0.26
    t1 = np.genfromtxt(EQ_PARAMS.format('_freen' * (not fixn) + '_cv' * correct_velocity),
                       dtype=EPDTYPE, names=True, delimiter=',')
    Mw1, fc1, sd1, n1 = t1['Mw'], t1['fc'], t1['stressdrop'], t1['n']
    ind8 = np.isin(t1['evid'], EVIDS8)

    fig = plt.figure(figsize=(10, 5))
    box1 = [0.45, 0.1, 0.45, 0.85]
    box4 = [0.45, 0.96, 0.45, 0.1]
    box3 = [0.05, 0.8+dxxx, 0.3, 0.23]
    box5 = [0.05, 0.45+dxxx, 0.3, 0.23]
    box6 = [0.05, 0.1, 0.3, 0.23]
    box7 = [0.05, 0.1+dxxx, 0.3, 0.23]
    akwargs = dict(xy=(0, 1), xytext=(5, -5),
                   xycoords='axes fraction', textcoords='offset points',
                   va='top', size='large')

    # main axes f) Mw vs fc
    ax1 = fig.add_axes(box1)
    ax1.plot(fc1[~ind8], Mw1[~ind8], 'x', color='C0')
    ax1.plot(fc1[ind8], Mw1[ind8], 'x', color='C6')
    m, b, _, _, m_stderr = linregress(np.log10(fc1), Mw1)
    print(f'L2 fit, fc independent variable: M0 ∝ fc^{1.5*m:.2f}+-{1.5*m_stderr:.2f}')
    m, b = _linear_fit_L1(Mw1, np.log10(fc1), m, b)
    print(f'L1 fit, fc independent variable: M0 ∝ fc^{1.5*m:.2f}')
    m2, b2 = _linear_fit_L1(np.log10(fc1), Mw1, m, b)
    m, b = 1 / m2, -b2 / m2
    print(f'L1 fit, Mw independent variable: M0 ∝ fc^{1.5*m:.2f}')

    kw = dict(rotation=-25, rotation_mode='anchor',
              xytext=(0, -15), textcoords='offset points')
    ax1.annotate('1 MPa', (1.7, fc2Mw(1e6, 1.7)), **kw)
    ax1.annotate('10 MPa', (1.7, fc2Mw(10e6, 1.7)), **kw)
    ax1.annotate('100 MPa', (7, fc2Mw(100e6, 7)), **kw)
    fcx = 10 ** (moment_magnitude(1e13) / m - b / m)
    sdx = fc2stress_drop(fcx, 1e13) / 1e6
    print(f'fc={fcx:.1f}Hz and stress drop {sdx:.1f} Mpa for M0=1e13 Nm')
    ax1.set_ylabel(f'moment magnitude {Mws}')
    ax1.set_xlabel(f'corner frequency {fcs} (Hz)')
    ax1.set_xscale('log')
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    ax1.xaxis.set_minor_formatter(ScalarFormatter())
    ax1.tick_params(top=False)
    _secondary_yaxis_seismic_moment(ax1)

    # b) corner frequency histogram
    ax4 = fig.add_axes(box4, sharex=ax1)
    bins = np.logspace(np.log10(1.7), np.log10(10), 31)
    ax4.hist(fc1, bins=bins, rwidth=0.9, zorder=10)
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.tick_params(bottom=False, labelbottom=False, right=False, top=False)
    ax4.axvline(np.median(fc1), color='C0', ls='--')

    # a) n histogram
    ax3 = fig.add_axes(box3)
    bins = np.arange(1.475, 3.525, 0.05)
    ax3.hist(n1, bins=bins, rwidth=0.9, zorder=10)
    ax3.axvline(np.median(n1), color='C0', ls='--')

    ax3.set_ylabel('counts')
    ax3.set_xlabel('high frequency fall-off $n$')
    ax1.annotate('f)', **akwargs)
    ax3.annotate('a)', **akwargs)
    if ax4:
        ax4.set_ylabel('counts')
        ax4.annotate('b)', **akwargs)

    # c) stress drop histogram
    bins = 10 ** np.arange(-0., 2.05, 0.1)
    ax5 = fig.add_axes(box5)
    ax5.axvline(np.median(sd1), color='C0', ls='--')
    ax5.hist(sd1, bins=bins, rwidth=0.9, zorder=10)

    ax5.set_xscale('log')
    ax5.set_xlabel('stress drop (MPa)')
    ax5.set_ylabel('counts')
    ax5.xaxis.set_major_formatter(ScalarFormatter())
    ax5.annotate('c)', **akwargs)

    # e) stress drop vs binned Mw
    ax6 = fig.add_axes(box6)
    Mwmid = np.arange(2.4, 4.6, 0.2)
    sdv1 = [sd1[np.logical_and(mwmidv-0.1 <= Mw1, Mw1 < mwmidv+0.1)]
             for mwmidv in  Mwmid]
    sdm1 = np.array(list(map(np.median, sdv1)))
    sderr1 = np.transpose(list(map(geomad, sdv1)))
    ax6.errorbar(Mwmid-0.01, sdm1, yerr=sderr1, marker='x', ls='', label='2018')
    ax6.set_yscale('log')
    ax6.yaxis.set_major_formatter(ScalarFormatter())
    # ax6.set_ylim(1/1.1, 100*1.1)

    ax6.set_ylabel('stress drop\n(MPa)')
    ax6.set_xlabel(f'moment magnitude {Mws}')
    ax6.annotate('e)', **akwargs)
    del akwargs['size']
    ax6.annotate('         median and MAD', **akwargs)
    # ax6.set_xticks(Mwmid)

    @matplotlib.ticker.FuncFormatter
    def myformatter(x, pos):
        # print(x)
        return f'{x:.1f}' if 0.25 <= x <= 0.35 else ''
    ax6.yaxis.set_minor_formatter(myformatter)

    # d) fc vs binned Mw
    fcv1 = [fc1[np.logical_and(mwmidv-0.1 <= Mw1, Mw1 < mwmidv+0.1)]
            for mwmidv in  Mwmid]
    fcm1 = np.array(list(map(np.median, fcv1)))
    fcerr1 = np.transpose(list(map(geomad, fcv1)))
    ax7 = fig.add_axes(box7, sharex=ax6)
    ax7.errorbar(Mwmid-0.01, fcm1, yerr=fcerr1, marker='x', ls='')
    ax7.set_yscale('log')
    ax7.set_yticks([2, 3, 4, 5, 7, 10])
    ax7.set_yticks([], minor=True)
    ax7.yaxis.set_major_formatter(ScalarFormatter())
    ax7.tick_params(labelbottom=False)
    ax7.set_ylabel('corner frequency\n(Hz)')
    ax7.annotate('d)', **akwargs)
    ax7.annotate('         median and MAD', **akwargs)
    # ax7.set_xlim(0.35, 2.05)

    # line plots for ax1 and ax7
    fclim = np.array([1.5, 5, 10])
    m, b, m_stderr, label = _logregr(Mw1, fc1)
    print('m, b = ', m, b)
    ax1.plot(fclim, np.log10(fclim)*m+b, zorder=-1, label=label, color='C0')

    mwlim = np.array((2, 5))
    m, b, m_stderr, label = _logregr(Mwmid, fcm1)
    ax7.plot(mwlim, 10**((mwlim-b)/m), '-.', zorder=-1, label=label, color='C0', alpha=0.5)
    ax1.plot(10**((mwlim-b)/m), mwlim, '-.', zorder=-1, label=label, color='C0', alpha=0.5)

    kw = dict(ls='--', color='0.5', zorder=-1)
    ax1.plot(fclim, fc2Mw(1e6, fclim), **kw)
    ax1.plot(fclim, fc2Mw(10e6, fclim), label=f'$M_0$ ∝ {fcs}$^{{-3}}$', **kw)
    ax1.plot(fclim, fc2Mw(100e6, fclim), **kw)
    ax1.legend()

    fname = '../figs/eqparams{}.pdf'.format('_freen' * (not fixn) + '_cv' * correct_velocity)
    fig.savefig(fname, bbox_inches='tight', pad_inches=0.1)


def plot_sds_8events(fixn=True):
    results = _load_json(QOPENEVENTRESULTS2 if fixn else QOPENEVENTRESULTS)
    plot_all_sds(results, plot_only_ids=EVIDS8, nx=5, figsize=(7, 4), fname='../figs/sds8.pdf')


if __name__ == '__main__':
    plot_sds_in_one()
    calc_stressdrops()
    calc_stressdrops(correct_velocity=True)
    calc_stressdrops(fixn=False)
    calc_stressdrops(fixn=False, correct_velocity=True)
    calc_stressdrops(only8=True)
    calc_stressdrops(only8=True, correct_velocity=True)
    calc_stressdrops(only8=True, fixn=False)
    calc_stressdrops(only8=True, fixn=False, correct_velocity=True)
    compare_mags()
    fc_vs_Mw_vs_n()
    fc_vs_Mw()
    fc_vs_Mw(correct_velocity=True)
    # fc_vs_Mw(fixn=False)
    # fc_vs_Mw(fixn=False, correct_velocity=True)
    plot_sds_8events()
    print_submission_55events()
