#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 14:03:11 2020

@author: Javier Gorroño
"""

import seaborn as sn
import pandas as pd
import numpy as np
import os
from scipy import stats

import matplotlib

matplotlib.use('Agg')  # this does not show the plot on the screen
import matplotlib.pyplot as pt
import matplotlib.colors as cl
import matplotlib.dates as mdates

RESULTS_FOLDER = os.path.join(os.getcwd(), 'L2Aunc_results')


def plot_spectralcorrelation(values, cbarlabel, variable, tagfile):
    res = sn.heatmap(np.round(values, 2), annot=True)
    cbar = res.collections[0].colorbar
    cbar.set_label(label=r'$r$', rotation=270, fontsize=18, labelpad=25)
    cbar.mappable.set_clim(vmin=0, vmax=1)
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_fontsize(16)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize=16)
    res.tick_params(axis='x', labelrotation=90)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize=16)
    res.tick_params(axis='y', labelrotation=0)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_speccorr' + variable + '.png'))
    pt.close()


def plot_functioncorrelation(values, bandname, corrname, tagfile):
    res = sn.heatmap(np.round(values, 2), annot=True)
    cbar = res.collections[0].colorbar
    cbar.set_label(label=r'$r$', rotation=270, fontsize=18, labelpad=25)
    cbar.mappable.set_clim(vmin=-1, vmax=1)
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_fontsize(16)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize=16)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize=16)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_' + corrname + '_' + bandname + '.png'))
    pt.close()


def plot_etoacorrelation(r_etoa, cbarlabel, tagfile, bandnames):
    # r(E_toa,atm. function) plot for all bands at once. Select r_etoa[0, 1:] (first row and skip Etoa vs Etoa)
    df = pd.DataFrame(r_etoa, columns=[r'$L_{path}$', r'$\tau_{up}$', r'$E_{g(dir)}$', r'$E_{g(diff)}$', r'$s_{atm}$'],
                      index=[r'$E_{toa}' + bandname + '$' for bandname in bandnames])
    res = sn.heatmap(np.round(df, 2), annot=True)
    cbar = res.collections[0].colorbar
    cbar.set_label(label=r'$r$[' + cbarlabel + ']', rotation=270, fontsize=18, labelpad=25)
    cbar.mappable.set_clim(vmin=0, vmax=1)
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_fontsize(16)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize=16)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize=16)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_speccorr_etoa.png'))
    pt.close()


def plot_histograms(values, bandname, cbarlabel, variable, tagfile):
    pt.hist(values, density=True, bins='auto',
            label=r'$\sigma = ' + str(np.round(100 * np.std(values) / np.mean(values), 2)) + '\%$\n$\mu = ' + str(
                np.round(np.mean(values), 3)) + '$')
    x_data = np.linspace(np.min(values), np.max(values), num=1000)
    y_data = stats.norm.pdf(x_data, np.mean(values), np.std(values))
    pt.plot(x_data, y_data, linewidth=3, color='red', label='$Normal fit$')
    ax = pt.gca()
    ax.xaxis.set_tick_params(labelsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    pt.xlabel(cbarlabel, fontsize=20)
    pt.ylabel('Normalised probability', fontsize=20)
    # pt.legend(loc='best', fontsize=14)
    print(r'$\sigma = ' + str(np.round(100 * np.std(values) / np.mean(values), 2)) + '\%$\n$\mu = ' + str(
                np.round(np.mean(values), 3)) + '$' + variable + '_' + bandname)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_hist_' + variable + '_' + bandname + '.png'))
    pt.close()


def plot_histograms_unctheory(values, bandname, cbarlabel, variable, tagfile, l2a_unc):
    pt.hist(values, density=True, bins='auto',
            label=r'$\sigma = ' + str(np.round(100 * np.std(values) / np.mean(values), 2)) + '\%$\n$\mu = ' + str(
                np.round(np.mean(values), 3)) + '$')
    x_data = np.linspace(np.min(values), np.max(values), num=1000)
    ax = pt.gca()
    if variable == 'L2Arho':
        y_data = stats.norm.pdf(x_data, np.mean(values), l2a_unc)
        pt.plot(x_data, y_data, linewidth=3, color='red',
                label=r'$U_{LPU}= ' + str(np.round(100 * l2a_unc / np.mean(values), 2)) + '\%$')
        print(r'$\sigma = ' + str(np.round(100 * np.std(values) / np.mean(values), 2)) + '\%$\n$\mu = ' + str(
                np.round(np.mean(values), 3)) + '$\n'+ r'$U_{LPU}= ' + str(np.round(100 * l2a_unc / np.mean(values), 2)) + '\%$' + variable + '_' + bandname)
    else:
        y_data = stats.norm.pdf(x_data, np.mean(values), np.std(values))
        pt.plot(x_data, y_data, linewidth=3, color='red', label='$Normal fit$')
        print(r'$\sigma = ' + str(np.round(100 * np.std(values) / np.mean(values), 2)) + '\%$\n$\mu = ' + str(
            np.round(np.mean(values), 3)) + '$' + variable + '_' + bandname)

    ax.xaxis.set_tick_params(labelsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    pt.xlabel(cbarlabel, fontsize=20)
    pt.ylabel('Normalised probability', fontsize=20)
    # pt.legend(loc='best', fontsize=18)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_hist_' + variable + '_' + bandname + '.png'))
    pt.close()

def plot_validationmap(bandname, aots, wvs, gum_unc, mcm_std, mcm_unc, gum_unc_percent, mcm_std_percent,
                       mcm_unc_percent, tagfile):
    vals = np.clip(gum_unc_percent.reshape((len(aots), len(wvs))),1.08,2)
    norm = pt.Normalize(vmin=1.1, vmax=2)
    h = pt.contourf(wvs, aots, vals, extend='max',norm=norm)
    cbar = pt.colorbar()
    cbar.mappable.set_clim(vmin=1.1, vmax=2)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('$U_\mathregular{LPU}$ [%]', rotation=270, fontsize=18)
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_fontsize(16)
    pt.xlabel('WV[cm]', fontsize=18)
    pt.ylabel('AOT', fontsize=18)
    ax = pt.gca()
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_gumpercentmap_' + bandname + '.png'), dpi=300)
    pt.close()

    vals = np.clip(mcm_unc_percent.reshape((len(aots), len(wvs))),1.21,1.92)
    norm = pt.Normalize(vmin=1.2, vmax=2)
    h = pt.contourf(wvs, aots, vals, extend='max',norm=norm)
    cbar = pt.colorbar()
    cbar.mappable.set_clim(vmin=1.2, vmax=2)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('$U_\mathregular{MCM}$ [%]', rotation=270, fontsize=18)
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_fontsize(16)
    pt.xlabel('WV[cm]', fontsize=18)
    pt.ylabel('AOT', fontsize=18)
    ax = pt.gca()
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_mcmuncpercentmap_' + bandname + '.png'), dpi=300)
    pt.close()

    vals = np.clip(mcm_std_percent.reshape((len(aots), len(wvs))),1.21,1.98)
    norm = pt.Normalize(vmin=1.2, vmax=2)
    h = pt.contourf(wvs, aots, vals, extend='max',norm=norm)
    cbar = pt.colorbar()
    cbar.mappable.set_clim(vmin=1.2, vmax=2)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('$\sigma_\mathregular{MCM}$ [%]', rotation=270, fontsize=18)
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_fontsize(16)
    pt.xlabel('WV[cm]', fontsize=18)
    pt.ylabel('AOT', fontsize=18)
    ax = pt.gca()
    ax.xaxis.set_tick_params(labelsize=18)
    ax.yaxis.set_tick_params(labelsize=18)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_mcmstdpercentmap_' + bandname + '.png'), dpi=300)
    pt.close()