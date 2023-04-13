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
import matplotlib.dates as mdates

RESULTS_FOLDER = os.path.join(os.getcwd(), 'L2Aunc_results')


def plot_spectralcorrelation(values, cbarlabel, variable, tagfile):
    res = sn.heatmap(np.round(values, 2), annot=True)
    cbar = res.collections[0].colorbar
    cbar.set_label(label=r'$r^{2}$[' + cbarlabel + ']', rotation=270, fontsize=14, labelpad=25)
    cbar.mappable.set_clim(vmin=0, vmax=1)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize=12)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize=12)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_speccorr' + variable + '.png'))
    pt.close()


def plot_functioncorrelation(values, bandname, corrname, tagfile):
    res = sn.heatmap(values, annot=True)
    cbar = res.collections[0].colorbar
    cbar.set_label(label=r'$r^{2}$', rotation=270, fontsize=14, labelpad=25)
    cbar.mappable.set_clim(vmin=-1, vmax=1)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize=15)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize=15)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_' + corrname + '_' + bandname + '.png'))
    pt.close()


def plot_etoacorrelation(r_etoa, cbarlabel, tagfile, bandnames):
    # r2(E_toa,atm. function) plot for all bands at once. Select r_etoa[0, 1:] (first row and skip Etoa vs Etoa)
    df = pd.DataFrame(r_etoa, columns=[r'$L_{path}$', r'$\tau_{up}$', r'$E_{g(dir)}$', r'$E_{g(diff)}$', r'$s_{atm}$'],
                      index=[r'$E_{toa}' + bandname + '$' for bandname in bandnames])
    res = sn.heatmap(np.round(df, 2), annot=True)
    cbar = res.collections[0].colorbar
    cbar.set_label(label=r'$r^{2}$[' + cbarlabel + ']', rotation=270, fontsize=18, labelpad=25)
    cbar.mappable.set_clim(vmin=0, vmax=1)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize=12)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize=12)
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
    ax.xaxis.set_tick_params(labelsize=12)
    pt.xlabel(cbarlabel, fontsize=16)
    pt.ylabel('Normalised probability', fontsize=14)
    pt.legend(loc='best', fontsize=14)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_hist_' + variable + '_' + bandname + '.png'))
    pt.close()


def plot_histograms_unctheory(values, bandname, cbarlabel, variable, tagfile, l2a_unc):
    pt.hist(values, density=True, bins='auto',
            label=r'$\sigma = ' + str(np.round(100 * np.std(values) / np.mean(values), 2)) + '\%$\n$\mu = ' + str(
                np.round(np.mean(values), 3)) + '$')
    x_data = np.linspace(np.min(values), np.max(values), num=1000)
    if variable == 'L2Arho':
        y_data = stats.norm.pdf(x_data, np.mean(values), l2a_unc)
        pt.plot(x_data, y_data, linewidth=3, color='red',
                label=r'$U_{L2A}= ' + str(np.round(100 * l2a_unc / np.mean(values), 2)) + '\%$')
    else:
        y_data = stats.norm.pdf(x_data, np.mean(values), np.std(values))
        pt.plot(x_data, y_data, linewidth=3, color='red', label='$Normal fit$')
    ax = pt.gca()
    ax.xaxis.set_tick_params(labelsize=12)
    pt.xlabel(cbarlabel, fontsize=16)
    pt.ylabel('Normalised probability', fontsize=14)
    pt.legend(loc='best', fontsize=14)
    pt.tight_layout()
    pt.savefig(os.path.join(RESULTS_FOLDER, tagfile + '_hist_' + variable + '_' + bandname + '.png'))
    pt.close()


def evi_ndvi_trenduncertainty(ndvi_samp, evi_samp, meas_dates, ndvi_sampcorr, evi_sampcorr, ndvi_sampuncorr,
                              evi_sampuncorr,filename):
    fig, axs = pt.subplots(3, 1,sharex=True)
    axs[0].plot(meas_dates, np.mean(np.array(evi_samp), axis=1), color='#089FFF')
    axs[0].fill_between(meas_dates, np.mean(np.array(evi_samp), axis=1) - np.std(np.array(evi_samp), axis=1),
                    np.mean(np.array(evi_samp), axis=1) + np.std(np.array(evi_samp), axis=1),
                    alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF', linewidth=0)
    axs[0].errorbar(meas_dates, np.mean(np.array(evi_samp), axis=1), yerr=np.std(np.array(evi_samp),axis=1), color='#089FFF',
                ecolor='#089FFF', capsize=6,markeredgewidth=1,label='EVI')
    axs[0].plot(meas_dates, np.mean(np.array(ndvi_samp), axis=1), color='#FF9848')
    axs[0].fill_between(meas_dates, np.mean(np.array(ndvi_samp), axis=1) - np.std(np.array(ndvi_samp), axis=1),
                    np.mean(np.array(ndvi_samp), axis=1) + np.std(np.array(ndvi_samp), axis=1),
                    alpha=0.2, edgecolor='#CC4F1B', facecolor='#FF9848', linewidth=0)
    axs[0].errorbar(meas_dates, np.mean(np.array(ndvi_samp), axis=1), yerr=np.std(np.array(ndvi_samp), axis=1),
                color='#FF9848', ecolor='#FF9848', capsize=6,markeredgewidth=1,label='NDVI')

    axs[1].plot(meas_dates, 100 * np.std(np.array(ndvi_samp), axis=1) / np.mean(np.array(ndvi_samp), axis=1),
                    color='#FF9848', marker = 'o', label='estimated')
    axs[1].plot(meas_dates, 100 * np.std(np.array(ndvi_sampcorr), axis=1) / np.mean(np.array(ndvi_sampcorr), axis=1),
                    color='#FF9848', marker = '^', label='correlated')
    axs[1].plot(meas_dates, 100 * np.std(np.array(ndvi_sampuncorr), axis=1) / np.mean(np.array(ndvi_sampuncorr), axis=1),
                    color='#FF9848', marker = 'v', label='uncorrelated')

    axs[2].plot(meas_dates, 100 * np.std(np.array(evi_samp), axis=1) / np.mean(np.array(evi_samp), axis=1),
                    color='#089FFF', marker = 'o', label='estimated')
    axs[2].plot(meas_dates, 100 * np.std(np.array(evi_sampcorr), axis=1) / np.mean(np.array(evi_sampcorr), axis=1),
                    color='#089FFF', marker = '^', label='correlated')
    axs[2].plot(meas_dates, 100 * np.std(np.array(evi_sampuncorr), axis=1) / np.mean(np.array(evi_sampuncorr), axis=1),
                    color='#089FFF', marker = 'v', label='uncorrelated')

    ratio = 0.5
    x_left, x_right = axs[0].get_xlim()
    y_low, y_high = axs[0].get_ylim()
    axs[0].set_aspect(abs((x_right - x_left) / (y_low - y_high)) * ratio)
    x_left, x_right = axs[1].get_xlim()
    y_low, y_high = axs[1].get_ylim()
    axs[1].set_aspect(abs((x_right - x_left) / (y_low - y_high)) * ratio)
    x_left, x_right = axs[2].get_xlim()
    y_low, y_high = axs[2].get_ylim()
    axs[2].set_aspect(abs((x_right - x_left) / (y_low - y_high)) * ratio)

    axs[0].yaxis.set_label_text('NDVI/EVI')
    axs[2].xaxis.set_label_text('Month-year')
    axs[1].yaxis.set_label_text('NDVI unc.[%]')
    axs[2].yaxis.set_label_text('EVI unc.[%]')
    myFmt = mdates.DateFormatter('%m-%y')
    axs[0].xaxis.set_major_formatter(myFmt)
    axs[0].xaxis.set_major_locator(mdates.MonthLocator())
    for tick in axs[2].xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    for tick in axs[0].yaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    for tick in axs[1].yaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    for tick in axs[2].yaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    # axs[0].grid(True, axis='both')
    axs[0].legend(loc='best', fontsize=6)
    axs[1].legend(loc='best', fontsize=6)
    axs[2].legend(loc='best', fontsize=6)
    pt.savefig(os.path.join(RESULTS_FOLDER, filename + '.png'),dpi=300)
    pt.close()
