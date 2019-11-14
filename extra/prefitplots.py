from __future__ import division
import matplotlib.pyplot as plt
import math,os,sys,random,settings
import ROOT, os, sys
from ROOT import gROOT, TPad, TProfile, TCanvas, gStyle, TH2D
from sklearn import svm
import numpy as np
#import random
from random import sample 
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from root_numpy import fill_hist
ROOT.gStyle.SetOptStat(0) #disables stat boxes

def prefit_plots(pt_array, ntrk_array, weight, gluon_pt_array, gluon_ntrk_array, gluon_weight_array, quark_pt_array, quark_ntrk_array, quark_weight_array):
    #diagnostic plots of quarks and gluons from bkgd+signal
    print 'total events: ' + str(len(pt_array))
    print 'gluons: ' + str(len(gluon_pt_array))
    print 'quarks: ' + str(len(quark_pt_array))
    bkg_ntrk_pt = plt.figure()
    plt.scatter(pt_array, ntrk_array)
    plt.title('prefit bkg ntrk vs pt')
    bkg_ntrk_pt.savefig(settings.outputdir + '/bkg_ntrk_pt.pdf')

    ntrk_pt_hist = plt.figure()
    plt.hist2d(pt_array, ntrk_array, bins=(50,max(ntrk_array)))
    plt.title('prefit bkg ntrk vs pt')
    ntrk_pt_hist.savefig(settings.outputdir + '/ntrk_pt_hist.pdf')

    gluon_ntrk_pt_hist = plt.figure()
    plt.hist2d(gluon_pt_array, gluon_ntrk_array, bins=(50,max(gluon_ntrk_array)))
    plt.title('prefit gluon ntrk vs pt')
    gluon_ntrk_pt_hist.savefig(settings.outputdir + '/gluon_ntrk_pt_hist.pdf')

    quark_ntrk_pt_hist = plt.figure()
    plt.title('prefit quark ntrk vs pt')
    plt.hist2d(quark_pt_array, quark_ntrk_array, bins=(50,max(quark_ntrk_array)))
    quark_ntrk_pt_hist.savefig(settings.outputdir + '/quark_ntrk_pt_hist.pdf')

def tag_eff(a, b, quark_pt_array, quark_ntrk_array, gluon_pt_array, gluon_ntrk_array):
    pass_quarks_pt = [x for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= a*x + b]
    pass_quarks_ntrk = [y for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= a*x + b]

    fail_gluons_pt = [x for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > a*x + b]
    fail_gluons_ntrk = [y for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > a*x + b]
    '''
    print 'nQuarks: ' + str(len(quark_pt_array))
    print 'nQuarks correctly classified: ' + str(len(pass_quarks_pt))
    '''
    q_tag_eff = len(quark_pt_array) > 0 and float(len(pass_quarks_pt)/len(quark_pt_array)) or 0
    #print 'q_tag_eff: ' + str(q_tag_eff)

    '''
    print 'nGluons: ' + str(len(gluon_pt_array))
    print 'nGluons correctly classified: ' + str(len(fail_gluons_pt))
    '''
    g_reject = len(gluon_pt_array) > 0 and float(len(fail_gluons_pt)/len(gluon_pt_array)) or 0
    #print 'g_rejection: ' + str(g_reject)
    return q_tag_eff, g_reject



def profile_plots(pt_array, ntrk_array, weight, gluon_pt_array, gluon_ntrk_array, gluon_weight_array, quark_pt_array, quark_ntrk_array, quark_weight_array, islogpt):
    #gStyle.SetOptFit()
    print 'total events: ' + str(len(pt_array))
    print 'gluons: ' + str(len(gluon_pt_array))
    print 'quarks: ' + str(len(quark_pt_array))
    c1 = TCanvas( 'c1', '', 200, 10, 700, 500 )
    fittype = 'pol1'
    if islogpt:
        ptmin = 4
        ptmax = 7
        nbins = 10
    else:
        ptmin = 0
        ptmax = 1300
        nbins = 10
    hquark  = TProfile( 'hquark', 'Quark Profile of nTrk vs ln(pt)', nbins, ptmin, ptmax, ntrkmin, ntrkmax)
    weighted_hquark  = TProfile( 'weighted_hquark', 'Quark Weighted Profile of nTrk vs ln(pt)', nbins, ptmin, ptmax, ntrkmin, ntrkmax)
    count = 0
    for i, j, k in zip(quark_pt_array, quark_ntrk_array, quark_weight_array):
        count +=1
        if count%100000==0:
            print 'processed ' + str(count)
        hquark.Fill(i,j)
        weighted_hquark.Fill(i,j,k)

    draw_hist_with_tag_eff(hquark, "quark", quark_pt_array, quark_ntrk_array, gluon_pt_array, gluon_ntrk_array)
    draw_hist_with_tag_eff(weighted_hquark, "weightedquark", quark_pt_array, quark_ntrk_array, gluon_pt_array, gluon_ntrk_array)


    hgluon  = TProfile( 'hgluon', 'Gluon Profile of nTrk vs ln(pt)', nbins, ptmin, ptmax, ntrkmin, ntrkmax)
    weighted_hgluon  = TProfile( 'weighted_hgluon', 'Gluon Weighted Profile of nTrk vs ln(pt)', nbins, ptmin, ptmax, ntrkmin, ntrkmax)
    for i, j, k in zip(gluon_pt_array, gluon_ntrk_array, gluon_weight_array):
        hgluon.Fill(i,j)
        weighted_hgluon.Fill(i,j,k)

    draw_hist_with_tag_eff(hgluon, "gluon", quark_pt_array, quark_ntrk_array, gluon_pt_array, gluon_ntrk_array)
    draw_hist_with_tag_eff(weighted_hgluon, "weightedgluon", quark_pt_array, quark_ntrk_array, gluon_pt_array, gluon_ntrk_array)

    c3 = TCanvas( 'c3', '', 200, 10, 700, 500 )
    weighted_hquark.SetMaximum(30)
    weighted_hquark.Draw("HIST")
    weighted_hgluon.Draw("HIST SAME")
    weighted_hquark.SetLineColor(2)

    c3.SaveAs(settings.outputdir + '/overlay.pdf')

def get_ntrk_cut(ptmin, ptmax, quark_pt_array, quark_ntrk_array, quark_weight_array, gluon_pt_array, gluon_ntrk_array, gluon_weight_array):
    c4 = TCanvas( 'c4', '', 200, 10, 700, 500 )
    b_quark_pt_array = [x for x in quark_pt_array if x >= ptmin and x < ptmax]
    b_quark_ntrk_array = [y for x,y in zip(quark_pt_array, quark_ntrk_array) if x >= ptmin and x < ptmax]
    b_quark_weight_array = [y for x,y in zip(quark_pt_array, quark_weight_array) if x >= ptmin and x < ptmax]
    b_gluon_pt_array = [x for x in gluon_pt_array if x >= ptmin and x < ptmax]
    b_gluon_ntrk_array = [y for x,y in zip(gluon_pt_array, gluon_ntrk_array) if x >= ptmin and x < ptmax]
    b_gluon_weight_array = [y for x,y in zip(gluon_pt_array, gluon_weight_array) if x >= ptmin and x < ptmax]

    qhist = ROOT.TH1F("qhist","",ntrkbins, ntrkmin, ntrkmax)
    ghist = ROOT.TH1F("ghist","",ntrkbins, ntrkmin, ntrkmax)
    for i,j in zip(b_quark_ntrk_array, b_quark_weight_array):
        qhist.Fill(i,j)
    for i,j in zip(b_gluon_ntrk_array,b_gluon_weight_array):
        ghist.Fill(i,j)
    ntrkcuts = range(10,22)
    print 'pT Range: ' + str(ptmin) + ' ' + str(ptmax)
    for i in ntrkcuts:
        qeff, geff = tag_eff(0, i, b_quark_pt_array, b_quark_ntrk_array, b_gluon_pt_array, b_gluon_ntrk_array)
        print 'nTrk: ' + str(i) + ' qeff: ' + str(qeff) + ' geff: ' + str(geff)

    histmax = max(qhist.GetMaximum(), ghist.GetMaximum())
    qhist.SetMaximum(1.1*histmax)
    qhist.GetXaxis().SetTitle("nTrk")
    qhist.SetTitle("nTrk for " + str(ptmin) + " < ln(pt) < " + str(ptmax))
    qhist.Draw("HIST")
    ghist.SetLineColor(2)
    ghist.Draw("HIST same")
    c4.Print(settings.outputdir + '/ntrk' + str(ptmin) + '_' + str(ptmax) + '.pdf')

def draw_hist_with_tag_eff(hist, histname, quark_pt_array, quark_ntrk_array, gluon_pt_array, gluon_ntrk_array, fittype='pol1'):
    c2 = TCanvas( 'c2', '', 200, 10, 700, 500 )
    hist.Draw()
    myfit = hist.Fit(fittype, 'S')
    inter = myfit.Parameter(0)
    slope = myfit.Parameter(1)
    tinter = '%.3f'%inter
    tslope = '%.3f'%slope
    quark_accept, gluon_rejection = tag_eff(slope, inter, quark_pt_array, quark_ntrk_array, gluon_pt_array, gluon_ntrk_array)
    qaccept = '%.3f'%quark_accept
    greject = '%.3f'%gluon_rejection
    mystring = "Best Fit (nTrk = " + str(tslope) + "x + " + str(tinter)  + "), qeff=" + str(qaccept) + " geff=" + str(greject)
    f2 = ROOT.TF1("f2", "[0] + [1]*x", 0, 10)
    f2.SetParameter(0, inter)
    f2.SetParameter(1, slope)
    leg = ROOT.TLegend(0.13,0.75,0.3,0.9)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(f2,mystring, "l")
    leg.Draw()
    c2.SaveAs(settings.outputdir + '/' + histname + ".pdf")
