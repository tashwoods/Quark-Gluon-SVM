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
ROOT.gROOT.SetBatch(True)

def main():
    #create output files
    os.system('rm -r ' + settings.outputdir)
    os.system('mkdir ' + settings.outputdir)
    myfile = open(settings.outputdir + '/linearpt_slope_intercept.txt','w')
    logmyfile = open(settings.outputdir + '/linearlogpt_slope_intercept.txt','w')
    rocmyfile = open(settings.outputdir + '/linearroc.txt','w')
    roclogmyfile = open(settings.outputdir + '/linearlogroc.txt','w')

    #load npz file
    npzfile=np.load(settings.bkgnpzfilelocation)
    classification=npzfile['classification']
    samples=npzfile['samples']

    #input file processing
    if settings.nevents!=-1:
      indices = np.random.choice(len(samples), settings.nevents, replace=False)
      random_samples = list()
      random_classification = list()
      for i in indices:
          random_samples.append(samples[i])
          random_classification.append(classification[i])
      samples = random_samples
      classification = random_classification

    if settings.dologpt:
        log_pt = [math.log(s[0]) for s in samples]
        ntrk = [s[1] for s in samples]
        samples = zip(log_pt, ntrk)
    weight=npzfile['weight']

    #standardize inputs
    if settings.dostandardization:
        scaler = StandardScaler()
        samples= scaler.fit_transform(samples)

    #filter samples by classification for diagnostic plotting
    pt_array = [s[0] for s in samples]
    ntrk_array = [s[1] for s in samples]

    gluon_pt_array = [a for a, b in zip(pt_array, classification) if b == 0]
    gluon_ntrk_array = [a for a, b in zip(ntrk_array, classification) if b == 0]
    gluon_weight_array = [a for a, b in zip(weight, classification) if b == 0]

    quark_pt_array = [a for a, b in zip(pt_array, classification) if b == 1]
    quark_ntrk_array = [a for a, b in zip(ntrk_array, classification) if b == 1]
    quark_weight_array = [a for a, b in zip(weight, classification) if b == 1]

    #svm
    if settings.dosvm:
        for lambdax in settings.lambdaarray:
            print '----------------------------'
            print 'C hyperparameter: ' + str(lambdax)
            clf = svm.LinearSVC(C=lambdax)#, loss="hinge")
            clf.fit(samples, classification)

            W=clf.coef_[0]
            I=clf.intercept_
            a = -W[0] / W[1]
            b = -I[0] / W[1]
            print 'slope: ' + str(a)
            print 'inter: ' + str(b)

            pass_quarks_pt = [x for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= a*x + b]
            pass_quarks_ntrk = [y for x,y in zip(quark_pt_array, quark_ntrk_array) if y <= a*x + b]

            fail_gluons_pt = [x for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > a*x + b]
            fail_gluons_ntrk = [y for x,y in zip(gluon_pt_array, gluon_ntrk_array) if y > a*x + b]

            #quark efficiency = nQuarks tagged quark / nQuarks
            q_tag_eff = float(len(pass_quarks_pt)/len(quark_pt_array))
            #gluon rejection = nGluons tagged gluon / nGluons
            g_rejection = float(len(fail_gluons_pt)/len(gluon_pt_array))

            print 'nQuarks: ' + str(len(quark_pt_array))
            print 'nQuarks correctly classified: ' + str(len(pass_quarks_pt))
            print 'q_tag_eff: ' + str(q_tag_eff)
            print 'nGluons: ' + str(len(gluon_pt_array))
            print 'nGluons correctly classified: ' + str(len(fail_gluons_pt))
            print 'g_rejection: ' + str(g_rejection)

            postfit_plots(pt_array, ntrk_array, classification, gluon_pt_array, gluon_ntrk_array, quark_pt_array, quark_ntrk_array, a, b, lambdax, q_tag_eff, g_rejection)

        if settings.doslowsvm:
            clf_no_weights = svm.SVC(gamma=1)
            clf_no_weights.fit(samples, classification)
            sample_weight_constant = np.ones(len(classification))
            svmfig = plt.figure()

            pt = [s[0] for s in samples]
            ntrk = [s[1] for s in samples]
            ptmin = min(pt)
            ptmax = max(pt)
            ntrkmin = min(ntrk)
            ntrkmax = max(ntrk)
            XX, YY = np.mgrid[ptmin:ptmax:200j, ntrkmin:ntrkmax:200j]
            Z = clf.decision_function(np.c_[XX.ravel(), YY.ravel()])

            plt.clf()
            #plt.hist2d(pt_array, ntrk_array, bins=(50,max(ntrk_array)))
            svmfig.savefig(settings.outputdir + '/svm_test.pdf')

        myfile.write(str(lambdax) + '\t' + str(a) + '\t' + str(b) + '\n')

def postfit_plots(pt_array, ntrk_array, classification, gluon_pt_array, gluon_ntrk_array, quark_pt_array, quark_ntrk_array, a, inter, lambdax, eff, rejec):
    mxx = np.linspace(0,max(pt_array))
    myy = a * mxx + inter

    #Create nTrk vs pT plots for all jets, quark only jets, gluon only jets
    allbkg = plt.figure()
    plt.scatter(pt_array, ntrk_array, c=classification,alpha=0.08)
    plt.title('PostFit All background pt vs ntrk (C=' + str(lambdax)+")") 
    plt.plot(mxx,myy,'y-', linewidth=3)
    allbkg.savefig(settings.outputdir + '/result' + str(lambdax) + 'linear".pdf')

    all_ntrk_pt_hist = plt.figure()
    plt.hist2d(pt_array, ntrk_array, bins=(50,max(ntrk_array)))
    plt.title('PostFit All pt vs ntrk (C=' + str(lambdax)+")") 
    plt.plot(mxx,myy,'y-', linewidth=3)
    all_ntrk_pt_hist.savefig(settings.outputdir + '/svm' + str(lambdax) + 'all_ntrk_pt_hist.pdf')

    gluon_ntrk_pt_hist = plt.figure()
    plt.hist2d(gluon_pt_array, gluon_ntrk_array, bins=(50,max(gluon_ntrk_array)))
    plt.title('PostFit Gluon pt vs ntrk (C=' + str(lambdax)+")") 
    plt.plot(mxx,myy,'y-', linewidth=3)
    gluon_ntrk_pt_hist.savefig(settings.outputdir + '/svm' + str(lambdax) + 'gluon_ntrk_pt_hist.pdf')

    quark_ntrk_pt_hist = plt.figure()
    quark_hist = plt.hist2d(quark_pt_array, quark_ntrk_array, bins=(50,max(quark_ntrk_array)))
    plt.title('PostFit Quark pt vs ntrk (C=' + str(lambdax)+")") 
    plt.plot(mxx,myy,'y-', linewidth=3)
    quark_ntrk_pt_hist.savefig(settings.outputdir + '/svm' + str(lambdax) + 'quark_ntrk_pt_hist.pdf')


    #Create heatmap of nTrk vs pT with fitted nTrk cut as a function of the percentage of quarks
    c1 = TCanvas( 'c1', 'c1')
    nptbins = int(max(pt_array)-min(pt_array)/settings.ptbinningfactor)
    histall = ROOT.TH2F("histall","", nptbins, min(pt_array), max(pt_array), settings.ntrkbins, settings.ntrkmin, settings.ntrkmax)
    fill_hist(histall,list(zip(pt_array,ntrk_array)))

    histquark = ROOT.TH2F("histquark","", nptbins, min(pt_array), max(pt_array), settings.ntrkbins, settings.ntrkmin, settings.ntrkmax)
    fill_hist(histquark,list(zip(quark_pt_array, quark_ntrk_array)))

    percentquark = histquark.Clone("percentquark")
    percentquark.Divide(histall)
    percentquark.GetXaxis().SetTitle("ln(p_{T})")
    percentquark.GetYaxis().SetTitle("Number of Tracks in Jet")
    titlestring = "Number of Tracks vs. ln(p_{T}), C = " + str(lambdax) + " ;ln(p_{T}); Number of Tracks in Jet; Percent of Quarks"
    percentquark.SetTitle(titlestring)

    fitted_cut = ROOT.TF1("fitted_cut","[0]*x + [1]",min(pt_array),max(pt_array));
    fitted_cut.SetParameters(a,inter)

    c1.SetRightMargin(0.15)
    c1.SetBottomMargin(0.15)
    percentquark.Draw("colz")
    fitted_cut.Draw("same")
    qeff = '%.2f'%eff
    grejec = '%.2f'%rejec
    form_slope = '%.2f'%a
    form_inter = '%.2f'%inter
    mystring = "Best Fit (nTrk = " + str(form_slope) + "x + " + str(form_inter)  + "), qeff=" + str(qeff) + " geff=" + str(grejec)
    leg = ROOT.TLegend(0.13,0.75,0.3,0.9)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(fitted_cut,mystring, "l")
    leg.Draw()
    c1.SaveAs(settings.outputdir + "/heatmap_allbkg_C" + str(lambdax) + ".pdf")


if __name__ == '__main__':
    main()

