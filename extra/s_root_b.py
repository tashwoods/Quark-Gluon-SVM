from __future__ import division
from array import array
import matplotlib.pyplot as plt
import math
import ROOT, os, sys
from ROOT import gROOT, TPad, TProfile, TCanvas, gStyle
#from sklearn.svm import LinearSVC
from sklearn import svm
import numpy as np
import random
from random import sample 
ROOT.gROOT.SetBatch(True)

#User Settings
nevents = -1
quicktest = 1
debug = 1
dologpt = 1
dontrk_perjet_optimization = 0
do_eventlevel_optimiziation = 1
ROOT.gStyle.SetOptStat(0)
ntrkmin = 0
ntrksearchmin = 1
ntrkmax = 25
ntrkbins = 10
ntrkbinmin = 0 
ntrkbinmax = 1000
bin_resmass_array = [0,300,330,370,410,450,490,540,590,640,700,760,820,890,960,1040,1130,1220,1320,1430,1540,1660,1790]
npzfiledirectory = "/gpfs/slac/atlas/fs1/d/woodsn/myscripts/training_minions/npzfiles/pass_any_signal_region_june6/"
bkgnpzfilelocation = npzfiledirectory + "bkg.npz"
signalfile = npzfiledirectory + 'hvtww300.npz'
#signalfilelocation = [npzfiledirectory + s for s in signalfiles]
def loadarrays(npzfile):
    #load bkg npz file
    npzfile=np.load(npzfile)
    ptarray1 = npzfile['ptarray1']
    ntrkarray1 = npzfile['ntrkarray1']
    qgtaggable1 = npzfile['qgtaggable1']

    ptarray2 = npzfile['ptarray2']
    ntrkarray2 = npzfile['ntrkarray2']
    qgtaggable2 = npzfile['qgtaggable2']

    weight = npzfile['weight']
    resmass = npzfile['resmass']

    if dologpt:
        ptarray1 = [math.log(s) for s in ptarray1]
        ptarray2 = [math.log(s) for s in ptarray2]
    weight = [w + 0 for w in weight]
    return ptarray1, ntrkarray1, qgtaggable1, ptarray2, ntrkarray2, qgtaggable2, weight, resmass

def main():
    #create outputfiles
    os.system('rm -r output')
    os.system('mkdir output')
    myfile = open('output/linearpt_slope_intercept.txt','w')

    #load signal and background arrays from npz file
    ptarray1, ntrkarray1, qgtaggable1, ptarray2, ntrkarray2, qgtaggable2, weight, resmass = loadarrays(bkgnpzfilelocation)
    sptarray1, sntrkarray1, sqgtaggable1, sptarray2, sntrkarray2, sqgtaggable2, sweight, sresmass = loadarrays(signalfile)

    if do_eventlevel_optimiziation:
        slope = 4
        inter = -5

    if dontrk_perjet_optimization: 
        ptbins = np.arange(4,7.25,0.25)
        print ptbins
        for i in range(len(ptbins)-1):
            s_root_b_array=[]
            print 'Background-----------------------------------------------'
            background_yields, background_event_count, background_ntrk_cuts = get_ntrk_cut(ptbins[i],ptbins[i+1], ptarray1, ntrkarray1, qgtaggable1, weight, resmass)
            print 'Signal-----------------------------------------------'
            signal_yields, signal_event_count, signal_ntrk_cuts = get_ntrk_cut(ptbins[i],ptbins[i+1], sptarray1, sntrkarray1, sqgtaggable1, sweight, sresmass)
            print 'background_yields:'
            print background_yields
            print 'signal_yields:'
            print signal_yields

            for i,j in zip(background_event_count, signal_event_count):
                if i != 0:
                    s_root_b = j/math.sqrt(i)
                else:
                    s_root_b = 0
                s_root_b_array.append(s_root_b)
            print 's_root_b_array'
            print s_root_b_array

            max_sig_index = s_root_b_array.index(max(s_root_b_array))
            maxsig = s_root_b_array[max_sig_index]
            maxsig_ntrk = background_ntrk_cuts[max_sig_index]
            sig_yield = signal_yields[max_sig_index]
            bkg_yield = background_yields[max_sig_index]

            print 'maxsig index: ' + str(max_sig_index) + ' maxsig value: ' + str(maxsig) + ' maxsig ntrk: ' + str(maxsig_ntrk) + ' bkgyield: ' + str(bkg_yield) + ' sigyield: ' + str(sig_yield)



def get_ntrk_cut(ptmin, ptmax, pt_array, ntrk_array, qgtaggable_array, weight_array, resmass_array):
    c4 = TCanvas( 'c4', '', 200, 10, 700, 500 )
    yield_array=[]
    event_count_array=[]
    ntrkcuts = range(ntrksearchmin, ntrkmax)
    cut_counts=[]

    pass_pt_array = [x for x,y in zip(pt_array, ntrk_array) if x >= ptmin and x < ptmax]
    pass_ntrk_array = [y for x,y in zip(pt_array, ntrk_array) if x >= ptmin and x < ptmax] 
    pass_weight_array = [y for x,y in zip(pt_array, weight_array) if x >= ptmin and x < ptmax]
    pass_resmass_array = [y for x,y in zip(pt_array, resmass_array) if x >= ptmin and x < ptmax]
    for thisntrack in ntrkcuts:
        print 'ntrk: ' + str(thisntrack)
        qgpass_pt_array = [x for x,y,z in zip(pass_pt_array, pass_ntrk_array, qgtaggable_array) if y >= ntrkmin and y <= thisntrack or z==0]
        qgpass_weight_array = [x for x,y,z in zip(pass_weight_array, pass_ntrk_array, qgtaggable_array) if y >= ntrkmin and y <= thisntrack or z==0]
        qgpass_resmass_array = [x for x,y,z in zip(pass_resmass_array, pass_ntrk_array, qgtaggable_array) if y >= ntrkmin and y <= thisntrack or z==0]
        qgpass_ntrk_array = [y for y,z in zip(pass_ntrk_array, qgtaggable_array) if y >= ntrkmin and y <= thisntrack or z==0] 

        whistall = ROOT.TH1F("whistall","",20, -10, 10)
        histall = ROOT.TH1F("histall", "", int(len(bin_resmass_array)-1), array('d', bin_resmass_array))

        passweightsum = 0
        for n,j in zip(pass_resmass_array, pass_weight_array):
            histall.Fill(n,j)
            whistall.Fill(j)
            passweightsum+=j

        dohist = 1
        domyweight = 1
        if dohist:
            c3 = ROOT.TCanvas("c3", "c3", 0,  0, 1500, 1000)
            whistpass = ROOT.TH1F("whistpass","",20, -10, 10)
            #histpass = ROOT.TH1F("histpass","",ntrkbins, ntrkbinmin, ntrkbinmax)
            histpass = ROOT.TH1F("histpass", "", int(len(bin_resmass_array)-1), array('d', bin_resmass_array))
            qgpassweightsum = 0
            for n,j in zip(qgpass_resmass_array, qgpass_weight_array):
                histpass.Fill(n,j)
                whistpass.Fill(j)
                qgpassweightsum+=j

            if histall.GetSumOfWeights() > 0:
                qgyield = histpass.GetSumOfWeights()/histall.GetSumOfWeights()
            if whistall.GetSumOfWeights() > 0:
                wqgyield = whistpass.Integral()/whistall.Integral()
            else:
                qgyield = 0

            histall.Draw("hist")
            histpass.SetLineColor(2)
            histpass.Draw("hist same")
            c3.Print("ntrk" + str(thisntrack) + '_pt_' + str(ptmin) + '_' + str(ptmax) + '.pdf')

            whistall.Draw("hist")
            whistpass.SetLineColor(2)
            whistpass.Draw("hist same")
            c3.Print("weight_ntrk" + str(thisntrack) + '_pt_' + str(ptmin) + '_' + str(ptmax) + '.pdf')

            event_count_array.append(histpass.GetSumOfWeights())
            yield_array.append(qgyield)


            if debug:
                if histall.Integral(0, histall.GetNbinsX()+1) != 0:
                    integralyield = histpass.Integral(0, histall.GetNbinsX()+1)/histall.Integral(0, histall.GetNbinsX()+1)
                else:
                    integralyield = 0
                if passweightsum != 0:
                    array_yield = qgpassweightsum/passweightsum
                else:
                    array_yield = 0
                print 'Integral: all: ' + str(histall.Integral(0, histall.GetNbinsX()+1)) + ' pass: ' + str(histpass.Integral(0, histall.GetNbinsX()+1)) + ' yield: ' + str(integralyield)
                print 'SumOfWeights: all: ' + str(histall.GetSumOfWeights()) + ' pass: ' + str(histpass.GetSumOfWeights()) + ' yield: ' + str(qgyield)
                print 'array weight all: ' + str(passweightsum) + ' pass: ' + str(qgpassweightsum) + ' yield: ' + str(array_yield)
                print 'len of QGpass_pt_array: ' + str(len(qgpass_pt_array))
                print 'len of QGpass_ntrk_array: ' + str(len(qgpass_ntrk_array))
                print 'len of QGpass_weight_array: ' + str(len(qgpass_weight_array))
                print 'len of QGpass_resmass_array: ' + str(len(qgpass_resmass_array))
            del histpass, histall, whistpass, whistall


    if debug:
        print 'ptmin: ' + str(ptmin) + ' ptmax: ' + str(ptmax)
        print 'len of pass_pt_array: ' + str(len(pass_pt_array))
        print 'len of pass_ntrk_array: ' + str(len(pass_ntrk_array))
        print 'len of pass_weight_array: ' + str(len(pass_weight_array))
        print 'len of pass_resmass_array: ' + str(len(pass_resmass_array))



    return yield_array, event_count_array, ntrkcuts


if __name__ == '__main__':
    main()

