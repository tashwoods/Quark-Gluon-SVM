#!/usr/bin/python
from __future__ import division
import ROOT, os, math, glob
import numpy as np
from ROOT import THStack
from ROOT import gROOT
from ROOT import TPad
import matplotlib.pyplot as plt



#-----------------Meta Settings--------------------------
debug = 0
nevents=-1
f = ROOT.TFile.Open("/gpfs/slac/atlas/fs1/d/woodsn/data/trees/nominal_rnn/nominal_rnn/all_onelep_trees/all.root")
#f = ROOT.TFile.Open("/gpfs/slac/atlas/fs1/d/woodsn/newrf/VVSemileptonicStats/SlimmedNtuples/onelep_mc16ade_rnn/Nominal/data/itest01/HVTWW_500_1lep.root")
t = f.Get("Nominal")

ptarray1=[]
ntrkarray1=[]
qgtaggable1=[]
ptarray2=[]
ntrkarray2=[]
qgtaggable2=[]
weight=[]
resmass=[]

#Save pdgid, ntrk, eta from background mc
count = 0
for event in t: 
    count+=1
    if count%1000==0:
        print 'processed ' + str(count)
    if count > nevents and nevents != -1:
        break
    if (t.Pass_Res_GGF_WW_SR or t.Pass_Res_GGF_WZ_Untag_SR or t.Pass_Res_VBF_WW_SR or t.Pass_Res_VBF_WZ_SR):
        ptarray1.append(t.sigWJ1_pt)
        ntrkarray1.append(t.sigWJ1_nTrk)
        qgtaggable1.append(t.sigWJ1_qgtaggable)
        ptarray2.append(t.sigWJ2_pt)
        ntrkarray2.append(t.sigWJ2_nTrk)
        qgtaggable2.append(t.sigWJ2_qgtaggable)
    weight.append(t.weight)
    resmass.append(t.X_resolved_WW_m)
print 'processed ' + str(count)

np.savez('july16_allbkg_eventtreetoarray.npz', ptarray1=ptarray1, ntrkarray1=ntrkarray1, qgtaggable1=qgtaggable1, ntrkarray2=ntrkarray2, ptarray2=ptarray2, qgtaggable2=qgtaggable2, weight=weight, resmass=resmass)


