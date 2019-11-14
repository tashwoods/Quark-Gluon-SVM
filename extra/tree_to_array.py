#!/usr/bin/python
from __future__ import division
import ROOT, os, math, glob, sys
import numpy as np
from ROOT import gROOT
from ROOT import TPad
import matplotlib.pyplot as plt
#-----------------Meta Settings--------------------------
#set debug = 1 if more verbose output is desired
debug = 0
#specify number of events to save, -1 will save all events
nevents=1000
ptmin = 50 #Minimum pT of jets that are saved
etamax = 2.1 #Maximum eta value of jets that are saved

#f = ROOT.TFile.Open("/gpfs/slac/atlas/fs1/d/woodsn/data/trees/nominal_rnn/nominal_rnn/all_onelep_trees/all.root")
inputfile = sys.argv[1]
savename = sys.argv[2]

f = ROOT.TFile.Open("/gpfs/slac/atlas/fs1/d/woodsn/data/qgthesis_files/files/" + inputfile)
t = f.Get("nominal")

#create arrays to store jet information for npz file
classification=[]
ptarray=[]
ntrkarray=[]
weight=[]
resmass=[]

#Save pdgid, ntrk, eta from input file
count = 0
for event in t: 
    count+=1
    if count%1000==0:
        print 'processed ' + str(count)
    if count > nevents and nevents != -1:
        break
    #only save events that pass signal region selection
    if (t.Pass_Res_GGF_WW_SR):
        if t.sigWJ1_pt > ptmin and t.sigWJ1_eta < etamax:
            if t.sigWJ1_pdgid == 21: #gluon selection
                ptarray.append(t.sigWJ1_pt)
                ntrkarray.append(t.sigWJ1_nTrk)
                classification.append(0)
                weight.append(t.weight)
                resmass.append(t.X_resolved_WW_m)
            elif t.sigWJ1_pdgid >= 0 and t.sigWJ1_pdgid <= 6: #quark selection
                ptarray.append(t.sigWJ1_pt)
                ntrkarray.append(t.sigWJ1_nTrk)
                classification.append(1)
                weight.append(t.weight)
                resmass.append(t.X_resolved_WW_m)
print 'processed ' + str(count)
samples = zip(ptarray, ntrkarray)

np.savez('npzfiles/' + savename + '.npz', samples=samples, classification=classification, weight=weight, resmass=resmass)


