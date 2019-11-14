#!/usr/bin/python
import ROOT, os, math, glob, sys
import numpy as np
from ROOT import THStack
from ROOT import gROOT
from ROOT import TPad
import matplotlib.pyplot as plt
thres = 0
def getrate(num,denom,t):
    print num
    print denom
    hist1 = ROOT.TH1F("hist1","", bins, mins, maxs)
    makehisto = variables + ">>hist1"
    t.Draw(makehisto,num,"goff")
    tagged = hist1.GetSumOfWeights()

    hist2 = ROOT.TH1F("hist2","", bins, mins, maxs)
    makehisto2 = variables + ">>hist2"
    t.Draw(makehisto2,denom,"goff")
    possible = hist2.GetSumOfWeights()

    if tagged!=0: yields = tagged/possible
    else: yields = 0
    print 'tagged: ' + str(tagged) + ' possible: ' + str(possible)
    print 'YIELD: ' + str(yields)

    return tagged, possible, yields

g = open('yields.txt', 'w')
g.write('low_pt_cut \t\t high_pt_cut \t slope \t\t intercept \t\t\t 300qtagged \t 300qpossible \t 300qeff \t 700qtagged \t\t\t 700qpossible \t\t 700qeff \t gtagged \t\t\t gpossible \t\t\t gtag eff \t\t\t rejection \n')
mymin = 0
#s_array = [0,1,2,3,4,5,6,7]
#i_array = [-100, -10,-5,0, 5, 10]
s_array=[4]
i_array=[-5]
thres_array = [10,11,12,13,14,15,16,17,18,19,20]
acceptances = []
rejections =[]
mass_array = [230,1500]
mass_string_array = []
array_s_over_b=[]
array_b=[]
array_s_rootb_complex=[]
array_s_over_root_b=[]
'''
for i in range(len(mass_array)-1):
    print i
    mass_string = "(X_resolved_WW_m > " + str(mass_array[i]) + " && X_resolved_WW_m < " + str(mass_array[i+1]) + ")"
    mass_string_array.append(mass_string)
print mass_string_array
'''
mass_string_array.append('(1==1)')
myselection = "weight*(Pass_Res_GGF_WW_SR)"
quark_selection = "(sigWJ1_pdgid > 0 && sigWJ1_pdgid < 7)"
gluon_selection = "(sigWJ1_pdgid==21)"

c1 = ROOT.TCanvas("c1", "c1", 0,  0, 1000, 700)

variables = "X_resolved_WW_m"
bins = 2
mins = 0
maxs = 5000

print 'using dy samples'

#x = ROOT.TFile.Open("/gpfs/slac/atlas/fs1/d/woodsn/ResonanceFinder/Run/test_run/WVlvqq_mc15_v28/data/itest01/Wjets.root")
x = ROOT.TFile.Open("/gpfs/slac/atlas/fs1/d/woodsn/data/qgthesis_files/files/allbkg-ade.root")
z = x.Get("nominal")
if thres == 0:
    for s in s_array:
        for i in i_array:
            for m_string in mass_string_array:
                selection = m_string + "*" + myselection
                print 'selection: ' + selection
                #~~~~~~~~~~quark efficiency~~~~~~~~~~~~~~~~~~~
                print "QUARK EFFICIENCY CALCULATION-------------------------------------"
                qgselection = "(sigWJ1_nTrk > " + str(mymin) + " && sigWJ1_nTrk <= (" + str(s) + "*TMath::Log(sigWJ1_pt) + " +str(i) +" ))"
                numselection = selection + " * " + qgselection + " * " + quark_selection
                denomselection = selection + " * " + quark_selection 
                low_pt_ntrk_cut = s*ROOT.TMath.Log(75)+i
                high_pt_ntrk_cut = s*ROOT.TMath.Log(175)+i
                yquark_selection = quark_selection + " * " + selection
                ygluon_selection = gluon_selection + "*" + selection

                print "W+jets --------------------------------------"
                print 'quark composition'
                wjets_q, total_wjets, qyield_wjets = getrate(yquark_selection, selection, z)
                print 'gluon composition'
                wjets_g, gtotal_wjets, gyield_wjets = getrate(ygluon_selection, selection, z)
                print 'quark tag efficiency'
                qtagged_wjets, qpossible_wjets, yields= getrate(numselection, denomselection, z)


                #~~~~~~~~~~~~~~gluon rejection~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                print "GLUON REJECTION CALCULATION-------------------------------------"
                gselection = "(sigWJ1_nTrk > (" + str(s) + "*TMath::Log(sigWJ1_pt) + " +str(i) +" ))"
                gnumselection = selection + "*" + gselection + " * " + gluon_selection
                gdenomselection = selection + " * " + gluon_selection 

                gtagged, gpossible, gyields = getrate(gnumselection, gdenomselection, z)
                if gyields != 1: rejection = 1/(1-gyields)
                else: rejection = 999

                acceptances.append(yields)
                rejections.append(rejection)

    


print 'b yields' 
print array_b
print 'acceptances:'
print acceptances
print 'rejections:'
print rejections
g.close()
fig = plt.figure()
ax1=fig.add_subplot(111)
ax1.scatter(acceptances,rejections,s=10, c='b', marker="s", label='300 GeV VBF Z')
#plt.legend(loc='upper right')
plt.xlabel("Quark Efficiency")
plt.ylabel("Gluon Rejection")
fig.savefig('roc.pdf')


print 'here'
finalacceptances=[0.16793917409493145, 0.693984157674405, 0.1003078491474742, 0.6236039025147464, 0.9229092832426578, 0.08119936098259366, 0.6009207869259742, 0.920655849504518, 0.988716035755272,  0.06706407057719604, 0.5723555946131719, 0.9154173761626906, 0.9888682308262193, 0.9985986969237992, 0.5403770132603188, 0.9081171285579256, 0.988651875088108, 0.9987112445392896, 0.9996524452896646, 0.8973840721015215, 0.9881839924820446, 0.998745349877317, 0.9996739742842943, 0.9997635007966162,  0.9873380735197838, 0.9987825460116031, 0.9996943309079295, 0.9997637139549789, 0.9997706416017658, 0.998773593360371, 0.9997026440840736, 0.9997660586969683, 0.9997706416017658, 0.9997729863437551]
finalrejections=[20.767928352628598, 2.450865683237959, 44.37553186870447, 3.081203013752214, 1.3122570801299442, 69.02566368195026, 3.4881398795430525, 1.333371565237154, 1.0565367091913802,  99.46625846478207, 3.970377605009699, 1.3631966794282286, 1.0563282123505575, 1.0081147396188015, 4.572835211543301, 1.4028014769080004, 1.0577010709572108, 1.007655765571756, 1.00102472743529, 1.4541485997578807, 1.0600928686444808, 1.0071518071638148, 1.0008396461988067, 1.000107268642809,  1.06436563428647, 1.0069253301922028, 1.000711157502614, 1.0000888310344154, 1.0000063411265006, 1.0068261821289495, 1.0006463559162717, 1.0000813733879093, 1.000001744653252, 1]
print 'finalacceptances length'
print len(finalacceptances)
print 'finalrejection length:'
print len(finalrejections)
finalacceptances.sort()
finalrejections.sort(reverse=True)
print finalacceptances
print finalrejections
fig = plt.figure()
ax1=fig.add_subplot(111)
ax1.scatter(finalacceptances,finalrejections,s=10, c='b', marker="s", label='300 GeV VBF Z')
#plt.legend(loc='upper right')
plt.ylim(-0.1,120)
plt.plot(finalacceptances,finalrejections, '-o')
plt.xlabel("Quark Efficiency")
plt.ylabel("Gluon Rejection")
fig.savefig('finalroc.pdf')


'''
ffig = plt.figure()
fax1=ffig.add_subplot(111)
fax1.scatter(aquarkfake300,agluonfake,s=10, c='b', marker="s", label='300 GeV VBF Z')
fax1.scatter(aquarkfake700,agluonfake,s=10, c='yellow', marker="s", label='700 GeV VBF Z')
plt.legend(loc='upper right')
plt.xlabel("Quark Fake Rate")
plt.ylabel("Gluon Fake Rate")
ffig.savefig('fakerate.pdf')
'''
