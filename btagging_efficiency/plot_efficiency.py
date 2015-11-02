##############################################################################################################
##										Imports and Other Header Stuff					 					##
##############################################################################################################
import os
import glob
import math
import numpy
import ROOT
from ROOT import *
import sys
from DataFormats.FWLite import Events, Handle
from optparse import OptionParser
from array import *


# Read efficiency files
# F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_AK5PF_CSVM_bTaggingEfficiencyMap.root'
# directory = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/codes/tagging/dir_jetCands/dir_efficiency/'
directory = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/codes/tagging/dir_jetCands/bTagging_efficiency_outputs/'
channel = 'W4Jets_no_cuts'
F_eff = directory+channel+'_efficiency.root'
# F_eff = ' /uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_AK5PF_CSVM_bTaggingEfficiencyMap.root '
file_tmp = TFile(F_eff)
print 'Getting this file: '+F_eff
efficiency_b = file_tmp.Get('efficiency_b')
efficiency_c = file_tmp.Get('efficiency_c')
efficiency_udsg = file_tmp.Get('efficiency_udsg') 

def main() :
	plot()

def get_btag_eff (pt,eta,jet_flavor):
	# x,y of TH2F of efficiency are pt and eta
	if jet_flavor == 5 :
		binx = efficiency_b.GetXaxis().FindBin(pt)
		biny = efficiency_b.GetYaxis().FindBin(eta)
		bin = efficiency_b.GetBin(binx,biny)
		return efficiency_b.GetBinContent(bin)
	elif jet_flavor == 4 :
		binx = efficiency_c.GetXaxis().FindBin(pt)
		biny = efficiency_c.GetYaxis().FindBin(eta)
		bin = efficiency_c.GetBin(binx,biny)
		return efficiency_c.GetBinContent(bin)	
	else :
		binx = efficiency_udsg.GetXaxis().FindBin(pt)
		biny = efficiency_udsg.GetYaxis().FindBin(eta)
		bin = efficiency_udsg.GetBin(binx,biny)
		return efficiency_udsg.GetBinContent(bin)

eta = [0.8,1.1,1.6,2.3]
pt = 75
n =40
pt_max = 1000
step = pt_max/n

def plot() :
	title =[]
	for i in range(len(eta)) :
		title.append(channel+' btagging efficiency eta = '+str(eta[i])) # +' flavor = '+str(jet_flavor)
	x_title = 'pt (GeV)'
	y_title = 'btag_eff'
	file_name = ['btag_eff_0p8_'+channel+'.eps','btag_eff_1p1_'+channel+'.eps','btag_eff_1p6_'+channel+'.eps','btag_eff_2p3_'+channel+'.eps']
	x = array('d')
	flavors = [5,4,3]
	flavor_names = ['b','c','udsg']
	color_list = [2,4,7]

	for i in range(n):
		pt_tmp = i*step
		x.append(pt_tmp)

	for eta_bin in range(len(eta)) :
		y = []
		for i in range(len(flavors)) :
			y_tmp = array('d')
			for j in range(n) :
				pt_tmp = j*step
				btag_eff = get_btag_eff(pt_tmp,eta[eta_bin],flavors[i])
				y_tmp.append(btag_eff)
			y.append(y_tmp)
		# print x
		# print y
		c1 = TCanvas("c1","A Simple Graph ",200,10,700,500)
		c1.SetGrid();
		mg = TMultiGraph();
		leg = TLegend(0.7,0.7,0.9,0.9)
		for i in range(len(flavors)):
			# i=2
			gr = TGraph(n,x,y[i]);
			gr.SetLineColor(color_list[i]);
			gr.SetLineWidth(4);
			gr.SetMarkerColor(1);
			gr.SetMarkerStyle(21);
			#debug
			# print title[eta_bin]
			gr.SetTitle(title[eta_bin]);
			gr.GetXaxis().SetTitle(x_title);
			gr.GetYaxis().SetTitle(y_title);
			mg.Add(gr)
			leg.AddEntry(gr,flavor_names[i],'LP')

		mg.SetTitle(title[eta_bin]);	
		mg.Draw("ACP");
		leg.Draw()

		c1.SaveAs(file_name[eta_bin])

	file_tmp.Close()

main()
