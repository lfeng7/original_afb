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
F_eff = 'Powheg_WnJets_efficiency.root'
file_tmp = TFile(F_eff)
efficiency_b = file_tmp.Get('efficiency_b')
efficiency_c = file_tmp.Get('efficiency_c')
efficiency_udsg = file_tmp.Get('efficiency_udsg') 

# Global vars
count_out_of_bounds = 0 
count_cut_off = 0
count_jets = 0


def main():
	# print getptbin_for_btag(35)
	# print get_eta_bin_jet(2.3)
	# print SFb_error
	# main_loop()
	# main_loop()
	make_plots()
	
def main_loop():

	# COMMAND LINE OPTIONS
	parser = OptionParser()
	parser.add_option('--dir', metavar='F', type='string', action='store',dest='directory',help='') ## Sets which files to run on
	parser.add_option('--nMaxEvents', metavar='F', type='int', action='store',default='-1',dest='nMaxEvents',help='') ## choose the maximum number of iterations to go through
	(options, args) = parser.parse_args()

	#Read in input files
	FILES = options.directory + '/*.root'
	print FILES
	files = glob.glob(FILES)
	print 'Getting these files : '
	for fname in files:
		print fname
	postfix = ""
	events = Events(files)

	#AK5 jet variables
	jetPtHandle         = Handle( "std::vector<float>" )
	jetPtLabel    = ( "pfShyftTupleJets" + postfix,   "pt" )
	jetEtaHandle         = Handle( "std::vector<float>" )
	jetEtaLabel    = ( "pfShyftTupleJets" + postfix,   "eta" )
	jetPhiHandle         = Handle( "std::vector<float>" )
	jetPhiLabel    = ( "pfShyftTupleJets" + postfix,   "phi" )
	jetMassHandle         = Handle( "std::vector<float>" )
	jetMassLabel    = ( "pfShyftTupleJets" + postfix,   "mass" )
	jetCSVHandle         = Handle( "std::vector<float>" )
	jetCSVLabel    = ( "pfShyftTupleJets" + postfix,   "csv" )
	jetFlavorHandle = Handle("vector<float>")
	jetFlavorLabel = ("pfShyftTupleJets"+postfix,"jetFlavor")

	# Hists 
	hist_btag_w = TH1F('btag_w','b tagging correction weight',40,0.8,1.2)

	count = 0
	count_cut_events = 0

	##############################################################################################################
	##											Main Loop Over Events							 				##
	##############################################################################################################
	for event in events:
		count = count + 1
		# Progress report
		if count%10000 == 1 :
			print 'Processing #'+str(count)+' events!'
		#Hard cut if nMaxEvents isn't -1 for "all data"
		if count == options.nMaxEvents + 1 :
			print 'ITERATION '+str(options.nMaxEvents)+' REACHED'
			break
		#jets
		event.getByLabel (jetPtLabel, jetPtHandle)
		if not jetPtHandle.isValid() :
			continue
		toomanyjetPts = jetPtHandle.product()
		if len(toomanyjetPts) < 4 :
			count_cut_events += 1
			continue	
		event.getByLabel (jetEtaLabel, jetEtaHandle)
		toomanyjetEtas = jetEtaHandle.product()
		event.getByLabel (jetPhiLabel, jetPhiHandle)
		toomanyjetPhis = jetPhiHandle.product()
		event.getByLabel (jetMassLabel, jetMassHandle)
		toomanyjetMasss = jetMassHandle.product()
		event.getByLabel (jetCSVLabel, jetCSVHandle)
		toomanyjetCSVs = jetCSVHandle.product()
		event.getByLabel (jetFlavorLabel,jetFlavorHandle)
		toomanyjetFlavors = jetFlavorHandle.product()
		# Get the event weight
		w_result = get_weight_btag(toomanyjetPts,toomanyjetEtas,toomanyjetFlavors,toomanyjetCSVs)
		b_tag_w = w_result[0]
		b_tag_w_err = w_result[1]
		# Fill hist
		hist_btag_w.Fill(b_tag_w)
		# print 'b tagging weight : '+str(b_tag_w)+' +/- '+str(b_tag_w_err)
		# print 'Event # '+str(count)

	print 'total number of jets passed cuts are :'+str(count_jets)		
	print 'total number of jets out of bounds : '+str(count_out_of_bounds)	
	print 'total number of events been cut off is :'+str(count_cut_events)
	c2 = TCanvas('c2','b tagging weight',200,10,700,500)
	hist_btag_w.Draw()
	c2.SaveAs('b_tagg_w.eps')

def make_plots():
	eta = 1.5
	jet_flavor = 5
	pt = 75
	n =40
	pt_max = 1000
	step = pt_max/n

	title = 'btag_eff eta = '+str(eta) # +' flavor = '+str(jet_flavor)
	x_title = 'pt (GeV)'
	y_title = 'btag_eff'
	file_name = 'btag_eff_1p5_WJets.eps'
	x = array('d')
	y = []
	flavors = [5,4,3]
	for i in range(n):
		pt_tmp = i*step
		x.append(pt_tmp)
	for i in range(len(flavors)) :
		y_tmp = array('d')
		for j in range(n) :
			pt_tmp = j*step
			btag_eff = get_btag_eff(pt_tmp,eta,flavors[i])
			y_tmp.append(btag_eff)
		y.append(y_tmp)
	# print x
	# print y
	c1 = TCanvas("c1","A Simple Graph ",200,10,700,500)
	c1.SetGrid();
	color_list = [2,4,7]
	mg = TMultiGraph();
	for i in range(len(flavors)):
		# i=2
		gr = TGraph(n,x,y[i]);
		gr.SetLineColor(color_list[i]);
		gr.SetLineWidth(4);
		gr.SetMarkerColor(1);
		gr.SetMarkerStyle(21);
		gr.SetTitle(title);
		gr.GetXaxis().SetTitle(x_title);
		gr.GetYaxis().SetTitle(y_title);
		mg.Add(gr)
	mg.Draw("ACP");
	
	c1.SaveAs(file_name)

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


def getptbin_for_btag(pt):	# checked
	if(pt<30) : pt_bin = 0;
	elif(pt<40) : pt_bin = 1;
	elif(pt<50) : pt_bin = 2;
	elif(pt<60) : pt_bin = 3;
	elif(pt<70) : pt_bin = 4;
	elif(pt<80) : pt_bin = 5;
	elif(pt<100) : pt_bin = 6;
	elif(pt<120) : pt_bin = 7;
	elif(pt<160) : pt_bin = 8;
	elif(pt<210) : pt_bin = 9;
	elif(pt<260) : pt_bin = 10;
	elif(pt<320) : pt_bin = 11;
	elif(pt<400) : pt_bin = 12;
	elif(pt<500) : pt_bin = 13;
	elif(pt<600) : pt_bin = 14;
	else : pt_bin = 15;
	return pt_bin


def get_eta_bin_jet(eta) : # checked
	eta = math.fabs(eta);
	# float etaNBins[5] = {0., 0.9, 1.2, 2.1, 2.4};
	if(eta<0.9) : return 0;
	elif(eta<1.2) : return 1;
	elif(eta<2.1) : return 2;
	elif(eta<2.4) : return 3;
	else : return -1;


# // for CSVM only 
# Checked
SFb_error = [
	0.0415694,
	0.023429,
	0.0261074,
	0.0239251,
	0.0232416,
	0.0197251,
	0.0217319,
	0.0198108,
	0.0193,
	0.0276144,
	0.0205839,
	0.026915,
	0.0312739,
	0.0415054,
	0.0740561,
	0.0598311 
]


# //only for CSVM point 
def get_SF_btag(ptJet,etaJet,flavJet):		# checked
	global count_out_of_bounds
	x = ptJet # the pt of the jet 
	eta = fabs(etaJet); # abs(eta) 

	result = []			# the first is SF, second is SF error

	if(eta>2.4) :
		print 'warning SF_btag_eta>2.4 ??  ' +str(eta); 
		result.append(1)
		result.append(0)
		return result
	
	if(x<20) : x=20; 
	if(x>800) : x= 800;

	if(abs(flavJet)==5 or abs(flavJet) == 4): # for b or c. flavJet[indj] refers to the MC-true flavor of the jet
		SF  = (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x)); 
		ptbin = getptbin_for_btag( ptJet ); #--> ptJet[indj] refers to the pt of this jet
		SFerr = SFb_error[ptbin];
		if(x>800 or x<20 ) : SFerr *= 2;
		if(abs(flavJet) == 4 ) : SFerr *= 2; 
	else : #SFlight
		if(eta<=0.8):
			SF = ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
			max = ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
			min = ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));		
		elif(eta<=1.6):
			SF = ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
			max = ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
			min = ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
		elif (eta>1.6 and eta<=2.4):
			SF = ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
			max = ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
			min = ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
		SFerr = math.fabs(max-SF) if math.fabs(max-SF)>math.fabs(min-SF) else math.fabs(min-SF);

	result.append(SF)
	result.append(SFerr)
	return result

# ///this is the functin to call to get the event-by-event b-tag weight
# ///as implemented, it only works for CSVM. 
# Input : given a list of selected jets (with pt, ID, etc)

def get_weight_btag(jet_Pts, jet_etas, jet_flavors,jet_csvs) :
	global count_cut_off,count_jets,count_out_of_bounds
	# Initiate probability and its uncertainty
	# Probability
	mcTag = 1.;
	mcNoTag = 1.;
	dataTag = 1.;
	dataNoTag = 1.;
	# uncertainty
	err1 = 0; 
	err2 = 0; 
	err3 = 0; 
	err4 = 0; 
	# Loop over all jet candidates
	for i in range(len(jet_Pts)):
		count_jets += 1
		jet_Pt = jet_Pts[i]
		jet_eta = math.fabs(jet_etas[i])
		jet_flavor = jet_flavors[i]
		jet_csv = jet_csvs[i]
		# Basic cuts to make sure things go right.
		if(jet_eta>2.4): 
			count_out_of_bounds += 1
			continue
		if(jet_flavor==0) :
			count_out_of_bounds += 1
			continue	#for jets with flavor 0, we ignore. 
		etabin = get_eta_bin_jet(jet_eta);
		# Get b-tagging efficiency using pt,eta and jet_flavor infor
		eff = get_btag_eff(jet_Pt,jet_eta,jet_flavor)
		# if eff == 0 or eff == 1:
		# 	print 'jet pt,eta,flavor,csv '+str(jet_Pt)+','+str(jet_eta)+','+str(jet_flavor)+','+str(jet_csv)
		# 	eff = 0.5
		# Get SF for this jet
		SF_result = get_SF_btag(jet_Pt,jet_eta,jet_flavor)
		jet_SF = SF_result[0]
		jet_SFerr = SF_result[1]
		# Get probability
		istag = jet_csv > 0.679 and math.fabs(jet_eta)<2.4 ;
		if istag :
			mcTag *= eff; 
			dataTag *= eff*jet_SF; 
			# debug
			if eff == 0 or eff*jet_SF == 0:
				print 'jet pt,eta,flavor,csv,eff,SF: '+str(jet_Pt)+','+str(jet_eta)+','+str(jet_flavor)+','+str(jet_csv)+','+str(eff)+','+str(jet_SF)
			# Get error
			if(jet_flavor==5 or jet_flavor ==4) : 
				err1 += jet_SFerr/jet_SF; # correlated for b/c
			else : 
				err3 += jet_SFerr/jet_SF; # correlated for light								
		else :
			mcNoTag *= (1- eff); 
			dataNoTag *= (1- eff*jet_SF); 
			#debug
			if (1-eff) == 0 or (1-eff*jet_SF)==0 :
				print 'jet pt,eta,flavor,csv,eff,SF: '+str(jet_Pt)+','+str(jet_eta)+','+str(jet_flavor)+','+str(jet_csv)+','+str(eff)+','+str(jet_SF)
			# Get error
			if(jet_flavor==5 or jet_flavor ==4 ) : 
				err2 += (-eff*jet_SFerr)/(1-eff*jet_SF); # correlated for b/c
			else : 
				err4 +=  (-eff*jet_SFerr)/(1-eff*jet_SF);  # correlated for light
	# Check if any of the probability is zero
	if dataTag*dataNoTag == 0  or mcTag*mcNoTag == 0:
		print 'one of the probability is zero!'
	# Get event weight for this event
	wtbtag = (dataNoTag * dataTag ) / ( mcNoTag * mcTag ); 
	wtbtagErr = math.sqrt( pow(err1+err2,2) + pow( err3 + err4,2)) * wtbtag; # un-correlated for b/c and light
	# Set return
	to_return = []
	to_return.append(wtbtag)
	to_return.append(wtbtagErr)
	return to_return

main()

