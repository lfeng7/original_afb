##############################################################################################################
##										Imports and Other Header Stuff					 					##
##############################################################################################################
# This code will read in selected events root file and generate b-tagging efficiency root files. The efficiency files are
# the input of tagger.
# This version is made on 10/9/2014
import os
import glob
import math
import ROOT
from ROOT import *
import sys
from DataFormats.FWLite import Events, Handle
from optparse import OptionParser
from array import array

cutoff = -1 
# Read efficiency files
dir_path = 'file_jetCands/'
sample_name = 'TTJets_FullLep_all'
input_file = dir_path+sample_name+'.root'
file_tmp = TFile(input_file)
t = file_tmp.Get('output')
print 'Get this input_file: '+input_file
#Set up output file
output_dir = './bTagging_efficiency_outputs/'
outputname = output_dir+sample_name+'_efficiency.root'
f = TFile( outputname, "Recreate" )
print 'The output file would be : '+outputname
tree1 = TTree("output", "output")
tree1.SetDirectory(0)
# constants
pt_min = 10
pt_max = 500
eta_min = 0
eta_max = 2.4
CSVM = 0.679
# Binning
# bins_pt = array('d',[0.5,30.0,40.0,50.0,60.0,70.0,80.0,100.0,120.0,160.0,210.0,260.0,320.0,400.0,500.0,600.0,4000.0]) # For signal
bins_pt_b = array('d',[0., 40., 60., 80., 100., 150., 200., 300., 400., 500., 1000.]) # For ttbar # the same as Sal's code
bins_eta_b = array('d',[0., 0.6, 1.2, 2.4])
bins_pt_c = array('d',[0., 40., 60., 80., 100., 150., 200., 300., 400., 1000.]) # For ttbar # the same as Sal's code
bins_eta_c = array('d',[0., 0.6, 1.2, 2.4])
bins_pt_udsg = array('d',[0., 40., 60., 80., 100., 150., 200., 300., 400., 1000.]) # For ttbar # the same as Sal's code
bins_eta_udsg = array('d',[0., 0.6, 1.2, 2.4])

pt_binning = [bins_pt_b,bins_pt_c,bins_pt_udsg]
eta_binning = [bins_eta_b,bins_eta_c,bins_eta_udsg]

pt_bin = []
# len(bins_pt)-1 # this is what should be
eta_bin = []
# len(bins_eta)-1
for i in range(len(pt_binning)) :
	pt_bin.append(len(pt_binning[i])-1)
	eta_bin.append(len(eta_binning[i])-1)
# Set up 2D histograms
h_name_num = ['b_num','c_num','udsg_num']
h_name_denom = ['b_denom','c_denom','udsg_denom']
h_name_eff = ['efficiency_b','efficiency_c','efficiency_udsg']
eff_num_list = []
eff_denom_list = []
eff_list = []
for i in range(len(h_name_num)):
	eff_num_list.append(TH2D(h_name_num[i],h_name_num[i],pt_bin[i],pt_binning[i],eta_bin[i],eta_binning[i]))
	eff_denom_list.append(TH2D(h_name_denom[i],h_name_denom[i],pt_bin[i],pt_binning[i],eta_bin[i],eta_binning[i]))

# reset counts
jet_count = evt_count = 0
# Loop over trees
nev = t.GetEntries()
print 'Sample size is '+str(nev)

for i in range(nev):

#	if i%10 != 1 :
#		continue
	evt_count += 1
	if evt_count == cutoff:
		break
	# Report progress
	if evt_count%10000 == 1 :
		print 'finishing event # '+str(i)
		print 'Progress ' + str(100*evt_count/nev)+'%'
	t.GetEntry(i)

	njets_denom = 0
	njets_num = 0

	if t.nValidJets == 4:
		num_jets = 4
	if t.nValidJets == 5:
		num_jets = 5

	valid_jets = 0
	for j in range(len(t.pt)):
		j_csv = t.CSV[j]
		j_PID = abs(t.PID[j])
		j_pt  = t.pt[j]
		j_eta = math.fabs(t.eta[j])
		# Skip empty events
		if not ( j_pt > 0 and j_PID > 0 ):
			continue
		# debug
		if j_pt == 0 :
			print 'pt and PID '+ str(j_pt)+','+str(j_PID)
		jet_count += 1
		valid_jets += 1
		# Fill numerator
		if j_PID==5 and j_csv>CSVM:
			eff_num_list[0].Fill(j_pt,j_eta)
		if j_PID==4 and j_csv>CSVM:
			eff_num_list[1].Fill(j_pt,j_eta)
		if j_PID!=4 and j_PID!=5 and j_csv>CSVM:
			eff_num_list[2].Fill(j_pt,j_eta)
		# Fill denominator
		if j_PID==5:
			eff_denom_list[0].Fill(j_pt,j_eta)
		if j_PID==4:
			eff_denom_list[1].Fill(j_pt,j_eta)
		if j_PID!=5 and j_PID!=4:
			eff_denom_list[2].Fill(j_pt,j_eta)
	# print 'Number of valid jets is : '+str(valid_jets) # For debug

print 'Done Getting numerators and denominators!'
# Get efficiency
for i in range(len(eff_num_list)):
	h_tmp = eff_num_list[i].Clone()
	h_tmp.SetName(h_name_eff[i])
	h_tmp.Divide(eff_denom_list[i])
	eff_list.append(h_tmp)
print 'All done!'
print 'Total number of events : '+str(evt_count)
print 'Total number of Jets   : '+str(jet_count)
print 'Average number of jets per event is : '+'%.2f'%(float(jet_count)/float(evt_count))
# Close up
f.Write()
file_tmp.Close()
f.Close()
