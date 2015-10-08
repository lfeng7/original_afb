##############################################################################################################
##	  Data and Monte Carlo Tagger/Matcher for Top Quark Forward/Backward Asymmetry Study at CMS				##
## 									Morris Swartz, Raymond Feng, Nick Eminizer								##
##											   		Summer 2014												##
##############################################################################################################

####### Editing Log ##########
### Edited on 5/4/14 by Raymond Feng 
### Add option Pdf_Version (cteq66 or CT10). Adding PDF_weights to ntuples as a branch called Pdf_weights
### Based on newer version of Nick's tagger. Get on 5/16/14
### Edited on 5/16/14 by Raymond
### Add GenEventInfo, with qScale and x1,x2,id1,id2 for PDF reweighting later.
### Edited on 6/27/14 by Nick
### Style overhaul, finalize algorithm, add pileup reweighting, btag efficiency, and top pT reweighting
### Edited on 7/15/14 by Nick
### Fix bug in plotting pileup distributions, added sideband option, split dileptonic and hadronic cases
### Synced with Nick's version on 9/2014
### Edited on 10/11/14 by Raymond
### Add the b-tagging efficiency input file for madgraph sample. The correct efficiency input file will be chosen
### based on the sample_type. Add an alternative CSV distribution input file from madgraph sample. 
### Now all ttbar sample (signal or dilep/had bck) have PDF weight information for later use
##############################################################################################################
##								Imports, user options, global constants						 				##
##############################################################################################################
import os
import pickle
import glob
import math
import ROOT
from ROOT import *
import sys
from DataFormats.FWLite import Events, Handle
from optparse import OptionParser
from array import array

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--out', metavar='F', type='string', action='store',dest='out',help='') ## Sets filename for outputfile
parser.add_option('--mc', metavar='F', type='string', action='store',dest='mc',help='') ## Is the file a MC
parser.add_option('--lep', metavar='F', type='string', action='store',default='mu',dest='lep',help='') ## which leptons do we use
parser.add_option('--dR', metavar='F', type='float', action='store',default='0.3',dest='dR',help='') ## choose the dR requirement for matching
parser.add_option('--tMotherCheck', metavar='F', type='string', action='store',default='none',dest='tMotherCheck',help='') 
							## 'none' for no restrictions, 
							# 'qq' for qq (ttbar mother IDs < 6 and flavors match), 'gg' for everything not qqbar.
parser.add_option('--decayType', metavar='F', type='string', action='store',default='none',dest='decayType',help='') 
							## 'semilep' or 'semileptonic', 'dilep' or 'dileptonic', or 'had' or 'hadronic'
							## use only for ttbar samples
parser.add_option('--nMaxEvents', metavar='F', type='int', action='store',default='-1',dest='nMaxEvents',help='') ## choose the maximum number of iterations to go through
parser.add_option('--printEvery', metavar='F', type='int', action='store',default='1000',dest='printEvery',help='') ## choose the number of events to go through between prints
parser.add_option('--nTotalJobs', metavar='F', type='int', action='store',default='1',dest='nTotalJobs',help='') # Total number of jobs that will be run on the grid
parser.add_option('--iJob', metavar='F', type='int', action='store',default='0',dest='iJob',help='') # Offset from 0 for this job on the grid
parser.add_option('--leptonCut', metavar='F', type='string', action='store',default='tight',dest='leptonCut',help='') # Use tight or loose leptons
parser.add_option('--sample_type', metavar='F', type='string', action='store',default='signal',dest='sample_type',help='') #Type of sample, needed to apply btagging efficiency
parser.add_option('--sideband', metavar='F', type='string', action='store',default='no',dest='sideband',help='') #Whether to look in the isolation sideband
parser.add_option('--PdfVersion', metavar='F', type='string', action='store',default='none',dest='PdfVersion',help='') # Options for different PDF set, only meaningful for signal MC
(options, args) = parser.parse_args()
if options.iJob < 0 or options.iJob >= options.nTotalJobs :
	print 'CHECK NUMBERING CONVENTION FOR GRID JOBS: SOMETHING IS VERY WRONG AND NOTHING IS GOOD!'

#Global constants
MW = 80.4 
MT = 173.3
if options.mc == 'yes' :
	MT = 172.5
QW = (MW*MW)/(MT*MT)
ZW = (2.0*2.0)/(MT*MT)
ZT = (1.4*1.4)/(MT*MT)
CDFT  = (math.acos(0.)+math.atan(1./sqrt(ZT)))/sqrt(ZT)
CDFW  = 0.5+2.*QW+(1.5*QW*QW-0.5*ZW*QW-1.5)*math.log(((1.-QW)*(1.-QW)+ZW*QW)/(QW*QW+QW*ZW))
CDFW += ((QW*QW*QW-3.*ZW*QW*QW-3.*QW+2.)/math.sqrt(ZW*QW))*(math.atan((1.-QW)/math.sqrt(ZW*QW))+math.atan(QW/math.sqrt(ZW*QW)))
SIGMAJ = 0.10
SIGMAL = 0.03
MMUON     = 0.105658
MELECTRON = 0.000511
if options.lep == 'mu' :
	MLEP = MMUON
elif options.lep == 'el' :
	MLEP = MELECTRON
# Global vars
count_out_of_bounds = 0 
count_jets = 0
# btagging efficiency constants for CSVM
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
	0.0598311]
#tracking efficiency constants and errors
tracking_eff_consts = [[-2.4,-2.1,0.9869,0.07],
					   [-2.1,-1.6,0.9948,0.02],
					   [-1.6,-1.2,0.9967,0.02],
					   [-1.2,-0.9,0.9974,0.02],
					   [-0.9,-0.6,0.9980,0.01],
					   [-0.6,-0.3,0.9980,0.01],
					   [-0.3,-0.2,0.9972,0.02],
					   [-0.2, 0.2,0.9963,0.01],
					   [ 0.2, 0.3,0.9978,0.02],
					   [ 0.3, 0.6,0.9977,0.01],
					   [ 0.6, 0.9,0.9976,0.01],
					   [ 0.9, 1.2,0.9968,0.02],
					   [ 1.2, 1.6,0.9959,0.03],
					   [ 1.6, 2.1,0.9970,0.02],
					   [ 2.1, 2.4,0.9836,0.08]]
#unpickle the lepton ID, isolation, and trigger files
#grid filenames
muon_id_filename  = './tardir/MuonEfficiencies_Run2012ReReco_53X.pkl'
muon_iso_filename = './tardir/MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl'
muon_trigger_filename = './tardir/SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.pkl'
input_files_list = open('./tardir/input_files_list.txt','r')
# #local filenames
# muon_id_filename  = '../MuonEfficiencies_Run2012ReReco_53X.pkl'
# muon_iso_filename = '../MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl'
# muon_trigger_filename = '../SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.pkl'
# input_files_list = open('./input_files_list.txt','r')

muon_id_file  = open(muon_id_filename)
muon_iso_file = open(muon_iso_filename)
muon_trigger_file = open(muon_trigger_filename)
muon_id_dict  = pickle.load(muon_id_file)
muon_iso_dict = pickle.load(muon_iso_file)
muon_trigger_dict = pickle.load(muon_trigger_file)

##############################################################################################################
##										Set up output file structure						 				##
##############################################################################################################

#kinematics histogram declarations
CUTFLOW = 0 #0=no cuts except that the event requires a lepton, met, and at least 4 jets
			#1=selected lepton
			#2=selected MET
			#3=passed kinematic jet cuts with 4 or 5 jets
			#4=4 hard enough jets
			#5=sufficient number of btags (ready to be reconstructed)
nSteps = 6
#leptons
lpts   = []
letas  = []
pfisos = []
#MET
metpts = []
#Jets
hardestjetpts = []
secondhardestjetpts = []
thirdhardestjetpts  = []
fourthhardestjetpts = []
fifthhardestjetpts  = []
hardestjetetas = []
secondhardestjetetas = []
thirdhardestjetetas  = []
fourthhardestjetetas = []
fifthhardestjetetas  = []
hardestjetmasss = []
secondhardestjetmasss = []
thirdhardestjetmasss  = []
fourthhardestjetmasss = []
fifthhardestjetmasss  = []
hardestjetcsvs = []
secondhardestjetcsvs = []
thirdhardestjetcsvs  = []
fourthhardestjetcsvs = []
fifthhardestjetcsvs  = []
combttbarmasss = []
for i in range(nSteps) :
	lpts.append(TH1F(  'lpt'+str(i),  'Lepton Pt at step '+str(i)+'; Pt (GeV)',35,0,350))
	letas.append(TH1F( 'leta'+str(i), 'Lepton Eta at step '+str(i)+'; #eta',   50,-2.5,2.5))
	pfisos.append(TH1F('pfiso'+str(i),'Lepton PFiso at step '+str(i)+'; pfiso',25,0.,0.14))
	metpts.append(TH1F('metpt'+str(i),'MET Pt at step '+str(i)+'; Pt (GeV)',   35,0,350))
	hardestjetpts.append(TH1F('hardestjetpt'+str(i),'Pt of 1st hardest jet at step '+str(i)+'; Pt (GeV)',35,0,350))
	secondhardestjetpts.append(TH1F('secondhardestjetpt'+str(i),'Pt of 2nd hardest jet at step '+str(i)+'; Pt (GeV)',35,0,350))
	thirdhardestjetpts.append(TH1F( 'thirdhardestjetpt'+str(i), 'Pt of 3rd hardest jet at step '+str(i)+'; Pt (GeV)',35,0,350))
	fourthhardestjetpts.append(TH1F('fourthhardestjetpt'+str(i),'Pt of 4th hardest jet at step '+str(i)+'; Pt (GeV)',35,0,350))
	fifthhardestjetpts.append(TH1F( 'fifthhardestjetpt'+str(i), 'Pt of 5th hardest jet at step '+str(i)+'; Pt (GeV)',35,0,350))
	hardestjetetas.append(TH1F('hardestjeteta'+str(i),'Eta of 1st hardest jet at step '+str(i)+'; #eta',50,-2.5,2.5))
	secondhardestjetetas.append(TH1F('secondhardestjeteta'+str(i),'Eta of 2nd hardest jet at step '+str(i)+'; #eta',64,-3.2,3.2))
	thirdhardestjetetas.append(TH1F( 'thirdhardestjeteta'+str(i), 'Eta of 3rd hardest jet at step '+str(i)+'; #eta',64,-3.2,3.2))
	fourthhardestjetetas.append(TH1F('fourthhardestjeteta'+str(i),'Eta of 4th hardest jet at step '+str(i)+'; #eta',64,-3.2,3.2))
	fifthhardestjetetas.append(TH1F( 'fifthhardestjeteta'+str(i), 'Eta of 5th hardest jet at step '+str(i)+'; #eta',64,-3.2,3.2))
	hardestjetmasss.append(TH1F('hardestjetmass'+str(i),'Mass of 1st hardest jet at step '+str(i)+'; Mass (GeV)',30,0.0,60.0))
	secondhardestjetmasss.append(TH1F('secondhardestjetmass'+str(i),'Mass of 2nd hardest jet at step '+str(i)+'; Mass (GeV)',30,0.0,60.0))
	thirdhardestjetmasss.append(TH1F( 'thirdhardestjetmass'+str(i), 'Mass of 3rd hardest jet at step '+str(i)+'; Mass (GeV)',30,0.0,60.0))
	fourthhardestjetmasss.append(TH1F('fourthhardestjetmass'+str(i),'Mass of 4th hardest jet at step '+str(i)+'; Mass (GeV)',30,0.0,60.0))
	fifthhardestjetmasss.append(TH1F( 'fifthhardestjetmass'+str(i), 'Mass of 5th hardest jet at step '+str(i)+'; Mass (GeV)',30,0.0,60.0))
	hardestjetcsvs.append(TH1F('hardestjetcsv'+str(i),'CSV of 1st hardest jet at step '+str(i)+'; CSV value',50,0.0,1.0))
	secondhardestjetcsvs.append(TH1F('secondhardestjetcsv'+str(i),'CSV of 2nd hardest jet at step '+str(i)+'; CSV value',50,0.0,1.0))
	thirdhardestjetcsvs.append(TH1F( 'thirdhardestjetcsv'+str(i), 'CSV of 3rd hardest jet at step '+str(i)+'; CSV value',50,0.0,1.0))
	fourthhardestjetcsvs.append(TH1F('fourthhardestjetcsv'+str(i),'CSV of 4th hardest jet at step '+str(i)+'; CSV value',50,0.0,1.0))
	fifthhardestjetcsvs.append(TH1F( 'fifthhardestjetcsv'+str(i), 'CSV of 5th hardest jet at step '+str(i)+'; CSV value',50,0.0,1.0))
	combttbarmasss.append(TH1F('combttbarmass'+str(i),'Combined t#bar{t} Mass at step '+str(i)+'; M_{t#bar{t}} (GeV)',50,300,1600))
leptonHistoLists = [lpts,letas,pfisos]
metHistoLists    = [metpts]
hardestjetHistoLists = [hardestjetpts,hardestjetetas,hardestjetmasss,hardestjetcsvs]
secondhardestjetHistoLists = [secondhardestjetpts,secondhardestjetetas,secondhardestjetmasss,secondhardestjetcsvs]
thirdhardestjetHistoLists  = [thirdhardestjetpts, thirdhardestjetetas, thirdhardestjetmasss, thirdhardestjetcsvs ]
fourthhardestjetHistoLists = [fourthhardestjetpts,fourthhardestjetetas,fourthhardestjetmasss,fourthhardestjetcsvs]
fifthhardestjetHistoLists  = [fifthhardestjetpts, fifthhardestjetetas, fifthhardestjetmasss, fifthhardestjetcsvs ]
#final chis and pileup
finalChiHist = TH1F('finalChiHist','Final #Chi^{2} (-2ln(L)) of fully reconstructed events; #Chi^{2}',150,-50.,100.)
pileupHist = TH1F('pileupHist','Pileup Events; npv',80,0,80)
allHistos = []
for i in range(nSteps) :
	for histolist in leptonHistoLists :
		allHistos.append(histolist[i])
	for histolist in metHistoLists :
		allHistos.append(histolist[i])
	for histolist in hardestjetHistoLists :
		allHistos.append(histolist[i])
	for histolist in secondhardestjetHistoLists :
		allHistos.append(histolist[i])
	for histolist in thirdhardestjetHistoLists :
		allHistos.append(histolist[i])
	for histolist in fourthhardestjetHistoLists :
		allHistos.append(histolist[i])
	for histolist in fifthhardestjetHistoLists :
		allHistos.append(histolist[i])
	allHistos.append(combttbarmasss[i])
allHistos.append(finalChiHist)
allHistos.append(pileupHist)
for histo in allHistos :
	histo.SetDirectory(0)

#output file branches
outputname = options.out+'.root'
f = TFile(outputname,'Recreate')
t = TTree('output','output')
t.SetDirectory(0)
#branches for reconstructed quantities
nParticles = 8 #lepton, neutrino, leptonic b, leptonic W, hadronic b, and hadronic W, plus the leptonic and hadronic t/tbar
pt   = array('f',nParticles*[0.])
eta  = array('f',nParticles*[0.])
phi  = array('f',nParticles*[0.])
mass = array('f',nParticles*[0.])
px   = array('f',nParticles*[0.])
py   = array('f',nParticles*[0.])
pz   = array('f',nParticles*[0.])
E    = array('f',nParticles*[0.])
PID  = array('f',nParticles*[0])
finalChi   = array('f',nParticles*[0.])
nbTags     = array('I',[0])
nValidJets = array('I',[0])
bestFitParValues = array('f',6*[0.])
motherParticles  = array('i',2*[0])
is_leptonic_side = array('I',nParticles*[0])
pileup_events = array('f',[0.])
t.Branch('pt',  pt,  'pt[8]/F')
t.Branch('eta', eta, 'eta[8]/F')
t.Branch('phi', phi, 'phi[8]/F')
t.Branch('mass',mass,'mass[8]/F')
t.Branch('px',  px,  'px[8]/F')
t.Branch('py',  py,  'py[8]/F')
t.Branch('pz',  pz,  'pz[8]/F')
t.Branch('E',   E,   'E[8]/F')
t.Branch('PID', PID, 'PID[8]/F')
t.Branch('finalChi',  finalChi,  'finalChi[8]/F')
t.Branch('nbTags',    nbTags,    'nbTags/i')
t.Branch('nValidJets',nValidJets,'nValidJets/i')
t.Branch('bestFitParValues',bestFitParValues,'bestFitParValues[6]/F')
t.Branch('motherParticles', motherParticles,'motherParticles[2]/I')
t.Branch('is_leptonic_side',is_leptonic_side,'is_leptonic_side[8]/i')
t.Branch('pileup_events',pileup_events,'pileup_events/F')
#branches for MC truth quantities				  
mc_pt   = array('f',(nParticles+2)*[0.])
mc_eta  = array('f',(nParticles+2)*[0.])
mc_phi  = array('f',(nParticles+2)*[0.])
mc_mass = array('f',(nParticles+2)*[0.])
mc_px   = array('f',(nParticles+2)*[0.])
mc_py   = array('f',(nParticles+2)*[0.])
mc_pz   = array('f',(nParticles+2)*[0.])
mc_E    = array('f',(nParticles+2)*[0.])
weight_GJR_scale     = array('f',[0.])
weight_CT10_scale    = array('f',[0.])
weight_cteq_scale    = array('f',[0.])
weight_top_pT        = array('f',[0.])
weight_btag_eff      = array('f',[0.])
weight_btag_eff_err  = array('f',[0.])
weight_pileup        = array('f',[0.])
weight_tracking      = array('f',[0.])
weight_tracking_low  = array('f',[0.])
weight_tracking_hi   = array('f',[0.])
weight_lep_ID        = array('f',[0.])
weight_lep_ID_low    = array('f',[0.])
weight_lep_ID_hi     = array('f',[0.])
weight_lep_iso       = array('f',[0.])
weight_lep_iso_low   = array('f',[0.])
weight_lep_iso_hi    = array('f',[0.])
weight_trig_eff      = array('f',[0.])
weight_trig_eff_low  = array('f',[0.])
weight_trig_eff_hi   = array('f',[0.])
weight_for_histos    = array('f',[0.])
mc_pileup_events = array('f',[0.])
mc_was_matched = array('I',nParticles*[0])
t.Branch('mc_pt',  mc_pt,  'mc_pt[10]/F')
t.Branch('mc_eta', mc_eta, 'mc_eta[10]/F')
t.Branch('mc_phi', mc_phi, 'mc_phi[10]/F')
t.Branch('mc_mass',mc_mass,'mc_mass[10]/F')
t.Branch('mc_px',  mc_px,  'mc_px[10]/F')
t.Branch('mc_py',  mc_py,  'mc_py[10]/F')
t.Branch('mc_pz',  mc_pz,  'mc_pz[10]/F')
t.Branch('mc_E',   mc_E,   'mc_E[10]/F')
t.Branch('weight_GJR_scale',    weight_GJR_scale,    'weight_GJR_scale/F')
t.Branch('weight_CT10_scale',   weight_CT10_scale,   'weight_CT10_scale/F')
t.Branch('weight_cteq_scale',   weight_cteq_scale,   'weight_cteq_scale/F')
t.Branch('weight_top_pT',	    weight_top_pT,       'weight_top_pT/F')
t.Branch('weight_btag_eff',     weight_btag_eff,     'weight_btag_eff/F')
t.Branch('weight_btag_eff_err', weight_btag_eff_err, 'weight_btag_eff_err/F')
t.Branch('weight_pileup',       weight_pileup,       'weight_pileup/F')
t.Branch('weight_tracking',     weight_tracking,     'weight_tracking/F')
t.Branch('weight_tracking_low', weight_tracking_low, 'weight_tracking_low/F')
t.Branch('weight_tracking_hi',  weight_tracking_hi,  'weight_tracking_hi/F')
t.Branch('weight_lep_ID',       weight_lep_ID,       'weight_lep_ID/F')
t.Branch('weight_lep_ID_low',   weight_lep_ID_low,   'weight_lep_ID_low/F')
t.Branch('weight_lep_ID_hi',    weight_lep_ID_hi,    'weight_lep_ID_hi/F')
t.Branch('weight_lep_iso',      weight_lep_iso,      'weight_lep_iso/F')
t.Branch('weight_lep_iso_low',  weight_lep_iso_low,  'weight_lep_iso_low/F')
t.Branch('weight_lep_iso_hi',   weight_lep_iso_hi,   'weight_lep_iso_hi/F')
t.Branch('weight_trig_eff',     weight_trig_eff,     'weight_trig_eff/F')
t.Branch('weight_trig_eff_low', weight_trig_eff_low, 'weight_trig_eff_low/F')
t.Branch('weight_trig_eff_hi',  weight_trig_eff_hi,  'weight_trig_eff_hi/F')
t.Branch('weight_for_histos',   weight_for_histos,   'weight_for_histos/F')
t.Branch('mc_pileup_events',mc_pileup_events,'mc_pileup_events/F')
t.Branch('mc_was_matched',mc_was_matched,'mc_was_matched[8]/i')
#generated PDF set
Pdf_weights = array('d',100*[0.])
t.Branch('Pdf_weights',Pdf_weights,'Pdf_weights[100]/D')
#generator info
qscale    = array('f',[0.])
parton_id = array('f',2*[0.])
parton_x  = array('f',2*[0.])
t.Branch('qScale',   qscale,   'qscale/F')
t.Branch('parton_id',parton_id,'parton_id[2]/F')
t.Branch('parton_x', parton_x, 'parton_x[2]/F')

f.cd()

##############################################################################################################
##							Handle input files for kinematic fit and data						 			##
##############################################################################################################

#Set up file and read histograms used for CSV info in kinematic fit.
f1 = TFile('/uscms_data/d3/eminizer/CMSSW_5_3_11/src/Analysis/EDSHyFT/test/tagging/dumped_Powheg_TT.root')  # The CSV distribution input files for Powheg signal
#f1 = TFile('/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/codes/tagging/tagged_files_9_14/dump_TTJets_SemiLep.root') # The CSV input for Madgraph signals
bdisc = f1.Get('bDisc')
Wsubdisc = f1.Get('WsubDisc')
otherdisc = f1.Get('otherDisc')
bdiscInt = bdisc.Integral()
WsubdiscInt = Wsubdisc.Integral()
otherdiscInt = otherdisc.Integral()

#Set up btag efficiency files
if options.sample_type == 'signal' or options.sample_type == 'Powheg_dilep_TT' or options.sample_type == 'Powheg_had_TT' or options.sample_type == 'Powheg_qq_semilep_TT' or options.sample_type == 'Powheg_gg_semilep_TT' :
	F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/TT_CT10_TuneZ2star_8TeV-powheg-tauola_AK5PF_CSVM_bTaggingEfficiencyMap.root'
elif options.sample_type =='TTJets_SemiLep' or options.sample_type == 'TTJets_Hadronic' or options.sample_type == 'TTJets_FullLep' :
	F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/codes/tagging/btagging_efficiency_files/TTJets_SemiLep_all_efficiency.root'
elif options.sample_type == 'WJets' :
	F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_AK5PF_CSVM_bTaggingEfficiencyMap.root'
elif options.sample_type == 'DYJets' :
	F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_AK5PF_CSVM_bTaggingEfficiencyMap.root'
elif options.sample_type == 'T_s-channel' :
	F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/T_s-channel_TuneZ2star_8TeV-powheg-tauola_AK5PF_CSVM_bTaggingEfficiencyMap.root'
elif options.sample_type == 'T_t-channel' :
	F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/T_t-channel_TuneZ2star_8TeV-powheg-tauola_AK5PF_CSVM_bTaggingEfficiencyMap.root'
elif options.sample_type == 'T_tW-channel' :
	F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_AK5PF_CSVM_bTaggingEfficiencyMap.root'
elif options.sample_type == 'Tbar_s-channel' :
	F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_AK5PF_CSVM_bTaggingEfficiencyMap.root'
elif options.sample_type == 'Tbar_t-channel' :
	F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_AK5PF_CSVM_bTaggingEfficiencyMap.root'
elif options.sample_type == 'Tbar_tW-channel' :
	F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_AK5PF_CSVM_bTaggingEfficiencyMap.root'
file_tmp = TFile(F_eff)
efficiency_b = file_tmp.Get('efficiency_b')
efficiency_c = file_tmp.Get('efficiency_c')
efficiency_udsg = file_tmp.Get('efficiency_udsg') 

#Set up pileup distribution files used in pileup reweighting
data_pufile = TFile('/uscms_data/d3/eminizer/CMSSW_5_3_11/src/Analysis/EDSHyFT/test/tagging/data_pileup_distribution.root')
data_pu_dist = data_pufile.Get('pileup')
MC_pu_dist   = f1.Get('pileup')
data_pu_dist.Scale(1.0/data_pu_dist.Integral())
MC_pu_dist.Scale(1.0/MC_pu_dist.Integral())
#Width and Central Value of Wmass peak
wmass  = 80.08
wwidth = 33.5

#Read in input data files
files = []
print 'Getting these files: '
for input_file in input_files_list :
	print input_file.rstrip()
	files.append(input_file.rstrip())
events = Events(files)

#Set up handles and labels for use
#muon variables
if options.lep == 'mu':
	if options.leptonCut == 'tight' :
		leptonPtHandle     = Handle('std::vector<float>')
		leptonEtaHandle    = Handle('std::vector<float>')
		leptonPhiHandle    = Handle('std::vector<float>')
		leptonPfisoHandle  = Handle('std::vector<float>')
		leptonChargeHandle = Handle('std::vector<float>')
		leptonPtLabel      = ('pfShyftTupleMuons','pt')
		leptonEtaLabel     = ('pfShyftTupleMuons','eta' )
		leptonPhiLabel     = ('pfShyftTupleMuons','phi' )
		leptonPfisoLabel   = ('pfShyftTupleMuons','pfisoPU' )
		leptonChargeLabel  = ('pfShyftTupleMuons','charge' )
		otherleptonPtHandle     = Handle('std::vector<float>')
		otherleptonEtaHandle    = Handle('std::vector<float>')
		otherleptonPtLabel      = ('pfShyftTupleElectrons','pt')
		otherleptonEtaLabel     = ('pfShyftTupleElectrons','eta')
	elif options.leptonCut == 'loose' :
		leptonPtHandle     = Handle('std::vector<float>')
		leptonEtaHandle	   = Handle('std::vector<float>')
		leptonPhiHandle    = Handle('std::vector<float>')
		leptonPfisoHandle  = Handle('std::vector<float>')
		leptonChargeHandle = Handle('std::vector<float>')
		leptonPtLabel	   = ('pfShyftTupleMuonsLoose','pt')
		leptonEtaLabel	   = ('pfShyftTupleMuonsLoose','eta')
		leptonPhiLabel	   = ('pfShyftTupleMuonsLoose','phi')
		leptonPfisoLabel   = ('pfShyftTupleMuonsLoose','pfisoPU')
		leptonChargeLabel  = ('pfShyftTupleMuonsLoose','charge')
		otherleptonPtHandle     = Handle('std::vector<float>')
		otherleptonEtaHandle    = Handle('std::vector<float>')
		otherleptonPtLabel      = ('pfShyftTupleElectronsLoose','pt')
		otherleptonEtaLabel     = ('pfShyftTupleElectronsLoose','eta')
	MLEP = MMUON
#electron variables
if options.lep == 'el':
	if options.leptonCut == 'tight' :
		leptonPtHandle     = Handle('std::vector<float>')
		leptonEtaHandle    = Handle('std::vector<float>')
		leptonPhiHandle    = Handle('std::vector<float>')
		leptonPfisoHandle  = Handle('std::vector<float>')
		leptonChargeHandle = Handle('std::vector<float>')
		leptonPtLabel      = ('pfShyftTupleElectrons','pt')
		leptonEtaLabel     = ('pfShyftTupleElectrons','eta')
		leptonPhiLabel     = ('pfShyftTupleElectrons','phi')
		leptonPfisoLabel   = ('pfShyftTupleElectrons','pfiso')
		leptonChargeLabel  = ('pfShyftTupleElectrons','charge')
		otherleptonPtHandle     = Handle('std::vector<float>')
		otherleptonEtaHandle    = Handle('std::vector<float>')
		otherleptonPtLabel      = ('pfShyftTupleMuons','pt')
		otherleptonEtaLabel     = ('pfShyftTupleMuons','eta' )
	elif options.leptonCut == 'loose' :
		leptonPtHandle     = Handle('std::vector<float>')
		leptonEtaHandle    = Handle('std::vector<float>')
		leptonPhiHandle    = Handle('std::vector<float>')
		leptonPfisoHandle  = Handle('std::vector<float>')
		leptonChargeHandle = Handle('std::vector<float>')
		leptonPtLabel      = ('pfShyftTupleElectronsLoose','pt')
		leptonEtaLabel     = ('pfShyftTupleElectronsLoose','eta')
		leptonPhiLabel     = ('pfShyftTupleElectronsLoose','phi')
		leptonPfisoLabel   = ('pfShyftTupleElectronsLoose','pfiso')
		leptonChargeLabel  = ('pfShyftTupleElectronsLoose','charge')
		otherleptonPtHandle     = Handle('std::vector<float>')
		otherleptonEtaHandle	= Handle('std::vector<float>')
		otherleptonPtLabel	    = ('pfShyftTupleMuonsLoose','pt')
		otherleptonEtaLabel	    = ('pfShyftTupleMuonsLoose','eta')
	MLEP = MELECTRON
#neutrino (MET) variables
if options.leptonCut == 'tight' :
	metPtHandle  = Handle('std::vector<float>')
	metPhiHandle = Handle('std::vector<float>')
	metPtLabel   = ('pfShyftTupleMET','pt')
	metPhiLabel  = ('pfShyftTupleMET','phi')
elif options.leptonCut == 'loose' :
	metPtHandle  = Handle('std::vector<float>')
	metPhiHandle = Handle('std::vector<float>')
	metPtLabel   = ('pfShyftTupleMETLoose','pt')
	metPhiLabel  = ('pfShyftTupleMETLoose','phi')
#AK5 jet variables
jetPtHandle   = Handle('std::vector<float>')
jetEtaHandle  = Handle('std::vector<float>')
jetPhiHandle  = Handle('std::vector<float>')
jetMassHandle = Handle('std::vector<float>')
jetCSVHandle  = Handle('std::vector<float>')
jetPtLabel    = ('pfShyftTupleJets','pt')
jetEtaLabel   = ('pfShyftTupleJets','eta')
jetPhiLabel   = ('pfShyftTupleJets','phi')
jetMassLabel  = ('pfShyftTupleJets','mass')
jetCSVLabel   = ('pfShyftTupleJets','csv')
#pileup for data
dataPileupHandle = Handle('unsigned int')
dataPileupLabel  = ('pileup','npv')
#MC GenParticle variables
if options.mc == 'yes':
	GenHandle         = Handle('vector<reco::GenParticle>')	#GenParticles strcture
	GenInfoHandle     = Handle('GenEventInfoProduct') 		# GenInfo, with x1 x2 and Q^2
	PdfHandle         = Handle('vector<double>')			#Pdf weights
	PdfHandle_CT10    = Handle('vector<double>')
	PdfHandle_cteq    = Handle('vector<double>')
	PdfHandle_GJR     = Handle('vector<double>')
	npvRealTrueHandle = Handle('int')
	jetFlavorHandle   = Handle('vector<float>')
	GenLabel          = ('prunedGenParticles','')
	GenInfoLabel      = ('generator','')
	PdfLabel_cteq     = ('pdfWeights','cteq66')
	PdfLabel_CT10     = ('pdfWeights','CT10')
	PdfLabel_GJR      = ('pdfWeights','GJR08VFnloE')
	npvRealTrueLabel  = ('pileup','npvRealTrue')
	jetFlavorLabel = ('pfShyftTupleJets','jetFlavor')


##############################################################################################################
##											Accessory Functions									 			##
##############################################################################################################

#btagging efficiency
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
# //only for CSVM point 
def get_SF_btag(ptJet,etaJet,flavJet):		# checked
	global count_out_of_bounds
	x = ptJet # the pt of the jet 
	eta = fabs(etaJet); # abs(eta) 
	result = []			# the first is SF, second is SF error
	if(eta>2.4) :
		print 'warning SF_btag_eta>2.4 ??  ' +str(eta)+'' 
		result.append(1)
		result.append(0)
		return result
	if(x<20) : x=20; 
	if(x>800) : x= 800;
	if(abs(flavJet)==5 or abs(flavJet) == 4) : # for b or c. flavJet[indj] refers to the MC-true flavor of the jet
		SF  = (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x)) 
		ptbin = getptbin_for_btag( ptJet ) #--> ptJet[indj] refers to the pt of this jet
		SFerr = SFb_error[ptbin]
		if(x>800 or x<20 ) : SFerr *= 2;
		if(abs(flavJet) == 4 ) : SFerr *= 2; 
	else : #SFlight
		if(eta<=0.8) :
			SF = ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)))
			max = ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)))
			min = ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)))		
		elif(eta<=1.6) :
			SF = ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)))
			max = ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)))
			min = ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)))
		elif (eta>1.6 and eta<=2.4) :
			SF = ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)))
			max = ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)))
			min = ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)))
		if math.fabs(max-SF)>math.fabs(min-SF) :
			SFerr = math.fabs(max-SF)
		else :
			SFerr = math.fabs(min-SF)
	result.append(SF)
	result.append(SFerr)
	return result
# ///this is the functin to call to get the event-by-event b-tag weight
# ///as implemented, it only works for CSVM. 
# Input : given a list of selected jets (with pt, ID, etc)
def get_weight_btag(jet_Pts, jet_etas, jet_flavors,jet_csvs) :
	global count_jets,count_out_of_bounds
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
				err4 += (-eff*jet_SFerr)/(1-eff*jet_SF);  # correlated for light
	# Check if any of the probability is zero. If so, then set weight to be 1, aka do nothing to this event
	if dataTag*dataNoTag == 0  or mcTag*mcNoTag == 0:
		print 'one of the probability is zero!'
		wtbtag = 1
		wtbtagErr = 0
	else :
		# Get event weight for this event
		wtbtag = (dataNoTag * dataTag ) / ( mcNoTag * mcTag ); 
		wtbtagErr = math.sqrt( pow(err1+err2,2) + pow( err3 + err4,2)) * wtbtag; # un-correlated for b/c and light
	# # debug
	# print 'wtbtag = '+str(wtbtag)+'+/-'+str(wtbtagErr)
	# Set return
	to_return = []
	to_return.append(wtbtag)
	to_return.append(wtbtagErr)
	return to_return

#four vector constrained rescaling function for use in kinematic fit
def pscale(scalefactor, vector) :
	p2 = scalefactor*scalefactor*(vector.Px()*vector.Px()+vector.Py()*vector.Py()+vector.Pz()*vector.Pz())
	m2 = vector.M2()
	newE = math.sqrt(p2+m2)
	return ROOT.TLorentzVector(scalefactor*vector.Px(),scalefactor*vector.Py(),scalefactor*vector.Pz(),newE)

#Minimization function for kinematic fitting
minf = 1000000000.
def fcn(npar, deriv, f, par, flag) :
	global minf
	lnL = 0.0
	l = pscale(par[1], lepton)
	bl = pscale(par[2], blep)
	bh = pscale(par[3], bhad)
	whs1 = pscale(par[4],Wsub1)
	whs2 = pscale(par[5],Wsub2)
	wh = whs1+whs2
	newmetx = met.Px()+(1.0-par[1])*lepton.Px()+(1.0-par[2])*blep.Px()+(1.0-par[3])*bhad.Px()+(1.0-par[4])*Wsub1.Px()+(1.0-par[5])*Wsub2.Px()
	newmety = met.Py()+(1.0-par[1])*lepton.Py()+(1.0-par[2])*blep.Py()+(1.0-par[3])*bhad.Py()+(1.0-par[4])*Wsub1.Py()+(1.0-par[5])*Wsub2.Py()
	v = pscale(1.0,met)
	v.SetPx(newmetx)
	v.SetPy(newmety)
	v.SetPz(par[0])
	v.SetE(math.sqrt(v.Px()*v.Px()+v.Py()*v.Py()+v.Pz()*v.Pz()))
	#print ' par0='+str(par[0])+' par1='+str(par[1])+' par2='+str(par[2])+' par3='+str(par[3])+' par4='+str(par[4])+' par5='+str(par[5])+''
	wl = v + l
	tl = wl + bl
	th = wh + bh
	mwl2 = wl.M2()
	mtl2 = tl.M2()
	mwh2 = wh.M2()
	mth2 = th.M2()
	ql = mwl2/(MT*MT)
	xl = mtl2/(MT*MT)
	qh = mwh2/(MT*MT)
	xh = mth2/(MT*MT)
	pdftl = 1./((xl - 1.)*(xl - 1.) + ZT)
	pdfth = 1./((xh - 1.)*(xh - 1.) + ZT)
	pdfwl = (1. - ql)*(1. - ql)*(2. + ql)/((ql-QW)*(ql-QW)+ZW*QW)
	pdfwh = (1. - qh)*(1. - qh)*(2. + qh)/((qh-QW)*(qh-QW)+ZW*QW)
	pdf = pdftl*pdfth*pdfwl*pdfwh/(CDFT*CDFT*CDFW*CDFW)
	if pdf > 0.0 :
		lnL += math.log(pdf)    #need positive f
	else :
		print('WARNING -- pdf is negative!!!')
		pdf = 1.e-50
		lnL += math.log(pdf)
	pdiscblep = bdisc.GetBinContent(bdisc.FindFixBin(blepCSV))/bdiscInt
	pdiscbhad = bdisc.GetBinContent(bdisc.FindFixBin(bhadCSV))/bdiscInt
	if WCSV1 == -1.0 :
		pdiscw1 = 1.0
	else :
		pdiscw1 = Wsubdisc.GetBinContent(Wsubdisc.FindFixBin(WCSV1))/WsubdiscInt
	if WCSV2 == -1.0 :
		pdiscw2 = 1.0
	else :
		pdiscw2 = Wsubdisc.GetBinContent(Wsubdisc.FindFixBin(WCSV2))/WsubdiscInt
	if extraCSV == -1.0 :
		pdiscextra = 1.0
	else :
		pdiscextra = otherdisc.GetBinContent(otherdisc.FindFixBin(extraCSV))/otherdiscInt
	if (pdiscblep*pdiscbhad*pdiscw1*pdiscw2*pdiscextra) <= 0 :
		print('WARNING -- problem with discriminant values!!!')
		print 'Values are: '+str(blepCSV)+' '+str(bhadCSV)+' '+str(WCSV1)+' '+str(WCSV2)+''
		pdf = 1.e-50
		f[0] += math.log(pdf)
	else :
		f[0] = -2.0*lnL + (par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL) + (par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
		f[0] = f[0] + (par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ) + (par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) + (par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ)
		f[0] = f[0] - 2.0*math.log(pdiscblep*pdiscbhad*pdiscw1*pdiscw2*pdiscextra)
	if f[0] < minf :
		minf = f[0]

#Event closeout function to fill histograms up to cutflow point and fill tree
def endEvent(cutflow) :
	weight_for_histos[0] = weight_pileup[0]*weight_btag_eff[0]*weight_top_pT[0]*weight_tracking[0]*weight_lep_ID[0]*weight_lep_iso[0]*weight_trig_eff[0]
	if options.sideband == 'yes' :
		weight_for_histos[0] = weight_pileup[0]*weight_btag_eff[0]*weight_top_pT[0]*weight_tracking[0]*weight_lep_ID[0]*weight_trig_eff[0]
	#put histogram variables in histograms up through cutflow number
	for i in range(cutflow+1) :
		lpts[i].Fill(lpt,weight_for_histos[0])
		letas[i].Fill(leta,weight_for_histos[0])
		pfisos[i].Fill(pfiso,weight_for_histos[0])
		metpts[i].Fill(metpt,weight_for_histos[0])
		hardestjetpts[i].Fill(hardestjetpt,weight_for_histos[0])
		secondhardestjetpts[i].Fill(secondhardestjetpt,weight_for_histos[0])
		thirdhardestjetpts[i].Fill(thirdhardestjetpt,weight_for_histos[0])
		fourthhardestjetpts[i].Fill(fourthhardestjetpt,weight_for_histos[0])
		fifthhardestjetpts[i].Fill(fifthhardestjetpt,weight_for_histos[0])
		hardestjetetas[i].Fill(hardestjeteta,weight_for_histos[0])
		secondhardestjetetas[i].Fill(secondhardestjeteta,weight_for_histos[0])
		thirdhardestjetetas[i].Fill(thirdhardestjeteta,weight_for_histos[0])
		fourthhardestjetetas[i].Fill(fourthhardestjeteta,weight_for_histos[0])
		fifthhardestjetetas[i].Fill(fifthhardestjeteta,weight_for_histos[0])
		hardestjetmasss[i].Fill(hardestjetmass,weight_for_histos[0])
		secondhardestjetmasss[i].Fill(secondhardestjetmass,weight_for_histos[0])
		thirdhardestjetmasss[i].Fill(thirdhardestjetmass,weight_for_histos[0])
		fourthhardestjetmasss[i].Fill(fourthhardestjetmass,weight_for_histos[0])
		fifthhardestjetmasss[i].Fill(fifthhardestjetmass,weight_for_histos[0])
		hardestjetcsvs[i].Fill(hardestjetcsv,weight_for_histos[0])
		secondhardestjetcsvs[i].Fill(secondhardestjetcsv,weight_for_histos[0])
		thirdhardestjetcsvs[i].Fill(thirdhardestjetcsv,weight_for_histos[0])
		fourthhardestjetcsvs[i].Fill(fourthhardestjetcsv,weight_for_histos[0])
		fifthhardestjetcsvs[i].Fill(fifthhardestjetcsv,weight_for_histos[0])
		combttbarmasss[i].Fill(combttbarmass,weight_for_histos[0])

#progress reporting count variables
realCount=0 #iterations
count=0     #events processed
#electrons
count_lepCands=0
count_lepPassedKin=0
count_lepPassedIso=0
count_lepPassedNoExtra=0
#neutrinos
count_neutrinoCands=0
count_neutrinoOneOption=0
#jets
count_jetEvents=0
count_jetCands=0
count_jetsPassedKin=0
count_lt4cands=0
count_mt5cands=0
count_fourCands=0
count_fiveCands=0
count_hardEnough=0
count_nobtag=0
count_onebtag=0
count_twobtags=0
count_threebtags=0
count_fourbtags=0
count_WPassed=0
count_jetsPassedAll=0
#MC matching
count_fullEvents = 0
count_mcmtleps = 0
count_mcmthads = 0
count_mcmleptons = 0
count_mcmneutrinos = 0
count_mcmbleps = 0
count_mcmlepWs = 0
count_mcmbhads = 0
count_mcmhadWs = 0
count_qq_events = 0
count_not_qq_events = 0
count_semileptonic = 0
count_not_semileptonic = 0
count_tau_events = 0
count_ele_events = 0
count_mu_events = 0
print 'Files opened, starting analysis. . . .'

##############################################################################################################
##											Main Loop Over Events							 				##
##############################################################################################################
for event in events:
	######################################################
	##	  Grid Split, check event type, check hard cut	##
	######################################################
	#options for running on the grid: split whole sample 'in parallel'
	realCount = realCount + 1
	if ((realCount-1)-options.iJob) % options.nTotalJobs != 0 :
		continue
	#set the original PIDs to be zero for now
	motherParticles[0] = 0
	motherParticles[1] = 0
	#Check to make sure the event type is correct
	is_qq = False
	is_semilep = False
	original_particles = []
	lep_ids = []
	if options.mc == 'yes' :
		event.getByLabel(GenLabel,GenHandle)
		if GenHandle.isValid() :
			GenParticles = GenHandle.product()
			for ig in GenParticles :
				if ig.pt()<0 :
					continue
				#Look through all the particles for protons; append their first daughters to the list
				if ig.pdgId() == 2212 :
					original_particles.append(ig.daughter(0).pdgId())
				#Look through particles for all ts
				if math.fabs(ig.pdgId()) == 6 and ig.status() == 3 :
					#look through all the daughters for Ws.
					for i in range(ig.numberOfDaughters()) :
						dau = ig.daughter(i)
						if math.fabs(dau.pdgId()) == 24 :
							#if the W doesn't have two daughters, I don't know what the hell happened.
							if dau.numberOfDaughters() != 2 :
								info = 'W without two daughters. PARTICLE: ' + getId(ig.pdgId()) + ', DAUGHTERS : '
								for j in range(ig.numberOfDaughters()) :
									info = info + getId(ig.daughter(i)) + ' '
								print info
								continue
							#check if the W decayed leptonically. 
							if (math.fabs(dau.daughter(0).pdgId()) == 11 or math.fabs(dau.daughter(0).pdgId()) == 13 or math.fabs(dau.daughter(0).pdgId()) == 15) :
								lep_ids.append(dau.daughter(0).pdgId())
							if (math.fabs(dau.daughter(1).pdgId()) == 11 or math.fabs(dau.daughter(1).pdgId()) == 13 or math.fabs(dau.daughter(1).pdgId()) == 15) :
								lep_ids.append(dau.daughter(1).pdgId())
			#cut on MC event type
			is_qq = len(original_particles) == 2 and math.fabs(original_particles[0]) < 6 and math.fabs(original_particles[1]) < 6 and original_particles[0] + original_particles[1] == 0
			is_semilep = len(lep_ids) == 1
			if options.tMotherCheck == 'qq' and not is_qq :
				continue
			if options.tMotherCheck == 'gg' and is_qq :
				continue
			if (options.decayType == 'semilep' or options.decayType == 'semileptonic') and not len(lep_ids) == 1 :
				continue
			if (options.decayType == 'dilep' or options.decayType == 'dileptonic') and not len(lep_ids) == 2 :
				continue
			if (options.decayType == 'had' or options.decayType == 'hadronic') and not len(lep_ids) == 0 :
				continue
			#Set the PIDs of the mothers as long as there are only two of them
			if len(original_particles) == 2 :
				motherParticles[0] = original_particles[0]
				motherParticles[1] = original_particles[1]

	count = count + 1
	#Hard cut if nMaxEvents isn't -1 for 'all data'
	if count == options.nMaxEvents + 1 :
		print 'ITERATION '+str(options.nMaxEvents)+' REACHED'
		break

	######################################################
	##				   Report Progress					##
	######################################################
	if count%options.printEvery == 0 :
		print 'Iteration: ' + str(count) + ''
		print 'Fully Reconstructed Events: ' + str(count_fullEvents) + ''
		if count_lepCands!=0:
			print '------------------Leptons-------------------------'
			print 'Lepton Candidates: ' + str(count_lepCands) + ''
			print '----Passed Kinematic Cuts: '+str(count_lepPassedKin)+' (%.4f%%)'%(100.0*count_lepPassedKin/count_lepCands)
			print '----Passed Isolation Cuts: '+str(count_lepPassedIso)+' (%.4f%%)'%(100.0*count_lepPassedIso/count_lepCands)
			print '----Passed Excess Muon Cut: '+str(count_lepPassedNoExtra)+' (%.4f%%)'%(100.0*count_lepPassedNoExtra/count_lepCands)
		if count_neutrinoCands!=0:
			print '------------------Neutrinos-----------------------'
			print 'Neutrino Candidates: ' + str(count_neutrinoCands) + ''
			print '----Neutrinos with only one Pz solution: '+str(count_neutrinoOneOption)+' (%.4f%%)'%(100.0*count_neutrinoOneOption/count_neutrinoCands)
		if count_jetCands!=0 and count_jetEvents!=0 :
			print '--------------------Jets--------------------------'
			print 'Jet Candidates: ' + str(count_jetCands) + ''
			print '----Passed Kinematic Cuts: '+str(count_jetsPassedKin)+' (%.4f%%)'%(100.0*count_jetsPassedKin/count_jetCands)
			print '----Passed All Cuts: '+str(count_jetsPassedAll)+' (%.4f%%)'%(100.0*count_jetsPassedAll/count_jetCands)
			print '----Events with 4 jet candidates: '+str(count_fourCands)+' (%.4f%%)'%(100.0*count_fourCands/count_jetEvents)
			print '----Events with 5 jet candidates: '+str(count_fiveCands)+' (%.4f%%)'%(100.0*count_fiveCands/count_jetEvents)
			print '----------Events with  <4 jet candidates: '+str(count_lt4cands)+' (%.4f%%)'%(100.0*count_lt4cands/count_jetEvents)
			print '----------Events with  >5 jet candidates: '+str(count_mt5cands)+' (%.4f%%)'%(100.0*count_mt5cands/count_jetEvents)
			print '----Events with hard enough jets: '+str(count_hardEnough)+' (%.4f%%)'%(100.0*count_hardEnough/count_jetEvents)
			print '----------Events with 0 btags (rejected): '+str(count_nobtag)+' (%.4f%%)'%(100.0*count_nobtag/count_jetEvents)
			print '----------Events with 1 btag : '+str(count_onebtag)+' (%.4f%%)'%(100.0*count_onebtag/count_jetEvents)
			print '----------Events with 2 btags: '+str(count_twobtags)+' (%.4f%%)'%(100.0*count_twobtags/count_jetEvents)
			print '----------Events with 3 btags: '+str(count_threebtags)+' (%.4f%%)'%(100.0*count_threebtags/count_jetEvents)
			print '----------Events with 4 btags: '+str(count_fourbtags)+' (%.4f%%)'%(100.0*count_fourbtags/count_jetEvents)
		if count_fullEvents!=0 and options.mc == 'yes' :
			print '-------------------MC Matching---------------------'
			print '---Matched leptons: '+str(count_mcmleptons)+' (%.4f%%)'%(100.0*count_mcmleptons/count_fullEvents)
			print '---Matched neutrinos: '+str(count_mcmneutrinos)+' (%.4f%%)'%(100.0*count_mcmneutrinos/count_fullEvents)
			print '---Matched leptonic bs: '+str(count_mcmbleps)+' (%.4f%%)'%(100.0*count_mcmbleps/count_fullEvents)
			print '---Matched leptonic Ws: '+str(count_mcmlepWs)+' (%.4f%%)'%(100.0*count_mcmlepWs/count_fullEvents)
			print '---Matched hadronic bs: '+str(count_mcmbhads)+' (%.4f%%)'%(100.0*count_mcmbhads/count_fullEvents)
			print '---Matched hadronic Ws: '+str(count_mcmhadWs)+' (%.4f%%)'%(100.0*count_mcmhadWs/count_fullEvents)
			print '---Matched leptonic ts: '+str(count_mcmtleps)+' (%.4f%%)'%(100.0*count_mcmtleps/count_fullEvents)
			print '---Matched hadronic ts: '+str(count_mcmthads)+' (%.4f%%)'%(100.0*count_mcmthads/count_fullEvents)
			print '# of        qq        events: '+str(count_qq_events)+' (%.4f%%)'%(100.0*count_qq_events/count_fullEvents)
			print '# of   not     qq     events: '+str(count_not_qq_events)+' (%.4f%%)'%(100.0*count_not_qq_events/count_fullEvents)
			print '# of   semileptonic   events: '+str(count_semileptonic)+' (%.4f%%)'%(100.0*count_semileptonic/count_fullEvents)
			print '# of not semileptonic events: '+str(count_not_semileptonic)+' (%.4f%%)'%(100.0*count_not_semileptonic/count_fullEvents)
			print '# of     electron     events: '+str(count_ele_events)+' (%.4f%%)'%(100.0*count_ele_events/count_fullEvents)			
			print '# of       muon       events: '+str(count_mu_events)+' (%.4f%%)'%(100.0*count_mu_events/count_fullEvents)			
			print '# of       tau        events: '+str(count_tau_events)+' (%.4f%%)'%(100.0*count_tau_events/count_fullEvents)			
		print '--------------------------------------------------'
		print ''
		print ''

	######################################################
	##		Check for valid event, set variables 		##
	##			    for selection histos				##
	######################################################

	#start all the weights off at one
	weight_GJR_scale[0]     = 1.0
	weight_CT10_scale[0]    = 1.0
	weight_cteq_scale[0]    = 1.0
	weight_top_pT[0]        = 1.0
	weight_btag_eff[0]      = 1.0
	weight_btag_eff_err[0]  = 1.0
	weight_pileup[0]        = 1.0
	weight_tracking[0]      = 1.0
	weight_tracking_low[0]  = 1.0
	weight_tracking_hi[0]   = 1.0
	weight_lep_ID[0]        = 1.0
	weight_lep_ID_low[0]    = 1.0
	weight_lep_ID_hi[0]     = 1.0
	weight_lep_iso[0]       = 1.0
	weight_lep_iso_low[0]   = 1.0
	weight_lep_iso_hi[0]    = 1.0
	weight_trig_eff[0]      = 1.0
	weight_trig_eff_low[0]  = 1.0
	weight_trig_eff_hi[0]   = 1.0
	weight_for_histos[0]    = 1.0
	
	#get variables from pat-tuple, make sure event is complete
	#leptons
	event.getByLabel(leptonPtLabel,leptonPtHandle)
	if not leptonPtHandle.isValid() :
		continue
	leptonPts = leptonPtHandle.product()
	if len(leptonPts) < 1 :
		continue
	lpt = leptonPts[0]
	event.getByLabel(leptonEtaLabel,leptonEtaHandle)
	leptonEtas = leptonEtaHandle.product()
	leta = leptonEtas[0]
	event.getByLabel(leptonPhiLabel,leptonPhiHandle)
	leptonPhis = leptonPhiHandle.product()
	event.getByLabel(leptonPfisoLabel,leptonPfisoHandle)
	leptonPfisos = leptonPfisoHandle.product()
	pfiso = leptonPfisos[0] / lpt
	events.getByLabel(leptonChargeLabel,leptonChargeHandle)
	leptonCharges = leptonChargeHandle.product()
	dummylep = ROOT.TLorentzVector()
	dummylep.SetPtEtaPhiM(lpt,leta,leptonPhis[0],MLEP)
	#other leptons
	event.getByLabel(otherleptonPtLabel,otherleptonPtHandle)
	if not otherleptonPtHandle.isValid() :
		continue
	otherleptonPts = otherleptonPtHandle.product()
	event.getByLabel(otherleptonEtaLabel,otherleptonEtaHandle)
	otherleptonEtas = otherleptonEtaHandle.product()
	#met
	metPt = 0
	metPhi = 0
	event.getByLabel(metPtLabel,metPtHandle)
	if not metPtHandle.isValid() :
		continue
	metPt = metPtHandle.product()
	metpt = metPt[0]
	event.getByLabel(metPhiLabel,metPhiHandle)
	metPhi = metPhiHandle.product()
	dummyMET = ROOT.TLorentzVector()
	dummyMET.SetPtEtaPhiM(metpt,0.0,metPhi[0],0.0)
	#jets
	event.getByLabel(jetPtLabel,jetPtHandle)
	if not jetPtHandle.isValid() :
		continue
	toomanyjetPts = jetPtHandle.product()
	if len(toomanyjetPts) < 4 :
		continue	
	event.getByLabel(jetEtaLabel,jetEtaHandle)
	toomanyjetEtas = jetEtaHandle.product()
	event.getByLabel(jetPhiLabel,jetPhiHandle)
	toomanyjetPhis = jetPhiHandle.product()
	event.getByLabel(jetMassLabel,jetMassHandle)
	toomanyjetMasss = jetMassHandle.product()
	event.getByLabel(jetCSVLabel,jetCSVHandle)
	toomanyjetCSVs = jetCSVHandle.product()
	# Get jet flavors
	if options.mc == 'yes' :
		event.getByLabel (jetFlavorLabel,jetFlavorHandle)
		toomanyjetFlavors = jetFlavorHandle.product()
	#order jets by pt
	dummypts = []
	for i in range(len(toomanyjetPts)) :
		pt_tuple = -1*toomanyjetPts[i],i
		dummypts.append(pt_tuple)
	dummypts.sort()
	hardestjetpt       = toomanyjetPts[dummypts[0][1]]
	secondhardestjetpt = toomanyjetPts[dummypts[1][1]]
	thirdhardestjetpt  = toomanyjetPts[dummypts[2][1]]
	fourthhardestjetpt = toomanyjetPts[dummypts[3][1]]
	fifthhardestjetpt  = -10.
	hardestjeteta       = toomanyjetEtas[dummypts[0][1]]
	secondhardestjeteta = toomanyjetEtas[dummypts[1][1]]
	thirdhardestjeteta  = toomanyjetEtas[dummypts[2][1]]
	fourthhardestjeteta = toomanyjetEtas[dummypts[3][1]]
	fifthhardestjeteta  = -10.
	hardestjetmass       = toomanyjetMasss[dummypts[0][1]]
	secondhardestjetmass = toomanyjetMasss[dummypts[1][1]]
	thirdhardestjetmass  = toomanyjetMasss[dummypts[2][1]]
	fourthhardestjetmass = toomanyjetMasss[dummypts[3][1]]
	fifthhardestjetmass  = -10.
	hardestjetcsv       = toomanyjetCSVs[dummypts[0][1]]
	secondhardestjetcsv = toomanyjetCSVs[dummypts[1][1]]
	thirdhardestjetcsv  = toomanyjetCSVs[dummypts[2][1]]
	fourthhardestjetcsv = toomanyjetCSVs[dummypts[3][1]]
	fifthhardestjetcsv  = -10.
	jet1 = ROOT.TLorentzVector()
	jet1.SetPtEtaPhiM(hardestjetpt,hardestjeteta,toomanyjetPhis[dummypts[0][1]],hardestjetmass)
	jet2 = ROOT.TLorentzVector()
	jet2.SetPtEtaPhiM(secondhardestjetpt,secondhardestjeteta,toomanyjetPhis[dummypts[1][1]],secondhardestjetmass)
	jet3 = ROOT.TLorentzVector()
	jet3.SetPtEtaPhiM(thirdhardestjetpt,thirdhardestjeteta,toomanyjetPhis[dummypts[2][1]],thirdhardestjetmass)
	jet4 = ROOT.TLorentzVector()
	jet4.SetPtEtaPhiM(fourthhardestjetpt,fourthhardestjeteta,toomanyjetPhis[dummypts[3][1]],fourthhardestjetmass)
	combttbar = jet1+jet2+jet3+jet4+dummylep+dummyMET
	if len(dummypts) > 4 :
		fifthhardestjetpt  = toomanyjetPts[dummypts[4][1]]
		fifthhardestjeteta  = toomanyjetEtas[dummypts[4][1]]
		fifthhardestjetmass  = toomanyjetMasss[dummypts[4][1]]
		fifthhardestjetcsv  = toomanyjetCSVs[dummypts[4][1]]
		jet5 = ROOT.TLorentzVector()
		jet5.SetPtEtaPhiM(fifthhardestjetpt,fifthhardestjeteta,toomanyjetPhis[dummypts[4][1]],fifthhardestjetmass)
		combttbar = combttbar + jet5
	combttbarmass = combttbar.M()


	lepton = ROOT.TLorentzVector()
	leptonCharge=0
	######################################################
	##				  lepton selection					##
	######################################################
	nLepVal = 0
	nLepExtra = 0
	for i in range(len(leptonPts)):
		count_lepCands += 1
		#see if the lepton is valid but for a lower pt: help exclude dilepton events
		if options.lep == 'mu' :
			#look for extra muons
			if (leptonPts[i] < 26.0 and leptonPts[i] > 10.0 and math.fabs(leptonEtas[i]) < 2.5) :
				nLepExtra += 1
				continue
			#pt and eta cut, want large muon within detector range
			if (leptonPts[i] < 26.0 or math.fabs(leptonEtas[i]) > 2.1) :
				continue
			count_lepPassedKin += 1
			#isolation cut with sideband switch
			if (options.sideband == 'no' and leptonPfisos[i] / leptonPts[i] > 0.12) or (options.sideband == 'yes' and (leptonPfisos[i] / leptonPts[i] < 0.13 or leptonPfisos[i] / leptonPts[i] > 0.20)) :
				continue
			lepindex = i
			nLepVal += 1
			count_lepPassedIso = count_lepPassedIso + 1
		elif options.lep == 'el' :
			#look for extra electrons
			if (leptonPts[i] < 30.0 and leptonPts[i] > 15.0 and math.fabs(leptonEtas[i]) < 2.5 and not (math.fabs(leptonEtas[i]) > 1.44 and math.fabs(leptonEtas[i]) < 1.56)) :
				nLepExtra += 1
				continue
			#pt and eta cut, want large electron within detector range (sensitivity degraded in transition region)
			if (leptonPts[i] < 30.0 or math.fabs(leptonEtas[i])>2.5 or (math.fabs(leptonEtas[i]) > 1.44 and math.fabs(leptonEtas[i]) < 1.56)) :
				continue
			count_lepPassedKin += 1
			#isolation cut with sideband switch
			if (options.sideband == 'no' and leptonPfisos[i] / leptonPts[i] > 0.10) or (options.sideband == 'yes' and (leptonPfisos[i] / leptonPts[i] < 0.11 or leptonPfisos[i] / leptonPts[i] > 0.15)) :
				continue
			lepindex = i
			nLepVal += 1
			count_lepPassedIso = count_lepPassedIso + 1
	#check veto on other lepton type
	for i in range(len(otherleptonPts)) :
		if options.lep == 'mu' :
			#look for extra electrons
			if (otherleptonPts[i] > 15.0 and math.fabs(otherleptonEtas[i]) < 2.5) :
				if not (math.fabs(otherleptonEtas[i]) > 1.44 and math.fabs(otherleptonEtas[i]) < 1.56) :
					nLepExtra += 1
		elif options.lep == 'el' :
			#look for extra muons
			if (otherleptonPts[i] > 10.0 and math.fabs(otherleptonEtas[i]) < 2.5) :
				nLepExtra += 1
	if nLepVal != 1 or nLepExtra > 0: #if there's more or less than 1 good lepton, or if there's an extra almost-valid lepton
		CUTFLOW = 0
		endEvent(CUTFLOW)
		continue
	count_lepPassedNoExtra = count_lepPassedNoExtra + 1
	#add the one valid lepton to the new root file
	pt[0]=leptonPts[lepindex]
	eta[0]=leptonEtas[lepindex]
	phi[0]=leptonPhis[lepindex]
	mass[0]=MLEP
	lepton.SetPtEtaPhiM(pt[0], eta[0], phi[0], mass[0])
	px[0]=lepton.Px()
	py[0]=lepton.Py()
	pz[0]=lepton.Pz()
	E[0]=lepton.E()
	if options.lep == 'mu' :
		PID[0]=-1*leptonCharges[lepindex]*13
	elif options.lep == 'el' :
		PID[0]=-1*leptonCharges[lepindex]*11
	#Reset the histogram variables to include the real lepton
	lpt = lepton.Pt()
	leta = lepton.Eta()
	pfiso = leptonPfisos[lepindex] / lpt
	combttbar = combttbar - dummylep + lepton
	combttbarmass = combttbar.M()
	#set the charge of the lepton
	leptonCharge = leptonCharges[lepindex]

	prewiggle_lpt  = leptonPts[lepindex]
	prewiggle_leta = leptonEtas[lepindex]

	met  = ROOT.TLorentzVector()
	met1 = ROOT.TLorentzVector()
	met2 = ROOT.TLorentzVector()
	######################################################
	##			  		MET setup for fit				##
	######################################################
	count_neutrinoCands += 1
	#find the best first guess of the Pz based on the W mass and the lepton we found
	met1.SetPtEtaPhiM(metPt[0],0.0,metPhi[0],0.0)
	met2   = met1
	pTv    = metPt[0]
	phivec = [math.cos(metPhi[0]),math.sin(metPhi[0])]
	Elep   = lepton.E()
	plep   = math.sqrt(lepton.Px()*lepton.Px()+lepton.Py()*lepton.Py()+lepton.Pz()*lepton.Pz())
	pZlep  = lepton.Pz()
	pPhi   = lepton.Px()*phivec[0]+lepton.Py()*phivec[1]
	arg0   = MW*MW+plep*plep-Elep*Elep+2.*pTv*pPhi
	arg    = Elep*Elep*(4.*pTv*pTv*(pZlep*pZlep-Elep*Elep)+arg0*arg0) #discriminant in the quadratic equation solution
	if not arg > 0 : #If discriminant is imaginary
		count_neutrinoOneOption += 1
		pzv1 = pZlep*arg0/(2.*(Elep*Elep-pZlep*pZlep))
		met1.SetPz(pzv1)
		met1.SetE(math.sqrt(met1.Px()*met1.Px()+met1.Py()*met1.Py()+met1.Pz()*met1.Pz()))
		met2 = met1
	else : #have two choices for the neutrino Pz from the quadratic equation
		pzv1 = (pZlep*arg0+sqrt(arg))/(2.*(Elep*Elep-pZlep*pZlep))
		met1.SetPz(pzv1)
		met1.SetE(math.sqrt(met1.Px()*met1.Px()+met1.Py()*met1.Py()+met1.Pz()*met1.Pz()))
		pzv2 = (pZlep*arg0-sqrt(arg))/(2.*(Elep*Elep-pZlep*pZlep))
		met2.SetPz(pzv2)
		met2.SetE(math.sqrt(met2.Px()*met2.Px()+met2.Py()*met2.Py()+met2.Pz()*met2.Pz()))
	#Reset the histogram variables to use a better met
	combttbar = combttbar - dummyMET + met1
	combttbarmass = combttbar.M()

	jetCands = []
	######################################################
	##   Jet Selection (AK5, check W helicity paper)    ##
	######################################################
	count_jetEvents=count_jetEvents+1
	#get the number of btags (just for interest, selection hasn't begun yet)
	numberOfbTags = 0
	for csv in toomanyjetCSVs :
		if csv > 0.679 :
			numberOfbTags = numberOfbTags + 1
	nbTags[0] = numberOfbTags	
	count_jetCands += len(toomanyjetPts)
	#pt and eta cuts
	jetCandCSVs = []
	jetCandFlavors = []
	for i in range(len(dummypts)) :
		j = dummypts[i][1]
		if toomanyjetPts[j] > 20. and math.fabs(toomanyjetEtas[j]) < 2.5 :
			jetCands.append(ROOT.TLorentzVector())
			jetCands[len(jetCands)-1].SetPtEtaPhiM(toomanyjetPts[j],toomanyjetEtas[j],toomanyjetPhis[j],toomanyjetMasss[j])
			jetCandCSVs.append(toomanyjetCSVs[j])
			if options.mc == 'yes' :
				jetCandFlavors.append(toomanyjetFlavors[j])
	count_jetsPassedKin += len(jetCands)
	#Reset the jets used for the histograms
	jets_reset = [jet1,jet2,jet3,jet4]
	jet_csvs_reset = [hardestjetcsv,secondhardestjetcsv,thirdhardestjetcsv,fourthhardestjetcsv]
	if len(dummypts) > 4 :
		jets_reset.append(jet5)
		jet_csvs_reset.append(fifthhardestjetcsv)
	for i in range(len(jetCands)) :
		if i >= len(jets_reset) :
			break
		jets_reset[i] = jetCands[i]
		jet_csvs_reset[i] = jetCandCSVs[i]
	hardestjetpt         = jets_reset[0].Pt()
	secondhardestjetpt   = jets_reset[1].Pt()
	thirdhardestjetpt    = jets_reset[2].Pt()
	fourthhardestjetpt   = jets_reset[3].Pt()
	hardestjeteta        = jets_reset[0].Eta()
	secondhardestjeteta  = jets_reset[1].Eta()
	thirdhardestjeteta   = jets_reset[2].Eta()
	fourthhardestjeteta  = jets_reset[3].Eta()
	hardestjetmass       = jets_reset[0].M()
	secondhardestjetmass = jets_reset[1].M()
	thirdhardestjetmass  = jets_reset[2].M()
	fourthhardestjetmass = jets_reset[3].M()
	hardestjetcsv        = jet_csvs_reset[0]
	secondhardestjetcsv  = jet_csvs_reset[1]
	thirdhardestjetcsv   = jet_csvs_reset[2]
	fourthhardestjetcsv  = jet_csvs_reset[3]
	combttbarmass = (jets_reset[0] + jets_reset[1] + jets_reset[2] + jets_reset[3] + lepton + met1).M()
	if len(dummypts) > 4 :
		fifthhardestjetpt   = jets_reset[4].Pt()
		fifthhardestjeteta  = jets_reset[4].Eta()
		fifthhardestjetmass = jets_reset[4].M()
		fifthhardestjetcsv  = jet_csvs_reset[4]
		combttbarmass = (jets_reset[0] + jets_reset[1] + jets_reset[2] + jets_reset[3] + jets_reset[4] + lepton + met1).M()
	#require at least 4 jet candidates
	if len(jetCands) < 4 :
		count_lt4cands += 1
		CUTFLOW = 2
		endEvent(CUTFLOW)
		continue
	#require no more than 5 jet candidates
	if len(jetCands) > 5 :
		count_mt5cands += 1
		CUTFLOW = 2
		endEvent(CUTFLOW)
		continue
	if len(jetCands) == 4 :
		#Reset the histogram variables to disregard the fifth jet
		fifthhardestjetpt   = -10.
		fifthhardestjeteta  = -10.
		fifthhardestjetmass = -10.
		fifthhardestjetcsv  = -10.
		combttbarmass = (jets_reset[0] + jets_reset[1] + jets_reset[2] + jets_reset[3] + lepton + met1).M()
		count_fourCands += 1
		#Record that the event had four jets
		for i in range(nParticles) :
			nValidJets[0] = 4
	elif len(jetCands) == 5 :
		#Record that the event had five jets
		for i in range(nParticles) :
			nValidJets[0] = 5
		count_fiveCands += 1
	#cut on pts of hardest 4 jets (cuts from A_C paper)
	if (jetCands[0].Pt()<45.0 or jetCands[1].Pt()<35.0 or jetCands[2].Pt()<20.0 or jetCands[3].Pt()<20.0) :
		CUTFLOW = 3
		endEvent(CUTFLOW)
		continue
	count_hardEnough += 1
	#Cut that we have at least the right number of btagged jets with pt > 20
	nbtags = 0
	for i in range(len(jetCands)) :
		if jetCands[i].Pt()>20.0 :
			if jetCandCSVs[i] > 0.679 :
				nbtags += 1
	nbTags[0] = nbtags
	if nbtags < 1 :
		count_nobtag += 1
		CUTFLOW = 4
		endEvent(CUTFLOW)
		continue
	if nbtags == 1 :
		count_onebtag += 1
	elif nbtags == 2 :
		count_twobtags += 1
	elif nbtags == 3 :
		count_threebtags += 1
	elif nbtags == 4 :
		count_fourbtags += 1
	count_jetsPassedAll += len(jetCands)
	#print str(len(jetCands))+' jet cands total, top four pts: '+str(jetCands[0].Pt())+' '+str(jetCands[1].Pt())+' '+str(jetCands[2].Pt())+' '+str(jetCands[3].Pt())+' '

	############################################################################################################
	#At this point, we have selected all the leptons, met, and jets for the event. The fourvector of the lepton#
	#is called 'lepton', the 2 cursory (unfitted) fourvectors of the met are in 'met1/2', and the list of jet  #
	#fourvectors is called 'jetCands'. Next we want to reconstruct the t and tbar from the event. 			   #
	############################################################################################################
	CUTFLOW = 5

	######################################################
	##				EVENT RECONSTRUCTION				##
	######################################################
	#make a list of all combinations of the jet candidates
	combos    = []
	comboCSVs = []
	icombos   = [] #this just checks the combinatorics. Print it if you want.
	extraJetCSVs = []
	if len(jetCands) == 4 :
		for i in range(len(jetCands)) :
			iblep = i
			iothers = []
			for j in range(len(jetCands)) :
				if j!=iblep :
					iothers.append(j)
			for j in range(len(iothers)) :
				combos.append([])
				combos[len(combos)-1].append(jetCands[iblep])
				comboCSVs.append([])
				comboCSVs[len(comboCSVs)-1].append(jetCandCSVs[iblep])
				icombos.append([])
				icombos[len(icombos)-1].append(iblep)
				extraJetCSVs.append(-1.0)
				for k in range(len(iothers)) :
					combos[len(combos)-1].append(jetCands[iothers[(k-j+len(iothers))%len(iothers)]])
					comboCSVs[len(comboCSVs)-1].append(jetCandCSVs[iothers[(k-j+len(iothers))%len(iothers)]])
					icombos[len(icombos)-1].append(iothers[(k-j+len(iothers))%len(iothers)])
	elif len(jetCands) == 5 :
		for i in range(len(jetCands)) :
			iblep = i
			for j in range(len(jetCands)) :
				if j == iblep :
					continue
				iextraJet = j
				iothers   = []
				for k in range(len(jetCands)) :
					if k != iextraJet and k != iblep :
						iothers.append(k)
				for k in range(len(iothers)) :
					combos.append([])
					combos[len(combos)-1].append(jetCands[iblep])
					comboCSVs.append([])
					comboCSVs[len(comboCSVs)-1].append(jetCandCSVs[iblep])
					icombos.append([])
					icombos[len(icombos)-1].append(str(iblep))
					extraJetCSVs.append(jetCandCSVs[iextraJet])
					for n in range(len(iothers)) :
						combos[len(combos)-1].append(jetCands[iothers[(n-k+len(iothers))%len(iothers)]])
						comboCSVs[len(comboCSVs)-1].append(jetCandCSVs[iothers[(n-k+len(iothers))%len(iothers)]])
						icombos[len(icombos)-1].append(str(iothers[(n-k+len(iothers))%len(iothers)]))
	#Combos now includes all possible orderings of the jets (order: blep, bhad, Wsub, Wsub)
	#For each combo, the extra jet goes in comboExtraJet
	#Remove any combos that have neither supposed b tagged as such, or any combo where the CSV value excludes one of the bs from being a b
	i = 0
	while True :
		twountagged = comboCSVs[i][0] < 0.679 and comboCSVs[i][1] < 0.679
		oneuntagged = comboCSVs[i][0] < 0.679 or comboCSVs[i][1] < 0.679
		if (nbtags == 1 and twountagged) or (nbtags >= 2 and oneuntagged) or comboCSVs[i][0] == -1.0 or comboCSVs[i][1] == -1.0 :
			combos.pop(i)
			comboCSVs.pop(i)
			icombos.pop(i)
			extraJetCSVs.pop(i)
		else :
			i += 1
		if i >= len(combos) :
			break
	i = 0
	if not len(combos) > 0 :
		print 'no valid combos'
		continue
	#Variables for the kinematic fits
	bestParValues1 = []
	bestParValues2 = []
	bestParValues = []
	Chis1 = []
	Chis2 = []
	parErrs = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
	parNames = ['pZv','scaleLep','scaleblep','scalebhad','scaleWsub1','scaleWsub2']
	for i in range(len(combos)) :
		bestParValues1.append([0.0,0.0,0.0,0.0,0.0,0.0])
		bestParValues2.append([0.0,0.0,0.0,0.0,0.0,0.0])
	#Now perform kinematic fits for each combination of jetsand both neutrino solutions.
	nFits = 1
	if met1.Pz() != met2.Pz() :
		nFits = 2
	for iFit in range(nFits) : #one for each neutrino solution
		#set which neutrino to use for this iteration
		if iFit == 0 :
			met = met1
		elif iFit == 1:
			met = met2
		#Stuff common to every fit
		minuit = ROOT.TMinuit(6)
		minuit.SetFCN(fcn)
		ierflag = ROOT.Long(1)
		arglist = array( 'd', 1*[0.] )
		arglist[0] = -1.
		minuit.mnexcm('SET PRINT', arglist, 1,ierflag)
		minuit.mnexcm('SET NOWARNINGS',arglist,1,ierflag)
		arglist[0] = 100000.
		#do fits for every comination in combos
		for i in range(len(combos)) :
			#set the jet particle variables
			blep    = combos[i][0]
			bhad    = combos[i][1]
			Wsub1   = combos[i][2]
			Wsub2   = combos[i][3]
			blepCSV = comboCSVs[i][0]
			bhadCSV = comboCSVs[i][1]
			WCSV1   = comboCSVs[i][2]
			WCSV2   = comboCSVs[i][3]
			extraCSV = extraJetCSVs[i]
			#set the parameters in minuit
			minuit.mnparm(0,parNames[0],ROOT.Double(met.Pz()),1.0,0,0,ierflag)
			for j in range(1,6) :
				minuit.mnparm(j,parNames[j],1.0,1.0,0,0,ierflag)
			#minimize
			minuit.mnexcm('MIGRAD', arglist, 1,ierflag)
			if ierflag != 0 :
				print 'PROBLEM IN FIT: ierflag = '+str(ierflag)+''
			#Set fit Chi of this particular combination
			if iFit == 0 :
				Chis1.append((minf,i))
			elif iFit == 1:
				Chis2.append((minf,i))
			minf = 1000000000.
			#Get the best parameters back from minuit
			for j in range(6) :
				tmp = ROOT.Double(1.0)
				minuit.GetParameter(j,tmp,ROOT.Double(parErrs[j]))
				if iFit == 0 :
					bestParValues1[i][j] = tmp
				elif iFit == 1:
					bestParValues2[i][j] = tmp
	#Find the best fit for this event and record it
	Chis1.sort()
	Chis2.sort()
	plot_final_chi = 1.0
	if len(Chis2) > 0 and Chis2[0][0] < Chis1[0][0] :
		j = Chis2[0][1]
		met = met2
		blep = combos[j][0]
		bhad = combos[j][1]
		Wsub1 = combos[j][2]
		Wsub2 = combos[j][3]
		for i in range(6) :
			bestParValues.append(bestParValues2[j][i])
		for i in range(nParticles) :
			finalChi[i] = Chis2[0][0]
		plot_final_chi = Chis2[0][0]
	else :
		j = Chis1[0][1]
		met = met1
		blep = combos[j][0]
		bhad = combos[j][1]
		Wsub1 = combos[j][2]
		Wsub2 = combos[j][3]
		for i in range(6) :
			bestParValues.append(bestParValues1[j][i])
		for i in range(nParticles) :
			finalChi[i] = Chis1[0][0]
		plot_final_chi = Chis1[0][0]
	#Rescale the particle fourvectors based on the optimal parameters
	lepton = pscale(bestParValues[1], lepton)
	blep = pscale(bestParValues[2], blep)
	bhad = pscale(bestParValues[3], bhad)
	Wtag = pscale(bestParValues[4],Wsub1)+pscale(bestParValues[5],Wsub2)
	newmetx = met.Px()+ (1.0-bestParValues[1])*lepton.Px()+(1.0-bestParValues[2])*blep.Px()+(1.0-bestParValues[3])*bhad.Px()
	newmetx += (1.0-bestParValues[4])*Wsub1.Px()+(1.0-bestParValues[5])*Wsub2.Px()
	newmety = met.Py()+ (1.0-bestParValues[1])*lepton.Py()+(1.0-bestParValues[2])*blep.Py()+(1.0-bestParValues[3])*bhad.Py()
	newmety += (1.0-bestParValues[4])*Wsub1.Py()+(1.0-bestParValues[5])*Wsub2.Py()
	met.SetPx(newmetx)
	met.SetPy(newmety)
	met.SetPz(bestParValues[0])
	met.SetE(math.sqrt(met.Px()*met.Px()+met.Py()*met.Py()+met.Pz()*met.Pz()))
	for i in range(6) :
		bestFitParValues[i] = bestParValues[i]

	######################################################
	##			WRITE PARTICLES TO BRANCHES				##
	######################################################
	#Add the rescaled particles to branches
	particles = [lepton,met,blep,(lepton + met),bhad,Wtag,(lepton + met + blep),(Wtag + bhad)]
	dummy_is_leptonic_side = [1,1,1,1,0,0,1,0]
	for i in range(len(particles)) :
		pt[i] = particles[i].Pt()
		eta[i]=particles[i].Eta()
		phi[i]=particles[i].Phi()
		mass[i]=particles[i].M()
		px[i]=particles[i].Px()
		py[i]=particles[i].Py()
		pz[i]=particles[i].Pz()
		E[i]=particles[i].E()
		is_leptonic_side[i]=dummy_is_leptonic_side[i]
	#The PIDs have to be done inidividually though :(
	if options.lep == 'el' :
		PID[1] = 12*leptonCharge
	if options.lep == 'mu' :
		PID[1] = 14*leptonCharge
	PID[2] = 5 * leptonCharge
	PID[3] = 24* leptonCharge
	PID[4] = 5 * leptonCharge*(-1)
	PID[5] = 24* leptonCharge*(-1)
	PID[6] = 6 * leptonCharge
	PID[7] = 6 * leptonCharge*(-1)
	#Reset the combined ttbar mass to the new and better value
	combttbarmass = (particles[6]+particles[7]).M()

	#If we made it here, we have a fully reconstructed event!
	count_fullEvents = count_fullEvents + 1

	############################################################
	## 	   			MC MATCHING FOR SIMULATED DATA		  	  ##
	############################################################

	if options.mc == 'yes' :
		#Monte Carlo ts, tbars, leptons, neutrinos, leptonic bs, hadronic bs, leptonic Ws, and hadronic Ws
		mc_particles = []
		mc_particles_found = []
		mc_dRs = []
		for i in range(8) :
			mc_particles.append(ROOT.TLorentzVector())
			mc_particles_found.append(False)
		#find out whether the t or tbar is the leptonic side
		t_is_leptonic_side = leptonCharge>0
		#start with an empty fourvector for the quark
		mother1 = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		mother2 = ROOT.TLorentzVector(1.0,0.0,0.0,1.0)
		#Open genparticles information from the file
		if GenHandle.isValid() :
			GenParticles = GenHandle.product()
			for ig in GenParticles :
				if ig.pt()<0 :
					continue
				#Look through particles for all of the protons to get the q and qbar if it's a qqbar event
				if ig.pdgId() == 2212 :
					if ig.daughter(0).pdgId() > 0 :
						mother1.SetPtEtaPhiM(ig.daughter(0).pt(),ig.daughter(0).eta(),ig.daughter(0).phi(),ig.daughter(0).mass())
					elif ig.daughter(0).pdgId() < 0 :
						mother2.SetPtEtaPhiM(ig.daughter(0).pt(),ig.daughter(0).eta(),ig.daughter(0).phi(),ig.daughter(0).mass())
				#Look through particles for all ts
				ist = ig.pdgId() == 6 and ig.status() == 3
				istbar = ig.pdgId() == -6 and ig.status() == 3
				#Add ts to list and then search through their daughters for bs and Ws (everything also done for charge conjugates)
				if ist or istbar :
					if (ist and t_is_leptonic_side) or (istbar and not t_is_leptonic_side) : #it's the leptonic t
						mc_particles[6].SetPtEtaPhiM(ig.pt(),ig.eta(),ig.phi(),ig.mass())
						mc_particles_found[6] = True
					elif (istbar and t_is_leptonic_side) or (ist and not t_is_leptonic_side) : #it's the hadronic t
						mc_particles[7].SetPtEtaPhiM(ig.pt(),ig.eta(),ig.phi(),ig.mass())
						mc_particles_found[7] = True
					for i in range(ig.numberOfDaughters()) :
						dau = ig.daughter(i)
						if dau.pt() < 0 or dau.status() != 3 :
							continue
						isb = ist and dau.pdgId() == 5
						isbbar = istbar and dau.pdgId() == -5
						isWplus = ist and dau.pdgId() == 24
						isWminus = istbar and dau.pdgId() == -24
						#Add any bs found to the approriate list depending on whether the t is the hadronic or leptonic side of the decay
						if (isb and t_is_leptonic_side) or (isbbar and not t_is_leptonic_side) : #it's the leptonic b
							mc_particles[2].SetPtEtaPhiM(dau.pt(),dau.eta(),dau.phi(),dau.mass())
							mc_particles_found[2] = True
						if (isb and not t_is_leptonic_side) or (isbbar and t_is_leptonic_side) : #it's the hadronic b
							mc_particles[4].SetPtEtaPhiM(dau.pt(),dau.eta(),dau.phi(),dau.mass())
							mc_particles_found[4] = True
						#set any Ws found to be the leptonic or hadronic Ws appropriately and look past the leptonic W for leptons and neutrinos
						if (isWplus and not t_is_leptonic_side) or (isWminus and t_is_leptonic_side) : #it's the hadronic W
							n_W_subjets = 0
							for j in range(dau.numberOfDaughters()) :
								if math.fabs(dau.daughter(j).pdgId()) < 6 :
									n_W_subjets += 1
							if n_W_subjets > 1 :
								mc_particles[5].SetPtEtaPhiM(dau.pt(),dau.eta(),dau.phi(),dau.mass())
								mc_particles_found[5] = True
						if (isWplus and t_is_leptonic_side) or (isWminus and not t_is_leptonic_side) : #it's the leptonic W
							mc_particles[3].SetPtEtaPhiM(dau.pt(),dau.eta(),dau.phi(),dau.mass())
							mc_particles_found[3] = True
							for j in range(dau.numberOfDaughters()) :
								subDau = dau.daughter(j)
								if subDau.pt() < 0 or subDau.status() != 3 :
									continue
								islep = math.fabs(subDau.pdgId()) == 11 or math.fabs(subDau.pdgId()) == 13 or math.fabs(subDau.pdgId()) == 15
								isnv  = math.fabs(subDau.pdgId()) == 12 or math.fabs(subDau.pdgId()) == 14 or math.fabs(subDau.pdgId()) == 16
								#Add any found leptons and neutrinos to the lists
								if islep :
									mc_particles[0].SetPtEtaPhiM(subDau.pt(),subDau.eta(),subDau.phi(),subDau.mass())
									mc_particles_found[0] = True
									if math.fabs(subDau.pdgId()) == 11 :
										count_ele_events += 1
									if math.fabs(subDau.pdgId()) == 13 :
										count_mu_events += 1
									if math.fabs(subDau.pdgId()) == 15 :
										count_tau_events += 1
								if isnv :
									mc_particles[1].SetPtEtaPhiM(subDau.pt(),subDau.eta(),subDau.phi(),subDau.mass())
									mc_particles_found[1] = True
		#find out the type of ttbar event
		if is_qq :
			count_qq_events = count_qq_events + 1
		else :
			count_not_qq_events = count_not_qq_events + 1
		if is_semilep :
			count_semileptonic = count_semileptonic + 1
		else :
			count_not_semileptonic = count_not_semileptonic + 1
		#the two mother particles go in the 9th and 10th spots in the array
		mc_pt[8]   = mother1.Pt()
		mc_eta[8]  = mother1.Eta()
		mc_phi[8]  = mother1.Phi()
		mc_mass[8] = mother1.M()
		mc_px[8]   = mother1.Px()
		mc_py[8]   = mother1.Py()
		mc_pz[8]   = mother1.Pz()
		mc_E[8]    = mother1.E()
		mc_pt[9]   = mother2.Pt()
		mc_eta[9]  = mother2.Eta()
		mc_phi[9]  = mother2.Phi()
		mc_mass[9] = mother2.M()
		mc_px[9]   = mother2.Px()
		mc_py[9]   = mother2.Py()
		mc_pz[9]   = mother2.Pz()
		mc_E[9]    = mother2.E()
		#try to match each particle separately
		counts = [count_mcmleptons,count_mcmneutrinos,count_mcmbleps,count_mcmlepWs,count_mcmbhads,count_mcmhadWs,count_mcmtleps,count_mcmthads]
		#fill out all of the particles
		for i in range(len(particles)) :
			if mc_particles_found[i] :
				if mc_particles[i].DeltaR(particles[i]) < options.dR and mc_particles[i].E()/particles[i].E()>0.2 :
					counts[i]  = counts[i] + 1
					mc_pt[i]   = mc_particles[i].Pt()
					mc_eta[i]  = mc_particles[i].Eta()
					mc_phi[i]  = mc_particles[i].Phi()
					mc_mass[i] = mc_particles[i].M()
					mc_px[i]   = mc_particles[i].Px()
					mc_py[i]   = mc_particles[i].Py()
					mc_pz[i]   = mc_particles[i].Pz()
					mc_E[i]    = mc_particles[i].E()
					mc_was_matched[i] = 1
				else : 
					mc_pt[i]   = mc_particles[i].Pt()
					mc_eta[i]  = mc_particles[i].Eta()
					mc_phi[i]  = mc_particles[i].Phi()
					mc_mass[i] = mc_particles[i].M()
					mc_px[i]   = mc_particles[i].Px()
					mc_py[i]   = mc_particles[i].Py()
					mc_pz[i]   = mc_particles[i].Pz()
					mc_E[i]    = mc_particles[i].E()
					mc_was_matched[i] = 0
		#set the counts to the new values
		count_mcmleptons = counts[0]
		count_mcmneutrinos = counts[1]
		count_mcmbleps = counts[2]
		count_mcmlepWs = counts[3]
		count_mcmbhads = counts[4]
		count_mcmhadWs = counts[5]
		count_mcmtleps = counts[6]
		count_mcmthads = counts[7]

	############################################################
	## 	   		Write out some reweighting factors		  	  ##
	############################################################

	event.getByLabel(dataPileupLabel,dataPileupHandle)
	if not dataPileupHandle.isValid() :
		continue
	data_pileup_number = dataPileupHandle.product()
	pileup_events[0] = 1.0*data_pileup_number[0]

	if options.mc == 'yes' :
		#PDF REWEIGHTING
		if options.tMotherCheck == 'qq' or options.tMotherCheck == 'gg' or options.PdfVersion != 'none' :
 #or options.sample_type == 'Powheg_dilep_TT' or options.sample_type == 'Powheg_had_TT' or options.sample_type == 'TTJets_FullLep' or options.sample_type == 'TTJets_Hadronic'  or options.sample_type == 'signal':
			#Get PDF weights
			if options.PdfVersion == 'CT10' :
				event.getByLabel(PdfLabel_CT10, PdfHandle_CT10)
				if not PdfHandle_CT10.isValid() :
					continue
				Pdf_ws = PdfHandle_CT10.product()
			elif options.PdfVersion == 'cteq66' :
				event.getByLabel(PdfLabel_cteq, PdfHandle_cteq)
				if not PdfHandle_cteq.isValid() :
					continue
				Pdf_ws = PdfHandle_cteq.product()
			if len(Pdf_ws) < 1 :
				continue
			Pdf_N = Pdf_ws[0]
			for i in range(len(Pdf_ws)) :
				Pdf_weights[i] = Pdf_ws[i]/Pdf_N
			# CT10 PDF reweighting
			event.getByLabel(PdfLabel_CT10, PdfHandle_CT10)
			if not PdfHandle_CT10.isValid() :
				continue
			Pdf_ws_CT10 = PdfHandle_CT10.product()
			if len(Pdf_ws_CT10) < 1 :
				continue
			weight_CT10_scale[0] = 1.0*Pdf_ws_CT10[0]/Pdf_ws[0]
			# cteq66 PDF reweighting
			event.getByLabel(PdfLabel_cteq, PdfHandle_cteq)
			if not PdfHandle_cteq.isValid() :
				continue
			Pdf_ws_cteq = PdfHandle_cteq.product()
			if len(Pdf_ws_cteq) < 1 :
				continue
			weight_cteq_scale[0] = 1.0*Pdf_ws_cteq[0]/Pdf_ws[0]
			# GJR08VFnlo PDF reweighting
			event.getByLabel(PdfLabel_GJR, PdfHandle_GJR)
			if not PdfHandle_GJR.isValid() :
				continue
			Pdf_ws_GJR = PdfHandle_GJR.product()
			if len(Pdf_ws_GJR) < 1 :
				continue
			weight_GJR_scale[0] = 1.0*Pdf_ws_GJR[0]/Pdf_ws[0]

		#qscale, etc. geninfo
		event.getByLabel(GenInfoLabel, GenInfoHandle)
		if not GenInfoHandle.isValid() :
			continue
		partons = GenInfoHandle.product()
		qscale[0] = partons.qScale()
		parton_id[0] = partons.pdf().id.first
		parton_id[1] = partons.pdf().id.second
		parton_x[0] = partons.pdf().x.first
		parton_x[1] = partons.pdf().x.second
		
		#pileup reweighting
		event.getByLabel(npvRealTrueLabel,npvRealTrueHandle)
		if not npvRealTrueHandle.isValid() :
			continue
		npvRealTrue = npvRealTrueHandle.product()
		weight_pileup[0] = data_pu_dist.GetBinContent(data_pu_dist.FindFixBin(1.0*npvRealTrue[0]))/MC_pu_dist.GetBinContent(MC_pu_dist.FindFixBin(1.0*npvRealTrue[0]))
		mc_pileup_events[0] = 1.0*npvRealTrue[0]
		
		#btag efficiency correction
		jetCandPts = []
		jetCandEtas = [] 
		for i in range(len(jetCands)):
			jetCandPts.append(jetCands[i].Pt())
			jetCandEtas.append(jetCands[i].Eta())
		w_result = get_weight_btag(jetCandPts,jetCandEtas,jetCandFlavors,jetCandCSVs)
		weight_btag_eff[0]     = w_result[0]
		weight_btag_eff_err[0] = w_result[1]

		#top pT reweighting
		if options.sample_type == 'signal' or options.sample_type == 'Powheg_qq_semilep_TT' or options.sample_type == 'Powheg_gg_semilep_TT' or options.sample_type == 'TTJets_SemiLep' :
			weight_top_pT[0] = math.sqrt(math.exp(0.159+(-0.00141)*mc_particles[6].Pt())*math.exp(0.159+(-0.00141)*mc_particles[7].Pt()))
		elif options.sample_type == 'Powheg_dilep_TT' or options.sample_type == 'TTJets_FullLep':
			weight_top_pT[0] = math.sqrt(math.exp(0.148+(-0.00129)*mc_particles[6].Pt())*math.exp(0.148+(-0.00129)*mc_particles[7].Pt()))

		#tracking efficiency SFs
		for i in range(len(tracking_eff_consts)) :
			if prewiggle_leta > tracking_eff_consts[i][0] and prewiggle_leta < tracking_eff_consts[i][1] :
				weight_tracking[0]     = tracking_eff_consts[i][2]
				weight_tracking_hi[0]  = tracking_eff_consts[i][3]
				weight_tracking_low[0] = tracking_eff_consts[i][3]

		#Lepton ID, isolation, and trigger scale factors
		if options.lep == 'el' :
			print 'WARNING: no electron ID, isolation, or trigger scale factors programmed yet!'
		if options.lep == 'mu' :
			#vertex corrections
			for vtx_key in muon_id_dict['Tight']['vtxpt20-500'].keys() :
				#print 'vtx_key = '+vtx_key+'comp 1 = '+vtx_key.split('_')[0]+' comp 2 = '+vtx_key.split('_')[1]+' pileup = '+str(pileup_events[0])
				if (pileup_events[0] > float(vtx_key.split('_')[0]) and pileup_events[0] < float(vtx_key.split('_')[1])) or (vtx_key == '28.5_30.5' and pileup_events[0] > 30) :
					muon_id_vtx_weight = muon_id_dict['Tight']['vtxpt20-500'][vtx_key]['data/mc']['efficiency_ratio']
					muon_id_vtx_weight_low = muon_id_dict['Tight']['vtxpt20-500'][vtx_key]['data/mc']['err_low']
					muon_id_vtx_weight_hi = muon_id_dict['Tight']['vtxpt20-500'][vtx_key]['data/mc']['err_hi']
					muon_iso_vtx_weight = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['vtxpt20-500'][vtx_key]['data/mc']['efficiency_ratio']
					muon_iso_vtx_weight_low = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['vtxpt20-500'][vtx_key]['data/mc']['err_low']
					muon_iso_vtx_weight_hi = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['vtxpt20-500'][vtx_key]['data/mc']['err_hi']
					muon_trigger_vtx_weight = muon_trigger_dict['IsoMu24']['TightID_IsodB']['VTX'][vtx_key]['data']['efficiency']
					muon_trigger_vtx_weight_low = muon_trigger_dict['IsoMu24']['TightID_IsodB']['VTX'][vtx_key]['data']['err_low']
					muon_trigger_vtx_weight_hi = muon_trigger_dict['IsoMu24']['TightID_IsodB']['VTX'][vtx_key]['data']['err_hi']
			#eta corrections
			for eta_key in muon_id_dict['Tight']['etapt20-500'].keys() :
				if prewiggle_leta > float(eta_key.split('_')[0]) and prewiggle_leta < float(eta_key.split('_')[1]) :
					muon_id_eta_weight = muon_id_dict['Tight']['etapt20-500'][eta_key]['data/mc']['efficiency_ratio']
					muon_id_eta_weight_low = muon_id_dict['Tight']['etapt20-500'][eta_key]['data/mc']['err_low']
					muon_id_eta_weight_hi = muon_id_dict['Tight']['etapt20-500'][eta_key]['data/mc']['err_hi']
					muon_iso_eta_weight = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['etapt20-500'][eta_key]['data/mc']['efficiency_ratio']
					muon_iso_eta_weight_low = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['etapt20-500'][eta_key]['data/mc']['err_low']
					muon_iso_eta_weight_hi = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['etapt20-500'][eta_key]['data/mc']['err_hi']
					muon_trigger_eta_weight = muon_trigger_dict['IsoMu24']['TightID_IsodB']['ETA'][eta_key]['data']['efficiency']
					muon_trigger_eta_weight_low = muon_trigger_dict['IsoMu24']['TightID_IsodB']['ETA'][eta_key]['data']['err_low']
					muon_trigger_eta_weight_hi = muon_trigger_dict['IsoMu24']['TightID_IsodB']['ETA'][eta_key]['data']['err_hi']
			#pt corrections
			nextKey = ''
			nextKey_trig = ''
			if abs(prewiggle_leta) < 0.9 :
				nextKey = 'ptabseta<0.9'
				nextKey_trig = 'PT_ABSETA_Barrel_0to0p9'
			elif abs(prewiggle_leta) < 1.2 :
				nextKey = 'ptabseta0.9-1.2'
				nextKey_trig = 'PT_ABSETA_Transition_0p9to1p2'
			elif abs(prewiggle_leta) < 2.1 :
				nextKey = 'ptabseta1.2-2.1'
				nextKey_trig = 'PT_ABSETA_Endcaps_1p2to2p1'
			elif abs(prewiggle_leta) < 2.4 :
				nextKey = 'ptabseta2.1-2.4'
				print 'WARNING: lepton found with abs(eta)>2.4! Cannot apply pt-based trigger SF!'
			else :
				print 'WARNING: lepton found with abs(eta) > 2.4! Cannot apply pt-based ID, iso, or trigger SFs!'
			for pt_key in muon_id_dict['Tight'][nextKey].keys() :
				if (prewiggle_lpt > float(pt_key.split('_')[0]) and prewiggle_lpt < float(pt_key.split('_')[1])) or (pt_key == '140_300' and prewiggle_lpt > 300) :
					muon_id_pt_weight = muon_id_dict['Tight'][nextKey][pt_key]['data/mc']['efficiency_ratio']
					muon_id_pt_weight_low = muon_id_dict['Tight'][nextKey][pt_key]['data/mc']['err_low']
					muon_id_pt_weight_hi = muon_id_dict['Tight'][nextKey][pt_key]['data/mc']['err_hi']
					muon_iso_pt_weight = muon_iso_dict['combRelIsoPF04dBeta<012_Tight'][nextKey][pt_key]['data/mc']['efficiency_ratio']
					muon_iso_pt_weight_low = muon_iso_dict['combRelIsoPF04dBeta<012_Tight'][nextKey][pt_key]['data/mc']['err_low']
					muon_iso_pt_weight_hi = muon_iso_dict['combRelIsoPF04dBeta<012_Tight'][nextKey][pt_key]['data/mc']['err_hi']
			for pt_key in muon_trigger_dict['IsoMu24']['TightID_IsodB'][nextKey_trig].keys() :
				if prewiggle_lpt > float(pt_key.split('_')[0]) and prewiggle_lpt < float(pt_key.split('_')[1]) :
					muon_trigger_pt_weight = muon_trigger_dict['IsoMu24']['TightID_IsodB'][nextKey_trig][pt_key]['data']['efficiency']
					muon_trigger_pt_weight_low = muon_trigger_dict['IsoMu24']['TightID_IsodB'][nextKey_trig][pt_key]['data']['err_low']
					muon_trigger_pt_weight_hi = muon_trigger_dict['IsoMu24']['TightID_IsodB'][nextKey_trig][pt_key]['data']['err_hi']
			weight_lep_ID[0]       = muon_id_vtx_weight*muon_id_eta_weight*muon_id_pt_weight
			weight_lep_ID_low[0]   = weight_lep_ID[0]*math.sqrt((muon_id_vtx_weight_low/muon_id_vtx_weight)**2+(muon_id_eta_weight_low/muon_id_eta_weight)**2+(muon_id_pt_weight_low/muon_id_pt_weight)**2)
			weight_lep_ID_hi[0]    = weight_lep_ID[0]*math.sqrt((muon_id_vtx_weight_hi/muon_id_vtx_weight)**2+(muon_id_eta_weight_hi/muon_id_eta_weight)**2+(muon_id_pt_weight_hi/muon_id_pt_weight)**2)
			weight_lep_iso[0]      = muon_iso_vtx_weight*muon_iso_eta_weight*muon_iso_pt_weight
			weight_lep_iso_low[0]  = weight_lep_iso[0]*math.sqrt((muon_iso_vtx_weight_low/muon_iso_vtx_weight)**2+(muon_iso_eta_weight_low/muon_iso_eta_weight)**2+(muon_iso_pt_weight_low/muon_iso_pt_weight)**2)
			weight_lep_iso_hi[0]   = weight_lep_iso[0]*math.sqrt((muon_iso_vtx_weight_hi/muon_iso_vtx_weight)**2+(muon_iso_eta_weight_hi/muon_iso_eta_weight)**2+(muon_iso_pt_weight_hi/muon_iso_pt_weight)**2)
			weight_trig_eff[0]     = muon_trigger_vtx_weight*muon_trigger_eta_weight*muon_trigger_pt_weight
			weight_trig_eff_low[0] = weight_trig_eff[0]*math.sqrt((muon_trigger_vtx_weight_low/muon_trigger_vtx_weight)**2+(muon_trigger_eta_weight_low/muon_trigger_eta_weight)**2+(muon_trigger_pt_weight_low/muon_trigger_pt_weight)**2)
			weight_trig_eff_hi[0]  = weight_trig_eff[0]*math.sqrt((muon_trigger_vtx_weight_hi/muon_trigger_vtx_weight)**2+(muon_trigger_eta_weight_hi/muon_trigger_eta_weight)**2+(muon_trigger_pt_weight_hi/muon_trigger_pt_weight)**2)


	weight_for_histos[0] = weight_pileup[0]*weight_btag_eff[0]*weight_top_pT[0]*weight_tracking[0]*weight_lep_ID[0]*weight_lep_iso[0]*weight_trig_eff[0]
	if options.sideband == 'yes' :
		weight_for_histos[0] = weight_pileup[0]*weight_btag_eff[0]*weight_top_pT[0]*weight_tracking[0]*weight_lep_ID[0]*weight_trig_eff[0]
	finalChiHist.Fill(plot_final_chi,weight_for_histos[0])
	pileupHist.Fill(pileup_events[0],weight_for_histos[0])
	#End the event
	endEvent(CUTFLOW)
	t.Fill()



##############################################################################################################
##												Closeout Processes						 					##
##############################################################################################################

#Plot plots, write branches, and close out the root file
f.cd()
t.Write()
for histo in allHistos :
	histo.Write()
f.Write()
f.Close()
	
