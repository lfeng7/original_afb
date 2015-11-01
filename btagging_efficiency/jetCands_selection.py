##############################################################################################################
##	  Data and Monte Carlo Tagger/Matcher for Top Quark/Antiquark Forward/Backward Asymmetry Study at CMS	##
## 									Morris Swartz, Raymond Feng, Nick Eminizer								##
##											   		Fall 2013												##
##############################################################################################################

####### Editing Log ##########
### Edited on 5/4/14 by Raymond Feng 
### Add option Pdf_Version (cteq66 or CT10). Adding PDF_weights to ntuples as a branch called Pdf_weights
### Based on newer version of Nick's tagger. Get on 5/16/14
### Edited on 5/16/14 by Raymond
### Add GenEventInfo, with qScale and x1,x2,id1,id2 for PDF reweighting later.
### Edited 6/15/14 by Raymond
### Add components to generate btagging SF correction given pt,eta,flavor and csv values
### btag SF weight is calculated based on candidate jets only, after all selection cuts are
### applied. Reference https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG

##############################################################################################################
##										Imports and Other Header Stuff					 					##
##############################################################################################################
import os
import glob
import math
import ROOT
from ROOT import *
import sys
from DataFormats.FWLite import Events, Handle
from optparse import OptionParser
from array import array


# Read efficiency files
F_eff = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/analysis_new/EDSHyFT/data/TT_CT10_TuneZ2star_8TeV-powheg-tauola_AK5PF_CSVM_bTaggingEfficiencyMap.root '
file_tmp = TFile(F_eff)
efficiency_b = file_tmp.Get('efficiency_b')
efficiency_c = file_tmp.Get('efficiency_c')
efficiency_udsg = file_tmp.Get('efficiency_udsg') 
# Global vars
count_out_of_bounds = 0 
count_jets = 0

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--dir', metavar='F', type='string', action='store',dest='directory',help='') ## Sets which files to run on
parser.add_option('--out', metavar='F', type='string', action='store',dest='out',help='') ## Sets filename for outputfile
parser.add_option('--mc', metavar='F', type='string', action='store',dest='mc',help='') ## Is the file a MC
parser.add_option('--lep', metavar='F', type='string', action='store',default='mu',dest='lep',help='') ## which leptons do we use
parser.add_option('--dR', metavar='F', type='float', action='store',default='0.3',dest='dR',help='') ## choose the dR requirement for matching
parser.add_option('--tMotherCheck', metavar='F', type='string', action='store',default='none',dest='tMotherCheck',help='') 
							## "none" for no restrictions, 
							# "qq" for qq (ttbar mother IDs < 6 and flavors match), "gg" for everything not qqbar.
parser.add_option('--semilep_only', metavar='F', type='string', action='store',default='neither',dest='semilep_only',help='') 
							## set to "yes" if sample contains more than semileptonic ttbar but that's all you want
							## set to "no" if sample contains more than semileptonic ttbar but all you want is background
parser.add_option('--lepton_truth', metavar='F', type='string', action='store',default='none',dest='lepton_truth',help='') 
							## set to "match" if you want only events with the relevant lepton from truth.
						   ## set to "sideband" if you want only events with the wrong lepton from truth.
						   ## leave at "none" if you want all lepton types in the sample.
						   ## ONLY USE FOR SEMILEPTONIC EVENTS OR THE CODE WILL CRASH
parser.add_option('--nMaxEvents', metavar='F', type='int', action='store',default='-1',dest='nMaxEvents',help='') ## choose the maximum number of iterations to go through
parser.add_option('--printEvery', metavar='F', type='int', action='store',default='1000',dest='printEvery',help='') ## choose the dR requirement for matching
parser.add_option('--nTotalJobs', metavar='F', type='int', action='store',default='1',dest='nTotalJobs',help='') # Total number of jobs that will be run on the grid
parser.add_option('--iJob', metavar='F', type='int', action='store',default='0',dest='iJob',help='') # Offset from 0 for this job on the grid
parser.add_option('--leptonCut', metavar='F', type='string', action='store',default='tight',dest='leptonCut',help='') # Use tight or loose leptons
parser.add_option('--sample_type', metavar='F', type='string', action='store',default='none',dest='sample_type',help='') #Type of sample, needed to apply background scale factors
# Options for different PDF set
parser.add_option('--PdfVersion', metavar='F', type='string', action='store',default='CT10',dest='PdfVersion',help='')
(options, args) = parser.parse_args()
if options.iJob < 0 or options.iJob >= options.nTotalJobs :
	print 'CHECK NUMBERING CONVENTION FOR GRID JOBS: SOMETHING IS VERY WRONG AND NOTHING IS GOOD!'


MMUON = 0.105658
MELECTRON = 0.000511

#Set up output file
outputname = options.out+'.root'
f = TFile( outputname, "Recreate" )
t = TTree("output", "output")
t.SetDirectory(0)
#branches for reconstructed quantities
nParticles = 5 # 4 or 5 jet candidates
pt = array('f',nParticles*[0.])
eta = array('f',nParticles*[0.])
phi = array('f',nParticles*[0.])
mass = array('f',nParticles*[0.])
px = array('f',nParticles*[0.])
py = array('f',nParticles*[0.])
pz = array('f',nParticles*[0.])
E  = array('f',nParticles*[0.])
PID = array('f',nParticles*[0.])
CSV = array('f',nParticles*[0.])
nbTags = array('I',[0])
nValidJets = array('I',[0])
t.Branch('pt', pt,'pt[5]/F')
t.Branch('eta', eta,'eta[5]/F')
t.Branch('phi',phi,'phi[5]/F')
t.Branch('mass',mass,'mass[5]/F')
t.Branch('px',px,'px[5]/F')
t.Branch('py',py,'py[5]/F')
t.Branch('pz',pz,'pz[5]/F')
t.Branch('E',E,'E[5]/F')
t.Branch('PID',PID,'PID[5]/F')
t.Branch('CSV',CSV,'CSV[5]/F')
t.Branch('nbTags',nbTags,'nbTags/i')
t.Branch('nValidJets',nValidJets,'nValidJets/i')

f.cd()

#Read in input files
FILES = options.directory + '/*.root'
print FILES
files = glob.glob(FILES)
print 'Getting these files : '
for fname in files:
	print fname
postfix = ""
events = Events(files)

#Width and Central Value of Wmass peak
wmass = 80.08
wwidth = 33.5

#Set up handles and labels for use
#muon variables
if options.lep == "mu":
	if options.leptonCut == 'tight' :
		leptonPtHandle         = Handle( "std::vector<float>" )
		leptonPtLabel    = ( "pfShyftTupleMuons" + postfix,   "pt" )
		leptonEtaHandle         = Handle( "std::vector<float>" )
		leptonEtaLabel    = ( "pfShyftTupleMuons" + postfix,   "eta" )
		leptonPhiHandle         = Handle( "std::vector<float>" )
		leptonPhiLabel    = ( "pfShyftTupleMuons" + postfix,   "phi" )
		leptonPfisoHandle         = Handle( "std::vector<float>" )
		leptonPfisoLabel    = ( "pfShyftTupleMuons" + postfix,   "pfisoPU" )
		leptonChargeHandle         = Handle( "std::vector<float>" )
		leptonChargeLabel    = ( "pfShyftTupleMuons" + postfix,   "charge" )
	elif options.leptonCut == 'loose' :
		leptonPtHandle         = Handle( "std::vector<float>" )
		leptonPtLabel    = ( "pfShyftTupleMuonsLoose" + postfix,   "pt" )
		leptonEtaHandle         = Handle( "std::vector<float>" )
		leptonEtaLabel    = ( "pfShyftTupleMuonsLoose" + postfix,   "eta" )
		leptonPhiHandle         = Handle( "std::vector<float>" )
		leptonPhiLabel    = ( "pfShyftTupleMuonsLoose" + postfix,   "phi" )
		leptonPfisoHandle         = Handle( "std::vector<float>" )
		leptonPfisoLabel    = ( "pfShyftTupleMuonsLoose" + postfix,   "pfisoPU" )
		leptonChargeHandle         = Handle( "std::vector<float>" )
		leptonChargeLabel    = ( "pfShyftTupleMuonsLoose" + postfix,   "charge" )
	MLEP = MMUON
#electron variables
if options.lep == "el":
	if options.leptonCut == 'tight' :
		leptonPtHandle         = Handle( "std::vector<float>" )
		leptonPtLabel    = ( "pfShyftTupleElectrons" + postfix,   "pt" )
		leptonEtaHandle         = Handle( "std::vector<float>" )
		leptonEtaLabel    = ( "pfShyftTupleElectrons" + postfix,   "eta" )
		leptonPhiHandle         = Handle( "std::vector<float>" )
		leptonPhiLabel    = ( "pfShyftTupleElectrons" + postfix,   "phi" )
		leptonPfisoHandle         = Handle( "std::vector<float>" )
		leptonPfisoLabel    = ( "pfShyftTupleElectrons" + postfix,   "pfiso" )
		leptonChargeHandle         = Handle( "std::vector<float>" )
		leptonChargeLabel    = ( "pfShyftTupleElectrons" + postfix,   "charge" )
	elif options.leptonCut == 'loose' :
		leptonPtHandle         = Handle( "std::vector<float>" )
		leptonPtLabel    = ( "pfShyftTupleElectronsLoose" + postfix,   "pt" )
		leptonEtaHandle         = Handle( "std::vector<float>" )
		leptonEtaLabel    = ( "pfShyftTupleElectronsLoose" + postfix,   "eta" )
		leptonPhiHandle         = Handle( "std::vector<float>" )
		leptonPhiLabel    = ( "pfShyftTupleElectronsLoose" + postfix,   "phi" )
		leptonPfisoHandle         = Handle( "std::vector<float>" )
		leptonPfisoLabel    = ( "pfShyftTupleElectronsLoose" + postfix,   "pfiso" )
		leptonChargeHandle         = Handle( "std::vector<float>" )
		leptonChargeLabel    = ( "pfShyftTupleElectronsLoose" + postfix,   "charge" )
	MLEP = MELECTRON
#neutrino (MET) variables
if options.leptonCut == 'tight' :
	metPtHandle = Handle( "std::vector<float>" )
	metPtLabel = ("pfShyftTupleMET",   "pt" )
	metPhiHandle         = Handle( "std::vector<float>" )
	metPhiLabel    = ( "pfShyftTupleMET",   "phi" )
elif options.leptonCut == 'loose' :
	metPtHandle = Handle( "std::vector<float>" )
	metPtLabel = ("pfShyftTupleMETLoose",   "pt" )
	metPhiHandle         = Handle( "std::vector<float>" )
	metPhiLabel    = ( "pfShyftTupleMETLoose",   "phi" )
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
# Jet Flavor handle
jetFlavorHandle = Handle("vector<float>")
jetFlavorLabel = ("pfShyftTupleJets"+postfix,"jetFlavor")	
#MC GenParticle variables
if options.mc == "yes":
	GenHandle = Handle( "vector<reco::GenParticle>" )
	GenLabel = ( "prunedGenParticles", "" )

print 'Files opened, starting analysis. . . .'

count = 0
count_all_cuts_but_btag_cuts = 0
count_all_cuts = 0
realCount = 0

count_nobtag=count_onebtag=count_twobtags=count_threebtags=count_fourbtags=0

##############################################################################################################
##											Main Loop Over Events							 				##
##############################################################################################################
for event in events:

	######################################################
	##	  Grid Split, check event type, check hard cut	##
	######################################################
	#options for running on the grid: split whole sample "in parallel"
	realCount = realCount + 1
	if ((realCount-1)-options.iJob) % options.nTotalJobs != 0 :
		continue

	count = count + 1

	#Hard cut if nMaxEvents isn't -1 for "all data"
	if count == options.nMaxEvents + 1 :
		print 'ITERATION '+str(options.nMaxEvents)+' REACHED'
		break

	#################### REPORT PROGRESS ####################
	if count % options.printEvery == 1 :
		print 'entry #'+str(count)
		print 'events passed all but btag cuts : '+str(count_all_cuts_but_btag_cuts)
		print 'events passed all cuts : '+str(count_all_cuts)
	#################### Hard cut off #######################

	#Check to make sure the event type is correct
	original_particles = []
	orig_partons = []
	lep_ids = []
	if options.mc == 'yes' :
		event.getByLabel( GenLabel, GenHandle )
		if GenHandle.isValid() :
			GenParticles = GenHandle.product()
			for ig in GenParticles :
				if ig.pt()<0 :
					continue
				#Look through all the particles for protons; append their first daughters to the list
				if ig.pdgId() == 2212 :
					original_particles.append(ig.daughter(0).pdgId())
					for i in range(ig.numberOfDaughters()):
						orig_partons.append(ig.daughter(i).pdgId())
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
							dau1_is_lep = math.fabs(dau.daughter(0).pdgId()) == 11 or math.fabs(dau.daughter(0).pdgId()) == 13 or math.fabs(dau.daughter(0).pdgId()) == 15
							dau2_is_lep = math.fabs(dau.daughter(1).pdgId()) == 11 or math.fabs(dau.daughter(1).pdgId()) == 13 or math.fabs(dau.daughter(1).pdgId()) == 15
							if dau1_is_lep :
								lep_ids.append(dau.daughter(0).pdgId())
							if dau2_is_lep :
								lep_ids.append(dau.daughter(1).pdgId())
			is_qq = len(original_particles) == 2 and original_particles[0] < 6 and original_particles[1] < 6 and original_particles[0] + original_particles[1] == 0
			is_semilep = len(lep_ids) == 1
			if options.semilep_only == 'yes' and not is_semilep :
				continue
			if options.semilep_only == 'no' and is_semilep :
				continue
			if options.lepton_truth == 'match' and ((options.lep == 'mu' and math.fabs(lep_ids[0]) != 13) or (options.lep == 'el' and math.fabs(lep_ids[0]) != 11)):
				continue
			if options.lepton_truth == 'sideband' and ((options.lep == 'mu' and math.fabs(lep_ids[0]) == 13) or (options.lep == 'el' and math.fabs(lep_ids[0]) == 11)):
				continue


	######################################################
	##		Check for valid event, set variables 		##
	##			    for selection histos				##
	######################################################
		
	#get variables from pat-tuple, make sure event is complete
	#leptons
	event.getByLabel (leptonPtLabel, leptonPtHandle)
	if not leptonPtHandle.isValid() :
		continue
	leptonPts = leptonPtHandle.product()
	if len(leptonPts) < 1 :
		continue
	lpt = leptonPts[0]
	event.getByLabel (leptonEtaLabel, leptonEtaHandle)
	leptonEtas = leptonEtaHandle.product()
	leta = leptonEtas[0]
	event.getByLabel (leptonPhiLabel, leptonPhiHandle)
	leptonPhis = leptonPhiHandle.product()
	event.getByLabel (leptonPfisoLabel, leptonPfisoHandle)
	leptonPfisos = leptonPfisoHandle.product()
	pfiso = leptonPfisos[0] / lpt
	events.getByLabel(leptonChargeLabel,leptonChargeHandle)
	leptonCharges = leptonChargeHandle.product()
	dummylep = ROOT.TLorentzVector()
	dummylep.SetPtEtaPhiM(lpt,leta,leptonPhis[0],MLEP)
	#met
	metPt = 0
	metPhi = 0
	event.getByLabel (metPtLabel, metPtHandle)
	if not metPtHandle.isValid() :
		continue
	metPt = metPtHandle.product()
	metpt = metPt[0]
	event.getByLabel (metPhiLabel, metPhiHandle)
	metPhi = metPhiHandle.product()
	dummyMET = ROOT.TLorentzVector()
	dummyMET.SetPtEtaPhiM(metpt,0.0,metPhi[0],0.0)
	#jets
	event.getByLabel (jetPtLabel, jetPtHandle)
	if not jetPtHandle.isValid() :
		continue
	toomanyjetPts = jetPtHandle.product()
	if len(toomanyjetPts) < 4 :
		continue	
	event.getByLabel (jetEtaLabel, jetEtaHandle)
	toomanyjetEtas = jetEtaHandle.product()
	event.getByLabel (jetPhiLabel, jetPhiHandle)
	toomanyjetPhis = jetPhiHandle.product()
	event.getByLabel (jetMassLabel, jetMassHandle)
	toomanyjetMasss = jetMassHandle.product()
	event.getByLabel (jetCSVLabel, jetCSVHandle)
	toomanyjetCSVs = jetCSVHandle.product()
	# Get jet flavors
	event.getByLabel (jetFlavorLabel,jetFlavorHandle)
	toomanyjetFlavors = jetFlavorHandle.product()		
	#order jets by pt
	dummypts = []
	for i in range(len(toomanyjetPts)) :
		pt_tuple = -1*toomanyjetPts[i],i
		dummypts.append(pt_tuple)
	dummypts.sort()

	lepton = ROOT.TLorentzVector()
	leptonCharge=0
	######################################################
	##				  lepton selection					##
	######################################################
	nLepVal = 0
	nLepExtra = 0
	for i in range(len(leptonPts)):
		#see if the lepton is valid but for a lower pt: help exclude dilepton events
		if options.lep == "mu" :
			if (leptonPts[i] < 26.0 and leptonPts[i] > 10.0 and math.fabs(leptonEtas[i]) < 2.5) :
				nLepExtra = nLepExtra + 1
				continue
			#pt and eta cut, want large muon within detector range
			if (leptonPts[i] < 26.0 or math.fabs(leptonEtas[i]) > 2.1) :
				continue
			#isolation cut
			if leptonPfisos[i] / leptonPts[i] > 0.12 :
				continue
			lepindex = i
			nLepVal = nLepVal + 1
		elif options.lep == "el" :
			if (leptonPts[i] < 30.0 and leptonPts[i] > 15.0 and math.fabs(leptonEtas[i]) < 2.5 and not (math.fabs(leptonEtas[i]) > 1.44 and math.fabs(leptonEtas[i]) < 1.56)) :
				nLepExtra = nLepExtra + 1
				continue
			#pt and eta cut, want large electron within detector range (sensitivity degraded in transition region)
			if (leptonPts[i] < 30.0 or math.fabs(leptonEtas[i])>2.5 or (math.fabs(leptonEtas[i]) > 1.44 and math.fabs(leptonEtas[i]) < 1.56)) :
				continue
			#isolation cut
			if leptonPfisos[i] / leptonPts[i] > 0.10 :
				continue
			lepindex = i
			nLepVal = nLepVal + 1
	if nLepVal != 1 or nLepExtra > 0: #if there's more or less than 1 good lepton, or if there's an extra almost-valid lepton
		continue

	jetCands = []
	######################################################
	##   Jet Selection (AK5, check W helicity paper)    ##
	######################################################
	#get the number of btags (just for interest, selection hasn't begun yet)
	numberOfbTags = 0
	for csv in toomanyjetCSVs :
		if csv > 0.679 :
			numberOfbTags = numberOfbTags + 1
	nbTags[0] = numberOfbTags	
	#pt and eta cuts
	jetCandCSVs = []
	jetCandFlavors = []
	for i in range(len(dummypts)) :
		j = dummypts[i][1]
		if toomanyjetPts[j] > 20. and math.fabs(toomanyjetEtas[j]) < 2.5 :
			jetCands.append(ROOT.TLorentzVector())
			jetCands[len(jetCands)-1].SetPtEtaPhiM(toomanyjetPts[j],toomanyjetEtas[j],toomanyjetPhis[j],toomanyjetMasss[j])
			jetCandCSVs.append(toomanyjetCSVs[j])
			jetCandFlavors.append(toomanyjetFlavors[j])
	#require at least 4 jet candidates
	if len(jetCands) < 4 :
		continue
	#require no more than 5 jet candidates
	if len(jetCands) > 5 :
		continue
	if len(jetCands) == 4 :
		#Record that the event had four jets
		nValidJets[0] = 4
	elif len(jetCands) == 5 :
		#Record that the event had five jets
		nValidJets[0] = 5
	#cut on pts of hardest 4 jets (cuts from A_Q paper)
	if (jetCands[0].Pt()<45.0 or jetCands[1].Pt()<35.0 or jetCands[2].Pt()<20.0 or jetCands[3].Pt()<20.0) :
		continue

	########################### All cuts but cut on b-tagging are done here. ###################################
	count_all_cuts_but_btag_cuts += 1

	# #fill the raw pts of the four hardest jets
	# print 'event # '+str(count)
	# print len(jetCands)
	# print nValidJets[0]

	##########################  reinitiate branch vars ###########################
	for i in range(len(pt)) :
		pt[i] = eta[i] = phi[i] = mass[i] = px[i] = py[i] = pz[i] = E[i] = PID[i] = CSV [i] = 0

	for i in range(len(jetCands)) :
		pt[i] = jetCands[i].Pt()
		eta[i] = jetCands[i].Eta()
		phi[i]=jetCands[i].Phi()
		mass[i]=jetCands[i].M()
		px[i]=jetCands[i].Px()
		py[i]=jetCands[i].Py()
		pz[i]=jetCands[i].Pz()
		E[i]=jetCands[i].E()
		PID[i]=jetCandFlavors[i]
		CSV[i]=jetCandCSVs[i]
		
	t.Fill()



##############################################################################################################
##												Closeout Processes						 					##
##############################################################################################################

#Plot plots, write branches, and close out the root file
f.cd()
t.Write()
f.Write()
f.Close()
file_tmp.Close()
	