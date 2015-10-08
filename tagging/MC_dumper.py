import os
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
parser.add_option('--file', metavar='F', type='string', action='store',dest='directory',help='') ## Sets which files to run on
parser.add_option('--out', metavar='F', type='string', action='store',dest='out',help='') ## Sets filename for outputfile
parser.add_option('--start', metavar='F', type='int', action='store',dest='start',default='0',help='') ## Sets starting event for run
parser.add_option('--end', metavar='F', type='int', action='store',dest='end',default='1000000000',help='') ## Sets ending event for run
(options, args) = parser.parse_args()

#Set up output file
outputname = options.out+".root"
f = TFile( outputname, "Recreate" )
t = TTree("output", "output")
t.SetDirectory(0)
#branches for reconstructed quantities
nParticles = 10 #lepton, neutrino, leptonic b, leptonic W, hadronic W subjet 1, hadronic W subjet 2, and hadronic b, and hadronic W
#plus the leptonic and hadronic t/tbar

#histograms of discriminant values
bDisc = TH1F('bDisc','CSV discriminant values for bs; CSV value',50,0.0,1.0)
WsubDisc = TH1F('WsubDisc','CSV discriminant values for jets from W decay; CSV value',50,0.0,1.0)
otherDisc = TH1F('otherDisc','CSV discriminant values for other jets; CSV value',50,0.0,1.0)
pileup = TH1F('pileup','Event Pileup; npvRealTrue',80,0,80)
bDisc.SetDirectory(0)
WsubDisc.SetDirectory(0)
otherDisc.SetDirectory(0)
pileup.SetDirectory(0)

#branches for MC truth quantities  
pt  = array('f',nParticles*[0.])
eta = array('f',nParticles*[0.])
phi = array('f',nParticles*[0.])
mass = array('f',nParticles*[0.])
px = array('f',nParticles*[0.])
py = array('f',nParticles*[0.])
pz = array('f',nParticles*[0.])
E  = array('f',nParticles*[0.])
PID  = array('f',nParticles*[0.])
n_btags = array('I',nParticles*[0])
n_jet_cands = array('I',nParticles*[0])
is_leptonic_side  = array('f',nParticles*[0.])
t.Branch('pt', pt,'pt[10]/F')
t.Branch('eta', eta,'eta[10]/F')
t.Branch('phi',phi,'phi[10]/F')
t.Branch('mass',mass,'mass[10]/F')
t.Branch('px',px,'px[10]/F')
t.Branch('py',py,'py[10]/F')
t.Branch('pz',pz,'pz[10]/F')
t.Branch('E',E,'E[10]/F')
t.Branch('PID',PID,'PID[10]/F')
t.Branch('is_leptonic_side',is_leptonic_side,'is_leptonic_side[10]/F')
t.Branch('n_btags',n_btags,'n_btags[10]/i')
t.Branch('n_jet_cands',n_jet_cands,'n_jet_cands[10]/i')
f.cd()

#Read in input files
FILES = options.directory + "/*.root"
files = glob.glob(FILES)
print 'Getting these files : '
for fname in files:
    print fname
postfix = ""
events = Events(files)

#Set up handles and labels for use
#jet variables
jetPtHandle       = Handle("std::vector<float>")
jetEtaHandle      = Handle("std::vector<float>")
jetPhiHandle      = Handle("std::vector<float>")
jetMassHandle     = Handle("std::vector<float>")
jetCSVHandle      = Handle("std::vector<float>")
npvRealTrueHandle = Handle('int')
jetPtLabel        = ("pfShyftTupleJets","pt")
jetEtaLabel       = ("pfShyftTupleJets","eta")
jetPhiLabel       = ("pfShyftTupleJets","phi")
jetMassLabel      = ("pfShyftTupleJets","mass")
jetCSVLabel       = ("pfShyftTupleJets","csv")
npvRealTrueLabel  = ('pileup','npvRealTrue')
#MC GenParticle variables
GenHandle = Handle( "vector<reco::GenParticle>" )
GenLabel = ( "prunedGenParticles", "" )

print 'Files opened'

count = 0
count_qq = 0
count_not_qq = 0
count_events = 0
for event in events :
    if count == 0 :
        print 'Loop started'
        print 'Total number of events get is: '+str(event.size())
    if count == options.end :
        break
    count = count + 1
    if count < options.start :
        continue
    run_event = options.end
    if count%(run_event/100) == 0 :
        print ''+str(100*count/run_event)+'% done'
        if count_events!=0 :
            print '# of qqbar  events: '+str(count_qq)+' (%.4f%%)'%(100.0*count_qq/count_events)
            print '# of not qq events: '+str(count_not_qq)+' (%.4f%%)'%(100.0*count_not_qq/count_events)

    #put pileup value in histogram
    event.getByLabel(npvRealTrueLabel,npvRealTrueHandle)
    if not npvRealTrueHandle.isValid() :
        continue
    pileup_value = npvRealTrueHandle.product()
    pileup.Fill(pileup_value[0])

    #loop through all jets
    event.getByLabel(jetPtLabel,jetPtHandle)
    event.getByLabel(jetEtaLabel,jetEtaHandle)
    event.getByLabel(jetPhiLabel,jetPhiHandle)
    event.getByLabel(jetMassLabel,jetMassHandle)
    event.getByLabel(jetCSVLabel,jetCSVHandle)
    event.getByLabel(GenLabel,GenHandle)
    if not jetPtHandle.isValid() or not GenHandle.isValid() :
        continue
    jetPts   = jetPtHandle.product()
    jetEtas  = jetEtaHandle.product()
    jetPhis  = jetPhiHandle.product()
    jetMasss = jetMassHandle.product()
    jetCSVs  = jetCSVHandle.product()
    GenParticles = GenHandle.product()
    for i in range(len(jetPts)) :
        #if jetPts[i] < 30 :
        #    print jetPts[i]
        if jetPts[i] < 20. or abs(jetEtas[i])>2.5 :
            continue
        thisjet = TLorentzVector()
        thisjet.SetPtEtaPhiM(jetPts[i],jetEtas[i],jetPhis[i],jetMasss[i])
        #look through all the MC jets to match it
        for particle in GenParticles :
            pid = abs(particle.pdgId())
            #only need u,d,s,c,g, and b jets
            wanted = [1,2,3,4,21,5]
            if not pid in wanted :
                continue
            #see if it matches
            thisparticle = TLorentzVector()
            thisparticle.SetPtEtaPhiM(particle.pt(),particle.eta(),particle.phi(),particle.mass())
            if thisparticle.DeltaR(thisjet) < 0.3 :
                if pid == 5 :
                    bDisc.Fill(jetCSVs[i])
                else :
                    isWsub = False
                    for j in range(particle.numberOfMothers()) :
                        if abs(particle.mother(j).pdgId()) == 24 :
                            isWsub = True
                            break
                    if isWsub :
                        WsubDisc.Fill(jetCSVs[i])
                    else :
                        otherDisc.Fill(jetCSVs[i])




#    ######################################################
#    ##                  Event Type Cut                  ##
#    ######################################################
#    #Monte Carlo leptons, neutrinos, leptonic bs, hadronic bs, leptonic Ws, and hadronic Ws
#    particles = []
#    particles_found = []
#    t_is_leptonic_side = False
#    for i in range(10) :
#        particles.append(ROOT.TLorentzVector())
#        particles_found.append(False)
#    wsubjets = []
#    #Open genparticles information from the file
#    event.getByLabel( GenLabel, GenHandle )
#    if GenHandle.isValid() :
#        GenParticles = GenHandle.product()
#        #Get the damn lepton to set its charge
#        leptonCharge = 0
#        for ig in GenParticles :
#            if ig.pt()<0 :
#                continue
#            ist = ig.pdgId() == 6 and ig.status() == 3
#            istbar = ig.pdgId() == -6 and ig.status() == 3
#            if ist or istbar :
#                for i in range(ig.numberOfDaughters()) :
#                    dau = ig.daughter(i)
#                    if dau.pt() < 0 :
#                        continue
#                    isW = ist and math.fabs(dau.pdgId()) == 24 and dau.status() == 3
#                    if isW :
#                        for j in range(dau.numberOfDaughters()) :
#                            subDau = dau.daughter(j)
#                            if subDau.pt() < 0 :
#                                continue
#                            islep = math.fabs(subDau.pdgId()) == 13 # or math.fabs(subDau.pdgId()) == 11 or math.fabs(subDau.pdgId()) == 15
#                            if islep :
#                                leptonCharge = -1*(subDau.pdgId()/math.fabs(subDau.pdgId())) #finally.
#        #find out whether the t or tbar is the leptonic side
#        if leptonCharge == 0 :
#            continue
#        t_is_leptonic_side = leptonCharge>0#

#    jetCands = []
#    ######################################################
#    ##                  Jet Selection                   ##
#    ######################################################
#    event.getByLabel (jetPtLabel, jetPtHandle)
#    if jetPtHandle.isValid() :
#        toomanyjetPts = jetPtHandle.product()    
#        event.getByLabel (jetEtaLabel, jetEtaHandle)
#        toomanyjetEtas = jetEtaHandle.product()
#        event.getByLabel (jetPhiLabel, jetPhiHandle)
#        toomanyjetPhis = jetPhiHandle.product()
#        event.getByLabel (jetMassLabel, jetMassHandle)
#        toomanyjetMasss = jetMassHandle.product()
#        event.getByLabel (jetCSVLabel, jetCSVHandle)
#        toomanyjetCSVs = jetCSVHandle.product() 
#        #begin selection: order jets by Pt
#        dummypts = []
#        for i in range(len(toomanyjetPts)) :
#            pt_tuple = -1*toomanyjetPts[i],i
#            dummypts.append(pt_tuple)
#        dummypts.sort()
#        #pt and eta cuts
#        jetCandCSVs = []
#        for i in range(len(dummypts)) :
#            j = dummypts[i][1]
#            if (toomanyjetPts[j] > 20.0 and math.fabs(toomanyjetEtas[j]) < 2.5) :
#                jetCands.append(ROOT.TLorentzVector())
#                jetCands[len(jetCands)-1].SetPtEtaPhiM(toomanyjetPts[j],toomanyjetEtas[j],toomanyjetPhis[j],toomanyjetMasss[j])
#                jetCandCSVs.append(toomanyjetCSVs[j])
#        #require more than 3 jet candidates
#        if len(jetCands) < 4 :
#            continue
#        #require fewer than 6 jet candidates
#        if len(jetCands) > 5 :
#            continue
#        #Cut that we have at least the right number of btagged jets with pt > 20
#        nbtags = 0
#        for i in range(len(jetCands)) :
#            if jetCands[i].Pt()>20.0 :
#                if jetCandCSVs[i] > 0.679 :
#                    nbtags = nbtags + 1
#        if nbtags < 1 :
#            continue
#        #make cuts on jetCand pt (orthogonal to selection in tagger)
#        if (jetCands[0].Pt()>45.0 and jetCands[1].Pt()>35.0 and jetCands[2].Pt()>20.0 and jetCands[3].Pt()>20.0) :
#            continue
#        for i in range(nParticles) :
#            n_btags[i] = nbtags
#            n_jet_cands[i] = len(jetCands)
#    else :
#        continue
#    count_events = count_events + 1#

#    ############################################################
#    ##            MC truth dumping for all particles          ##
#    ############################################################
#    #count the number of ts/tbars coming from gluons and from quarks
#    original_particles = []
#    #loop again
#    for ig in GenParticles :
#        if ig.pt()<0 :
#            continue
#        #Look for all protons and add their first daughters to the list of original particles
#        if ig.pdgId() == 2212 :
#            original_particles.append(ig.daughter(0).pdgId())
#        #Look through particles for all ts
#        ist = ig.pdgId() == 6 and ig.status() == 3
#        istbar = ig.pdgId() == -6 and ig.status() == 3
#        #Add ts to list and then search through their daughters for bs and Ws (everything also done for charge conjugates)
#        if ist or istbar:
#            #Find out which type of semileptonic ttbar event is going on
#            if (ist and t_is_leptonic_side) or (istbar and not t_is_leptonic_side) : #it's the leptonic t
#                particles[8].SetPtEtaPhiM(ig.pt(),ig.eta(),ig.phi(),ig.mass())
#                particles_found[8] = True
#                PID[8]=ig.pdgId()
#            elif (istbar and t_is_leptonic_side) or (ist and not t_is_leptonic_side) : #it's the hadronic t
#                particles[9].SetPtEtaPhiM(ig.pt(),ig.eta(),ig.phi(),ig.mass())
#                particles_found[9] = True
#                PID[9]=ig.pdgId()
#            for i in range(ig.numberOfDaughters()) :
#                dau = ig.daughter(i)
#                if dau.pt() < 0 :
#                    continue
#                isb = ist and dau.pdgId() == 5 and dau.status() == 3
#                isbbar = istbar and dau.pdgId() == -5 and dau.status() == 3
#                isWplus = ist and dau.pdgId() == 24 and dau.status() == 3
#                isWminus = istbar and dau.pdgId() == -24 and dau.status() == 3
#                #Add any bs found to the approriate list depending on whether the t is the hadronic or leptonic side of the decay
#                if (isb and t_is_leptonic_side) or (isbbar and not t_is_leptonic_side) : #it's the leptonic b
#                    particles[2].SetPtEtaPhiM(dau.pt(),dau.eta(),dau.phi(),dau.mass())
#                    particles_found[2] = True
#                    PID[2]=dau.pdgId()
#                if (isb and not t_is_leptonic_side) or (isbbar and t_is_leptonic_side) : #it's the hadronic b
#                    particles[6].SetPtEtaPhiM(dau.pt(),dau.eta(),dau.phi(),dau.mass())
#                    particles_found[6] = True
#                    PID[6]=dau.pdgId()
#                #set any Ws found to be the leptonic or hadronic Ws appropriately and look past the leptonic W for leptons and neutrinos
#                #and past the hadronic W for its subjets
#                if (isWplus and not t_is_leptonic_side) or (isWminus and t_is_leptonic_side) : #it's the hadronic W
#                    particles[7].SetPtEtaPhiM(dau.pt(),dau.eta(),dau.phi(),dau.mass())
#                    particles_found[7] = True
#                    PID[7]=dau.pdgId()
#                    for j in range(dau.numberOfDaughters()) :
#                        subDau = dau.daughter(j)
#                        if subDau.pt() < 0 :
#                            continue
#                        if subDau.status() == 3 : #Add any valid daughter of the W to the list of W daughters
#                            wsubjets.append((subDau.pt(),ROOT.TLorentzVector()))
#                            wsubjets[len(wsubjets)-1][1].SetPtEtaPhiM(subDau.pt(),subDau.eta(),subDau.phi(),subDau.mass())
#                if (isWplus and t_is_leptonic_side) or (isWminus and not t_is_leptonic_side) : #it's the leptonic W
#                    particles[3].SetPtEtaPhiM(dau.pt(),dau.eta(),dau.phi(),dau.mass())
#                    particles_found[3] = True
#                    PID[3]=dau.pdgId()
#                    for j in range(dau.numberOfDaughters()) :
#                        subDau = dau.daughter(j)
#                        if subDau.pt() < 0 :
#                            continue
#                        islep = (math.fabs(subDau.pdgId()) == 11 or math.fabs(subDau.pdgId()) == 13 or math.fabs(subDau.pdgId()) == 17) and subDau.status() == 3
#                        isnv = (math.fabs(subDau.pdgId()) == 12 or math.fabs(subDau.pdgId()) == 14 or math.fabs(subDau.pdgId()) == 18) and subDau.status() == 3
#                        #Add any found leptons and neutrinos to the lists
#                        if islep :
#                            particles[0].SetPtEtaPhiM(subDau.pt(),subDau.eta(),subDau.phi(),subDau.mass())
#                            particles_found[0] = True
#                            PID[0]=subDau.pdgId()
#                        if isnv :
#                            particles[1].SetPtEtaPhiM(subDau.pt(),subDau.eta(),subDau.phi(),subDau.mass())
#                            particles_found[1] = True
#                            PID[1]=subDau.pdgId()
#    #find out the type of semileptonic ttbar event
#    is_qq = len(original_particles) == 2 and original_particles[0] < 6 and original_particles[1] < 6 and original_particles[0] + original_particles[1] == 0
#    if is_qq :
#        count_qq += 1
#    else :
#        count_not_qq += 1
#    #order the W subjets based on pt
#    wsubjets.sort()
#    particles[4] = wsubjets[0][1]
#    particles_found[4] = True
#    PID[4] = 0
#    PID[5] = 0
#    particles[5] = wsubjets[1][1]
#    particles_found[5] = True
#    dummy_is_leptonic_side = [1,1,1,1,0,0,0,0,1,0]
#    #fill out all of the particles
#    for i in range(len(particles)) :
#        if particles_found[i] :
#            pt[i] = particles[i].Pt()
#            eta[i] = particles[i].Eta()
#            phi[i] = particles[i].Phi()
#            mass[i] = particles[i].M()
#            px[i] = particles[i].Px()
#            py[i] = particles[i].Py()
#            pz[i] = particles[i].Pz()
#            E[i] = particles[i].E()
#            is_leptonic_side[i] = dummy_is_leptonic_side[i]
#    #Match and fill histograms for leptonic, hadronic bs and w subjets
#    for i in range(len(jetCands)) :
#        if jetCands[i].DeltaR(particles[2])<0.3 :
#            blepDisc.Fill(jetCandCSVs[i])
#        if jetCands[i].DeltaR(particles[6])<0.3 :
#            bhadDisc.Fill(jetCandCSVs[i])
#        if jetCands[i].DeltaR(particles[4])<0.3 or jetCands[i].DeltaR(particles[5])<0.3 :
#            WhadDisc.Fill(jetCandCSVs[i])#

#    t.Fill()

f.cd()
#t.Write()
bDisc.Write()
WsubDisc.Write()
otherDisc.Write()
pileup.Write()
f.Write()
f.Close()
