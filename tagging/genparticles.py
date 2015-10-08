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
parser.add_option('--file', metavar='F', type='string', action='store',
				  dest='directory',
				  help='') ## Sets which files to run on
(options, args) = parser.parse_args()

#Read in input files
FILES = options.directory + "/*.root"
files = glob.glob(FILES)
#print 'Getting these files : '
#for fname in files:
#	print fname
postfix = ""
events = Events(files)

#Set up handles and labels for use
#MC GenParticle variables
GenEventHandle = Handle( "GenEventInfoProduct" )
GenEventLabel = ( "generator", "" )
GenHandle = Handle( "vector<reco::GenParticle>" )
GenLabel = ( "prunedGenParticles", "" )

ids = [(2212,'p'),(1,'d'),(2,'u'),(3,'s'),(4,'c'),(5,'b'),(6,'t'),(11,'e'),(12,'ve'),(13,'mu'),(14,'vmu'),(15,'tau'),(16,'vtau'),(21,'g'),(22,'photon'),(23,'Z'),(24,'W')]
def getId(pdgid) :
	name = ''
	if pdgid < 0 :
		name = name+'~'
		pdgid = -1*pdgid
	for tup in ids :
		if tup[0] == pdgid :
			name = name+tup[1]
	if name == '~' or name == '' :
		name = name + str(pdgid)
	return name

count = 0
count_qq = 0
count_not_qq = 0
count_semilep = 0
count_not_semilep = 0
count_W_noJets = 0
count_W_bJets = 0
count_W_cJets = 0
count_W_LFJets = 0
for event in events :
	if count == 1000000 :
		break
	if count!=0 and count%50000 == 0 :
		print '-----------------------------------------------------------------------'
		print 'iteration: '+str(count)+''
		print 'semileptonic fraction : '+str(count_semilep)+' (%.4f%%)'%(100.0*count_semilep/count)
		print 'not semilep  fraction : '+str(count_not_semilep)+' (%.4f%%)'%(100.0*count_not_semilep/count)
		print 'qqbar     fraction : '+str(count_qq)+' (%.4f%%)'%(100.0*count_qq/count)
		print 'not qqbar fraction : '+str(count_not_qq)+' (%.4f%%)'%(100.0*count_not_qq/count)
		print 'total count from semilep check: '+str((count_semilep+count_not_semilep))+' (%.4f%%)'%(100.0*(count_semilep+count_not_semilep)/count)
		print 'total count from mother  check: '+str((count_qq+count_not_qq))+' (%.4f%%)'%(100.0*(count_qq+count_not_qq)/count)
		print '----------------W+Jets event types-----------------'
		print '---W + 0 Jets: '+str(count_W_noJets)+' (%.4f%%)'%(100.0*count_W_noJets/count)
		print '---W + b-Jets: '+str(count_W_bJets)+' (%.4f%%)'%(100.0*count_W_bJets/count)
		print '---W + c-Jets: '+str(count_W_cJets)+' (%.4f%%)'%(100.0*count_W_cJets/count)
		print '---W + LF-Jets: '+str(count_W_LFJets)+' (%.4f%%)'%(100.0*count_W_LFJets/count)
	count = count + 1

#    print '----------------------------------------------------------------'
#    print 'NEW EVENT'
	#open genEvent Information
	event.getByLabel(GenEventLabel,GenEventHandle)
	if GenEventHandle.isValid() :
		GenEvent = GenEventHandle.product()
		if GenEvent.weight() != 1.0 :
			print 'EVENT WEIGHT = '+str(GenEvent.weight())

#	#find out whether the event was qqbar or not
#	mother_ids = []
#	event.getByLabel( GenLabel, GenHandle )
#	lep_charge = 0
#	is_qq = False
#	is_semilep = False
#	if GenHandle.isValid() :
#		GenParticles = GenHandle.product()
#		for ig in GenParticles :
#			#check if it's a proton
#			if ig.pdgId() == 2212 :
#				mother_ids.append(ig.daughter(0).pdgId())
#			#is it a t or a tbar?
#			if math.fabs(ig.pdgId()) == 6 and ig.status() == 3 :
#				#look through all the daughters for Ws.
#				for i in range(ig.numberOfDaughters()) :
#					dau = ig.daughter(i)
#					if math.fabs(dau.pdgId()) == 24 :
#						#if the W doesn't have two daughters, I don't know what the hell happened.
#						if dau.numberOfDaughters() != 2 :
#							info = 'W without two daughters. PARTICLE: ' + getId(ig.pdgId()) + ', DAUGHTERS : '
#							for j in range(ig.numberOfDaughters()) :
#								info = info + getId(ig.daughter(i)) + ' '
#							print info
#							continue
#						#check if the W decayed leptonically. If it did, add the charge of the lepton to the running total.
#						dau1_is_lep = math.fabs(dau.daughter(0).pdgId()) == 11 or math.fabs(dau.daughter(0).pdgId()) == 13 or math.fabs(dau.daughter(0).pdgId()) == 15
#						dau2_is_lep = math.fabs(dau.daughter(1).pdgId()) == 11 or math.fabs(dau.daughter(1).pdgId()) == 13 or math.fabs(dau.daughter(1).pdgId()) == 15
#						if dau1_is_lep or dau2_is_lep :
#							if not is_semilep :
#								is_semilep = True
#							else :
#								is_semilep = False
#	#add event to counts of stuff
#	if is_semilep : 
#		count_semilep+=1
#	if not is_semilep :
#		count_not_semilep+=1
#	if len(mother_ids) != 2 :
#		print 'did not find two protons. daughters of protons: '
#		print mother_ids
#	else :
#		if math.fabs(mother_ids[0]) < 6 and math.fabs(mother_ids[1]) < 6 and mother_ids[0] + mother_ids[1] == 0 :
#			is_qq = True
#			count_qq+=1
#		else :
#			count_not_qq+=1

	event.getByLabel( GenLabel, GenHandle )
	if GenHandle.isValid() :
		GenParticles = GenHandle.product()
		extraJets = []
		#get the list of mothers of the W
		Wmothers = []
		for ig in GenParticles :
			if ig.pt()<0 or ig.status() != 3 :
				continue
			if math.fabs(ig.pdgId()) == 24 :
				Wmothers.append(ig.mother(0))
				Wmothers.append(ig.mother(1))
		#loop again to find extra jet particles
		for ig in GenParticles :
			if ig.pt()<0 or ig.status() != 3 :
				continue
			if math.fabs(ig.pdgId()) == 24 or (math.fabs(ig.pdgId()) > 10 and math.fabs(ig.pdgId()) < 19):
				continue
			if ig.numberOfDaughters() == 0 :
				extraJets.append(getId(ig.pdgId()))
			else :
				checkDaughters = []
				for i in range(ig.numberOfDaughters()) :
					if math.fabs(ig.daughter(i).pdgId()) < 11 or math.fabs(ig.daughter(i).pdgId()) > 18 or math.fabs(ig.daughter(i).pdgId()) == 24 :
						checkDaughters.append(getId(ig.daughter(i).pdgId()))
				if len(checkDaughters) == 0 :
					extraJets.append(getId(ig.pdgId()))

#	#Open genparticles information from the file
#	event.getByLabel( GenLabel, GenHandle )
#	if GenHandle.isValid() :
#		GenParticles = GenHandle.product()
#		print 'NEW EVENT'
#		for ig in GenParticles :
#			if ig.pt()<0 or ig.status() != 3 :
#				continue
#			ist = ig.pdgId() == 6 and ig.status() == 3
#			istbar = ig.pdgId() == -6 and ig.status() == 3
#			#if not (ist or istbar) :
#			#    continue
#			#print out a whole bunch of shit about this event
#			nMothers = ig.numberOfMothers()
#			nDaughters = ig.numberOfDaughters()
#			s = '[ '
#			for i in range(nMothers) :
#				s = s + getId(ig.mother(i).pdgId()) + ' '
#			s = s + '] '
#			if nMothers != 0 :
#				s = s + '-> '
#			s = s + '(' + getId(ig.pdgId()) + ') '
#			if nDaughters != 0 :
#				s = s + '-> { '
#			for i in range(nDaughters) :
#				s = s + getId(ig.daughter(i).pdgId()) + ' '
#			if nDaughters != 0 :
#				s = s + '}'
#			print s
##            for i in range(nDaughters) :
##                s = '   [ '
##                nOtherMothers   = ig.daughter(i).numberOfMothers()
##                nOtherDaughters = ig.daughter(i).numberOfDaughters()
##                for j in range(nOtherMothers) :
##                    s = s + getId(ig.daughter(i).mother(j).pdgId()) + ' '
##                s = s + '] '
##                if nOtherMothers != 0 :
##                    s = s + '-> '
##                s = s + '(' + getId(ig.daughter(i).pdgId()) + ') '
##                if nOtherDaughters != 0 :
##                    s = s + '-> { '
##                for j in range(nOtherDaughters) :
##                    s = s + getId(ig.daughter(i).daughter(j).pdgId()) + ' '
##                if nOtherDaughters != 0 :
##                    s = s + '}'
##                print s
#	
##	print 'is_semilep = '+str(is_semilep)
##	print 'is_qq = '+str(is_qq)
		
		#print 'EXTRA JETS: ' + str(extraJets)
		if 'b' in extraJets or '~b' in extraJets :
			count_W_bJets += 1
			#print 'event classified as W+b-Jets'
		elif 'c' in extraJets or '~c' in extraJets :
			count_W_cJets += 1
			#print 'event classified as W+c-Jets'
		elif len(extraJets) > 0 :
			count_W_LFJets += 1
			#print 'event classified as W+LF-Jets'
		else :
			count_W_noJets += 1
			#print 'event classified as W, no jets'
		#print ''





