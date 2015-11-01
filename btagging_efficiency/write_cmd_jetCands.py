#dir_input = '/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/training_set_ntuples/W4Jets'
#output = 'WJets_JetCands_before_btag_cuts_'
dir_input = '/eos/uscms/store/user/eminizer/nTuples/muons/TTJets_SemiLep'
output = 'TTJets_SemiLep_'
num_jobs = 200 
for i in range(num_jobs):
#	cmd = 'python ./tardir/jetCands_selection.py --out '+ output +str(i)+' --nTotalJobs '+str(num_jobs)+' --iJob '+str(i)+'  --mc yes --dir '+dir_input
	cmd = 'python ./tardir/jetCands_selection.py --out '+ output +str(i)+' --nTotalJobs '+str(num_jobs)+' --iJob '+str(i)+' --mc yes --semilep_only yes  --dir '+ dir_input
	print cmd
