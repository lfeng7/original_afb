#! /bin/sh
tar czvf tarball.tgz ../ttbar_tagger.py ../dumped_Powheg_TT.root ../MuonEfficiencies_Run2012ReReco_53X.pkl ../MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl ../SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.pkl ./input_files_list.txt
voms-proxy-init --voms cms

/uscms_data/d3/eminizer/runManySections/runManySections.py --createCommandFile --cmssw --addLog --setTarball=tarball.tgz \ana.listOfJobs commands.cmd
/uscms_data/d3/eminizer/runManySections/runManySections.py --submitCondor commands.cmd
condor_q eminizer
