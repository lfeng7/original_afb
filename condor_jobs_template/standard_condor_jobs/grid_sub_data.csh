#! /bin/sh
tar czvf tarball_data.tgz ../angles_data.C ../angles_data.h ../angles.C ../angles.h ../main_data.C ttbar.C ../ttbar.h all_MC.txt input_files/*.txt 

~cplager/bin/development/runManySections.py --createCommandFile --cmssw --addLog --setTarball=tarball_data.tgz \ana_data.listOfJobs commands_data.cmd
~cplager/bin/runManySections.py --submitCondor commands_data.cmd
condor_q lfeng7 
