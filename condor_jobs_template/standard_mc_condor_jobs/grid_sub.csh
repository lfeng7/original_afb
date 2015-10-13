#! /bin/sh
tar czvf tarball.tgz ../*.C ../*.h input_files/*.txt

~cplager/bin/development/runManySections.py --createCommandFile --cmssw --addLog --setTarball=tarball.tgz \ana.listOfJobs commands.cmd
~cplager/bin/runManySections.py --submitCondor commands.cmd
condor_q lfeng7 
