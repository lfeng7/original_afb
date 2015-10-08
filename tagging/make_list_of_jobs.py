from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option('--run_name', metavar='F', type='string', action='store',dest='run_name',help='')
parser.add_option('--decayType', metavar='F', type='string', action='store',dest='decayType',default='none',help='')
parser.add_option('--mother_check', metavar='F', type='string', action='store',dest='mother_check',default='none',help='')
parser.add_option('--sample_type', metavar='F', type='string', action='store',dest='sample_type',default='none',help='')
parser.add_option('--PdfVersion', metavar='F', type='string', action='store',dest='PdfVersion',default='none',help='')
parser.add_option('--is_mc', metavar='F', type='string', action='store',dest='is_mc',help='')
parser.add_option('--is_sb', metavar='F', type='string', action='store',dest='is_sb',help='')
parser.add_option('--nJobs', metavar='F', type='int', action='store',dest='nJobs',help='')

(options, args) = parser.parse_args()

for i in range(options.nJobs) :
	cmd = 'echo "python ./tardir/ttbar_tagger.py --out '+options.run_name+'_'+str(i)+' --nTotalJobs '
	cmd = cmd+str(options.nJobs)+' --iJob '+str(i)+' --mc '+options.is_mc+' --sideband '+options.is_sb+' '
	if not options.decayType == 'none' :
		cmd = cmd + '--decayType '+options.decayType+' '
	if not options.mother_check == 'none' :
		cmd = cmd + '--tMotherCheck '+options.mother_check+' '
	if not options.sample_type == 'none' :
		cmd = cmd + '--sample_type '+options.sample_type+' '
	if not options.PdfVersion == 'none' :
		cmd = cmd + '--PdfVersion '+options.PdfVersion+' '
	cmd = cmd + '" >> ana.listOfJobs'
	print cmd
	os.system(cmd)

print 'done. Completed file: '
os.system('cat ana.listOfJobs')
