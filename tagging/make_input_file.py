import glob
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--dir', metavar='F', type='string', action='store',dest='directory',help='') ## Sets which files to run on
(options, args) = parser.parse_args()

FILES = '/eos/uscms/store/user/eminizer/nTuples/muons/' + options.directory + '/*.root'
print FILES
files = glob.glob(FILES)
for fname in files :
    cmd = 'echo '+fname+' >> input_files_list.txt'
    print cmd
    os.system(cmd)
