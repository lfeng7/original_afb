#imports
from ROOT import *
#consts
prepend ='total_tagged_samples/'
# file names to read
filenames = []
samples = []
filelist = []
# samples.append('Powheg_qq_semilep_TT_all')
samples.append('Powheg_qq_semilep_TT_all_unscaled')
samples.append('W4Jets_all_unscaled')
#samples.append('SingleMu_Run2012A_all')
# samples.append('T_s-channel_all')
# samples.append('Powheg_Dilep_had_ttbar')
# samples.append('DYnJets')

color = [2,3,4,6,7,9]
for sp in samples:
	filenames = prepend+sp+'.root'
	filelist.append(TFile(filenames))
# Create histograms for plots
hist_list = [] # A list of list of histograms, each list is the histograms of a var for all samples
vars_list = ['finalChiHist','combttbarmass5']
#vars_list.extend('combttbarmass5')
#vars_list.extend(('hardestjetpt5','hardestjeteta5','hardestjetmass5','hardestjetcsv5'))
#vars_list.extend(('secondhardestjetpt5','secondhardestjeteta5','secondhardestjetmass5','secondhardestjetcsv5'))
# # vars_list.extend(('thirdhardestjetpt5','thirdhardestjeteta5','thirdhardestjetmass5','thirdhardestjetcsv5'))
# vars_list.extend(('fourthhardestjetpt5','fourthhardestjeteta5','fourthhardestjetmass5','fourthhardestjetcsv5'))
# vars_list.extend(('fifthhardestjetpt5','fifthhardestjeteta5','fifthhardestjetmass5','fifthhardestjetcsv5'))
plots_name = []
for ibr in range(len(vars_list)):
	plots_name.append('tagged_var_plots/candidate_jets/'+vars_list[ibr]+'.png')

#lepton, neutrino, leptonic b, leptonic W, hadronic b, and hadronic W, plus the leptonic and hadronic t/tbar
#0		1 			2 			3			4				5 						6 			7

# Loop over sample files
for i in range(len(samples)):
	print 'Sample type is: '+samples[i]
	# Open sample file
	tree = filelist[i].Get('output')
	# Loop over entries
	nev = tree.GetEntries()
	# Normalize to one is a good choice
	# weight = 1/float(nev)
	# Get premade histograms and rescale and color them
	tmp_list = []
	# Loop over interested vars
	for ibr in range(len(vars_list)):	
		tmp_hist = filelist[i].Get(vars_list[ibr])
		weight = 1/float(tmp_hist.Integral())
		tmp_hist.Scale(weight)
		tmp_hist.SetLineColor(color[i])
		tmp_list.append(tmp_hist)
	# Put the histograms from current sample in the list	
	hist_list.append(tmp_list)
	# Report
	print 'number of events: '+str(nev)

#make a legend
leg = TLegend(0.7,0.9,0.9,0.7)
for isample in range(len(samples)):
	tmp = hist_list[isample][0]
	leg.AddEntry(tmp,samples[isample],"L")

# Set a proper y-max and x-axis tile
for ibr in range(len(vars_list)):
	# Find the y max among all samples
	ymax = 0
	for isample in range(len(samples)):
		 tmp_max = hist_list[isample][ibr].GetMaximum()
		 if tmp_max > ymax :
		 	ymax = tmp_max
	# Set the y-axis max a bit larger than the y max
	for isample in range(len(samples)):
		hist_list[isample][ibr].SetMaximum(ymax*1.1)
		hist_list[isample][ibr].GetXaxis().SetTitle(vars_list[ibr])

# Plotting
gStyle.SetOptStat(0)
for ibr in range(len(vars_list)):
	canv = TCanvas(vars_list[ibr],vars_list[ibr],1200,900) # May need to use ""
	canv.cd()
	hist_list[0][ibr].Draw()
	for isample in range(len(samples)):
		hist_list[isample][ibr].Draw('same')
	leg.Draw()
	canv.SaveAs(plots_name[ibr])

# Close files
for isample in range(len(samples)):
	filelist[isample].Close()

		
