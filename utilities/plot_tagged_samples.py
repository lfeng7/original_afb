#imports
from ROOT import *
#consts
prepend ='small_tagged_samples/'
# file names to read
filenames = []
samples = []
filelist = []
# samples.append('Powheg_qq_semilep_TT_all')
samples.append('Powheg_qq_semilep_TT_all_unscaled')
samples.append('W4Jets_all_unscaled')
# samples.append('Single_top')
# samples.append('Powheg_Dilep_had_ttbar')
# samples.append('DYnJets')

color = [2,3,4,6,7,9]
for sp in samples:
	filenames = prepend+sp+'.root'
	filelist.append(TFile(filenames))
# Create histograms for plots
hist_list = [] # A list of list of histograms, each list is the histograms of a var for all samples
vars_list = ['nbTags']
vars_list.extend(('M_W_lep','M_W_had','M_b_lep','M_b_had','M_t_lep','M_t_had')) # Mass, 6
vars_list.extend(('Pt_W_had','Eta_W_had','Phi_W_had','E_W_had')) # W_had stuff,4
#vars_list.extend(('pZv'))#,'scale_Lep','scale_blep','scale_bh','scale_h1','scale_h2')) # p scale , 6 of them
# brach_list = [tree.finalChi,tree.nbTags]
plots_name = []
# var_hist = [] # each histogram will be the histograms of a var from different samples
num_bins =[4,50,50,50,50,50 ,50 ,100,50, 100,100,100,100,100,100,100,100]
bin_min = [1,70,70,2 ,2 ,168,168,0  ,-4,-3.3,0.0,0,0,0,0,0,0]
bin_max = [5,90,90,50,50,180,180,300, 4, 3.3,600,2,2,2,2,2,2]
for i in range(len(vars_list)):
	tmp = []
	for isample in range(len(samples)):
		tmp.append(TH1F(vars_list[i]+'_'+str(isample),samples[isample],num_bins[i],bin_min[i],bin_max[i]))
	hist_list.append(tmp) #len(var_hist) = len(samples)
	plots_name.append('tagged_var_plots/'+vars_list[i]+'.png')

#lepton, neutrino, leptonic b, leptonic W, hadronic b, and hadronic W, plus the leptonic and hadronic t/tbar
#0		1 			2 			3			4				5 						6 			7
number_signals = 0
# Fill the histograms
for i in range(len(samples)):
	print 'Sample type is: '+samples[i]
	# Open sample file
	tree = filelist[i].Get('output')
	# Loop over entries
	nev = tree.GetEntries()
	# Get signal qq number of events for normalization
	if samples[i] == 'Powheg_qq_semilep_TT_all' :
		number_signals = nev
	# Fill histograms of all branches
	# Loop over events
	for iev in range(nev):
		tree.GetEntry(iev)
		hist_list[0][i].Fill(tree.nbTags)	
		# Mass of b W and top				
		hist_list[1][i].Fill(tree.mass[3]) # 5 is W_had
		hist_list[2][i].Fill(tree.mass[5])
		hist_list[3][i].Fill(tree.mass[2])
		hist_list[4][i].Fill(tree.mass[4])
		hist_list[5][i].Fill(tree.mass[6])
		hist_list[6][i].Fill(tree.mass[7])
		# W_had stuff
		hist_list[7][i].Fill(tree.pt[5])
		hist_list[8][i].Fill(tree.eta[5])
		hist_list[9][i].Fill(tree.phi[5])
		hist_list[10][i].Fill(tree.E[5])
		# Scale stuff
#		hist_list[11][i].Fill(tree.bestFitParValues[0])
#		hist_list[12][i].Fill(tree.bestFitParValues[1])
#		hist_list[13][i].Fill(tree.bestFitParValues[2])
#		hist_list[14][i].Fill(tree.bestFitParValues[3])
#		hist_list[15][i].Fill(tree.bestFitParValues[4])
#		hist_list[16][i].Fill(tree.bestFitParValues[5])

	# End loop over events

	# Normalize to one is a good choice
	weight = 1/float(nev)
	# Rescale each histogram and assign colors
	for ibr in range(len(vars_list)):
		sc = 1/hist_list[ibr][i].Integral()
		hist_list[ibr][i].Scale(sc)
		hist_list[ibr][i].SetLineColor(color[i])
	# Report
	print 'number of events: '+str(nev)
	# print 'weights: '+str(weight)
# print 'Signal qqbar->ttbar has events: '+str(number_signals)

#make a legend
leg = TLegend(0.7,0.9,0.9,0.7)
for isample in range(len(samples)):
	tmp = hist_list[0][isample]
	leg.AddEntry(tmp,samples[isample],"L")

# Set a proper y-max and x-axis tile
for ibr in range(len(vars_list)):
	# Find the y max among all samples
	ymax = 0
	for isample in range(len(samples)):
		 tmp_max = hist_list[ibr][isample].GetMaximum()
		 if tmp_max > ymax :
		 	ymax = tmp_max
	# Set the y-axis max a bit larger than the y max
	for isample in range(len(samples)):
		hist_list[ibr][isample].SetMaximum(ymax*1.1)
		hist_list[ibr][isample].GetXaxis().SetTitle(vars_list[ibr])

# Plotting
gStyle.SetOptStat(0)
for ibr in range(len(vars_list)):
	canv = TCanvas(vars_list[ibr],vars_list[ibr],1200,900) # May need to use ""
	canv.cd()
	hist_list[ibr][0].Draw()
	for isample in range(len(samples)):
		hist_list[ibr][isample].Draw('same')
	leg.Draw()
	canv.SaveAs(plots_name[ibr])

# Close files
for isample in range(len(samples)):
	filelist[isample].Close()

		



