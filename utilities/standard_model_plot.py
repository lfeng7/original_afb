#Takes in output files from the tagger and puts together the plots and stacks them, then reoutputs them.
#imports
from ROOT import *

# If want to print fraction, set it to be 1
more_info = 0
powheg = 1
#lists of filenames, fill colors, and pointers
filenames = []
fillcolors = []
eventsGenerated = []
cross_section = []
prepend = 'MC_angles_files/'
fr = []
hist_MC = []
hist_data = []
nev_MC = []
nev_scaled = []
data_file_name = 'angles_data_SingleMu_Run2012_all.root'
hist_title = []
bin_cs = 20
bin_xf = 30
bin_mtt = 40
powheg_ratio = 1.0-(33./204.)
# Set title names
hist_title.append("cos_cs_stack_powheg.eps")
hist_title.append("xf_stack_powheg.eps")
hist_title.append("mtt_stack_powheg.eps")
#QCD samples
filenames.append(prepend+'angles_QCD_Pt-170to300.root')
eventsGenerated.append(5814398)
cross_section.append(34138.15)
filenames.append(prepend+'angles_QCD_Pt-300to470.root')
eventsGenerated.append(5978500)
cross_section.append(1759.549)
filenames.append(prepend+'angles_QCD_Pt-470to600.root')
eventsGenerated.append(3994848)
cross_section.append(113.8791)
filenames.append(prepend+'angles_QCD_Pt-600to800.root')
eventsGenerated.append(3996864)
cross_section.append(26.9921)
filenames.append(prepend+'angles_QCD_Pt-800to1000.root')
eventsGenerated.append(3998563)
cross_section.append(3.550036)
filenames.append(prepend+'angles_QCD_Pt-1000to1400.root')
eventsGenerated.append(1964088)
cross_section.append(0.737844)
filenames.append(prepend+'angles_QCD_Pt-1400to1800.root')
eventsGenerated.append(2000062)
cross_section.append(0.03352235)
filenames.append(prepend+'angles_QCD_Pt-1800.root')
eventsGenerated.append(971974)
cross_section.append(0.001829005)
for i in range(8) :
	fillcolors.append(kYellow)
# Single top samples (X-Sections from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopSigma8TeV)
filenames.append(prepend+'angles_T_s-channel.root')
eventsGenerated.append(259961)
cross_section.append(3.79)
filenames.append(prepend+'angles_T_t-channel.root')
eventsGenerated.append(3758227)
cross_section.append(56.4)
filenames.append(prepend+'angles_T_tW-channel.root')
eventsGenerated.append(497658)
cross_section.append(11.1)
filenames.append(prepend+'angles_TBar_s-channel.root')
eventsGenerated.append(139974)
cross_section.append(1.76)
filenames.append(prepend+'angles_TBar_t-channel.root')
eventsGenerated.append(1935072)
cross_section.append(30.7)
filenames.append(prepend+'angles_TBar_tW-channel.root')
eventsGenerated.append(493460)
cross_section.append(11.1)
for i in range(6) :
	fillcolors.append(kMagenta)
#ZJets sample
filenames.append(prepend+'angles_DYJets_all.root')
eventsGenerated.append(30459503)
cross_section.append(3503.71)
fillcolors.append(kAzure-2)
#WJets sample
filenames.append(prepend+'angles_WJets_all.root')
eventsGenerated.append(57709905)
cross_section.append(37509)
fillcolors.append(kGreen-3)
#hadronic sample
if powheg == 1 :
	filenames.append(prepend+'angles_hadronic_all.root')	# Powheg
	eventsGenerated.append(10537444*powheg_ratio)
else :
	filenames.append(prepend+'angles_hadronic_all.root') # Madgraph
	eventsGenerated.append(10537444)
cross_section.append(112.33)
fillcolors.append(kRed-7)
#fully leptonic sample
if powheg == 1 :
	filenames.append(prepend+'angles_Powheg_dilep_had_TT_all.root') # Powheg
	eventsGenerated.append(12119013*powheg_ratio)
else :
	filenames.append(prepend+'angles_dileptonic_all.root') # Madgraph
	eventsGenerated.append(12119013)
cross_section.append(25.9)
fillcolors.append(kRed-7)
#semileptonic sample
filenames.append(prepend+'angles_mg_semilep_all.root') # Madgraph
# filenames.append(prepend+'angles_mg_semilep_all.root') # Powheg
eventsGenerated.append(25424818)
cross_section.append(107.67)
fillcolors.append(kRed+1)

filelist = []
for filename in filenames :
	filelist.append(TFile(filename))
#lists of histogram stacks and histogram names
stacklist = []
stacklist.append(THStack('cs','cos(#theta^{*})'))
stacklist.append(THStack('xf','x_{F}'))
stacklist.append(THStack('mtt','M_{t#bar{t}'))

# Add data to compare with MC
file_data = TFile(data_file_name)
tree_data = file_data.Get('angles_data')
h_cs_data = TH1F('cs_data','cos#theta^{*}',bin_cs,-1,1)
h_xf_data = TH1F('xf_data','x_{F}',bin_xf,0,0.6)
h_mtt_data = TH1F('mtt_data','M_{t#bar{t}}',bin_mtt,350,1700)

nev_data = tree_data.GetEntries()
for iev in range(nev_data):
	tree_data.GetEntry(iev)
	h_cs_data.Fill(tree_data.cos_theta_cs)
	h_xf_data.Fill(tree_data.Feynman_x)
	h_mtt_data.Fill(tree_data.ttbar_mass)

# # Debug
# print 'Integral of hist_cs_data = ' + str(h_cs_data.Integral())
# print 'Integral of hist_xf_data = ' + str(h_xf_data.Integral())
# print 'Integral of hist_mtt_data = ' + str(h_mtt_data.Integral())

# Set data histogram format
ymax = h_cs_data.GetMaximum()
print '# of Events in data is : '+str(nev_data)
h_cs_data.GetXaxis().SetTitle("cos#theta^{*}")
h_cs_data.SetMaximum(ymax*1.2)
h_cs_data.SetMinimum(0)
hist_data.append(h_cs_data)

ymax = h_xf_data.GetMaximum()
h_xf_data.GetXaxis().SetTitle("x_{F}")
h_xf_data.SetMaximum(ymax*1.2)
h_xf_data.SetMinimum(0)
hist_data.append(h_xf_data)

ymax = h_mtt_data.GetMaximum()
h_mtt_data.GetXaxis().SetTitle("M_{t#bar{t}}")
h_mtt_data.SetMaximum(ymax*1.2)
h_mtt_data.SetMinimum(0)
hist_data.append(h_mtt_data)
	
# Make histogram for each file
h_list_cs = []
h_list_xf = []
h_list_mtt = []
for j in range(len(filenames)) :
	tmp = filelist[j].Get('angles')
	h_cs_temp = TH1F('cs_'+str(j),filenames[j],bin_cs,-1,1)
	h_xf_temp = TH1F('xf_'+str(j),filenames[j],bin_xf,0,0.6)
	h_mtt_temp = TH1F('mtt_'+str(j),filenames[j],bin_mtt,350,1700)
	nev = tmp.GetEntries()
	nev_MC.append(nev)
	for iev in range(nev):
		tmp.GetEntry(iev)
		h_cs_temp.Fill(tmp.cos_theta_cs)
		h_xf_temp.Fill(tmp.Feynman_x)
		h_mtt_temp.Fill(tmp.ttbar_mass)
	h_list_cs.append(h_cs_temp)
	h_list_xf.append(h_xf_temp)
	h_list_mtt.append(h_mtt_temp)
	# ######################## debugging #########################
	# print filenames[j]
	# print 'nev_MC = '+str(nev_MC[j])
	# print 'h_list_cs Integral = '+ str(h_cs_temp.Integral())
	# print 'h_list_xf Integral = '+ str(h_xf_temp.Integral())
	# print 'h_list_mtt Integral = '+str(h_mtt_temp.Integral()) 
	# ######################## debugging #########################	
hist_MC.append(h_list_cs)
hist_MC.append(h_list_xf)
hist_MC.append(h_list_mtt)
print 'Make histograms from MC templates done!'

# Make stacks by scaling each MC samples and get actual fractions
nev_total_scaled = []
for i in range(len(hist_MC)):
	nev_total = 0
	for j in range(len(filenames)) :
		sc = float((1.0*eventsGenerated[len(filenames)-1]/cross_section[len(filenames)-1])*cross_section[j]/eventsGenerated[j])
		hist_MC[i][j].Scale(sc)
		#Get actual event numbers and total number of events
		nev_tmp =  hist_MC[i][j].Integral()
		nev_total += nev_tmp
	nev_total_scaled.append(nev_total)
	# print 'nev_total_scaled = '+str(nev_total_scaled[i])
print 'make a stack histogram done!'

# Get actual fractions of each process and calculate actual scale of each process (to compare to data)
fr_actual = []
for i in range(len(hist_MC)) :
	tmp_list = []	
	for j in range(len(filenames)) :
		tmp = float(hist_MC[i][j].Integral()/nev_total_scaled[i])
		tmp_list.append(tmp)
	fr_actual.append(tmp_list)
print 'Calculate actual fractions and scales done!'

# Scale the stack plot to accommodate data
for i in range(len(hist_MC)) :
	nev_total_actual = 0
	int_data = hist_data[i].Integral()
	sc_data = float(int_data/nev_total_scaled[i])
	for j in range(len(filenames)):
		tmp = hist_MC[i][j]
		tmp.Scale(sc_data)
		tmp.SetFillColor(fillcolors[j])
		tmp.SetLineColor(fillcolors[j])
		tmp.SetMarkerStyle(21)
		stacklist[i].Add(tmp)
		#Get actual event numbers and total number of events
		nev_tmp =  tmp.Integral()
		nev_total_actual += nev_tmp
	if more_info == 1 :
		print 'sc_data_'+str(i)+'= '+str(sc_data)
		print 'number of events of MC after scaling is : '+str(nev_total_actual)
print 'Make a scaled stack histogram done!'

# Print actual fractions of each process
if more_info == 1 :
	for i in range(len(hist_MC)):
		print str(i)
		for j in range(len(filenames)):
			print str(fr_actual[i][j])

# #print out event numbers
# bck_events = 0
# for i in range(len(filenames)) :
# 	# print '# of events in '+filenames[i]+' is '+str(nev_scaled[i])+' actual fraction is : '+str(fr_actual[i])
# 	if i < len(filenames)-1 :
# 		bck_events += nev_scaled[i]
# # print 'total background events: '+str(bck_events)
# # print 'total signal events: '+str(nev_scaled[len(filenames)-1])

#make a legend
leg = TLegend(0.3,0.6,0.6,0.8)
tmp = h_list_cs[0]
leg.AddEntry(tmp,"QCD","F")
tmp = h_list_cs[8]
leg.AddEntry(tmp,"Single Top","F")
tmp = h_list_cs[14]
leg.AddEntry(tmp,"Z+Jets","F")
tmp = h_list_cs[15]
leg.AddEntry(tmp,"W+Jets","F")
tmp = h_list_cs[16]
leg.AddEntry(tmp,"Dileptonic/Hadronic t#bar{t}","F")
tmp = h_list_cs[18]
leg.AddEntry(tmp,"Semileptonic t#bar{t}","F")
# Data
tmp = h_cs_data
leg.AddEntry(tmp,"Data","F")

# Output
#Make a new root file and save the histogram stacks to it

f = TFile( 'stacked_plots.root', 'Recreate' )
f.cd()
for i in range(len(hist_data)) :
	canv = TCanvas("canv"+str(i),"plot_canvas",1200,900)
	canv.cd()
	hist_data[i].Draw('elp')
	stacklist[i].Draw('same')
	hist_data[i].Draw('elp same')
	leg.Draw()
	canv.SaveAs(hist_title[i])
	canv.Write()
	stacklist[i].Write()

f.Write()
f.Close()
for filep in filelist :
	filep.Close()
