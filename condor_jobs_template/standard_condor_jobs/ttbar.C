#define ttbar_cxx
#include "ttbar.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

// 5j + 1b

//Lorentz Declartions

const double sqrt_s=8000;
const double beam_energy=sqrt_s/2;
const int ppbar=1;
const double m_top_mc=172.5;
const double m_top_data=173.3;
const int num_events = 10000; 
TLorentzRotation S;

//variable for ouput file name
const char* output_filename = "angles.root";
//variable for output tree name (should be default except for total data file)
const char* output_treename = "angles";

void ttbar::Loop(int data_or_mc)
{
	//Histograms for tests

	int kept=0;
	int count=0;
    //declarations of output variables
    double  ttbar_mass,cos_theta_cs,x_f,w_a,w_a_opposite_cos_theta,Qt,cos_theta_mc;
    int Q_l;
	//open output file and create branches
    TFile * file=new TFile(output_filename,"Recreate");
    TTree * output=new TTree(output_treename,"Recreate");
    output->Branch("ttbar_mass",&ttbar_mass);
    output->Branch("Qt",&Qt);
    output->Branch("cos_theta_cs",&cos_theta_cs);
    output->Branch("Feynman_x",&x_f);
    output->Branch("w_a",&w_a);
    output->Branch("w_a_opposite_cos_theta",&w_a_opposite_cos_theta);
    output->Branch("Q_l",&Q_l);
    output->Branch("cos_theta_mc",&cos_theta_mc);
	
    //Setup loop
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	int forward = 0;
	int backward = 0;
	// Control over totoal number of events in the loop
	if ( nentries > num_events ){ nentries = num_events; }
	//MAIN LOOP
	for (Long64_t jentry=0; jentry< nentries ;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		count+=1;
		//intialize the rotations to the identity
		TLorentzRotation R_data;
		TLorentzRotation R_mc;
		
		//Supercool Load bar
		if (nentries > 100) {
			if ( count % (nentries/100) == 0 ){
				float ratio = count/(float)nentries;
				int   c     = ratio * 50;
				printf("%3d%% [", (int)(ratio*100) );
				for (int x=0; x<c; x++) {
					if (x+1>=c)
						printf(">");
					else
						printf("=");
				}
				for (int x=c; x<50; x++)
					printf(" ");
				printf("]\n\033[F\033[J");
			}//end supercool loadbar
		}		
		//Vectors of top and antitop
		TLorentzVector Top_MC;
		TLorentzVector ATop_MC;
		TLorentzVector Top_Data;
		TLorentzVector ATop_Data;
		TLorentzVector quark;
		TLorentzVector antiquark;
		
		//loop over particles in event and get the top and antitop fourvectors
        for (int j=6;j<8;j++) {
			if (PID[j] == 6) {
				Top_Data = *(new TLorentzVector(px[j],py[j],pz[j],E[j]));
				if (data_or_mc==1) 
					Top_MC = *(new TLorentzVector(mc_px[j],mc_py[j],mc_pz[j],mc_E[j]));
				else 
					Top_MC = *(new TLorentzVector(px[j],py[j],pz[j],E[j]));
			}
			if (PID[j] == -6) {
				ATop_Data = *(new TLorentzVector(px[j],py[j],pz[j],E[j]));
				if (data_or_mc==1) 
					ATop_MC = *(new TLorentzVector(mc_px[j],mc_py[j],mc_pz[j],mc_E[j]));
				else 
					ATop_MC = *(new TLorentzVector(px[j],py[j],pz[j],E[j]));
			}
		}//end of loop over particles
		//assign the quark and antiquark fourvectors from the last two particles in the list of mc_ particles
		if (data_or_mc==1) {
			quark = *(new TLorentzVector(mc_px[8],mc_py[8],mc_pz[8],mc_E[8]));
			antiquark = *(new TLorentzVector(mc_px[9],mc_py[9],mc_pz[9],mc_E[9]));
		}
		else {
			quark = *(new TLorentzVector(0,0,sqrt(beam_energy*beam_energy -1*1),beam_energy));
			antiquark = *(new TLorentzVector(0,0,-1*quark.Pz(),beam_energy));
		}
		//set the lepton charge
		if (PID[0] > 0)
			Q_l = -1;
		else if (PID[0] < 0)
			Q_l = 1;
		else
			printf("Something went wrong in getting the lepton charge, not sure what.\n");

		int weight_is_valid = 1;

		//If the event was one of the background ones that didn't actually have a t and/or tbar associated with it
		//Then set the fourvector of the MC top and antitop to be the same as the data ones, just make sure we don't
		//Get a weight for this event (won't matter anyway since there's no asymmetric distributions built from these)
		if (Top_MC.Mag() == 0.0 || ATop_MC.Mag() == 0.0) {
			Top_MC = Top_Data;
			ATop_MC = ATop_Data;
			weight_is_valid = 0;
		}
		//Likewise if there wasn't actually a quark or antiquark associated with the event
		if (quark.Mag() == 0.0 || antiquark.Mag() == 0.0) {
			quark = *(new TLorentzVector(0,0,sqrt(beam_energy*beam_energy -1*1),beam_energy));
			antiquark = *(new TLorentzVector(0,0,-1*quark.Pz(),beam_energy));
			weight_is_valid = 0;
		}
		
		//Make the 4-vector of the ttbar
		TLorentzVector Q_Data = Top_Data + ATop_Data;
		double ttbar_mass_data=Q_Data.Mag();
		TLorentzVector Q_MC = Top_MC + ATop_MC;
		double ttbar_mass_mc=Q_MC.Mag();

		//defining the Px, Py,and Pz, and energies to boost into the ttbar rest frame
		double Bx_data = -1*Q_Data.Px()/Q_Data.E();  
		double By_data = -1*Q_Data.Py()/Q_Data.E();  
		double Bz_data = -1*Q_Data.Pz()/Q_Data.E();
		Qt = sqrt(Q_Data.Px()*Q_Data.Px()+Q_Data.Py()*Q_Data.Py());
		double Bx_mc = -1*Q_MC.Px()/Q_MC.E();  
		double By_mc = -1*Q_MC.Py()/Q_MC.E();  
		double Bz_mc = -1*Q_MC.Pz()/Q_MC.E();
		
		double beta_mc, beta_data;
		if (data_or_mc == 1) {
			beta_mc = sqrt(1-4*(m_top_mc/ttbar_mass_mc)*(m_top_mc/ttbar_mass_mc));
			beta_data = sqrt(1-4*(m_top_mc/ttbar_mass_data)*(m_top_mc/ttbar_mass_data));
		}
		else {
			beta_mc = sqrt(1-4*(m_top_data/ttbar_mass_mc)*(m_top_data/ttbar_mass_mc));
			beta_data = sqrt(1-4*(m_top_data/ttbar_mass_data)*(m_top_data/ttbar_mass_data));
		}
		
		//Need beta to be real 
		if (TMath::IsNaN(beta_data) || TMath::IsNaN(beta_mc))// || nValidJets[0] != 5 || nbTags != 1)
			continue;

		//Feynman x (x_f) (only needed for reconstructed events, not MC truth)
		x_f = 2*Q_Data.Pz()/sqrt_s;
		
		//Creating the Lorentz 4 vectors of the protons
		TLorentzVector Proton_data = *(new TLorentzVector(0,0,sqrt(beam_energy*beam_energy -1*1),beam_energy));
		TLorentzVector Proton2_data = *(new TLorentzVector(0,0,-1*Proton_data.Pz(),beam_energy));
		
		//Doing the boost
		R_data = R_data.Boost(Bx_data,By_data,Bz_data);
		Top_Data = R_data*Top_Data;
		ATop_Data = R_data*ATop_Data;
		Proton_data = R_data*Proton_data;
		Proton2_data = R_data*Proton2_data;
		R_mc = R_mc.Boost(Bx_mc,By_mc,Bz_mc);
		Top_MC = R_mc*Top_MC;
		ATop_MC = R_mc*ATop_MC;
		quark = R_mc*quark;
		antiquark = R_mc*antiquark;
		//Reset the boost
		R_data=S;
		R_mc=S;
		//Define three vectors for P,Pbar,top, quark and antiquark in ttbar c.m frame
		TVector3 top_data = *(new TVector3(Top_Data.Px(),Top_Data.Py(),Top_Data.Pz()));
		TVector3 proton_data = *(new TVector3(Proton_data.Px(),Proton_data.Py(),Proton_data.Pz()));
		TVector3 proton2_data = *(new TVector3(Proton2_data.Px(),Proton2_data.Py(),Proton2_data.Pz()));
		TVector3 top_mc = *(new TVector3(Top_MC.Px(),Top_MC.Py(),Top_MC.Pz()));
		TVector3 true_quark_direction = *(new TVector3(quark.Px(),quark.Py(),quark.Pz()));
		TVector3 true_antiquark_direction = *(new TVector3(antiquark.Px(),antiquark.Py(),antiquark.Pz()));
		
		//Flip the larger one between proton and proton2, and flip the antiquark direction
		if(proton_data.Mag()>proton2_data.Mag()){proton_data=-1.0*proton_data;}else{proton2_data=-1.0*proton2_data;}
		true_antiquark_direction = -1.0*true_antiquark_direction;
		
		//Normalize vectors
		top_data = top_data*(1.0/top_data.Mag());
		proton_data = proton_data*(1.0/proton_data.Mag());
		proton2_data = proton2_data*(1.0/proton2_data.Mag());
		top_mc = top_mc*(1.0/top_mc.Mag());
		true_quark_direction = true_quark_direction*(1.0/true_quark_direction.Mag());
		true_antiquark_direction = true_antiquark_direction*(1.0/true_antiquark_direction.Mag());
		//find the unit bisectors
		TVector3 bisector_data = (proton_data+proton2_data)*(1.0/(proton_data+proton2_data).Mag());
		TVector3 bisector_mc = (true_quark_direction+true_antiquark_direction)*(1.0/(true_quark_direction+true_antiquark_direction).Mag());
		//find the CS angle
		double cos_theta_cs_data=top_data*bisector_data;
		double cos_theta_cs_mc=top_mc*bisector_mc;

		//fill histograms
		cos_theta_mc = cos_theta_cs_mc;

		//count the forwards and backwards to get the asymmetry from MC truth
		if (data_or_mc==1) {
			if (cos_theta_cs_mc > 0.0)
				++forward;
			else if (cos_theta_cs_mc < 0.0)
				++backward;
		}
		
		//weights (cf equation 5 in note)
		if (weight_is_valid != 1) {
			w_a = 0.0;
			w_a_opposite_cos_theta = 0.0;
		}
		else {
			w_a=(2.0*((1.0+(1.0/3)*beta_mc*beta_mc+(1-beta_mc*beta_mc)))*cos_theta_cs_mc)
						/(1.0+beta_mc*beta_mc*cos_theta_cs_mc*cos_theta_cs_mc+(1.0-beta_mc*beta_mc));
			w_a_opposite_cos_theta=(2.0*((1.0+(1.0/3)*beta_mc*beta_mc+(1-beta_mc*beta_mc)))*(-1.0*cos_theta_cs_mc))
						/(1.0+beta_mc*beta_mc*(-1.0*cos_theta_cs_mc)*(-1.0*cos_theta_cs_mc)+(1.0-beta_mc*beta_mc));
		}
		//Set the other output variables to the correct values
		ttbar_mass=ttbar_mass_data;
		cos_theta_cs=cos_theta_cs_data;
		//Write all data into leaves
		output->Fill();
		kept +=1;

	} //End of loop over entries

	cout << "		Done." << endl;
	
    //write the output file
    output->Write();
	
	//close the file
	file->Close();
	if (count != 0) 
		printf("# OF EVENTS CUT FROM SAMPLE = %d (%f) %%\n",count-kept,100.*(count-kept)/count);
	if (data_or_mc == 1 && (forward+backward)!=0) {
		printf("Asymmetry from MC truth: %f = (%d - %d)/%d \n", 1.0*(forward-backward)/(forward+backward), forward, backward, (forward+backward));
	}
	
}

void ttbar::SetOutputFilename(const char* outname) {
	output_filename = outname;
}

void ttbar::SetOutputTreeName(const char* treename) {
	output_treename = treename;
}
