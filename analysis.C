#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void analysis::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L analysis.C
//      Root > analysis t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
double xs = 15.618;
TH1D* top_hadronic = new TH1D("top_from_jets","Top Quark from Jets Mass",60,0,300);
TH1D* w_from_jets = new TH1D("w_from_jets","W Boson from Jets Mass",60,0,200); 
TH1D* higgs = new TH1D("higgs","Higgs lep/had invariant mass",100,0,250);
TH1D* electrons = new TH1D("electrons","Leading electrons",60,0,200);
TH1D* leading_electron = new TH1D("leading_electron","Leading electron",750,0,200);
TH1D* leading_muon = new TH1D("leading_muon","leading muon",60,0,200);
TH1D* missingET = new TH1D("misingET","Missing Energy",60,0,200);
TH1D* M_b12 = new TH1D("M_b12","Invariant Mass of b+1+2 Jets",60,0,300); 
TH1D* M_b13 = new TH1D("M_b13","Invariant Mass of b+1+3 Jets",60,0,300); 
TH1D* M_b23 = new TH1D("M_23","Invariant Mass of b+2+3 Jets",60,0,300); 
TH1D* M_eta = new TH1D("M_eta","Invariant Mass of smallest eta",60,0,300);
TH1D* leading_lepton = new TH1D("leading_lepton","Leading Lepton PT",60,0,300);
TH1D* tau_hadron_PT = new TH1D("tau_hadron_PT","Hadron from tau Pt",60,0,300);
TH1D* higgs_leptonic = new TH1D("higgs_leptonic","Higgs invariant mass from e+mu opposite sign",60,0,300);
TH1D* higgs_hadronic = new TH1D("higgs_hadronic","Higgs mass from 2 hadronic taus",60,0,300);
TH1D* higgs_electron = new TH1D("higgs_electron","Higgs mass form 2 electrons",60,0,300);
TH1D* higgs_muon = new TH1D("higgs_muon","Higgs from 2 muons",60,0,300);
TH1D* Jet1_Pt = new TH1D("Jet1_Pt","Jet1 Pt",60,0,300);
TH1D* Jet2_Pt = new TH1D("Jet2_Pt","Jet2 Pt",60,0,300);
TH1D* Jet3_Pt = new TH1D("Jet3_Pt","Jet3 Pt",60,0,300);
TH1D* Jet4_Pt = new TH1D("Jet4_Pt","Jet4 Pt",60,0,300);
TH1D* M_12_postSelection = new TH1D("M_12_postSelection","Invariant Mass of j1+j2",60,0,300);
TH1D* M_b12_postSelection = new TH1D("M_b12_postSelection","Invariant Mass of b+j1+j2",60,0,300);
TH1D* M_3emu_postSelection = new TH1D("M_3h_postSelection","Invariant Mass of j3+tau+tau",60,0,300);
TH1D* M_12_col = new TH1D("M_12_col","Invariant mass of jets near b jet",60,0,300);
TH1D* bJet_Pt = new TH1D("bJet_Pt","bJet Pt",60,0,300);

int acceptance_hadronic=0;
int acceptance_lep=0;
int acceptance=0;
int acceptance_electron=0;
int acceptance_muon=0;
int events=0;
int count=0;
double M_b12_accuracy=0;
double M_b13_accuracy=0;
double M_b23_accuracy=0;
double M_eta_accuracy=0;
static double top_mass = 173.2;
static double w_mass=80.0;
 Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
TLorentzVector v_bJet;
TLorentzVector v_j1;
TLorentzVector v_j2;
TLorentzVector v_j3; 
TLorentzVector v_top_from_jets;
TLorentzVector v_w;
TLorentzVector v_tau1;
TLorentzVector v_tau2;
TLorentzVector v_higgs;
TLorentzVector v_lep;
TLorentzVector v_electronSum;
TLorentzVector v_M_b12;
TLorentzVector v_M_b13;
TLorentzVector v_M_b23;
TLorentzVector v_M_eta;
TLorentzVector v_missing;
TLorentzVector v_muon;
TLorentzVector v_electron;
TLorentzVector v_m1;
TLorentzVector v_m2;
TLorentzVector v_e1;
TLorentzVector v_e2;
TLorentzVector v_M_12;
TLorentzVector v_M_13;
TLorentzVector v_M_23;
TLorentzVector v_M_12_postSelection;
TLorentzVector v_M_b12_postSelection;
TLorentzVector v_M_3h_postSelection;

int numbJets=0;
int numLightJets=0;
int numTauJets=0;
int trashEvent=0;
//cout << "Event Number: " << jentry << endl;
//if(Jet_size>=3 && Jet_PT[0]>Jet_PT[1] && Jet_PT[1]>Jet_PT[2]){
//cout << Jet_PT[0] << endl;
//cout << Jet_PT[1] << endl;
//cout << Jet_PT[2] << endl;
//}
	for(int i=0;i<Jet_size;i++)
		{	
		
			if( numbJets==0 && Jet_BTag[i]==1){
				v_bJet.SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i],Jet_Phi[i],Jet_Mass[i]);	
				numbJets++;
				bJet_Pt->Fill(Jet_PT[i]);
				
				//cout << jentry << endl;
					}
			else if (numbJets==1 && Jet_BTag[i]==1){
				trashEvent=1;
					}
			else if(numLightJets==0 && Jet_BTag[i]==0 && Jet_TauTag[i]==0){	
				v_j1.SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i],Jet_Phi[i],Jet_Mass[i]);
				numLightJets++;
				
					}
			else if(numLightJets==1 && Jet_BTag[i]==0 && Jet_TauTag[i]==0){
				v_j2.SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i],Jet_Phi[i],Jet_Mass[i]);
				numLightJets++;
				
					}
			else if(numLightJets==2 && Jet_BTag[i]==0 && Jet_TauTag[i]==0){
                                v_j3.SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i],Jet_Phi[i],Jet_Mass[i]);
                                numLightJets++;
				
                                        }
		
			if(numTauJets==0 && Jet_TauTag[i]==1){
					v_tau1.SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i],Jet_Phi[i],Jet_Mass[i]);
				numTauJets++;
					}
			 else if(numTauJets==1 && Jet_TauTag[i]==1){
                       		v_tau2.SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i],Jet_Phi[i],Jet_Mass[i]);
                        	numTauJets++;
                                        }		
		}//end of Jet Loop
		//if(trashEvent==1){continue;}	
	
		//cout << "Jet PT's " << v_j1.Pt()<<" "<<v_j2.Pt()<<" "<<v_j3.Pt() <<endl;
		/*
		double  delta_eta_b12 = fabs(v_bJet.Eta()-v_j1.Eta())+fabs(v_bJet.Eta()-v_j2.Eta());
		double  delta_eta_b13 = fabs(v_bJet.Eta()-v_j1.Eta())+fabs(v_bJet.Eta()-v_j3.Eta());
		double  delta_eta_b23 = fabs(v_bJet.Eta()-v_j2.Eta())+fabs(v_bJet.Eta()-v_j3.Eta());
		
		double dRb1 = sqrt((v_bJet.Phi()-v_j1.Phi())*(v_bJet.Phi()-v_j1.Phi())+(v_bJet.Eta()-v_j1.Eta())*(v_bJet.Eta()-v_j1.Eta()));
		double dRb2 = sqrt((v_bJet.Phi()-v_j2.Phi())*(v_bJet.Phi()-v_j2.Phi())+(v_bJet.Eta()-v_j2.Eta())*(v_bJet.Eta()-v_j2.Eta()));
		double dRb3 = sqrt((v_bJet.Phi()-v_j3.Phi())*(v_bJet.Phi()-v_j3.Phi())+(v_bJet.Eta()-v_j3.Eta())*(v_bJet.Eta()-v_j3.Eta()));
		cout<<"debug!"<<endl;
		int etaFill=0;
		if(dRb1 < dRb3 && dRb2 < dRb3){
			v_M_eta=v_bJet+v_j1+v_j2;
			M_12_col->Fill((v_j2+v_j2).M());
			//cout << etaFill<<endl;
			//etaFill++;
		}
		  if(dRb1 < dRb2 && dRb3 < dRb2){
                        v_M_eta=v_bJet+v_j1+v_j3;
			M_12_col->Fill((v_j1+v_j3).M());
			//cout << etaFill<<endl;
			//etaFill++;
                }
		  if(dRb2 < dRb1 && dRb3 < dRb1){
                        v_M_eta=v_bJet+v_j2+v_j3;
			M_12_col->Fill((v_j2+v_j3).M());
			//cout << etaFill<<endl;
			//etaFill++;
                }

*/
		v_M_b12 = v_bJet+v_j1+v_j2;
		v_M_b13 = v_bJet+v_j1+v_j3;
		v_M_b23 = v_bJet+v_j2+v_j3;
		v_M_12 = v_j1 + v_j2;
		v_M_13 = v_j1 + v_j3;
		v_M_23 = v_j2 + v_j3;


//cout << " Masses " << v_M_b12.M() << " " << v_M_b13.M() <<" "<< v_M_b23.M() << endl;
//cout << numTauJets << endl;

//if(numTauJets==2)
//{
//v_higgs = v_tau1 + v_tau2;
//cout << v_higgs.M();
//higgs->Fill(v_higgs.M());
//}

if (Muon_size>1){
v_m1.SetPtEtaPhiM(Muon_PT[0],Muon_Eta[0],Muon_Phi[0],0);
v_m2.SetPtEtaPhiM(Muon_PT[1],Muon_Eta[1],Muon_Phi[1],0);
}
if(Electron_size>1){
v_e1.SetPtEtaPhiM(Electron_PT[0],Electron_Eta[0],Electron_Phi[0],0);
v_e2.SetPtEtaPhiM(Electron_PT[1],Electron_Eta[1],Electron_Phi[1],0);
}



if(Muon_size>0 && Electron_size>0 && Muon_Charge[0]!=Electron_Charge[0]){
v_electron.SetPtEtaPhiM(Electron_PT[0],Electron_Eta[0],Electron_Phi[0],0);
v_muon.SetPtEtaPhiM(Muon_PT[0],Muon_Eta[0],Muon_Phi[0],0);
}


if(Electron_size>0 && Muon_size==0){
leading_lepton->Fill(Electron_PT[0]);
v_lep.SetPtEtaPhiM(Electron_PT[0],Electron_Eta[0],Electron_Phi[0],0);
}
else if(Muon_size>0 && Electron_size==0){
leading_lepton->Fill(Muon_PT[0]);
v_lep.SetPtEtaPhiM(Muon_PT[0],Muon_Eta[0],Muon_Phi[0],0);

}
else if(Muon_PT[0]>Electron_PT[0]){
leading_lepton->Fill(Muon_PT[0]);
v_lep.SetPtEtaPhiM(Muon_PT[0],Muon_Eta[0],Muon_Phi[0],0);

}
else if(Electron_PT[0]>Muon_PT[0]){
leading_lepton->Fill(Electron_PT[0]);

v_lep.SetPtEtaPhiM(Electron_PT[0],Electron_Eta[0],Electron_Phi[0],0);
}
else
cout << "NO LL" <<endl;


M_b12->Fill(v_M_b12.M());
M_b13->Fill(v_M_b13.M());
M_b23->Fill(v_M_b23.M());
M_eta->Fill(v_M_eta.M());


electrons->Fill(v_electronSum.M());
missingET->Fill(MissingET_MET[0]);
if(v_tau1.Pt()>2.0)
tau_hadron_PT->Fill(v_tau1.Pt());


//v_top_from_jets = v_bJet+v_j1+v_j2;
//v_w = v_j1+v_j2;
//if(numbJets==1 && numLightJets==2)
//top_hadronic->Fill(v_top_from_jets.M());

//if(numLightJets==2)
//w_from_jets->Fill(v_w.M());
 
Jet1_Pt->Fill(Jet_PT[0]);
Jet2_Pt->Fill(Jet_PT[1]);
Jet3_Pt->Fill(Jet_PT[2]);
Jet4_Pt->Fill(Jet_PT[3]);

if(Electron_size>0)
{
leading_electron->Fill(Electron_PT[0]);
//cout << count <<endl;
//count++;
}
if(Muon_size>0)
{
leading_muon->Fill(Muon_PT[0]);
//cout << count <<endl;
//count++;
}


double DM12 = fabs(v_M_b12.M()-top_mass)+fabs(v_M_12.M()-w_mass);
double DM13 = fabs(v_M_b13.M()-top_mass)+fabs(v_M_13.M()-w_mass);
double DM23 = fabs(v_M_b23.M()-top_mass)+fabs(v_M_23.M()-w_mass);

M_b12_accuracy = M_b12_accuracy + DM12;
M_b13_accuracy = M_b13_accuracy + DM13;
M_b23_accuracy = M_b23_accuracy + DM23;
M_eta_accuracy = M_eta_accuracy + fabs(v_M_eta.M()-top_mass);
TLorentzVector v_swap;

if(DM13 < DM12)
{
	v_swap = v_j3;
	v_j3 = v_j2;
	v_j2 = v_swap;
}
if(DM23 < DM12)
{
	v_swap = v_j1;
	v_j1 = v_j3;
	v_j3 = v_swap;
}

M_12_postSelection->Fill((v_j1+v_j2).M());
M_b12_postSelection->Fill((v_bJet+v_j1+v_j2).M());

//This will be the computational part
//First is the lep/had decay mode
if(v_tau1.Pt()>2 && MissingET_MET[0]>2 && v_lep.Pt()>2){
v_missing.SetPtEtaPhiM(MissingET_MET[0],0,MissingET_Phi[0],0);
double theta_l=TMath::ASin(v_lep.Py()/v_lep.Pt());
double theta_h=TMath::ASin(v_tau1.Py()/v_tau1.Pt());
double theta_m=TMath::ASin(v_missing.Py()/v_missing.Pt());

double xl=(v_lep.Pt()*sin(theta_h-theta_l))/(v_lep.Pt()*sin(theta_h-theta_l)+v_missing.Pt()*sin(theta_h-theta_m));
double xh=(v_tau1.Pt()*sin(theta_h-theta_l))/(v_tau1.Pt()*sin(theta_h-theta_l)-v_missing.Pt()*sin(theta_l-theta_m));
//cout << "xl = " << xl << " , xh = " << xh <<endl;

TLorentzVector v_lep_scaled;
TLorentzVector v_had_scaled;
TLorentzVector v_higgs;

v_lep_scaled = v_lep*(1/xl);
v_had_scaled = v_tau1*(1/xh);

v_higgs = v_lep_scaled+v_had_scaled;

higgs->Fill(v_higgs.M());

if(v_higgs.M()>110 && v_higgs.M()<140)
acceptance++;
}
//Here begins the fully hadronic tau reconstruction
if(v_tau1.Pt()>2 && MissingET_MET[0]>2 && v_tau2.Pt()>2){

v_missing.SetPtEtaPhiM(MissingET_MET[0],0,MissingET_Phi[0],0);
double theta_l=TMath::ASin(v_tau2.Py()/v_tau2.Pt());
double theta_h=TMath::ASin(v_tau1.Py()/v_tau1.Pt());
double theta_m=TMath::ASin(v_missing.Py()/v_missing.Pt());



double xl=(v_tau2.Pt()*sin(theta_h-theta_l))/(v_tau2.Pt()*sin(theta_h-theta_l)+v_missing.Pt()*sin(theta_h-theta_m));
double xh=(v_tau1.Pt()*sin(theta_h-theta_l))/(v_tau1.Pt()*sin(theta_h-theta_l)-v_missing.Pt()*sin(theta_l-theta_m));
//cout << "xl = " << xl << " , xh = " << xh <<endl;


TLorentzVector v_had2_scaled;
TLorentzVector v_had_scaled;
TLorentzVector v_higgs;

v_had2_scaled = v_tau2*(1/xl);
v_had_scaled = v_tau1*(1/xh);

v_higgs = v_had2_scaled+v_had_scaled;

higgs_hadronic->Fill(v_higgs.M());

if(v_higgs.M()>110 && v_higgs.M()<140)
acceptance_hadronic++;
//cout << "exit if"<<endl;
}
//This begins the 2 electron channel
if(v_e1.Pt()>2 && MissingET_MET[0]>2 && v_e2.Pt()>2){
v_missing.SetPtEtaPhiM(MissingET_MET[0],0,MissingET_Phi[0],0);
double theta_l=TMath::ASin(v_e1.Py()/v_e1.Pt());
double theta_h=TMath::ASin(v_e2.Py()/v_e2.Pt());
double theta_m=TMath::ASin(v_missing.Py()/v_missing.Pt());

double xl=(v_e1.Pt()*sin(theta_h-theta_l))/(v_e1.Pt()*sin(theta_h-theta_l)+v_missing.Pt()*sin(theta_h-theta_m));
double xh=(v_e2.Pt()*sin(theta_h-theta_l))/(v_e2.Pt()*sin(theta_h-theta_l)-v_missing.Pt()*sin(theta_l-theta_m));
//cout << "xl = " << xl << " , xh = " << xh <<endl;

TLorentzVector v_lep_scaled;
TLorentzVector v_had_scaled;
TLorentzVector v_higgs;

v_lep_scaled = v_e1*(1/xl);
v_had_scaled = v_e2*(1/xh);

v_higgs = v_lep_scaled+v_had_scaled;

higgs_electron->Fill(v_higgs.M());
    if(v_higgs.M()>110 && v_higgs.M()<140)
    acceptance_electron++;

}
//This begins the 2 muon channel
if(v_m1.Pt()>2 && MissingET_MET[0]>2 && v_m2.Pt()>2){
v_missing.SetPtEtaPhiM(MissingET_MET[0],0,MissingET_Phi[0],0);
double theta_l=TMath::ASin(v_m1.Py()/v_m1.Pt());
double theta_h=TMath::ASin(v_m2.Py()/v_m2.Pt());
double theta_m=TMath::ASin(v_missing.Py()/v_missing.Pt());

double xl=(v_m1.Pt()*sin(theta_h-theta_l))/(v_m1.Pt()*sin(theta_h-theta_l)+v_missing.Pt()*sin(theta_h-theta_m));
double xh=(v_m2.Pt()*sin(theta_h-theta_l))/(v_m2.Pt()*sin(theta_h-theta_l)-v_missing.Pt()*sin(theta_l-theta_m));
//cout << "xl = " << xl << " , xh = " << xh <<endl;

TLorentzVector v_lep_scaled;
TLorentzVector v_had_scaled;
TLorentzVector v_higgs;

v_lep_scaled = v_m1*(1/xl);
v_had_scaled = v_m2*(1/xh);

v_higgs = v_lep_scaled+v_had_scaled;

higgs_muon->Fill(v_higgs.M());
    if(v_higgs.M()>110 && v_higgs.M()<140)
    acceptance_muon++;

}
//------------------------------------------------------------------------
//The below is for e+mu opposite sign
if(v_electron.Pt()>2 && MissingET_MET[0]>2 && v_muon.Pt()>2 && Muon_size==1 && Electron_size==1){
v_missing.SetPtEtaPhiM(MissingET_MET[0],0,MissingET_Phi[0],0);
double theta_e=TMath::ASin(v_electron.Py()/v_electron.Pt());
double theta_mu=TMath::ASin(v_muon.Py()/v_muon.Pt());
double theta_m=TMath::ASin(v_missing.Py()/v_missing.Pt());

double xe=(v_electron.Pt()*sin(theta_mu-theta_e))/(v_electron.Pt()*sin(theta_mu-theta_e)+v_missing.Pt()*sin(theta_mu-theta_m));
double xmu=(v_muon.Pt()*sin(theta_mu-theta_e))/(v_muon.Pt()*sin(theta_mu-theta_e)-v_missing.Pt()*sin(theta_e-theta_m));
//cout << "xl = " << xe << " , xh = " << xmu <<endl;

TLorentzVector v_electron_scaled;
TLorentzVector v_muon_scaled;
TLorentzVector v_higgs;

v_electron_scaled = v_electron*(1/xe);
v_muon_scaled = v_muon*(1/xmu);

v_higgs = v_electron_scaled+v_muon_scaled;

higgs_leptonic->Fill(v_higgs.M());

M_3emu_postSelection->Fill((v_higgs+v_j3).M());




if(v_higgs.M()>105 && v_higgs.M()<145)
acceptance_lep++;
}
//Ends the Compuations
//
      events++;
//if (jentry>10000){break;}
if(jentry % 100 == 0 )
cout << jentry << " events completed" <<endl;




  }//Event Loop
Jet1_Pt->Sumw2();
Jet2_Pt->Sumw2();
Jet3_Pt->Sumw2();
Jet4_Pt->Sumw2();
bJet_Pt->Sumw2();
leading_electron->Sumw2();
leading_muon->Sumw2();
missingET->Sumw2();
M_12_postSelection->Sumw2();
M_b12_postSelection->Sumw2();
M_3emu_postSelection->Sumw2();
M_b12->Sumw2();
M_b13->Sumw2();
M_b23->Sumw2();
higgs->Sumw2();
higgs_leptonic->Sumw2();
higgs_hadronic->Sumw2();
higgs_electron->Sumw2();
higgs_muon->Sumw2();


bJet_Pt->Scale(xs/ecents);
Jet1_Pt->Scale(xs/events);
Jet2_Pt->Scale(xs/events);
Jet3_Pt->Scale(xs/events);
Jet4_Pt->Scale(xs/events);
leading_electron->Scale(xs/events);
leading_muon->Scale(xs/events);
missingET->Scale(xs/events);
M_12_postSelection->Scale(xs/events);
M_b12_postSelection->Scale(xs/events);
M_3emu_postSelection->Scale(xs/events);
//M_12_col->Scale(xs/M_12_col->GetEntries());

TFile* output = new TFile("tch.root","RECREATE");
output->cd();
//TCanvas* c1 = new TCanvas("c1","c1",800,800);
//leading_lepton->Draw("e");
M_b12->Scale(xs/events);
M_b12->Write("M_b12");

//TCanvas* c2 = new TCanvas("c2","c2",800,800);
M_b13->Scale(xs/events);
M_b13->Write("M_b13");
//TCanvas* c3 = new TCanvas("c3","c3",800,800);
//tau_hadron_PT->Draw("e");
M_b23->Scale(xs/events);
M_b23->Write("M_b23");
//TCanvas* ce = new TCanvas("ce","ce",800,800);
//M_eta->Scale(xs/M_eta->GetEntries());
//M_eta->Write("M_eta");
//M_12_col->Write("M_12_col");
M_12_postSelection->Write("M_12_postSelection");
M_b12_postSelection->Write("M_b12_postSelection");
M_3emu_postSelection->Write("M_3emu_postSelection");

//TCanvas* c4 = new TCanvas("c4","c4",800,800);
higgs->Scale(xs/events);
higgs->Write("higgs_from_l/h");
//TCanvas* c5 = new TCanvas("c5","c5",800,800);
higgs_leptonic->Scale(xs/events);
higgs_leptonic->Write("higgs_leptonic");
//TCanvas* c6 = new TCanvas("c6","c6",800,800);
higgs_hadronic->Scale(xs/events);
higgs_hadronic->Write("higgs_hadronic");
//TCanvas* c7 = new TCanvas("c7","c7",800,800);
higgs_electron->Scale(xs/events);
higgs_electron->Write("higgs_electron");
//TCanvas* c8 = new TCanvas("c8","c8",800,800);
higgs_muon->Scale(xs/events);
higgs_muon->Write("higgs_muonic");
Jet1_Pt->Write("Jet1 Pt");
Jet2_Pt->Write("Jet2 Pt");
Jet3_Pt->Write("Jet3 Pt");
Jet4_Pt->Write("Jet4 Pt");
bJet_Pt->Write("bJet Pt");
leading_electron->Write("Leading Electron Pt");
leading_muon->Write("Leading Muon Pt");
missingET->Write("Missing Et");
cout << "M_b12_accuracy = " << M_b12_accuracy << endl;
cout << "M_b13_accuracy = " << M_b13_accuracy << endl;
cout << "M_b23_accuracy = " << M_b23_accuracy << endl;
cout << "M_eta_accuracy = " << M_eta_accuracy << endl;
cout << "total events that are in lep/had higgs window = " << acceptance <<endl << (float)(acceptance)/(float)(events) << " = Acceptance"<<endl;
cout << "total events that are in hadronic higgs window = " << acceptance_hadronic <<endl << (float)(acceptance_hadronic)/(float)(events) << " = Acceptance"<<endl;
cout << "total events that are in emu+- higgs window = " << acceptance_lep << endl<< (float)(acceptance_lep)/(float)(events) << " = Acceptance"<<endl;
cout << "total events that are in electron higgs window = " << acceptance_electron << endl<< (float)(acceptance_electron)/(float)(events) << " = Acceptance"<<endl;
cout << "total events that are in muon higgs window = " << acceptance_muon << endl<< (float)(acceptance_muon)/(float)(events) << " = Acceptance"<<endl;
output->Write();
}
