#define fnch_cxx
#include "fnch.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void fnch::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L fnch.C
//      Root > fnch t
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

TH1D* top_hadronic = new TH1D("top_from_jets","Top Quark from Jets Mass",50,0,300);
TH1D* w_from_jets = new TH1D("w_from_jets","W Boson from Jets Mass",50,0,200); 
TH1D* higgs = new TH1D("higgs","Higgs lep/had invariant mass",50,0,250);
TH1D* electrons = new TH1D("electrons","Leading electrons",30,0,200);
TH1D* leading_electron = new TH1D("leading_electron","Leading electron",30,0,200);
TH1D* missingET = new TH1D("misingET","Missing Energy",30,0,200);
TH1D* M_b12 = new TH1D("M_b12","Invariant Mass of b+1+2 Jets",35,0,300); 
TH1D* M_b13 = new TH1D("M_b13","Invariant Mass of b+1+3 Jets",35,0,300); 
TH1D* M_b23 = new TH1D("M_23","Invariant Mass of b+2+3 Jets",35,0,300); 
TH1D* M_eta = new TH1D("M_eta","Invariant Mass of smallest eta",35,0,300);
TH1D* leading_lepton = new TH1D("leading_lepton","Leading Lepton PT",35,0,300);
TH1D* tau_hadron_PT = new TH1D("tau_hadron_PT","Hadron from tau Pt",30,0,300);
TH1D* higgs_leptonic = new TH1D("higgs_leptonic","Higgs invariant mass from e+mu opposite sign",30,0,300);
int acceptance_lep=0;
int acceptance=0;
double M_b12_accuracy=0;
double M_b13_accuracy=0;
double M_b23_accuracy=0;
double M_eta_accuracy=0;
static double top_mass = 173.2;
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

int numbJets=0;
int numLightJets=0;
int numTauJets=0;
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
				
				//cout << jentry << endl;
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
			 if(numTauJets==1 && Jet_TauTag[i]==1){
                       		v_tau2.SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i],Jet_Phi[i],Jet_Mass[i]);
                        	numTauJets++;
                                        }		
		}//end of Jet Loop
		//cout << "Jet PT's " << v_j1.Pt()<<" "<<v_j2.Pt()<<" "<<v_j3.Pt() <<endl;
		
		double  delta_eta_b12 = fabs(v_bJet.Eta()-v_j1.Eta())+fabs(v_bJet.Eta()-v_j2.Eta());
		double  delta_eta_b13 = fabs(v_bJet.Eta()-v_j1.Eta())+fabs(v_bJet.Eta()-v_j3.Eta());
		double  delta_eta_b23 = fabs(v_bJet.Eta()-v_j2.Eta())+fabs(v_bJet.Eta()-v_j3.Eta());
		//cout<<"debug!"<<endl;
		if(delta_eta_b12 < delta_eta_b13 && delta_eta_b12 < delta_eta_b23){
			v_M_eta=v_bJet+v_j1+v_j2;
		}
		  if(delta_eta_b13 < delta_eta_b12 && delta_eta_b13 < delta_eta_b23){
                        v_M_eta=v_bJet+v_j1+v_j3;
                }
		  if(delta_eta_b23 < delta_eta_b12 && delta_eta_b23 < delta_eta_b13){
                        v_M_eta=v_bJet+v_j2+v_j3;
                }


		v_M_b12 = v_bJet+v_j1+v_j2;
		v_M_b13 = v_bJet+v_j1+v_j3;
		v_M_b23 = v_bJet+v_j2+v_j3;
//cout << " Masses " << v_M_b12.M() << " " << v_M_b13.M() <<" "<< v_M_b23.M() << endl;
//cout << numTauJets << endl;

//if(numTauJets==2)
//{
//v_higgs = v_tau1 + v_tau2;
//cout << v_higgs.M();
//higgs->Fill(v_higgs.M());
//}

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

//This will be the computational part
//
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

if(v_higgs.M()>105 && v_higgs.M()<145)
acceptance++;
}
//------------------------------------------------------------------------
//The below is for e+mu opposite sign
if(v_electron.Pt()>2 && MissingET_MET[0]>2 && v_muon.Pt()>2){
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


M_b12->Fill(v_M_b12.M());
M_b13->Fill(v_M_b13.M());
M_b23->Fill(v_M_b23.M());
M_eta->Fill(v_M_eta.M());
M_b12_accuracy = M_b12_accuracy + fabs(v_M_b12.M()-top_mass);
M_b13_accuracy = M_b13_accuracy + fabs(v_M_b13.M()-top_mass);
M_b23_accuracy = M_b23_accuracy + fabs(v_M_b23.M()-top_mass);
M_eta_accuracy = M_eta_accuracy + fabs(v_M_eta.M()-top_mass);



if(v_higgs.M()>105 && v_higgs.M()<145)
acceptance_lep++;
}
//Ends the Compuations
//

electrons->Fill(v_electronSum.M());
missingET->Fill(MissingET_MET[0]);
if(v_tau1.Pt()>2.0)
tau_hadron_PT->Fill(v_tau1.Pt());

//M_b12->Fill(v_M_b12.M());
//M_b13->Fill(v_M_b13.M());
//M_b23->Fill(v_M_b23.M());
//M_eta->Fill(v_M_eta.M());
//M_b12_accuracy = M_b12_accuracy + fabs(v_M_b12.M()-top_mass);
//M_b13_accuracy = M_b13_accuracy + fabs(v_M_b13.M()-top_mass);
//M_b23_accuracy = M_b23_accuracy + fabs(v_M_b23.M()-top_mass);
//M_eta_accuracy = M_eta_accuracy + fabs(v_M_eta.M()-top_mass);

//v_top_from_jets = v_bJet+v_j1+v_j2;
//v_w = v_j1+v_j2;
//if(numbJets==1 && numLightJets==2)
//top_hadronic->Fill(v_top_from_jets.M());

//if(numLightJets==2)
//w_from_jets->Fill(v_w.M());
if (jentry>5000)
break;

  }//Event Loop

TCanvas* c1 = new TCanvas("c1","c1",800,800);
//leading_lepton->Draw("e");
M_b12->Draw("e");
TCanvas* c2 = new TCanvas("c2","c2",800,800);
//missingET->Draw("e");
M_b13->Draw("e");
TCanvas* c3 = new TCanvas("c3","c3",800,800);
//tau_hadron_PT->Draw("e");
M_b23->Draw("e");
TCanvas* c6 = new TCanvas("c6","c6",800,800);
M_eta->Draw("e");
TCanvas* c4 = new TCanvas("c4","c4",800,800);
higgs->Draw("e");
TCanvas* c5 = new TCanvas("c5","c5",800,800);
higgs_leptonic->Draw("e");

cout << "M_b12_accuracy = " << M_b12_accuracy << endl;
cout << "M_b13_accuracy = " << M_b13_accuracy << endl;
cout << "M_b23_accuracy = " << M_b23_accuracy << endl;
cout << "M_eta_accuracy = " << M_eta_accuracy << endl;
cout << "total events that are in higgs window = " << acceptance <<endl;
cout << "total events that are in higgs_leptonic window = " << acceptance_lep << endl;
}
