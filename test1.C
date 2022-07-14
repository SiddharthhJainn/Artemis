#include "fastjet/ClusterSequence.hh"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
using namespace fastjet;

#include "Math/IFunction.h"
using namespace TMath;

#include <iostream>
using namespace std;

void test1(){

//------------------------------------------------------------------------------------------------------------
Double_t hpx[5000], hpy[5000] , hpz[5000] , he[5000];
Double_t epx[10] , epy[10] , epz[10] , ee[10];

Double_t jpx[50], jpy[50], jpz[50], je[50], jeta[50], jphi[50], jtheta[50], jpt[50], jee[50] , dist[50];
Double_t Asj_px[50], Asj_py[50], Asj_pz[50], Asj_e[50], Asj_eta[50], Asj_phi[50], Asj_theta[50], Asj_pt[50];
Int_t ne , nh , jnum, Asj;
//------------------------------------------------------------------------------------------------------------
double R = 0.8;
JetDefinition jet_def(kt_algorithm, R);
vector<PseudoJet> particles;
double EtMin = 10.;
double ycut = 0.0005;
//-------------------------------------------------------------------------------------------------------------
TFile *f = new TFile("input.root");
TTree *t1 = (TTree*)f->Get("etree");
TTree *t2 = (TTree*)f->Get("htree");

t1->SetBranchAddress("ne",&ne);
t1->SetBranchAddress("epx",epx);
t1->SetBranchAddress("epy",epy);
t1->SetBranchAddress("epz",epz);
t1->SetBranchAddress("ee",ee);

t2->SetBranchAddress("nh",&nh);
t2->SetBranchAddress("hpx",hpx);
t2->SetBranchAddress("hpy",hpy);
t2->SetBranchAddress("hpz",hpz);
t2->SetBranchAddress("he",he);
//-------------------------------------------------------------------------------------------------------------
  // Set up the ROOT TFile and TTree.

  TFile *file = new TFile("output.root","recreate");

  TTree *T = new TTree("Trree","ev1 Tree");                 // trree is the name of the tree. Be cautious!

  T->Branch("jnum",&jnum,"jnum/I");
  T->Branch("jpx",jpx,"jpx[jnum]/D");
  T->Branch("jpy",jpy,"jpy[jnum]/D");
  T->Branch("jpz",jpz,"jpz[jnum]/D");
  T->Branch("je",je,"je[jnum]/D");
  T->Branch("jeta",jeta,"jeta[jnum]/D");
  T->Branch("jphi",jphi,"jphi[jnum]/D");
  //T->Branch("jtheta",jtheta,"jtheta[jnum]/D");
  T->Branch("jpt",jpt,"jpt[jnum]/D");
  T->Branch("jee",jee,"jee[jnum]/D");
  //T->Branch("dist",dist,"dist[jnum]/D");

  T->Branch("Asj",&Asj,"Asj/I");
  T->Branch("Asj_px",Asj_px,"Asj_px[Asj]/D");
  T->Branch("Asj_py",Asj_py,"Asj_py[Asj]/D");
  T->Branch("Asj_pz",Asj_pz,"Asj_pz[Asj]/D");
  T->Branch("Asj_e",Asj_e,"Asj_e[Asj]/D");
  T->Branch("Asj_eta",Asj_eta,"Asj_eta[Asj]/D");
  T->Branch("Asj_phi",Asj_phi,"Asj_phi[Asj]/D");
  //T->Branch("Asj_theta",Asj_theta,"Asj_theta[Asj]/D");
  T->Branch("Asj_pt",Asj_pt,"Asj_pt[Asj]/D");
//-------------------------------------------------------------------------------------------------------------
Int_t e_entries = (Int_t)t1->GetEntries();
Int_t h_entries = (Int_t)t2->GetEntries();

cout << h_entries << endl;
cout << e_entries << endl;

TLorentzVector v;

Double_t etheta[100000], ephi[100000], eeta[100000], ye[100000], energy[100000];
double eElectron = 27.5;
Int_t b = 0;

// event loop for positrons
for(Int_t i = 0; i < e_entries; i++){
t1->GetEntry(i);

//positrons
 for(Int_t j = 0; j < ne; j++){
  
  v.SetPxPyPzE(epx[j],epy[j], epz[j],ee[j]);
  
  energy[b] = ee[j];
  etheta[b] = v.Theta();
  //cout << "next event" << endl;
  ephi[b] = v.Phi();
  eeta[b] = v.PseudoRapidity();

  ye[b] = 1-energy[b]*(1-cos(etheta[b]))/(2*eElectron);

  //cout << eeta[b] << endl;
  b++;
 }
 //cout << "positrons are: " << e << endl;
}

 cout<< "b is " << b << endl;

//-------------------------hadrons-----------------------------

for(Int_t i = 0; i < h_entries; i++){
t2->GetEntry(i);
particles.resize(0);

//e+ cuts
  if(energy[i]<10.) continue;
 //if(ee[i]<10.) continue;
  if(ye[i]>0.95) continue;
  if(etheta[i]>2.443) continue;

 for(Int_t j = 0; j < nh; j++){
 particles.push_back( PseudoJet( hpx[j],  hpy[j],  hpz[j], he[j]) );
 }
 
 //....................Jets........................................

 fastjet::ClusterSequence clustSeq1( particles , jet_def );
 Selector jet_selector = SelectorEtMin(EtMin);
 vector <fastjet::PseudoJet> inclusive_jets = sorted_by_pt( jet_selector(clustSeq1.inclusive_jets()) );
   
   int z = 0;
   int c = 0;
   int a = 0;

   for (unsigned int k = 0; k < inclusive_jets.size(); k++){
  
   jpx[z] = inclusive_jets[k].px();            // momentum and energy
   jpy[z] = inclusive_jets[k].py();
   jpz[z] = inclusive_jets[k].pz();
   je[z] = inclusive_jets[k].e();
   jeta[z] = inclusive_jets[k].eta();          // returns the pseudo-rapidity
   jphi[z] = inclusive_jets[k].phi();          //returns the azimuthal angle in range 0 . . . 2π
   //jtheta[z] = inclusive_jets[i].theta();
   jpt[z] = inclusive_jets[k].pt();
   jee[z] = inclusive_jets[k].Et();            // returns the transerve energy
   
   //Jet cuts
    dist[z] = sqrt(pow((jeta[z]-eeta[i]),2) + pow((jphi[z]-ephi[i]),2));
    if(dist[z]<1) continue;
    if(inclusive_jets.size() == 0) continue;

   if(jeta[z] > -1. && jeta[z] < 2. && jee[z]> 15.){
   
   double dcut =  ycut*pow(jee[z],2);

   vector<PseudoJet> AsubJets = sorted_by_pt(inclusive_jets[k].exclusive_subjets(dcut));
   //cout<< "Subjets:  " << AsubJets.size() << endl;

   for(unsigned int j = 0; j < AsubJets.size(); j++){             //subjet loop
  
   Asj_px[c] = AsubJets[j].px();            // momentum and energy
   Asj_py[c] = AsubJets[j].py();
   Asj_pz[c] = AsubJets[j].pz();
   Asj_e[c] = AsubJets[j].e();
   Asj_eta[c] = AsubJets[j].eta();          // returns the pseudo-rapidity
   Asj_phi[c] = AsubJets[j].phi();          //returns the azimuthal angle in range 0 . . . 2π
   //Asj_theta[c] = AsubJets[j].theta(); 
   Asj_pt[c] = AsubJets[j].pt();

   c++;
     }       // for loop subjet
     a++;  
    }         // if condition

   z++;
  }
   if(a == 0) continue;

   jnum = z;
   Asj = c;
   
    T->Fill();  
 }

cout << "Clustering with " << jet_def.description() << endl;

//TCanvas *c1 = new TCanvas();
//hist->Draw();
//h1->Draw();

f->Close();

file->Write();
file->Close();
//delete file;
}