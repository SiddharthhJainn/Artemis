#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "fastjet/ClusterSequence.hh"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
using namespace fastjet;

#include <TLorentzVector.h>
#include <vector>

#include "Math/IFunction.h"
using namespace TMath;

#include <iostream>
using namespace std;

void MyClass::Loop()
{

  // Set up FastJet jet finder.

  //-------------------------------------------------------------------------------------------
   double ycut = 0.0005;
   double R = 0.8;
   JetDefinition jet_def(kt_algorithm, R);
   vector<PseudoJet> particles;
   double EtMin = 10.;
   Double_t eElectron = 27.5;
//---------------------------------------------------------------------------------------------

  std::vector <fastjet::PseudoJet> fjInputs;

//--------------------- Variables -------------------------------------------------------------
  
 Double_t hpx[500],hpy[500],hpz[500],he[500];

 Double_t jpx[50], jpy[50], jpz[50], je[50], jeta[50], jphi[50], jtheta[50], jpt[50], jee[50] , dist[50];

 //Double_t jjpx[50]; // jjpy[50], jjpz[50], jje[50], jjeta[50], jjphi[50], jjtheta[50], jjpt[50], jjee[50];

 Double_t Asj_px[50], Asj_py[50], Asj_pz[50], Asj_e[50], Asj_eta[50], Asj_phi[50], Asj_theta[50], Asj_pt[50];

 Double_t epx , epy , epz , ee, ept, etheta, ye, eeta, ephi;

 Int_t nh, jnum, Asj;
//----------------------------------------------------------------------------------------------

 // Set up the ROOT TFile and TTree.
 
  TFile *file = TFile::Open("008mymain08.root","recreate");

  TTree *T = new TTree("Tree","ev1 Tree");
  
  T->Branch("epx",&epx);
  T->Branch("epy",&epy);
  T->Branch("epz",&epz);
  T->Branch("ee",&ee);
  T->Branch("ept",&ept);
  T->Branch("etheta",&etheta);
  T->Branch("ye",&ye);
  
  T->Branch("nh",&nh,"nh/I");
  T->Branch("hpx",hpx,"hpx[nh]/D");
  T->Branch("hpy",hpy,"hpy[nh]/D");
  T->Branch("hpz",hpz,"hpz[nh]/D");
  T->Branch("he",he,"he[nh]/D");

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
  T->Branch("dist",dist,"dist[jnum]/D");

  T->Branch("Asj",&Asj,"Asj/I");
  T->Branch("Asj_px",Asj_px,"Asj_px[Asj]/D");
  T->Branch("Asj_py",Asj_py,"Asj_py[Asj]/D");
  T->Branch("Asj_pz",Asj_pz,"Asj_pz[Asj]/D");
  T->Branch("Asj_e",Asj_e,"Asj_e[Asj]/D");
  T->Branch("Asj_eta",Asj_eta,"Asj_eta[Asj]/D");
  T->Branch("Asj_phi",Asj_phi,"Asj_phi[Asj]/D");
  //T->Branch("Asj_theta",Asj_theta,"Asj_theta[Asj]/D");
  T->Branch("Asj_pt",Asj_pt,"Asj_pt[Asj]/D");

   int b = 0;
   int del;
   int e;

   TLorentzVector v;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      /*
      int a = 0;
      cout<<"Total particles: " << particles_ << endl;

      for(int i = 0; i < particles_ ; i++){
         if (particles_pid[i]== -11 && particles_status[i] == 1) a++;

         cout << i << "  " << particles_pid[i] << "   status  " << particles_status[i] << "   prodV.  " << links1[i] << "   decayV.  " << links2[i] << endl;
      }
      if(a > 1)
      cout << "  " << endl;
      cout << jentry << " with a =  " << a << endl;
      cout << "  " << endl; 
       */
       
       // /*
       
       for(int i = 0; i < particles_ ; i++){

       if(particles_pid[i]== -11 && particles_status[i] == 4){
            e = links2[i];
       }

       //cout<< d << endl; 
       if(links1[i]== e && particles_pid[i] == -11 && particles_status[i] == 1){
         del = i;
         //cout << del << endl;
         }
       }
         
         epx = particles_momentum_m_v1[del];
         epy = particles_momentum_m_v2[del];
         epz = particles_momentum_m_v3[del];
         ee = particles_momentum_m_v4[del];

         v.SetPxPyPzE(epx,epy, epz,ee);
          //energy = v.E();
          etheta = v.Theta();
          ephi = v.Phi();
          eeta = v.PseudoRapidity();
          ye = 1-ee*(1-cos(etheta))/(2*eElectron);

          if(ee<10.) continue;
          if(ye>0.95) continue;
          if(etheta>2.443) continue;
        

        fjInputs.resize(0);
        int y = 0;
        // hadrons
        for(int i = 0; i < particles_ ; i++){
        if(i != del && particles_status[i]==1){

        hpx[y] = particles_momentum_m_v1[i];
        hpy[y] = particles_momentum_m_v2[i];
        hpz[y] = particles_momentum_m_v3[i];
        he[y] = particles_momentum_m_v4[i];

        fjInputs.push_back( PseudoJet( hpx[y], hpy[y], hpz[y], he[y]) );
        
        y++;
        }
        nh = y;
       }
    
    
    fastjet::ClusterSequence clustSeq1( fjInputs, jet_def );
    Selector jet_selector = SelectorEtMin(EtMin);
    vector <fastjet::PseudoJet> inclusive_jets = sorted_by_pt( jet_selector(clustSeq1.inclusive_jets()) );
   
    cout << inclusive_jets.size() << endl;
       if(inclusive_jets.size()==0) continue;
   
   int z = 0;
   int c = 0;
   int a = 0;

   for (unsigned int i = 0; i < inclusive_jets.size(); i++){
  
   jpx[z] = inclusive_jets[i].px();            // momentum and energy
   jpy[z] = inclusive_jets[i].py();
   jpz[z] = inclusive_jets[i].pz();
   je[z] = inclusive_jets[i].e();
   jeta[z] = inclusive_jets[i].eta();          // returns the pseudo-rapidity
   jphi[z] = inclusive_jets[i].phi();          //returns the azimuthal angle in range 0 . . . 2π
   //jtheta[z] = inclusive_jets[i].theta();
   jpt[z] = inclusive_jets[i].pt();
   jee[z] = inclusive_jets[i].Et();            // returns the transerve energy

   dist[z] = sqrt(pow((jeta[z]-eeta),2) + pow((jphi[z]-ephi),2));
   if(dist[z]<1) break;

   if(jeta[z] > -1. && jeta[z] < 2. && jee[z]> 15.){
      double dcut =  ycut*pow(jee[z],2);
    vector<PseudoJet> AsubJets = sorted_by_pt(inclusive_jets[i].exclusive_subjets(dcut));
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
     }       // subjet loop
     a++;

   }
   z++;
   }

   if(dist[z]<1) continue;
   if(a == 0) continue;
   
    jnum = z;
    Asj = c;

    T->Fill();
    }  // end of event

    //Write tree
    T->Print();
    T->Write();
    
    delete file;
    
   //  file->Write();
   //  file->Close();
}