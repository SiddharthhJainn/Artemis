// Basic setup for Deeply Inelastic Scattering at HERA.

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include <ROOT/TDataFrame.hxx>

// ROOT, for saving Pythia events as trees in a file.
#include "TTree.h"
#include "TFile.h"
#include <TLorentzVector.h>
#include <vector>
#include "Pythia8Plugins/FastJet3.h"

using namespace fastjet;

using namespace Pythia8; 

int main() {

// Beam energies, minimal Q2, number of events to generate

  double eProton   = 820.;                            // proton energy
  double eElectron = 27.5;                            // positron energy
  double Q2min     = 125.;                            // Q2 minimum
  int    nEvent    = 100000;                          // number of events
  double EtMin = 10.;
  double R = 0.8;                                     // Jet radius

  TFile *file = TFile::Open("008mymain08.root","recreate");
  double ycut = 0.0005;

  // Generator.
  Pythia pythia;
  
  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Set up incoming beams, for frame with unequal beam energies.
  pythia.readString("Beams:frameType = 2");
  // BeamA = proton.
  pythia.readString("Beams:idA = 2212");
  pythia.settings.parm("Beams:eA", eProton);
  // BeamB = electron.
  pythia.readString("Beams:idB = -11");
  pythia.settings.parm("Beams:eB", eElectron);

  // Set up DIS process within some phase space.
  // Neutral current (with gamma/Z interference).
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  
  // Phase-space cut: minimal Q2 of process.
  pythia.settings.parm("PhaseSpace:Q2Min", Q2min);

  // Set dipole recoil on. Necessary for DIS + shower.
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // QED radiation off lepton not handled yet by the new procedure.
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");


 //Specifying critical lifetime of particles
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 3");

  // Initialize.
  pythia.init();

  // Set up FastJet jet finder.

  //JetDefinition jet_def(kt_algorithm, R,E_scheme);
  JetDefinition jet_def(kt_algorithm, R,E_scheme);

  std::vector <fastjet::PseudoJet> fjInputs;

//--------------------- Variables -------------------------------------------------------------
  
 Double_t hpx[500],hpy[500],hpz[500],he[500];

 Double_t jpx[50], jpy[50], jpz[50], je[50], jeta[50], jphi[50], jtheta[50], jpt[50], jee[50] , dist[50];

 //Double_t jjpx[50]; // jjpy[50], jjpz[50], jje[50], jjeta[50], jjphi[50], jjtheta[50], jjpt[50], jjee[50];

 Double_t Asj_px[50], Asj_py[50], Asj_pz[50], Asj_e[50], Asj_eta[50], Asj_phi[50], Asj_theta[50], Asj_pt[50];

 Double_t epx , epy , epz , ee, ept, etheta, ye, eeta, ephi;

 Int_t nb, jnum, jjnum, Asj;
//----------------------------------------------------------------------------------------------

  // Set up the ROOT TFile and TTree.

  //TFile *file = TFile::Open("008mymain08.root","recreate");
  
  //double ycut = 0.0005;

  TTree *T = new TTree("Tree","ev1 Tree");
  
  T->Branch("epx",&epx);
  T->Branch("epy",&epy);
  T->Branch("epz",&epz);
  T->Branch("ee",&ee);
  T->Branch("ept",&ept);
  T->Branch("etheta",&etheta);
  T->Branch("ye",&ye);
  
  T->Branch("nb",&nb,"nb/I");
  T->Branch("hpx",hpx,"hpx[nb]/D");
  T->Branch("hpy",hpy,"hpy[nb]/D");
  T->Branch("hpz",hpz,"hpz[nb]/D");
  T->Branch("he",he,"he[nb]/D");

  T->Branch("jnum",&jnum,"jnum/I");
  T->Branch("jpx",jpx,"jpx[jnum]/D");
  T->Branch("jpy",jpy,"jpy[jnum]/D");
  T->Branch("jpz",jpz,"jpz[jnum]/D");
  T->Branch("je",je,"je[jnum]/D");
  T->Branch("jeta",jeta,"jeta[jnum]/D");
  T->Branch("jphi",jphi,"jphi[jnum]/D");
  T->Branch("jtheta",jtheta,"jtheta[jnum]/D");
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
  T->Branch("Asj_theta",Asj_theta,"Asj_theta[Asj]/D");
  T->Branch("Asj_pt",Asj_pt,"Asj_pt[Asj]/D");

  int d, d1 , d2; 

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;

 //...............................electrons.......................................................................

  for (int i=0; i < event.size(); ++i) {
  if(event[i].id()==-11 && event[i].mother1() == 0 && event[i].mother2() == 0) {
  
     d = i;
          
  for(int j=0; j < 100; ++j) {
  d1 = event[d].daughter1();
  
  if(event[d1].id()==-11){
  d = d1;
  
     if(event[d1].isFinal()){
        break;   
  }
 }
  
   d2 = event[d].daughter2();
  
  if(event[d2].id()==-11){
  d = d2;
  
     if(event[d2].isFinal()){
        break;
      }
     }
  }                                                                                 //for decay electron
 }                                                                                // condition for initial electron 
}                                                                               //end of electron loop
   
    epx = event[d].px();
    epy = event[d].py();
    epz = event[d].pz();
    ee = event[d].e();
    ept = event[d].pT();
    etheta = event[d].theta();

    eeta = event[d].eta();
    ephi = event[d].phi();

    ye = 1-ee*(1-cos(etheta))/(2*eElectron);
    
    if(ee<10) continue;
    if(ye>0.95) continue;
    if(etheta > 2.443) continue;

//--------------------------------------------------Hadrons--------------------------------------------------------------------------        

  //hadrons for fastjet
 // Begin FastJet analysis: extract particles from event record.
    fjInputs.resize(0);            
    Vec4   pTemp;                  // for four momentum
    double mTemp;		              // for storing mass
    int y = 0;        
 
   for (int l = 0; l < event.size(); ++l){
   //if(event[l].isFinal() && l != d && event[l].isHadron()){
   if(event[l].isFinal() && l != d){
   
   hpx[y] = event[l].px();
   hpy[y] = event[l].py();
   hpz[y] = event[l].pz();
   he[y] = event[l].e();

 // Create a PseudoJet from the complete Pythia particle.
      fastjet::PseudoJet particleTemp = event[l];    //All jets, as well as input particles to the clustering are PseudoJet objects
  
      // Store acceptable particles as input to Fastjet.
      // Conversion to PseudoJet is performed automatically
      // with the help of the code in FastJet3.h.
      fjInputs.push_back( particleTemp);
  
    }                                       // only hadron loop
   y++;
   }                                                             // all particle loop

 nb = y;
  
//------------------------------------Jets--------------------------------------------------------------------
   
    // Run Fastjet anti-kT algorithm and sort jets in pT order.
   //  ClusterSequence clust_seq(fjInputs, jet_def);
   // double ptmin = 5.0;
   //vector<PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin)); 

    fastjet::ClusterSequence clustSeq1( fjInputs, jet_def );
    Selector jet_selector = SelectorEtMin(EtMin);
    vector <fastjet::PseudoJet> inclusive_jets = sorted_by_pt( jet_selector(clustSeq1.inclusive_jets()) );
  
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
   jtheta[z] = inclusive_jets[i].theta();
   jpt[z] = inclusive_jets[i].pt();
   jee[z] = inclusive_jets[i].Et();            // returns the transerve energy
  
   dist[z] = sqrt(pow((jeta[z]-eeta),2) + pow((jphi[z]-ephi),2));
   
   if(dist[z]<1) break;

   if(jeta[z] > -1. && jeta[z] < 2. && jee[z]> 15.){   
   
   double dcut =  ycut*pow(jee[z],2);
   //Rs = R*pow(dcut,0.5);
   //JetDefinition jetDefCA( kt_algorithm, Rs,E_scheme);

   vector<PseudoJet> AsubJets = sorted_by_pt(inclusive_jets[i].exclusive_subjets(dcut));

   for(unsigned int j = 0; j < AsubJets.size(); j++){             //subjet loop
  
   Asj_px[c] = AsubJets[j].px();            // momentum and energy
   Asj_py[c] = AsubJets[j].py();
   Asj_pz[c] = AsubJets[j].pz();
   Asj_e[c] = AsubJets[j].e();
   Asj_eta[c] = AsubJets[j].eta();          // returns the pseudo-rapidity
   Asj_phi[c] = AsubJets[j].phi();          //returns the azimuthal angle in range 0 . . . 2π
   Asj_theta[c] = AsubJets[j].theta(); 
   Asj_pt[c] = AsubJets[j].pt();

   c++;
     }       // for loop subjet
     a++;  
    }         // if condition

   z++;
  }
   if(dist[z]<1) continue;

   if(a == 0) continue;

 jnum = z;
 Asj = c;
   
    T->Fill();  
 }    // End of event loop.

  cout<< "Jet Info : " << jet_def.description() << endl;          // prints out jet definitation descriptions
  
  pythia.stat();
  
  //Write tree
  T->Print();
  T->Write();

  delete file;

  // Done.
  return 0; 
}
