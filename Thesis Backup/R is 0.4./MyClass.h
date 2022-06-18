//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 13 16:45:29 2022 by ROOT version 6.24/02
// from TTree Tree/ev1 Tree
// found on file: 0005.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        epx;
   Double_t        epy;
   Double_t        epz;
   Double_t        ee;
   Double_t        ept;
   Double_t        etheta;
   Double_t        ye;
   Int_t           nb;
   Double_t        hpx[52];   //[nb]
   Double_t        hpy[52];   //[nb]
   Double_t        hpz[52];   //[nb]
   Double_t        he[52];   //[nb]
   Int_t           jnum;
   Double_t        jpx[4];   //[jnum]
   Double_t        jpy[4];   //[jnum]
   Double_t        jpz[4];   //[jnum]
   Double_t        je[4];   //[jnum]
   Double_t        jeta[4];   //[jnum]
   Double_t        jphi[4];   //[jnum]
   Double_t        jtheta[4];   //[jnum]
   Double_t        jpt[4];   //[jnum]
   Double_t        jee[4];   //[jnum]
   Double_t        dist[4];   //[jnum]
   Int_t           Asj;
   Double_t        Asj_px[13];   //[Asj]
   Double_t        Asj_py[13];   //[Asj]
   Double_t        Asj_pz[13];   //[Asj]
   Double_t        Asj_e[13];   //[Asj]
   Double_t        Asj_eta[13];   //[Asj]
   Double_t        Asj_phi[13];   //[Asj]
   Double_t        Asj_theta[13];   //[Asj]
   Double_t        Asj_pt[13];   //[Asj]

   // List of branches
   TBranch        *b_epx;   //!
   TBranch        *b_epy;   //!
   TBranch        *b_epz;   //!
   TBranch        *b_ee;   //!
   TBranch        *b_ept;   //!
   TBranch        *b_etheta;   //!
   TBranch        *b_ye;   //!
   TBranch        *b_nb;   //!
   TBranch        *b_hpx;   //!
   TBranch        *b_hpy;   //!
   TBranch        *b_hpz;   //!
   TBranch        *b_he;   //!
   TBranch        *b_jnum;   //!
   TBranch        *b_jpx;   //!
   TBranch        *b_jpy;   //!
   TBranch        *b_jpz;   //!
   TBranch        *b_je;   //!
   TBranch        *b_jeta;   //!
   TBranch        *b_jphi;   //!
   TBranch        *b_jtheta;   //!
   TBranch        *b_jpt;   //!
   TBranch        *b_jee;   //!
   TBranch        *b_dist;   //!
   TBranch        *b_Asj;   //!
   TBranch        *b_Asj_px;   //!
   TBranch        *b_Asj_py;   //!
   TBranch        *b_Asj_pz;   //!
   TBranch        *b_Asj_e;   //!
   TBranch        *b_Asj_eta;   //!
   TBranch        *b_Asj_phi;   //!
   TBranch        *b_Asj_theta;   //!
   TBranch        *b_Asj_pt;   //!

   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("1.root");
      }
      f->GetObject("Tree",tree);

   }
   Init(tree);
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("epx", &epx, &b_epx);
   fChain->SetBranchAddress("epy", &epy, &b_epy);
   fChain->SetBranchAddress("epz", &epz, &b_epz);
   fChain->SetBranchAddress("ee", &ee, &b_ee);
   fChain->SetBranchAddress("ept", &ept, &b_ept);
   fChain->SetBranchAddress("etheta", &etheta, &b_etheta);
   fChain->SetBranchAddress("ye", &ye, &b_ye);
   fChain->SetBranchAddress("nb", &nb, &b_nb);
   fChain->SetBranchAddress("hpx", hpx, &b_hpx);
   fChain->SetBranchAddress("hpy", hpy, &b_hpy);
   fChain->SetBranchAddress("hpz", hpz, &b_hpz);
   fChain->SetBranchAddress("he", he, &b_he);
   fChain->SetBranchAddress("jnum", &jnum, &b_jnum);
   fChain->SetBranchAddress("jpx", jpx, &b_jpx);
   fChain->SetBranchAddress("jpy", jpy, &b_jpy);
   fChain->SetBranchAddress("jpz", jpz, &b_jpz);
   fChain->SetBranchAddress("je", je, &b_je);
   fChain->SetBranchAddress("jeta", jeta, &b_jeta);
   fChain->SetBranchAddress("jphi", jphi, &b_jphi);
   fChain->SetBranchAddress("jtheta", jtheta, &b_jtheta);
   fChain->SetBranchAddress("jpt", jpt, &b_jpt);
   fChain->SetBranchAddress("jee", jee, &b_jee);
   fChain->SetBranchAddress("dist", dist, &b_dist);
   fChain->SetBranchAddress("Asj", &Asj, &b_Asj);
   fChain->SetBranchAddress("Asj_px", Asj_px, &b_Asj_px);
   fChain->SetBranchAddress("Asj_py", Asj_py, &b_Asj_py);
   fChain->SetBranchAddress("Asj_pz", Asj_pz, &b_Asj_pz);
   fChain->SetBranchAddress("Asj_e", Asj_e, &b_Asj_e);
   fChain->SetBranchAddress("Asj_eta", Asj_eta, &b_Asj_eta);
   fChain->SetBranchAddress("Asj_phi", Asj_phi, &b_Asj_phi);
   fChain->SetBranchAddress("Asj_theta", Asj_theta, &b_Asj_theta);
   fChain->SetBranchAddress("Asj_pt", Asj_pt, &b_Asj_pt);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
