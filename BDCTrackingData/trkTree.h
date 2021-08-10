//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug  9 10:01:57 2021 by ROOT version 6.22/06
// from TTree trkTree/BDC Track Tree
// found on file: bdcAnaTrack_Data_SJ_Run_520_selecttrack_20210808_v3.root
//////////////////////////////////////////////////////////

#ifndef trkTree_h
#define trkTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class trkTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event;
   Int_t           trckNumX;
   Int_t           trckNumY;
   Double_t        Xgrad[1];   //[trkNumX]
   Double_t        Ygrad[1];   //[trkNumY]
   Double_t        Xc[1];   //[trkNumX]
   Double_t        Yc[1];   //[trkNumY]
   Int_t           EvtTime;
   Double_t        dur_sec;
   Int_t           EvtTimeDif;
   Double_t        dur_secDif;
   Int_t           EvtTag;

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_trkNumX;   //!
   TBranch        *b_trkNumY;   //!
   TBranch        *b_Xgrad;   //!
   TBranch        *b_Ygrad;   //!
   TBranch        *b_Xc;   //!
   TBranch        *b_Yc;   //!
   TBranch        *b_EvtTime;   //!
   TBranch        *b_dur_sec;   //!
   TBranch        *b_EvtTimeDif;   //!
   TBranch        *b_dur_secDif;   //!
   TBranch        *b_EvtTag;   //!

   trkTree(TTree *tree=0);
   virtual ~trkTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef trkTree_cxx
trkTree::trkTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("bdcAnaTrack_Data_SJ_Run_520_selecttrack_20210808_v3.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("bdcAnaTrack_Data_SJ_Run_520_selecttrack_20210808_v3.root");
      }
      f->GetObject("trkTree",tree);

   }
   Init(tree);
}

trkTree::~trkTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t trkTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t trkTree::LoadTree(Long64_t entry)
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

void trkTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("trckNumX", &trckNumX, &b_trkNumX);
   fChain->SetBranchAddress("trckNumY", &trckNumY, &b_trkNumY);
   fChain->SetBranchAddress("Xgrad", Xgrad, &b_Xgrad);
   fChain->SetBranchAddress("Ygrad", Ygrad, &b_Ygrad);
   fChain->SetBranchAddress("Xc", Xc, &b_Xc);
   fChain->SetBranchAddress("Yc", Yc, &b_Yc);
   fChain->SetBranchAddress("EvtTime", &EvtTime, &b_EvtTime);
   fChain->SetBranchAddress("dur_sec", &dur_sec, &b_dur_sec);
   fChain->SetBranchAddress("EvtTimeDif", &EvtTimeDif, &b_EvtTimeDif);
   fChain->SetBranchAddress("dur_secDif", &dur_secDif, &b_dur_secDif);
   fChain->SetBranchAddress("EvtTag", &EvtTag, &b_EvtTag);
   Notify();
}

Bool_t trkTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void trkTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t trkTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef trkTree_cxx
