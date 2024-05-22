//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 16 18:22:34 2018 by ROOT version 6.06/01
// from TTree hits/HGC rechits
// found on file: muon_v10.root
//////////////////////////////////////////////////////////

#ifndef HGCNtupleVariables_h
#define HGCNtupleVariables_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class HGCNtupleVariables {
public:
  HGCNtupleVariables(TTree * /*tree*/ = 0) : fChain(0) {}
  ~HGCNtupleVariables() {}
  // void    Init(TTree *tree);
  void Init(TTree *tree, TTree *tree2);
  Bool_t Notify();
  Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
    return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
  }

  TTree *fChain;   //! pointer to the analyzed TTree or TChain
  TTree *fChain2;  //! pointer to the analyzed TTree or TChain
  Int_t fCurrent;  //! current Tree number in a TChain
  Int_t fCurrent2; //! current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  //

   vector<vector<float> > *iEta;
   vector<vector<float> > *iPhi;
   vector<vector <float> >   *Hit_ES_Eta;
   vector<vector <float> >   *Hit_ES_Phi;
   vector<vector <float> >   *Hit_ES_X;
   vector<vector <float> >   *Hit_ES_Y;
   vector<vector <float> >   *Hit_ES_Z;
   vector<vector <float> >   *ES_RecHitEn;
   vector<float>   *recoPhotonIdx;
   vector<vector<float> > *Hit_Eta;
   vector<vector<float> > *Hit_Phi;
   vector<vector<float> > *Hit_X;
   vector<vector<float> > *Hit_Y;
   vector<vector<float> > *Hit_Z;
   vector<vector<float> > *RecHitEn;
   vector<vector<float> > *RecHitFrac;
   vector<vector<int> > *RecHitGain1;
   vector<bool>    *RecHitQuality1;
   vector<vector<float> > *HitNoise;
   vector<bool>    *RecHitFlag_kGood;
   vector<bool>    *RecHitFlag_Endcap;
   vector<bool>    *RecHitFlag_kPoorReco;
   vector<bool>    *RecHitFlag_kOutOfTime;
   vector<bool>    *RecHitFlag_kFaultyHardware;
   vector<bool>    *RecHitFlag_kNoisy;
   vector<bool>    *RecHitFlag_kPoorCalib;
   vector<bool>    *RecHitFlag_kSaturated;
   vector<bool>    *RecHitFlag_kLeadingEdgeRecovered;
   vector<bool>    *RecHitFlag_kNeighboursRecovered;
   vector<bool>    *RecHitFlag_kTowerRecovered;
   vector<bool>    *RecHitFlag_kDead;
   vector<bool>    *RecHitFlag_kKilled;
   vector<bool>    *RecHitFlag_kTPSaturated;
   vector<bool>    *RecHitFlag_kL1SpikeFlag;
   vector<bool>    *RecHitFlag_kWeird;
   vector<bool>    *RecHitFlag_kDiWeird;
   vector<bool>    *RecHitFlag_kHasSwitchToGain6;
   vector<bool>    *RecHitFlag_kHasSwitchToGain1;
   vector<bool>    *RecHitFlag_kESGood;
   vector<bool>    *RecHitFlag_kESDead;
   vector<bool>    *RecHitFlag_kESHot;
   vector<bool>    *RecHitFlag_kESPassBX;
   vector<bool>    *RecHitFlag_kESTwoGoodRatios;
   vector<bool>    *RecHitFlag_kESBadRatioFor12;
   vector<bool>    *RecHitFlag_kESBadRatioFor23Upper;
   vector<bool>    *RecHitFlag_kESBadRatioFor23Lower;
   vector<bool>    *RecHitFlag_kESTS1Largest;
   vector<bool>    *RecHitFlag_kESTS3Largest;
   vector<bool>    *RecHitFlag_kESTS3Negative;
   vector<double>  *DrminLead;
   vector<double>  *DrminSubLead;
   vector<double>  *DrminOther;
   Int_t           nPhotons;
   vector<float>   *pt;
   vector<float>   *eta;
   vector<float>   *phi;
   vector<float>   *energy;
   vector<float>   *energy_ecal_mustache;
   vector<int>     *passMediumId;
   vector<int>     *passTightId;
   vector<int>     *passMVAMediumId;
   vector<float>   *Pho_R9;
   vector<float>   *Pho_S4;
   vector<float>   *Pho_SigIEIE;
   vector<float>   *Pho_SigIPhiIPhi;
   vector<float>   *Pho_SCEtaW;
   vector<float>   *Pho_SCPhiW;
   vector<float>   *Pho_CovIEtaIEta;
   vector<float>   *Pho_CovIEtaIPhi;
   vector<float>   *Pho_ESSigRR;
   vector<float>   *Pho_SCRawE;
   vector<float>   *Pho_SC_ESEnByRawE;
   vector<float>   *Pho_HadOverEm;
   vector<float>   *Pho_Gen_Pt;
   vector<float>   *Pho_Gen_Eta;
   vector<float>   *Pho_Gen_Phi;
   vector<float>   *Pho_Gen_E;
   vector<float>   *A_Gen_Pt;
   vector<float>   *A_Gen_Eta;
   vector<float>   *A_Gen_Phi;
   vector<float>   *A_Gen_E;
   vector<int>     *Pho_Gen_matchedIndex;
   Float_t         rho;
   Int_t           run;
   Int_t           event;
   Int_t           lumi;

   // List of branches
   TBranch        *b_iEta;   //!
   TBranch        *b_iPhi;   //!
   TBranch        *b_recoPhotonIdx;   //!
				      //
   TBranch        *b_Hit_ES_Eta;   //!
   TBranch        *b_Hit_ES_Phi;   //!
   TBranch        *b_Hit_ES_X;   //!
   TBranch        *b_Hit_ES_Y;   //!
   TBranch        *b_Hit_ES_Z;   //!
   TBranch        *b_ES_RecHitEn;   //!
				 //
   TBranch        *b_Hit_Eta;   //!
   TBranch        *b_Hit_Phi;   //!
   TBranch        *b_Hit_X;   //!
   TBranch        *b_Hit_Y;   //!
   TBranch        *b_Hit_Z;   //!
   TBranch        *b_RecHitEn;   //!
   TBranch        *b_RecHitFrac;   //!
   TBranch        *b_RecHitGain1;   //!
   TBranch        *b_RecHitQuality1;   //!
   TBranch        *b_HitNoise;   //!
   TBranch        *b_RecHitFlag_kGood;   //!
   TBranch        *b_RecHitFlag_Endcap;   //!
   TBranch        *b_RecHitFlag_kPoorReco;   //!
   TBranch        *b_RecHitFlag_kOutOfTime;   //!
   TBranch        *b_RecHitFlag_kFaultyHardware;   //!
   TBranch        *b_RecHitFlag_kNoisy;   //!
   TBranch        *b_RecHitFlag_kPoorCalib;   //!
   TBranch        *b_RecHitFlag_kSaturated;   //!
   TBranch        *b_RecHitFlag_kLeadingEdgeRecovered;   //!
   TBranch        *b_RecHitFlag_kNeighboursRecovered;   //!
   TBranch        *b_RecHitFlag_kTowerRecovered;   //!
   TBranch        *b_RecHitFlag_kDead;   //!
   TBranch        *b_RecHitFlag_kKilled;   //!
   TBranch        *b_RecHitFlag_kTPSaturated;   //!
   TBranch        *b_RecHitFlag_kL1SpikeFlag;   //!
   TBranch        *b_RecHitFlag_kWeird;   //!
   TBranch        *b_RecHitFlag_kDiWeird;   //!
   TBranch        *b_RecHitFlag_kHasSwitchToGain6;   //!
   TBranch        *b_RecHitFlag_kHasSwitchToGain1;   //!
   TBranch        *b_RecHitFlag_kESGood;   //!
   TBranch        *b_RecHitFlag_kESDead;   //!
   TBranch        *b_RecHitFlag_kESHot;   //!
   TBranch        *b_RecHitFlag_kESPassBX;   //!
   TBranch        *b_RecHitFlag_kESTwoGoodRatios;   //!
   TBranch        *b_RecHitFlag_kESBadRatioFor12;   //!
   TBranch        *b_RecHitFlag_kESBadRatioFor23Upper;   //!
   TBranch        *b_RecHitFlag_kESBadRatioFor23Lower;   //!
   TBranch        *b_RecHitFlag_kESTS1Largest;   //!
   TBranch        *b_RecHitFlag_kESTS3Largest;   //!
   TBranch        *b_RecHitFlag_kESTS3Negative;   //!
   TBranch        *b_DrminLead;   //!
   TBranch        *b_DrminSubLead;   //!
   TBranch        *b_DrminOther;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_energy_ecal_mustache;   //!
   TBranch        *b_passMediumId;   //!
   TBranch        *b_passTightId;   //!
   TBranch        *b_passMVAMediumId;   //!
   TBranch        *b_Pho_R9;   //!
   TBranch        *b_Pho_S4;   //!
   TBranch        *b_Pho_SigIEIE;   //!
   TBranch        *b_Pho_SigIPhiIPhi;   //!
   TBranch        *b_Pho_SCEtaW;   //!
   TBranch        *b_Pho_SCPhiW;   //!
   TBranch        *b_Pho_CovIEtaIEta;   //!
   TBranch        *b_Pho_CovIEtaIPhi;   //!
   TBranch        *b_Pho_ESSigRR;   //!
   TBranch        *b_Pho_SCRawE;   //!
   TBranch        *b_Pho_SC_ESEnByRawE;   //!
   TBranch        *b_Pho_HadOverEm;   //!
   TBranch        *b_Pho_Gen_Pt;   //!
   TBranch        *b_Pho_Gen_Eta;   //!
   TBranch        *b_Pho_Gen_Phi;   //!
   TBranch        *b_Pho_Gen_E;   //!
   TBranch        *b_A_Gen_Pt;   //!
   TBranch        *b_A_Gen_Eta;   //!
   TBranch        *b_A_Gen_Phi;   //!
   TBranch        *b_A_Gen_E;   //!
   TBranch        *b_Pho_Gen_matchedIndex;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
  // End of modification by Chirayu
};

#endif

#ifdef HGCNtupleVariables_cxx

void HGCNtupleVariables::Init(TTree *tree, TTree *tree2) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  // rechit_detid = 0;
  // rechit_module = 0;
  // rechit_layer = 0;
  // rechit_chip = 0;
  // rechit_channel = 0;
  // rechit_type = 0;
  // rechit_x = 0;
  // rechit_y = 0;
  // rechit_z = 0;
  // rechit_iu = 0;
  // rechit_iv = 0;
  // rechit_energy = 0;
  // rechit_energy_noHG = 0;
  // rechit_amplitudeHigh = 0;
  // rechit_amplitudeLow = 0;
  // rechit_hg_goodFit = 0;
  // rechit_lg_goodFit = 0;
  // rechit_hg_saturated = 0;
  // rechit_lg_saturated = 0;
  // rechit_fully_calibrated = 0;
  // rechit_TS2High = 0;
  // rechit_TS2Low = 0;
  // rechit_TS3High = 0;
  // rechit_TS3Low = 0;
  // rechit_Tot = 0;
  // rechit_time = 0;
  // rechit_timeMaxHG = 0;
  // rechit_timeMaxLG = 0;
  // rechit_toaRise = 0;
  // rechit_toaFall = 0;
  // Set branch addresses and branch pointers
  // Modified By Chirayu
   iEta = 0;
   iPhi = 0;
   recoPhotonIdx = 0;
   Hit_Eta = 0;
   Hit_Phi = 0;
   Hit_X = 0;
   Hit_Y = 0;
   Hit_Z = 0;
   RecHitEn = 0;

   Hit_ES_Eta = 0;
   Hit_ES_Phi = 0;
   Hit_ES_X = 0;
   Hit_ES_Y = 0;
   Hit_ES_Z = 0;
   ES_RecHitEn = 0;

   RecHitFrac = 0;
   RecHitGain1 = 0;
   RecHitQuality1 = 0;
   HitNoise = 0;
   RecHitFlag_kGood = 0;
   RecHitFlag_Endcap = 0;
   RecHitFlag_kPoorReco = 0;
   RecHitFlag_kOutOfTime = 0;
   RecHitFlag_kFaultyHardware = 0;
   RecHitFlag_kNoisy = 0;
   RecHitFlag_kPoorCalib = 0;
   RecHitFlag_kSaturated = 0;
   RecHitFlag_kLeadingEdgeRecovered = 0;
   RecHitFlag_kNeighboursRecovered = 0;
   RecHitFlag_kTowerRecovered = 0;
   RecHitFlag_kDead = 0;
   RecHitFlag_kKilled = 0;
   RecHitFlag_kTPSaturated = 0;
   RecHitFlag_kL1SpikeFlag = 0;
   RecHitFlag_kWeird = 0;
   RecHitFlag_kDiWeird = 0;
   RecHitFlag_kHasSwitchToGain6 = 0;
   RecHitFlag_kHasSwitchToGain1 = 0;
   RecHitFlag_kESGood = 0;
   RecHitFlag_kESDead = 0;
   RecHitFlag_kESHot = 0;
   RecHitFlag_kESPassBX = 0;
   RecHitFlag_kESTwoGoodRatios = 0;
   RecHitFlag_kESBadRatioFor12 = 0;
   RecHitFlag_kESBadRatioFor23Upper = 0;
   RecHitFlag_kESBadRatioFor23Lower = 0;
   RecHitFlag_kESTS1Largest = 0;
   RecHitFlag_kESTS3Largest = 0;
   RecHitFlag_kESTS3Negative = 0;
   DrminLead = 0;
   DrminSubLead = 0;
   DrminOther = 0;
   pt = 0;
   eta = 0;
   phi = 0;
   energy = 0;
   energy_ecal_mustache = 0;
   passMediumId = 0;
   passTightId = 0;
   passMVAMediumId = 0;
   Pho_R9 = 0;
   Pho_S4 = 0;
   Pho_SigIEIE = 0;
   Pho_SigIPhiIPhi = 0;
   Pho_SCEtaW = 0;
   Pho_SCPhiW = 0;
   Pho_CovIEtaIEta = 0;
   Pho_CovIEtaIPhi = 0;
   Pho_ESSigRR = 0;
   Pho_SCRawE = 0;
   Pho_SC_ESEnByRawE = 0;
   Pho_HadOverEm = 0;
   Pho_Gen_Pt = 0;
   Pho_Gen_Eta = 0;
   Pho_Gen_Phi = 0;
   Pho_Gen_E = 0;
   A_Gen_Pt = 0;
   A_Gen_Eta = 0;
   A_Gen_Phi = 0;
   A_Gen_E = 0;
   Pho_Gen_matchedIndex = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("iEta", &iEta, &b_iEta);
   fChain->SetBranchAddress("iPhi", &iPhi, &b_iPhi);
   fChain->SetBranchAddress("recoPhotonIdx", &recoPhotonIdx, &b_recoPhotonIdx);
   fChain->SetBranchAddress("Hit_Eta", &Hit_Eta, &b_Hit_Eta);
   fChain->SetBranchAddress("Hit_Phi", &Hit_Phi, &b_Hit_Phi);
   fChain->SetBranchAddress("Hit_X", &Hit_X, &b_Hit_X);
   fChain->SetBranchAddress("Hit_Y", &Hit_Y, &b_Hit_Y);
   fChain->SetBranchAddress("Hit_Z", &Hit_Z, &b_Hit_Z);
   fChain->SetBranchAddress("RecHitEn", &RecHitEn, &b_RecHitEn);

   fChain->SetBranchAddress("Hit_ES_Eta", &Hit_ES_Eta, &b_Hit_ES_Eta);
   fChain->SetBranchAddress("Hit_ES_Phi", &Hit_ES_Phi, &b_Hit_ES_Phi);
   fChain->SetBranchAddress("Hit_ES_X", &Hit_ES_X, &b_Hit_ES_X);
   fChain->SetBranchAddress("Hit_ES_Y", &Hit_ES_Y, &b_Hit_ES_Y);
   fChain->SetBranchAddress("Hit_ES_Z", &Hit_ES_Z, &b_Hit_ES_Z);
   fChain->SetBranchAddress("ES_RecHitEn", &ES_RecHitEn, &b_ES_RecHitEn);

   fChain->SetBranchAddress("RecHitFrac", &RecHitFrac, &b_RecHitFrac);
   fChain->SetBranchAddress("RecHitGain1", &RecHitGain1, &b_RecHitGain1);
   fChain->SetBranchAddress("RecHitQuality1", &RecHitQuality1, &b_RecHitQuality1);
   fChain->SetBranchAddress("HitNoise", &HitNoise, &b_HitNoise);
   fChain->SetBranchAddress("RecHitFlag_kGood", &RecHitFlag_kGood, &b_RecHitFlag_kGood);
   fChain->SetBranchAddress("RecHitFlag_Endcap", &RecHitFlag_Endcap, &b_RecHitFlag_Endcap);
   fChain->SetBranchAddress("RecHitFlag_kPoorReco", &RecHitFlag_kPoorReco, &b_RecHitFlag_kPoorReco);
   fChain->SetBranchAddress("RecHitFlag_kOutOfTime", &RecHitFlag_kOutOfTime, &b_RecHitFlag_kOutOfTime);
   fChain->SetBranchAddress("RecHitFlag_kFaultyHardware", &RecHitFlag_kFaultyHardware, &b_RecHitFlag_kFaultyHardware);
   fChain->SetBranchAddress("RecHitFlag_kNoisy", &RecHitFlag_kNoisy, &b_RecHitFlag_kNoisy);
   fChain->SetBranchAddress("RecHitFlag_kPoorCalib", &RecHitFlag_kPoorCalib, &b_RecHitFlag_kPoorCalib);
   fChain->SetBranchAddress("RecHitFlag_kSaturated", &RecHitFlag_kSaturated, &b_RecHitFlag_kSaturated);
   fChain->SetBranchAddress("RecHitFlag_kLeadingEdgeRecovered", &RecHitFlag_kLeadingEdgeRecovered, &b_RecHitFlag_kLeadingEdgeRecovered);
   fChain->SetBranchAddress("RecHitFlag_kNeighboursRecovered", &RecHitFlag_kNeighboursRecovered, &b_RecHitFlag_kNeighboursRecovered);
   fChain->SetBranchAddress("RecHitFlag_kTowerRecovered", &RecHitFlag_kTowerRecovered, &b_RecHitFlag_kTowerRecovered);
   fChain->SetBranchAddress("RecHitFlag_kDead", &RecHitFlag_kDead, &b_RecHitFlag_kDead);
   fChain->SetBranchAddress("RecHitFlag_kKilled", &RecHitFlag_kKilled, &b_RecHitFlag_kKilled);
   fChain->SetBranchAddress("RecHitFlag_kTPSaturated", &RecHitFlag_kTPSaturated, &b_RecHitFlag_kTPSaturated);
   fChain->SetBranchAddress("RecHitFlag_kL1SpikeFlag", &RecHitFlag_kL1SpikeFlag, &b_RecHitFlag_kL1SpikeFlag);
   fChain->SetBranchAddress("RecHitFlag_kWeird", &RecHitFlag_kWeird, &b_RecHitFlag_kWeird);
   fChain->SetBranchAddress("RecHitFlag_kDiWeird", &RecHitFlag_kDiWeird, &b_RecHitFlag_kDiWeird);
   fChain->SetBranchAddress("RecHitFlag_kHasSwitchToGain6", &RecHitFlag_kHasSwitchToGain6, &b_RecHitFlag_kHasSwitchToGain6);
   fChain->SetBranchAddress("RecHitFlag_kHasSwitchToGain1", &RecHitFlag_kHasSwitchToGain1, &b_RecHitFlag_kHasSwitchToGain1);
   fChain->SetBranchAddress("RecHitFlag_kESGood", &RecHitFlag_kESGood, &b_RecHitFlag_kESGood);
   fChain->SetBranchAddress("RecHitFlag_kESDead", &RecHitFlag_kESDead, &b_RecHitFlag_kESDead);
   fChain->SetBranchAddress("RecHitFlag_kESHot", &RecHitFlag_kESHot, &b_RecHitFlag_kESHot);
   fChain->SetBranchAddress("RecHitFlag_kESPassBX", &RecHitFlag_kESPassBX, &b_RecHitFlag_kESPassBX);
   fChain->SetBranchAddress("RecHitFlag_kESTwoGoodRatios", &RecHitFlag_kESTwoGoodRatios, &b_RecHitFlag_kESTwoGoodRatios);
   fChain->SetBranchAddress("RecHitFlag_kESBadRatioFor12", &RecHitFlag_kESBadRatioFor12, &b_RecHitFlag_kESBadRatioFor12);
   fChain->SetBranchAddress("RecHitFlag_kESBadRatioFor23Upper", &RecHitFlag_kESBadRatioFor23Upper, &b_RecHitFlag_kESBadRatioFor23Upper);
   fChain->SetBranchAddress("RecHitFlag_kESBadRatioFor23Lower", &RecHitFlag_kESBadRatioFor23Lower, &b_RecHitFlag_kESBadRatioFor23Lower);
   fChain->SetBranchAddress("RecHitFlag_kESTS1Largest", &RecHitFlag_kESTS1Largest, &b_RecHitFlag_kESTS1Largest);
   fChain->SetBranchAddress("RecHitFlag_kESTS3Largest", &RecHitFlag_kESTS3Largest, &b_RecHitFlag_kESTS3Largest);
   fChain->SetBranchAddress("RecHitFlag_kESTS3Negative", &RecHitFlag_kESTS3Negative, &b_RecHitFlag_kESTS3Negative);
   fChain->SetBranchAddress("DrminLead", &DrminLead, &b_DrminLead);
   fChain->SetBranchAddress("DrminSubLead", &DrminSubLead, &b_DrminSubLead);
   fChain->SetBranchAddress("DrminOther", &DrminOther, &b_DrminOther);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPho);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("energy", &energy, &b_energy);
   fChain->SetBranchAddress("energy_ecal_mustache", &energy_ecal_mustache, &b_energy_ecal_mustache);
   fChain->SetBranchAddress("passMediumId", &passMediumId, &b_passMediumId);
   fChain->SetBranchAddress("passTightId", &passTightId, &b_passTightId);
   fChain->SetBranchAddress("passMVAMediumId", &passMVAMediumId, &b_passMVAMediumId);
   fChain->SetBranchAddress("Pho_R9", &Pho_R9, &b_Pho_R9);
   fChain->SetBranchAddress("Pho_S4", &Pho_S4, &b_Pho_S4);
   fChain->SetBranchAddress("Pho_SigIEIE", &Pho_SigIEIE, &b_Pho_SigIEIE);
   fChain->SetBranchAddress("Pho_SigIPhiIPhi", &Pho_SigIPhiIPhi, &b_Pho_SigIPhiIPhi);
   fChain->SetBranchAddress("Pho_SCEtaW", &Pho_SCEtaW, &b_Pho_SCEtaW);
   fChain->SetBranchAddress("Pho_SCPhiW", &Pho_SCPhiW, &b_Pho_SCPhiW);
   fChain->SetBranchAddress("Pho_CovIEtaIEta", &Pho_CovIEtaIEta, &b_Pho_CovIEtaIEta);
   fChain->SetBranchAddress("Pho_CovIEtaIPhi", &Pho_CovIEtaIPhi, &b_Pho_CovIEtaIPhi);
   fChain->SetBranchAddress("Pho_ESSigRR", &Pho_ESSigRR, &b_Pho_ESSigRR);
   fChain->SetBranchAddress("Pho_SCRawE", &Pho_SCRawE, &b_Pho_SCRawE);
   fChain->SetBranchAddress("Pho_SC_ESEnByRawE", &Pho_SC_ESEnByRawE, &b_Pho_SC_ESEnByRawE);
   fChain->SetBranchAddress("Pho_HadOverEm", &Pho_HadOverEm, &b_Pho_HadOverEm);
   fChain->SetBranchAddress("Pho_Gen_Pt", &Pho_Gen_Pt, &b_Pho_Gen_Pt);
   fChain->SetBranchAddress("Pho_Gen_Eta", &Pho_Gen_Eta, &b_Pho_Gen_Eta);
   fChain->SetBranchAddress("Pho_Gen_Phi", &Pho_Gen_Phi, &b_Pho_Gen_Phi);
   fChain->SetBranchAddress("Pho_Gen_E", &Pho_Gen_E, &b_Pho_Gen_E);
   fChain->SetBranchAddress("A_Gen_Pt", &A_Gen_Pt, &b_A_Gen_Pt);
   fChain->SetBranchAddress("A_Gen_Eta", &A_Gen_Eta, &b_A_Gen_Eta);
   fChain->SetBranchAddress("A_Gen_Phi", &A_Gen_Phi, &b_A_Gen_Phi);
   fChain->SetBranchAddress("A_Gen_E", &A_Gen_E, &b_A_Gen_E);
   fChain->SetBranchAddress("Pho_Gen_matchedIndex", &Pho_Gen_matchedIndex, &b_Pho_Gen_matchedIndex);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
  if (!tree)
    return;
  Notify();
  return;
  // End of modification by Chirayu
  //  Set branch addresses and branch pointers
  if (!tree2)
    return;
  fChain2 = tree2;
  fCurrent2 = -1;
  fChain2->SetMakeClass(1);

  Notify();
}

Bool_t HGCNtupleVariables::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef HGCNtupleVariables_cxx
