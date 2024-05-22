#ifndef ANALYZEHGCMuons_H
#define ANALYZEHGCMuons_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HGCNtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TText.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TFrame.h"
#include "TFile.h"
#include "TKey.h"
#include "TLorentzVector.h"
#include "TDirectory.h"

void PaintOverflow(TH1 *h)
{
  // This function paint the histogram h with an extra bin for overflows

  const char *name = h->GetName();
  const char *title = h->GetTitle();
  Int_t nx = h->GetNbinsX() + 1;
  Double_t x1 = h->GetBinLowEdge(1);
  Double_t bw = h->GetBinWidth(nx);
  Double_t x2 = h->GetBinLowEdge(nx) + bw;

  // Book a temporary histogram having ab extra bin for overflows
  TH1F *htmp = new TH1F("rebinned", title, nx, x1, x2);

  // Fill the new hitogram including the extra bin for overflows
  for (Int_t i = 1; i <= nx; i++)
  {
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
  }

  // Fill the underflows
  htmp->Fill(x1 - 1, h->GetBinContent(0));

  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());

  // Draw the temporary histogram
  htmp->Draw();
  TText *t = new TText(x2 - bw / 2, h->GetBinContent(nx), "Overflow");
  t->SetTextAngle(90);
  t->SetTextAlign(12);
  t->SetTextSize(0.03);
  ;
  t->Draw();
}
class AnalyzeHGCMuons : public HGCNtupleVariables
{

public:
  AnalyzeHGCMuons(const TString &inputFileList = "foo.txt", const char *outFileName = "histo.root", const char *dataset = "data", const char *massP="0.1");
  ~AnalyzeHGCMuons();
  Bool_t FillChain(TChain *chain, TChain *chain2, const TString &inputFileList);
  void graphHist(TH2F *hist, char title[], char name[]);
  Long64_t LoadTree(Long64_t entry);
  void EventLoop(const char *);
  void BookHistogram(const char *);
  void graphSide(TH1F *hist1, TH1F *hist2, string name);
  void graphSide(TH2F *hist1, TH2F *hist2, string name);
  void graphOverlay(TH2F *hist1, TH2F *hist2, string name);
  void graphOverlay(TH1F *hist1, TH1F *hist2, string name);
  int fillhist = 1;
  int ctr_twogen = 0;
  int ctr_mcut = 0;
  int ctr_tworeco = 0;
  int ctr_onereco = 0;
  int ctr_onedr = 0;
  int ctr_twodr = 0;
//  string mass = "1.8";
  string gamma = "";
  string mass = "1.8";

  TFile *oFile;
  TDirectory *d_Layer1;

  TTree *treeout;
  double M_m_gen;
  double M_m_reco;
  double M_p_gen;
  double M_pt_gen;
  double M_eta_gen;
  double M_phi_gen;
  vector<float> M_pt_reco;
  vector<float> M_eta_reco;
  vector<float> M_phi_reco;
  vector<float> M_E_reco;

  vector<float> *M_HitX1;
  bool M_EEFlag1;
  bool M_EEFlag2;
  vector<float> *M_HitY1;
  vector<float> *M_HitZ1;
  vector<float> *M_HitEta1;
  vector<float> *M_HitPhi1;

  vector<float> *M_HitiEta1;
  vector<float> *M_HitiPhi1;
  vector<float> *M_HitEn1;
  vector<float> *M_HitET1;
  vector<float> *M_HitEZ1;
  vector<float> *M_HitX2;
  vector<float> *M_HitY2;
  vector<float> *M_HitZ2;
  vector<float> *M_HitEta2;
  vector<float> *M_HitPhi2;
  vector<float> *M_HitiEta2;
  vector<float> *M_HitiPhi2;
  vector<float> *M_HitEn2;
  vector<float> *M_HitET2;
  vector<float> *M_HitEZ2;
  vector<float> M_SigIEIE1;
  vector<float> M_SigIEIE2;
  vector<float> M_SigIPhiIPhi1;
  vector<float> M_SigIPhiIPhi2;
  vector<float> M_R9_Pho1;
  vector<float> M_R9_Pho2;
  
  vector<float> *M_ES_HitX1;
  vector<float> *M_ES_HitY1;
  vector<float> *M_ES_HitZ1;
  vector<float> *M_ES_HitEta1;
  vector<float> *M_ES_HitPhi1;
  vector<float> *M_ES_HitEn1;
  vector<float> *M_ES_HitX2;
  vector<float> *M_ES_HitY2;
  vector<float> *M_ES_HitZ2;
  vector<float> *M_ES_HitEta2;
  vector<float> *M_ES_HitPhi2;
  vector<float> *M_ES_HitEn2;

  double M_EBM;
  double M_angle_gen;
  int M_event;
  int M_lumi;
  int M_run;
  //  float m_gen;
  TH1F *P_gen_PhoLead_noDR;
  //  TH2F *Pt_gen_LeadSublead;
  TH1F *Pt_gen_PhoLead_noDR;
  TH1F *Pt_gen_PhoSublead_noDR;
  TH1F *Pt_gen_photons;
  TH2F *Pt_gen_LeadSublead_noDR;
  TH2F *M_Pt_gen;
  TH2F *Eta_Phi_gen;
  TH2F *Pt_Eta_gen;
  TH1F *Pt_gen;
  TH1F *Eta_gen;
  TH1F *Phi_gen;
  TH1F *M_gen;
  TH1F *Pt_gen_diphoton;
  TH1F *Eta_gen_diphoton;
  TH1F *Phi_gen_diphoton;
  TH1F *M_gen_diphoton;
  //  TH1F *Pt_reco_PhoLead ;
  //  TH1F *Pt_reco_PhoSublead ;
  //  TH2F *Pt_reco_LeadSublead;
  //  TH1F *mass_gen;
  TH1F *mass_gen_noDR;
  //  TH1F *mass_reco;
  TH1F *nGen;
  //  TH1F *angle_gen;
  TH1F *angle_gen_noDR;
  TH1F *mod_target;
  TH1F *nreco;
  TH1F *nRechits;
  TH1F *nRechits_ecut;
  TH1F *mod_mass;
  TH1F *mod_mass2;
  TH1F *mod_angle;
  TH1F *mod_angle2;
  TH1F *gamma_L;
  TH1F *gamma_L_noDR;

  TH1F *hit_eta;
  TH1F *hit_phi;
  TH1F *hit_E;
  TH1F *hit_logE;
  TH1F *total_rechitE;
  TH1F *total_rechitE_cuts;

  TH2F *hit_etaphiE;
  TH2F *hit_etaphi;

  TH1F *hitL_X;
  TH1F *hitL_Y;
  TH1F *hitL_Z;
  TH1F *hitL_E;
  TH1F *hitL_E2;
  TH1F *hitS_X;
  TH1F *hitS_Y;
  TH1F *hitS_Z;
  TH1F *hitS_E;
  TH1F *hitS_E2;


  //  TH1F *angle_reco;
  // TH2F  *deta_dphi_gen ;
  TH2F *deta_dphi_gen_noDR;
  TH2F  *deta_dphi_reco ;
  TH2F  *eta_phi_gen ;
  TH2F *eta_phi_gen_noDR;
  TH2F  *eta_phi_reco ;
  TH1F *h_ADChg[128];

  TH2F *P_Pt_gen_PhoLead;
  TH2F *P_eta_gen_PhoLead;
  TH2F *Pt_eta_gen_PhoLead;
  TH2F *P_Pt_reco_PhoLead;
  TH2F *P_eta_reco_PhoLead;
  TH2F *Pt_eta_reco_PhoLead;

  TH1F *P_gen_PhoLead;
  TH1F *Pt_gen_PhoLead;
  TH1F *Pt_gen_PhoSublead ;
  TH1F *Pt_reco_PhoSublead ;
  TH1F *Eta_gen_PhoLead;
  TH1F *P_reco_PhoLead;
  TH1F *Pt_reco_PhoLead;
  TH1F *Eta_reco_PhoLead;

  TH1F *Pt_genuponreco_PhoLead;
  TH1F *Pt_genuponreco_PhoSublead;
  TH1F *Pt_genuponreco_Total;

  TH1F *nHitsPho1;
  TH1F *nHitsPho2;

  TH1F *nHitsPhoT;
  
  TH1F *drminPhoLead;
  TH1F *drminPhoSubLead;
  TH1F *drminPhoOther;
  
  TH1F *E_gen;
  TH1F *Et_gen;
};
#endif

#ifdef ANALYZEHGCMuons_cxx

void AnalyzeHGCMuons::BookHistogram(const char *outFileName)
{

  //  char hname[200], htit[200];
  //  double xlow = 0.0,  xhigh = 2000.0;
  //  int nbins = 2000;

  oFile = new TFile(outFileName, "recreate");
  oFile->cd();
  if (fillhist)
  {
    /* 	  d_Layer1 = oFile->mkdir("Validation Plots");
        char hname[200];
        d_Layer1-> cd(); */
    //	 eta_phi_gen= new TH2F("Eta_Phi_gen" , " eta  Vs. phi of the two electron", 276,-2.4,2.4,369,-3.2,3.2);
    // eta_phi_gen= new TH2F("Eta_Phi_gen" , " #eta  Vs. #phi of the two photons in a single a->#gamma #gamma decay at ; #eta ; #phi", 32,-16,16,32,-16,16);
    E_gen = new TH1F("E_gen","gen energy",100,0,100);
    Et_gen = new TH1F("Et_gen","gen transverse energy",100,0,100);
    eta_phi_gen_noDR = new TH2F("Eta_Phi_gen_noDR", " #eta  Vs. #phi of the two gen photons without dR matching; #eta ; #phi", 32, -16, 16, 32, -16, 16);

    eta_phi_reco= new TH2F("Eta_Phi_reco" , " #eta  Vs. #phi of the two photons ; #eta ; #phi", 32,-16,16,32,-16,16);
    //  deta_dphi_gen= new TH2F("dEta_dPhi_gen" , " Fraction of events in #Delta #eta vs #Delta #phi between 2 photons ;#Delta #eta;#Delta #phi;", 50,0,50,50,0,50);
    deta_dphi_gen_noDR = new TH2F("dEta_dPhi_gen_noDR", "Fraction of events in #Delta #eta vs #Delta #phi between 2 gen level photons without dR matching;#Delta #eta;#Delta #phi;", 10, 0, 10, 10, 0, 10);
    deta_dphi_reco= new TH2F("dEta_dPhi_reco" , " Fraction of events in #Delta #eta vs #Delta #phi between 2 photons ;#Delta #eta;#Delta #phi;", 50,0,50,50,0,50);
    P_gen_PhoLead_noDR = new TH1F("P_gen_PhoLead_noDR", "P of Lead Gen Photon without dR matching;P", 130, 0, 130);
    Pt_gen_PhoLead_noDR = new TH1F("Pt_gen_PhoLead_noDR", "Pt of Lead Gen Photon without dR matching;Pt", 100, 0, 100);
    //  Pt_reco_PhoLead = new TH1F(("Pt_reco_PhoLead"+string("_Gamma")+gamma+"_M"+mass).c_str(),("Frequency distribution of Pt of Lead reco Photon ;Pt"),100,0,100);
    //  Pt_gen_PhoSublead = new TH1F(("Pt_gen_PhoSublead"+string("_Gamma")+gamma+"_M"+m2ass).c_str(),("Frequency distribution of Pt of SubLead Gen Photon ;Pt"),100,0,100);
    Pt_gen_PhoSublead_noDR = new TH1F("Pt_gen_PhoSublead_noDR", "Pt of SubLead Gen Photon without dR matching;Pt", 100, 0, 100);
    //  Pt_reco_PhoSublead = new TH1F(("Pt_reco_PhoSublead"+string("_Gamma")+gamma+"_M"+mass).c_str(),("Frequency distribution of Pt of SubLead reco Photon ;Pt"),100,0,100);
    //  Pt_gen_LeadSublead = new TH2F(("Pt_gen_LeadSublead"+string("_Gamma")+gamma+"_M"+mass).c_str(),("Pt of Lead vs Pt of Sublead Photon;Pt(Lead) ;Pt(Sublead Photon)"),100,0,100,100,0,100);
    Pt_gen_LeadSublead_noDR = new TH2F("Pt_gen_LeadSublead_noDR", "Pt of Lead vs Pt of Sublead Gen Photons without dR matching;Pt(Lead) ;Pt(Sublead Photon)", 100, 0, 100, 100, 0, 100);
    //  Pt_reco_LeadSublead = new TH2F(("Pt_reco_LeadSublead"+string("_Gamma")+gamma+"_M"+mass).c_str(),("Pt of Lead vs Pt of Sublead Photon;Pt(Lead) ;Pt(Sublead Photon)"),100,0,100,100,0,100);

    //  mass_gen = new TH1F("mass_gen","Invariant Mass of Gen DiPhotons",100,-1,2);
    mass_gen_noDR = new TH1F("mass_gen_noDR", "Invariant Mass of Gen DiPhotons without dR matching;mass", 100, 0.5, 0.7);
    //  mass_reco = new TH1F("mass_reco","Invariant Mass of reco DiPhotons",100,-1,2);

    angle_gen_noDR = new TH1F("angle_gen_noDR", "Angle between 2 gen photons without dR matching;#theta(Radians)", 50, 0, 0.2);
    //  angle_gen = new TH1F("angle_gen","Angle between 2 photons",50,0,0.2);
    //  angle_reco = new TH1F("angle_reco","Angle between 2 photons",50,0,0.2);

    P_Pt_gen_PhoLead = new TH2F("P_Pt_gen_PhoLead", "P vs Pt of dR matched lead gen Photons;p;pt", 130, 0, 130, 100, 0, 100);
    P_eta_gen_PhoLead = new TH2F("P_eta_gen_PhoLead", "P vs Eta of dR matched lead gen Photons;p;#eta", 130, 0, 130, 100, -1.5, 1.5);
    Pt_eta_gen_PhoLead = new TH2F("Pt_eta_gen_PhoLead", "Pt vs Eta of dR matched lead gen Photons;Pt;#eta", 100, 0, 100, 100, -1.5, 1.5);
    P_Pt_reco_PhoLead = new TH2F("P_Pt_reco_PhoLead", "P vs Pt of dR matched lead reco Photons;p;pt", 130, 0, 130, 100, 0, 100);
    P_eta_reco_PhoLead = new TH2F("P_eta_reco_PhoLead", "P vs Eta of dR matched lead reco Photons;p;#eta", 130, 0, 130, 100, -1.5, 1.5);
    Pt_eta_reco_PhoLead = new TH2F("Pt_eta_reco_PhoLead", "Pt vs Eta of dR matched lead reco Photons;Pt;#eta", 100, 0, 100, 100, -1.5, 1.5);

    Pt_gen_PhoLead = new TH1F("Pt_gen_PhoLead", "Pt of dR matched lead gen photons;Pt", 100, 0, 100);
    Pt_gen_PhoSublead = new TH1F("Pt_gen_PhoSublead", "Pt of dR matched sublead gen photons;Pt", 100, 0, 100);
    Pt_genuponreco_PhoLead = new TH1F("Pt_genuponreco_PhoLead", "Pt gen/Pt reco of dR matched lead photons;Pt gen/Pt reco", 100, 0, 2.5);
    Pt_genuponreco_PhoSublead = new TH1F("Pt_genuponreco_PhoSublead", "Pt gen/Pt reco of dR sublead matched photons;Pt gen/Pt reco", 100, 0, 2.5);
    Pt_genuponreco_Total = new TH1F("Pt_genuponreco_Total", "Total Pt gen/Total Pt reco of dR matched photons;Pt gen/Pt reco", 100, 0, 2.5);
    P_gen_PhoLead = new TH1F("P_gen_PhoLead", "P of dR matched lead gen photons;P", 130, 0, 130);
    Eta_gen_PhoLead = new TH1F("Eta_gen_PhoLead", "Eta of dR matched lead gen photons;#eta", 100, -1.5, 1.5);

    Pt_reco_PhoLead = new TH1F("Pt_reco_PhoLead", "Pt of dR matched lead reco photons;Pt", 100, 0, 100);
    Pt_reco_PhoSublead = new TH1F("Pt_reco_PhoSublead", "Pt of dR matched sublead reco photons;Pt", 100, 0, 100);
    P_reco_PhoLead = new TH1F("P_reco_PhoLead", "P of dR matched reco photons;P", 130, 0, 130);
    Eta_reco_PhoLead = new TH1F("Eta_reco_PhoLead", "Eta of dR matched reco photons;#eta", 100, -1.5, 1.5);
    nHitsPho1 = new TH1F((string(mass)+"nHitsPho1"s).c_str(),"Number of Rechits of Photon 1",150,0,150);
    nHitsPho2 = new TH1F((string(mass)+"nHitsPho2"s).c_str(),"Number of Rechits of Photon 2",120,0,120);
    nHitsPhoT = new TH1F((string(mass)+"nHitsPhoT"s).c_str(),"Total Number of rechits",150,0,150);

    mod_target = new TH1F("mod_target","",100,0,0.7);
    mod_mass = new TH1F("mod_mass","",100,0,1);
    mod_mass2 = new TH1F("mod_mass2","",100,0.4,2);
    mod_angle = new TH1F("mod_angle","",50,0,0.2);
    mod_angle2 = new TH1F("mod_angle2","",50,0,0.2);

    nreco = new TH1F("nReco","Number of reco photons without selections",4,0,4);
    nGen = new TH1F("nGen","Number of Gen Photons without selections",4,0,4);

    hit_eta = new TH1F("hit_eta","Eta;Eta",100,-1.5,1.5);
    hit_phi = new TH1F("hit_phi","phi;phi",100,-3.15,3.15);
    hit_E = new TH1F("hit_E","E;E",100,0,100);
    hit_logE = new TH1F("hit_logE","E;E",100,0,100);
    hit_etaphi = new TH2F("hit_etaphi","Eta vs Phi;Eta;Phi",100,-1.5,1.5,100,-3.15,3.15);
    hit_etaphiE = new TH2F("hit_etaphiE","Eta vs Phi weighted with energyj;Eta;Phi",100,-1.5,1.5,100,-3.15,3.15);

    drminPhoLead = new TH1F("drmin_lead","drmin of Lead Gen Photon(underflow indicates not matched);drmin",200,0,0.2);
    drminPhoSubLead = new TH1F("drmin_sublead","drmin of Sublead Gen Photon(underflow indicates not matched);drmin",200,0,0.2);
    drminPhoOther = new TH1F("drmin_other","drmin of non lead non sublead gen photons;drmin",200,0,0.2);

    hitL_X = new TH1F("hitL_X","Rechit distribution of X for training of lead photon",100,-200,200);
    hitL_Y = new TH1F("hitL_Y","Rechit distribution of Y for training of lead photon",100,-200,200);
    hitL_Z = new TH1F("hitL_Z","Rechit distribution of Z for training of lead photon",100,-330,330);
    hitL_E = new TH1F("hitL_E","Rechit distribution of E for training of lead photon",400,0,1);
    hitL_E2 = new TH1F("hitL_E2","Rechit distribution of E for training of lead photon",200,0,100);
    hitS_X = new TH1F("hitS_X","Rechit distribution of X for training of Sublead photon",100,-200,200);
    hitS_Y = new TH1F("hitS_Y","Rechit distribution of Y for training of Sublead photon",100,-200,200);
    hitS_Z = new TH1F("hitS_Z","Rechit distribution of Z for training of Sublead photon",100,-330,330);
    hitS_E = new TH1F("hitS_E","Rechit distribution of E for training of Sublead photon",400,0,1);
    hitS_E2 = new TH1F("hitS_E2","Rechit distribution of E for training of Sublead photon",200,0,100);

    nRechits = new TH1F("nRechits","Number of rechits before Ecut",100,0,100);
    nRechits_ecut = new TH1F("nRechits_ecut","Number of rechits after Ecut",100,0,100);

    total_rechitE = new TH1F("total_rechitE","Total rechit energy before ecut",100,0,100);
    total_rechitE_cuts = new TH1F("total_rechitE_cuts","Total rechit energy after ecut",100,0,100);
    
    gamma_L = new TH1F("gamma_L","Longitudinal boost",100,0,500);
    gamma_L_noDR = new TH1F("gamma_L_noDR","Longitudinal boost before DR matching",100,0,500);

    Pt_gen_photons = new TH1F("Pt_gen_photons", "Pt of dR matched gen photons;Pt", 100, 0, 100);

    M_Pt_gen = new TH2F("M_Pt_gen", ";True M_{a};True p_{T,a}",100,0.5,2.1, 100, 10, 110);
    Pt_Eta_gen = new TH2F("Pt_Eta_gen", "Pt vs Eta of A;Pt;Eta",100,10,110, 100, -3, 3);

    Eta_Phi_gen = new TH2F("Eta_Phi_gen",";True #eta_{a};True #phi_{a}",100,-1.6,1.6,100, -3.25,3.25);

    M_gen = new TH1F("M_gen","Mass of gen A;M",100,0,2.1);
    Pt_gen = new TH1F("Pt_gen", "Pt of gen A;Pt", 100, 10, 110);
    Eta_gen = new TH1F("Eta_gen", "Eta of gen A; Eta",100,-3.1,3.1);
    Phi_gen = new TH1F("Phi_gen", "Phi of gen A; Phi",100, -3.16,3.16);
    M_gen_diphoton = new TH1F("M_gen_diphoton","Mass of gen diphotons;M",100,0,2.1);
    Pt_gen_diphoton = new TH1F("Pt_gen_diphoton", "Pt of gen diphotons;Pt", 100, 10, 110);
    Eta_gen_diphoton = new TH1F("Eta_gen_diphoton", "Eta of gen diphotons; Eta",100,-3.1,3.1);
    Phi_gen_diphoton = new TH1F("Phi_gen_diphoton", "Phi of gen diphotons; Phi",100, -3.16,3.16);
  }
  //  for(int ii=0; ii<128; ii++) {
  //    sprintf(hname, "ADC_HG_%d", ii+1);
  //    h_ADChg[ii] = new TH1F(hname, hname, 100, 0, 400);
  //  }
}

AnalyzeHGCMuons::AnalyzeHGCMuons(const TString &inputFileList, const char *outFileName, const char *dataset, const char *massP)
{

  TChain *tree = new TChain("nTuplelize/T");
  TChain *tree2;
  mass = string(massP);

  if (!FillChain(tree, tree2, inputFileList))
  {
    std::cerr << "Cannot get the tree " << std::endl;
  }
  else
  {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  HGCNtupleVariables::Init(tree, tree2);

  BookHistogram(outFileName);
//cout<<massP<<endl;
  treeout = new TTree("fordrn", "ForDRN");
  treeout->Branch("Hit_X_Pho1", &M_HitX1);
  treeout->Branch("EEFlag1", &M_EEFlag1);
  treeout->Branch("EEFlag2", &M_EEFlag2);
  treeout->Branch("Hit_Y_Pho1", &M_HitY1);
  treeout->Branch("Hit_Z_Pho1", &M_HitZ1);
  treeout->Branch("Hit_Eta_Pho1", &M_HitEta1);
  treeout->Branch("Hit_Phi_Pho1", &M_HitPhi1);
  treeout->Branch("Hit_iEta_Pho1", &M_HitiEta1);
  treeout->Branch("Hit_iPhi_Pho1", &M_HitiPhi1);
  treeout->Branch("Hit_X_Pho2", &M_HitX2);
  treeout->Branch("Hit_Y_Pho2", &M_HitY2);
  treeout->Branch("Hit_Z_Pho2", &M_HitZ2);
  treeout->Branch("Hit_Eta_Pho2", &M_HitEta2);
  treeout->Branch("Hit_Phi_Pho2", &M_HitPhi2);
  treeout->Branch("Hit_iEta_Pho2", &M_HitiEta2);
  treeout->Branch("Hit_iPhi_Pho2", &M_HitiPhi2);
  treeout->Branch("RecHitEnPho1", &M_HitEn1);
  treeout->Branch("RecHitETPho1", &M_HitET1);
  treeout->Branch("RecHitEZPho1", &M_HitEZ1);
  treeout->Branch("RecHitEnPho2", &M_HitEn2);
  treeout->Branch("RecHitETPho2", &M_HitET2);
  treeout->Branch("RecHitEZPho2", &M_HitEZ2);
  treeout->Branch("m_gen", &M_m_gen);
  treeout->Branch("p_gen", &M_p_gen);
  treeout->Branch("pt_gen", &M_pt_gen);
  treeout->Branch("eta_gen", &M_eta_gen);
  treeout->Branch("phi_gen", &M_phi_gen);
  treeout->Branch("SigIEIEPho1",&M_SigIEIE1);
  treeout->Branch("SigIEIEPho2",&M_SigIEIE2);
  treeout->Branch("R9_Pho1",&M_R9_Pho1);
  treeout->Branch("R9_Pho2",&M_R9_Pho2);
  treeout->Branch("SigIPhiIPhiPho1",&M_SigIPhiIPhi1);
  treeout->Branch("SigIPhiIPhiPho2",&M_SigIPhiIPhi2);

  treeout->Branch("Hit_ES_X_Pho1", &M_ES_HitX1);
  treeout->Branch("Hit_ES_Y_Pho1", &M_ES_HitY1);
  treeout->Branch("Hit_ES_Z_Pho1", &M_ES_HitZ1);
  treeout->Branch("Hit_ES_Eta_Pho1", &M_ES_HitEta1);
  treeout->Branch("Hit_ES_Phi_Pho1", &M_ES_HitPhi1);
  treeout->Branch("ES_RecHitEnPho1", &M_ES_HitEn1);
  treeout->Branch("Hit_ES_X_Pho2", &M_ES_HitX2);
  treeout->Branch("Hit_ES_Y_Pho2", &M_ES_HitY2);
  treeout->Branch("Hit_ES_Z_Pho2", &M_ES_HitZ2);
  treeout->Branch("Hit_ES_Eta_Pho2", &M_ES_HitEta2);
  treeout->Branch("Hit_ES_Phi_Pho2", &M_ES_HitPhi2);
  treeout->Branch("ES_RecHitEnPho2", &M_ES_HitEn2);
}

Bool_t AnalyzeHGCMuons::FillChain(TChain *chain, TChain *chain2, const TString &inputFileList)
{

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if (!infile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while (1)
  {
    infile >> buffer;
    if (!infile.good())
      break;
    // std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    chain->Add(buffer.c_str());
    // chain2->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in chain  : " << chain->GetEntries() << std::endl;
  // std::cout << "No. of Entries in chain2 : " << chain2->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t AnalyzeHGCMuons::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (!fChain->InheritsFrom(TChain::Class()))
    return centry;
  TChain *chain = (TChain *)fChain;
  if (chain->GetTreeNumber() != fCurrent)
  {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  // Modified By Chirayu
  return centry;
  // End Modified

  if (!fChain2)
    return -5;
  Long64_t centry2 = fChain2->LoadTree(entry);
  if (centry2 < 0)
    return centry2;
  if (!fChain2->InheritsFrom(TChain::Class()))
    return centry2;
  TChain *chain2 = (TChain *)fChain2;
  if (chain2->GetTreeNumber() != fCurrent)
  {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }

  if (centry == centry2)
    return centry;
  else
    return -1;
}
void AnalyzeHGCMuons::graphSide(TH1F *hist1, TH1F *hist2, string name)
{
  TCanvas *c1 = new TCanvas(name.c_str(), "", 200, 10, 700, 500);
  c1->Divide(2, 1);
  c1->cd(1);
  hist1->Draw("colz");
  c1->cd(2);
  hist2->Draw("colz");

  c1->Update();
  c1->GetFrame()->SetBorderSize(12);
  c1->Modified();
  c1->Write();
  c1->SaveAs((name + string(".png")).c_str());
}
void AnalyzeHGCMuons::graphSide(TH2F *hist1, TH2F *hist2, string name)
{
  TCanvas *c1 = new TCanvas(name.c_str(), "", 200, 10, 700, 500);
  c1->Divide(2, 1);
  c1->cd(1);
  hist1->Draw("colz");
  c1->cd(2);
  hist2->Draw("colz");

  c1->Update();
  c1->GetFrame()->SetBorderSize(12);
  c1->Modified();
  c1->Write();
  c1->SaveAs((name + string(".png")).c_str());
}
void AnalyzeHGCMuons::graphOverlay(TH1F *hist1, TH1F *hist2, string name)
{
  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
  c3->cd();
  THStack *hist3 = new THStack(name.c_str(), "");
  hist1->SetLineColor(kRed);
  hist2->SetLineColor(kBlue);
  hist3->Add(hist1);
  hist3->Add(hist2);
  hist3->Draw("NOSTACKB");
  hist3->Write();
  gPad->BuildLegend(0.75, 0.75, 0.95, 0.95, "");
  c3->Update();
  c3->Modified();
  c3->Write();
  c3->SaveAs((name + string(".png")).c_str());
}
void AnalyzeHGCMuons::graphOverlay(TH2F *hist1, TH2F *hist2, string name)
{
  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
  c3->cd();
  THStack *hist3 = new THStack(name.c_str(), "");
  hist1->SetLineColor(kRed);
  hist2->SetLineColor(kBlue);
  hist3->Add(hist1);
  hist3->Add(hist2);
  hist3->Draw("NOSTACKB");
  hist3->Write();
  gPad->BuildLegend(0.75, 0.75, 0.95, 0.95, "");
  c3->Update();
  c3->Modified();
  c3->Write();
  c3->SaveAs((name + string(".png")).c_str());
}
void AnalyzeHGCMuons::graphHist(TH2F *hist, char title[], char name[])
{

  TCanvas *c1 = new TCanvas(name, title, 200, 10, 700, 500);
  TGraph *gr1 = new TGraph();
  for (int i = 0; i < hist->GetNbinsX(); i++)
  {
    for (int j = 0; j < hist->GetNbinsY(); j++)
    {
      double x = hist->GetXaxis()->GetBinCenter(i);
      double y = hist->GetYaxis()->GetBinCenter(j);
      gr1->SetPoint(gr1->GetN(), x, y);
    }
  }
  gr1->SetLineColor(2);
  gr1->SetLineWidth(4);
  gr1->SetMarkerStyle(0);
  gr1->SetTitle(title);
  gr1->GetXaxis()->SetTitle("Eta");
  gr1->GetYaxis()->SetTitle("Pt");
  gr1->Draw("ACP");
  gr1->Write(name);

  // TCanvas::Update() draws the frame, after which one can change it
  c1->Update();
  c1->GetFrame()->SetBorderSize(12);
  c1->Modified();
}
void DivideHistogramZValues(TH2F *hist, Double_t constant)
{
  Int_t nBinsX = hist->GetNbinsX();
  Int_t nBinsY = hist->GetNbinsY();

  // Iterate over each bin and divide the z-value by the constant
  for (Int_t i = 1; i <= nBinsX; i++)
  {
    for (Int_t j = 1; j <= nBinsY; j++)
    {
      Double_t binContent = hist->GetBinContent(i, j);
      Double_t newContent = binContent / constant;
      hist->SetBinContent(i, j, newContent);
    }
  }
}
AnalyzeHGCMuons::~AnalyzeHGCMuons()
{

  cout<<"Std of bin content in Pt vs Eta: "<<endl;
  cout<<"Pt: "<<Pt_Eta_gen->GetRMS(1)<<endl;
  cout<<"Eta: "<<Pt_Eta_gen->GetRMS(2)<<endl;
  cout<<"Std of bin content in M vs Pt: "<<endl;
  cout<<"M: "<<M_Pt_gen->GetRMS(1)<<endl;
  cout<<"Pt: "<<M_Pt_gen->GetRMS(2)<<endl;
  fillhist = 1;
  if (fillhist)
  {

    TH1F *plotsToPrintTH1F[] = {
//        mass_gen_noDR,
//        P_gen_PhoLead_noDR,
//        Pt_gen_PhoLead_noDR,
//        Pt_gen_PhoSublead_noDR,
	angle_gen_noDR,
//        P_gen_PhoLead,
//        Pt_gen_PhoLead,
//        Eta_gen_PhoLead,
//        P_reco_PhoLead,
//        Pt_reco_PhoLead,
//        Eta_reco_PhoLead,
//	nHitsPho1,
//	nHitsPho2,
//	nHitsPhoT,
//	mod_target,
//	mod_mass,
//	mod_mass2,
//	mod_angle,
//	mod_angle2,
//	nreco,
  hit_eta,
  hit_phi,
  hit_E,
  hit_logE,
  drminPhoLead,
  drminPhoSubLead,
  drminPhoOther,
  Pt_gen_PhoLead,
  Pt_gen_PhoSublead,
  Pt_genuponreco_PhoLead,
  Pt_genuponreco_PhoSublead,
  Pt_genuponreco_Total,
  Pt_reco_PhoLead,
  Pt_reco_PhoSublead,
  hitL_X,
  hitL_Y,
  hitL_Z,
  hitS_X,
  hitS_Y,
  hitS_Z,
  nRechits,
  nRechits_ecut,
  total_rechitE,
  total_rechitE_cuts,
  gamma_L,
  gamma_L_noDR,
  nGen,
  mass_gen_noDR,
  nreco,
  hitL_E,
  hitS_E,
  hitL_E2,
  hitS_E2,
  E_gen,
  Et_gen,
    };
    //TH1F *plotsToPrintTH1Flog[] = {};
    TH1F *plotsToPrintTH1Flog[] = {nGen,nreco,hitL_E,hitS_E,hitL_E2,hitS_E2 };

    Int_t b_max = angle_gen_noDR->GetMaximumBin();
   Double_t x_max = angle_gen_noDR->GetBinCenter(b_max);
    cout<<"ANGLE IS: "<<x_max<<endl;
    TH1F *plotsSideBySideTH1F[][2] = {};
    TH1F *plotsOverlayTH1F[][2] = {};
    // auto plotsToPrintTH2F = {Pt_gen_LeadSublead,deta_dphi_gen,deta_dphi_reco,eta_phi_gen,eta_phi_reco,Pt_reco_LeadSublead,Pt_gen_LeadSublead_noDR,deta_dphi_gen_noDR,eta_phi_gen_noDR};
    //TH2F *plotsToPrintTH2F[] = {  hit_etaphi,hit_etaphiE,eta_phi_gen_noDR, deta_dphi_gen_noDR, Pt_gen_LeadSublead_noDR, P_eta_gen_PhoLead, P_Pt_gen_PhoLead, Pt_eta_gen_PhoLead, P_eta_reco_PhoLead, P_Pt_reco_PhoLead, Pt_eta_reco_PhoLead,  eta_phi_reco, deta_dphi_reco};
	TH2F *plotsToPrintTH2F[] = {deta_dphi_gen_noDR, Pt_Eta_gen, M_Pt_gen , Eta_Phi_gen} ;
    // TH2F* plotsSideBySideTH2F[][2] = {{deta_dphi_gen,eta_phi_gen}};
    TH2F *plotsSideBySideTH2F[][2] = {};
    // TH2F* plotsOverlayTH2F[][2] = {{deta_dphi_gen,eta_phi_gen}};
    TH2F *plotsOverlayTH2F[][2] = {};
    // DivideHistogramZValues(deta_dphi_gen,deta_dphi_gen->GetEntries());
    DivideHistogramZValues(deta_dphi_gen_noDR, deta_dphi_gen_noDR->GetEntries());
    // DivideHistogramZValues(deta_dphi_reco,deta_dphi_reco->GetEntries());
    // deta_dphi_gen->SetMinimum(1e-6);
    TStyle *st1 = new TStyle("st1", "my style");
    st1->cd();
    st1->SetPaintTextFormat(".2f");
    st1->SetOptStat(111111);
    st1->SetFrameFillColor(10);
    st1->SetFillColor(10);
    st1->SetCanvasColor(10);

    for (auto i : plotsToPrintTH2F)
    {
      TCanvas *c3 = new TCanvas("c3", "c3", 1600, 1400);
      c3->SetLeftMargin(0.125);
      c3->cd();
      i->Draw("colz");
      i->Draw("colz");
      i->GetXaxis()->SetTitleSize(0.045);
      i->GetYaxis()->SetTitleSize(0.045);
      i->SetStats(0);
      c3->Update();
      c3->Modified();
      c3->Write(i->GetName());
      c3->SaveAs((string(i->GetName()) + string(".pdf")).c_str());
      c3->SaveAs((string(i->GetName()) + string(".png")).c_str());
    }
    for (auto i : plotsSideBySideTH2F)
    {
      auto p1 = i[0];
      auto p2 = i[1];
      string name = string(p1->GetName()) + "," + string(p2->GetName());
      graphSide(p1, p2, name);
    }
    for (auto i : plotsOverlayTH2F)
    {
      auto p1 = i[0];
      auto p2 = i[1];
      string name = string(p1->GetName()) + "+" + string(p2->GetName());
      graphOverlay(p1, p2, name);
    }
    for (auto i : plotsToPrintTH1F)
    {
      TCanvas *c3 = new TCanvas("c3", "c3", 1600, 1400);
      c3->cd();
  //c3->SetLogy();
      i->SetTitle((i->GetTitle()+string(" mass(MeV)=")+mass).c_str());
      cout<<i->GetTitle()<<endl;
      i->Draw("colz");
      c3->Update();
      c3->Modified();
      c3->Write(i->GetName());
      c3->SaveAs((string("0")+mass+string(i->GetName()) + string(".pdf")).c_str());
      c3->SaveAs((string("0")+mass+string(i->GetName()) + string(".png")).c_str());
    }
    for (auto i : plotsToPrintTH1Flog)
    {
      TCanvas *c3 = new TCanvas("c3", "c3", 1600, 1400);
      c3->cd();
  c3->SetLogy();
      i->SetTitle((i->GetTitle()+string(" mass(MeV)=")+mass).c_str());
      i->Draw("colz");
      c3->Update();
      c3->Modified();
      c3->Write(i->GetName());
      c3->SaveAs((string("0log")+mass+string(i->GetName()) + string(".pdf")).c_str());
      c3->SaveAs((string("0log")+mass+string(i->GetName()) + string(".png")).c_str());
    }
    for (auto i : plotsSideBySideTH1F)
    {
      auto p1 = i[0];
      auto p2 = i[1];
      string name = string(p1->GetName()) + "," + string(p2->GetName());
      graphSide(p1, p2, name);
    }
    for (auto i : plotsOverlayTH1F)
    {
      auto p1 = i[0];
      auto p2 = i[1];
      string name = string(p1->GetName()) + "+" + string(p2->GetName());
      graphOverlay(p1, p2, name);
    }
    cout << "Number of events with 2 gen Photons:" << ctr_twogen << endl;
    cout << "Number of events with 1 reco Photons:" << ctr_onereco << endl;
    cout << "Number of events with 2 reco Photons:" << ctr_tworeco << endl;
    cout << "Number of events with mass of A < 1 GeV:" << ctr_mcut << endl;
    // graphHist(Eta_Pt_gen,"Eta vs Pt Gen level","eta_pt_gen");
    // graphHist(angle_p_gen,"Angle(e- e+) vs P(Z) Gen level","angle_p_gen");
    // graphHist(angle_pt_gen,"Angle(e- e+) vs Pt(Z) Gen level","angle_pt_gen");
  }
  // if (!fChain || !fChain2) return;
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
  // delete fChain2->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  // fChain->CloneTree(-1,"");
  oFile->Close();
  // TFile myFile( "file.root", "RECREATE" );
  // cout<<"Writing Tree";
  // myFile.Close();
}

#endif
// TODO: Pt condition in dR matching (later), Plots vs mass using jupyter notebook, Publish , Write explainations( tomorrow)
