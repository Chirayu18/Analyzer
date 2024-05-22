#define ANALYZEHGCMuons_cxx

#include "AnalyzeHGCMuons.h"
#include "TLorentzVector.h"
#include <cstring>
#include <iostream>
#include <vector>

using namespace std;
const double drthreshold = 0.05;
int ctr = 0;

int main(int argc, char *argv[]) {

  if (argc < 2) {
    cerr << "Please give 3 arguments "
         << "runList "
         << " "
         << "outputFileName"
         << " "
         << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName = argv[2];
  const char *data = argv[3];
  const char *massP = argv[4];

  cout<<massP<<endl;
  AnalyzeHGCMuons hgcmuons(inputFileList, outFileName, data,massP);
  // cout << "dataset " << data << " " << endl;

  hgcmuons.EventLoop(data);

  return 0;
}

double DeltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  while (result > M_PI)
    result -= 2 * M_PI;
  while (result <= -M_PI)
    result += 2 * M_PI;
  return result;
}

double DeltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta * deta + dphi * dphi);
}
void AnalyzeHGCMuons::EventLoop(const char *data) {
  if (fChain == 0)
    return;

  TLorentzVector p1gen, p2gen, pAgen;
  TLorentzVector p1reco, p2reco,pTreco;
  Long64_t nentries = fChain->GetEntriesFast();
  // cout << "nentries " << nentries << endl;
  // cout << "Analyzing dataset " << data << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  int decade = 0;
  double genMass = 0;
  double genAngle = 0;
  int wentry=0;

  for (Long64_t jentry = 0; jentry < nentries; jentry++) {

    
  M_HitX1=0;
  M_HitY1=0;
  M_HitZ1=0;
  M_HitEn1=0;
  M_HitET1=0;
  M_HitEZ1=0;
  M_HitEta1=0;
  M_HitEta2=0;
  M_HitiEta1=0;
  M_HitiEta2=0;
  M_HitPhi1=0;
  M_HitPhi2=0;
  M_HitiPhi1=0;
  M_HitiPhi2=0;
  M_HitX2=0;
  M_HitY2=0;
  M_HitZ2=0;
  M_EEFlag1=0;
  M_EEFlag2=0;
  M_HitEn2=0;
  M_HitET2=0;
  M_HitEZ2=0;
  M_SigIEIE1.clear();
  M_SigIEIE2.clear();
  M_SigIPhiIPhi1.clear();
  M_SigIPhiIPhi2.clear();
  M_R9_Pho1.clear();
  M_R9_Pho2.clear();
  M_m_gen=-100;
  M_p_gen=-100;
  M_pt_gen=-100;
  M_eta_gen=-100;
  M_phi_gen=-100;

  M_ES_HitX1=0;
  M_ES_HitY1=0;
  M_ES_HitZ1=0;
  M_ES_HitEn1=0;
  M_ES_HitEta1=0;
  M_ES_HitPhi1=0;
  M_ES_HitX2=0;
  M_ES_HitY2=0;
  M_ES_HitZ2=0;
  M_ES_HitEn2=0;
  M_ES_HitEta2=0;
  M_ES_HitPhi2=0;
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int(progress);
    if (k > decade)
      // cout << 10 * k << " %" << endl;
      decade = k;

    // ===============read this entry == == == == == == == == == == ==
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    // Store if given reco electron is good
    int kgoodgen = 0;
    fillhist = 0;
    
	if(Pho_Gen_Pt->size() <2)
	  {
		continue;
	  }
        ctr_twogen++;
        p1gen.SetPtEtaPhiE(Pho_Gen_Pt->at(0), Pho_Gen_Eta->at(0),
                         Pho_Gen_Phi->at(0), Pho_Gen_E->at(0));
	p2gen.SetPtEtaPhiE(Pho_Gen_Pt->at(1), Pho_Gen_Eta->at(1),
	Pho_Gen_Phi->at(1), Pho_Gen_E->at(1));
        pAgen.SetPtEtaPhiE(A_Gen_Pt->at(0), A_Gen_Eta->at(0),
                         A_Gen_Phi->at(0), A_Gen_E->at(0));
        //genMass = (p1gen + p2gen).M();
	genMass = pAgen.M();
	/**
	if(genMass < 0.6 )
	{
		continue;
	}
	**/
	ctr_mcut++;
	nGen->Fill(Pho_Gen_Pt->size());
	nreco->Fill(pt->size());

  // if((p1gen+p2gen).M()<0.55 || (p1gen+p2gen).M() > 0.65 )
  //   {
  //     continue;
  //   }
    if(DrminLead->size()>0)
    {
	 drminPhoLead->Fill(DrminLead->at(0));
    }
    if(DrminSubLead->size()>0)
      {
    drminPhoSubLead->Fill(DrminSubLead->at(0));
    }
    if(DrminOther->size()>0)
    {
    drminPhoOther->Fill(DrminOther->at(0));
    }

      // Sort in E
      //using namespace std;
      //for (int i = 2; i < Pho_Gen_E->size(); i++) {
      //  if (Pho_Gen_Pt->at(i) > p1gen.Pt()) {
      //    p2gen = p1gen;
      //    p1gen.SetPtEtaPhiE(Pho_Gen_Pt->at(i), Pho_Gen_Eta->at(i),
      //                       Pho_Gen_Phi->at(i), Pho_Gen_E->at(i));
      //  }
      //}
      //if (p1gen.Pt() < p2gen.Pt()) {
      //  auto swap = p1gen;
      //  p1gen = p2gen;
      //  p2gen = swap;
      //}
      //DO dr matching
      int matchedLead=0;
      int matchedSublead=0;
//	vector<float>* hitX[3] = {Hit_X_Pho1,Hit_X_Pho2,Hit_X_Pho3};
//	vector<float>* hitY[3] = {Hit_Y_Pho1,Hit_Y_Pho2,Hit_Y_Pho3};
//	vector<float>* hitZ[3] = {Hit_Z_Pho1,Hit_Z_Pho2,Hit_Z_Pho3};
//	vector<float>* hitE[3] = {RecHitEnPho1,RecHitEnPho2,RecHitEnPho3};
  TLorentzVector pgen=p1gen+p2gen;
  double Ez = sqrt(pAgen.E()*pAgen.E() - pAgen.Et()*pAgen.Et());
  gamma_L_noDR->Fill(Ez/pgen.M());
	int nLead=0,necutLead=0;
	int nSublead=0,necutSublead=0;
	double etotal=0;
	double etotal_cuts=0;
    vector<float> ET2{};
    vector<float> EZ2{};
    vector<float> ET1{};
    vector<float> EZ1{};
	if(DrminSubLead->size()>0 && DrminSubLead->at(0)>0 && DrminSubLead->at(0)<0.1)
	{
        	int ridx = Pho_Gen_matchedIndex->at(1);
		p2reco.SetPtEtaPhiE(pt->at(ridx),eta->at(ridx), phi->at(ridx), energy->at(ridx));
		Pt_gen_PhoSublead->Fill(p2gen.Pt());
		Pt_reco_PhoSublead->Fill(p2reco.Pt());
		Pt_gen_photons->Fill(p2gen.Pt());
		Pt_genuponreco_PhoSublead->Fill(p2gen.Pt()/p2reco.Pt());
		matchedSublead=1;
		M_HitX2 = &Hit_X->at(ridx);
		M_HitY2 = &Hit_Y->at(ridx);
		M_HitZ2 = &Hit_Z->at(ridx);
    M_HitEta2 = &Hit_Eta->at(ridx);
    M_HitPhi2 = &Hit_Phi->at(ridx);
    M_HitiEta2 = &iEta->at(ridx);
    M_HitiPhi2 = &iPhi->at(ridx);
		M_HitEn2 = &RecHitEn->at(ridx);
    M_EEFlag2 = RecHitFlag_Endcap->at(ridx) ;

		M_ES_HitX2 = &Hit_ES_X->at(ridx);
		M_ES_HitY2 = &Hit_ES_Y->at(ridx);
		M_ES_HitZ2 = &Hit_ES_Z->at(ridx);
    M_ES_HitEta2 = &Hit_ES_Eta->at(ridx);
    M_ES_HitPhi2 = &Hit_ES_Phi->at(ridx);
		M_ES_HitEn2 = &ES_RecHitEn->at(ridx);
      M_SigIEIE2.push_back(Pho_SigIEIE->at(ridx));
      M_SigIPhiIPhi2.push_back(Pho_SigIPhiIPhi->at(ridx));
      M_R9_Pho2.push_back(Pho_R9->at(ridx));
    
		for(int i = 0 ; i < M_HitX2->size();i++)
		{
			//cout<<"I is "<<i<<endl;
			hitS_X->Fill(M_HitX2->at(i));
			hitS_Y->Fill(M_HitY2->at(i));

			hitS_Z->Fill(M_HitZ2->at(i));
			hitS_E->Fill(M_HitEn2->at(i));
			hitS_E2->Fill(M_HitEn2->at(i));
      ET2.push_back(M_HitEn2->at(i)*1/cosh(M_HitEta2->at(i)));
      EZ2.push_back(M_HitEn2->at(i)*tanh(M_HitEta2->at(i)));
			nSublead++;
			etotal+=M_HitEn2->at(i);
			if(M_HitEn2->at(i)>0.1)
			{
				necutSublead++;
				etotal_cuts+=M_HitEn2->at(i);
			}
		}
		M_HitET2 = &ET2;
		M_HitEZ2 = &EZ2;

	}
    if (DrminLead->size()>0 && DrminLead->at(0)>0 && DrminLead->at(0)<0.1) {
    int ridx = Pho_Gen_matchedIndex->at(0);
    p1reco.SetPtEtaPhiE(pt->at(ridx),eta->at(ridx), phi->at(ridx), energy->at(ridx));
    Pt_gen_PhoLead->Fill(p1gen.Pt());
    Pt_reco_PhoLead->Fill(p1reco.Pt());
    Pt_genuponreco_PhoLead->Fill(p1gen.Pt()/p1reco.Pt());
	Pt_gen_photons->Fill(p1gen.Pt());
    matchedLead=1;
		M_HitX1 = &Hit_X->at(ridx);
		M_HitY1 = &Hit_Y->at(ridx);
		M_HitZ1 = &Hit_Z->at(ridx);
    M_HitEta1 = &Hit_Eta->at(ridx);
    M_HitPhi1 = &Hit_Phi->at(ridx);
    M_HitiEta1 = &iEta->at(ridx);
    M_HitiPhi1 = &iPhi->at(ridx);
		M_HitEn1 = &RecHitEn->at(ridx);
    M_EEFlag1 = RecHitFlag_Endcap->at(ridx) ;

		M_ES_HitX1 = &Hit_ES_X->at(ridx);
		M_ES_HitY1 = &Hit_ES_Y->at(ridx);
		M_ES_HitZ1 = &Hit_ES_Z->at(ridx);
    M_ES_HitEta1 = &Hit_ES_Eta->at(ridx);
    M_ES_HitPhi1 = &Hit_ES_Phi->at(ridx);
		M_ES_HitEn1 = &ES_RecHitEn->at(ridx);
      M_SigIEIE1.push_back(Pho_SigIEIE->at(ridx));
      M_SigIPhiIPhi1.push_back(Pho_SigIPhiIPhi->at(ridx));
      M_R9_Pho1.push_back(Pho_R9->at(ridx));
		for(int i = 0 ;  i < M_HitX1->size();i++)
		{
			hitL_X->Fill(M_HitX1->at(i));
			hitL_Y->Fill(M_HitY1->at(i));
			hitL_Z->Fill(M_HitZ1->at(i));
			hitL_E->Fill(M_HitEn1->at(i));
			hitL_E2->Fill(M_HitEn1->at(i));
      ET1.push_back(M_HitEn1->at(i)*1/cosh(M_HitEta1->at(i)));
      EZ1.push_back(M_HitEn1->at(i)*tanh(M_HitEta1->at(i)));
			nLead++;
			etotal+=M_HitEn1->at(i);
			if(M_HitEn1->at(i)>0.1)
			{
				necutLead++;
				etotal_cuts+=M_HitEn1->at(i);
			}
		}
        M_HitET1 = &ET1;
        M_HitEZ1 = &EZ1;
    }
    int nofRechits=0;
    if(matchedLead)
    {
      if(M_HitX1->size() == M_HitY1->size() && M_HitX1->size() == M_HitZ1->size()
    && M_HitX1->size() == M_HitEn1->size() )
        {
          nRechits->Fill(M_HitX1->size());
        }
        else{
          throw 205;
        }
    } 

        genAngle = p1gen.Angle(p2gen.Vect());
        angle_gen_noDR->Fill(genAngle);
	if( matchedLead && matchedSublead)
	{
		ctr_tworeco++;
	}
    if(matchedLead || matchedSublead)
    {
	    ctr_onereco++;
  E_gen->Fill(pAgen.E());
  Et_gen->Fill(pAgen.Et());

  Pt_gen->Fill(pAgen.Pt());
  M_gen->Fill(pAgen.M());
  Eta_gen->Fill(pAgen.Eta());
  Phi_gen->Fill(pAgen.Phi());

  Pt_gen_diphoton->Fill(pgen.Pt());
  M_gen_diphoton->Fill(pgen.M());
  Eta_gen_diphoton->Fill(pgen.Eta());
  Phi_gen_diphoton->Fill(pgen.Phi());

  M_Pt_gen->Fill(pAgen.M(),pAgen.Pt());
  Pt_Eta_gen->Fill(pAgen.Pt(),pAgen.Eta());
  Eta_Phi_gen->Fill(pAgen.Eta(),pAgen.Phi());
      nRechits->Fill(nofRechits);
	    if(matchedLead && matchedSublead)
	    {
		    pTreco = p1reco+p2reco;
	    }
	    else if(matchedLead)
	    {
		    pTreco = p1reco;
	    }
	    else
	    {
		    pTreco = p2reco;
	    }
	//nRechits->Fill(nLead+nSublead);
	total_rechitE->Fill(etotal);
	total_rechitE_cuts->Fill(etotal_cuts);
	nRechits_ecut->Fill(necutSublead+necutLead);
	Pt_genuponreco_Total->Fill((p1gen+p2gen).Pt()/pTreco.Pt());
  TLorentzVector pgen=p1gen+p2gen;
  double Ez = sqrt(pgen.E()*pgen.E() - pgen.Et()*pgen.Et());
  gamma_L->Fill(Ez/pgen.M());
	
        // cout << "Photon 1(Eta,Phi,E): " << p1gen.Eta() / 0.0174 << "," <<
        // p1gen.Phi() / 0.0174 << "," << p1gen.E() << endl; cout << "Photon
        // 2(Eta,Phi,E): " << p2gen.Eta() / 0.0174 << "," << p2gen.Phi() /
        // 0.0174 << "," << p2gen.E() << endl; cout << endl;

        mass_gen_noDR->Fill(genMass);
        //cout<<"SEGMENTATION ERROR"<<endl;
        M_m_gen=genMass;
	//M_p_gen=(p1gen+p2gen).P();
	M_p_gen = pAgen.P();
	M_pt_gen = pAgen.Pt();
	M_eta_gen = pAgen.Eta();
	M_phi_gen = pAgen.Phi();
        treeout->Fill();
        float deta_gen = p1gen.Eta() - p2gen.Eta();
        float dphi_gen = DeltaPhi(p1gen.Phi(), p2gen.Phi());
        //float E = (p1gen + p2gen).E();
	float E = pAgen.E();
        eta_phi_gen_noDR->Fill(p1gen.Eta() / 0.0174, p1gen.Phi() / 0.0174,
                               p1gen.E());
        eta_phi_gen_noDR->Fill(p2gen.Eta() / 0.0174, p2gen.Phi() / 0.0174,
                               p2gen.E());
//        deta_dphi_gen_noDR->Fill(sqrt(deta_gen*deta_gen + dphi_gen*dphi_gen) / 0.0174,pgen.Gamma());
        deta_dphi_gen_noDR->Fill(deta_gen/0.0174,dphi_gen/0.0174);
        Pt_gen_LeadSublead_noDR->Fill(p1gen.Pt(), p2gen.Pt());
        P_gen_PhoLead_noDR->Fill(p1gen.P());
//	cout<<Hit_X_Pho1->size()<<endl;
//	nHitsPho1->Fill(Hit_X_Pho1->size());
//	nHitsPho2->Fill(Hit_X_Pho2->size());
//
//	nHitsPhoT->Fill(Hit_X_Pho1->size()+Hit_X_Pho2->size());
//      float totalrechiten=0;
//      for(int i =0;i<RecHitEnPho1->size();i++)
//      {
//	      totalrechiten+=RecHitEnPho1->at(i);
//      hit_eta->Fill(Hit_Eta_Pho1->at(i));
//      hit_E->Fill(RecHitEnPho1->at(i));
//      hit_logE->Fill(RecHitEnPho1->at(i));
//      hit_phi->Fill(Hit_Phi_Pho1->at(i));
//      hit_etaphi->Fill(Hit_Eta_Pho1->at(i),Hit_Phi_Pho1->at(i));
//      hit_etaphiE->Fill(Hit_Eta_Pho1->at(i),Hit_Phi_Pho1->at(i),RecHitEnPho1->at(i));
//
//      }
//      for(int i =0;i<RecHitEnPho2->size();i++)
//      {
//	     totalrechiten+=RecHitEnPho2->at(i);
//      hit_eta->Fill(Hit_Eta_Pho2->at(i));
//      hit_E->Fill(RecHitEnPho2->at(i));
//      hit_logE->Fill(RecHitEnPho2->at(i));
//      hit_phi->Fill(Hit_Phi_Pho2->at(i));
//      hit_etaphi->Fill(Hit_Eta_Pho2->at(i),Hit_Phi_Pho2->at(i));
//      hit_etaphiE->Fill(Hit_Eta_Pho2->at(i),Hit_Phi_Pho2->at(i),RecHitEnPho2->at(i));
//      }
//      
//
//	      
//
//      
//	double modT = sqrt((1-cos(genAngle))/(1+cos(genAngle)));
//	double modM = genAngle*(p1gen.P()+p2gen.P())/2;
//	double modM2 = genAngle*totalrechiten/(2);
//      cout<<"Total Rechit Energy:"<<totalrechiten<<endl;
//      cout<<"Corresponding Angle:"<<genAngle<<endl;
//      cout<<"LHS:"<<cos(genAngle)<<endl;
//      cout<<"RHS:"<<(totalrechiten*totalrechiten - 2*genMass*genMass)/(totalrechiten*totalrechiten + 2*genMass*genMass)<<endl;
//      cout<<"Mass:"<<totalrechiten*totalrechiten/2*(1-cos(genAngle))/(cos(genAngle)+1);
//      mod_target->Fill(modT);
//      mod_mass->Fill(modM);
//      mod_mass2->Fill(modM2);
//      mod_angle->Fill(2*genMass/totalrechiten);
//      mod_angle2->Fill(2*genMass/(p1gen.E()+p2gen.E()));
//      int kmatchedPho1 = 0, kmatchedPho2 = 0;
//      double drminPho1 = drthreshold;
//      double drminPho2 = drthreshold;
//

 //     if (energy->size() >= 1) {
 //       ctr_onereco++;
 //       p1reco.SetPtEtaPhiE(pt->at(0), eta->at(0), phi->at(0), energy->at(0));
 //       if (energy->size() >= 2) {
 //         ctr_tworeco++;
 //         p2reco.SetPtEtaPhiE(pt->at(1), eta->at(1), phi->at(1), energy->at(1));
 //         for (int i = 2; i < energy->size(); i++) {
 //           if (energy->at(i) > p1reco.Pt()) {
 //             p2reco = p1reco;
 //             p1reco.SetPtEtaPhiE(pt->at(i), eta->at(i), phi->at(i),
 //                                 energy->at(i));
 //           }
 //         }
 //         if (p1reco.Pt() < p2reco.Pt()) {
 //           auto swap = p1reco;
 //           p1reco = p2reco;
 //           p2reco = swap;
 //         }
 //         for (int i = 0; i < Pho_Gen_E->size(); i++) {
 //           double drPho1 = DeltaR(p1reco.Eta(), p1reco.Phi(),
 //                                  Pho_Gen_Eta->at(i), Pho_Gen_Phi->at(i));
 //           double drPho2 = DeltaR(p2reco.Eta(), p2reco.Phi(),
 //                                  Pho_Gen_Eta->at(i), Pho_Gen_Phi->at(i));

 //           // If both reco Photons have same dr with the reco Photon then
 //           // not matched
 //           if (drPho1 < drminPho1 && drPho1 <= drPho2) {
 //             kmatchedPho1 = 1;

 //             p1gen.SetPtEtaPhiE(Pho_Gen_Pt->at(i), Pho_Gen_Eta->at(i),
 //                                Pho_Gen_Phi->at(i), Pho_Gen_E->at(i));
 //             drminPho1 = drPho1;
 //           }
 //           if (drPho2 < drminPho2 && drPho2 < drPho1) { kmatchedPho2 = 1;
 //             p2gen.SetPtEtaPhiE(Pho_Gen_Pt->at(i), Pho_Gen_Eta->at(i),
 //                                Pho_Gen_Phi->at(i), Pho_Gen_E->at(i));
 //             drminPho2 = drPho2;
 //           }
 //         }
 //       } else {
 //         for (int i = 0; i < Pho_Gen_E->size(); i++) {
 //           double drPho1 = DeltaR(p1reco.Eta(), p1reco.Phi(),
 //                                  Pho_Gen_Eta->at(i), Pho_Gen_Phi->at(i));
 //           // If both reco Phoctrons have same dr with the reco Phoctron then
 //           // not matched
 //           if (drPho1 < drminPho1) {
 //             kmatchedPho1 = 1;

 //             p1gen.SetPtEtaPhiE(Pho_Gen_Pt->at(i), Pho_Gen_Eta->at(i),
 //                                Pho_Gen_Phi->at(i), Pho_Gen_E->at(i));
 //             drminPho1 = drPho1;
 //           }
 //         }
 //       }
 //       using namespace std;
 //     }
 //     if (kmatchedPho1) {
 //       ctr_onedr++;
 //     }
 //     if (kmatchedPho2) {
 //       ctr_twodr++;
 //     }
 //     if (fillhist && kmatchedPho1) {
 //       // jcout << "Photon 2(Eta,Phi,E): " << p2reco.Eta() / 0.0174 << "," <<
 //       // p2reco.Phi() / 0.0174 << "," << p2reco.E() << endl;
 //       /*         if (0)
 //               {

 //                 genMass = (p1gen + p2gen).M();
 //                 mass_gen->Fill(genMass);
 //                 double recoMass = (p1reco + p2reco).M();
 //                 mass_reco->Fill(recoMass);

 //                 double angle = p1gen.Angle(p2gen.Vect());
 //                 angle_gen->Fill(angle);
 //                 float deta_gen = p1gen.Eta() - p2gen.Eta();
 //                 float dphi_gen = DeltaPhi(p1gen.Phi(), p2gen.Phi());
 //                 float E = (p1gen + p2gen).E();
 //                 eta_phi_gen->Fill(p1gen.Eta() / 0.0174, p1gen.Phi() / 0.0174,
 //          p1gen.E()); eta_phi_gen->Fill(p2gen.Eta() / 0.0174, p2gen.Phi() /
 //          0.0174, p2gen.E()); deta_dphi_gen->Fill(abs(deta_gen) / 0.0174,
 //          abs(dphi_gen) / 0.0174); Pt_gen_PhoSublead->Fill(p2gen.Pt());
 //                 Pt_gen_LeadSublead->Fill(p1gen.Pt(), p2gen.Pt());

 //                 angle = p1reco.Angle(p2reco.Vect());
 //                 angle_reco->Fill(angle);
 //                 eta_phi_reco->Fill(p1reco.Eta() / 0.0174, p1reco.Phi() /
 //          0.0174, p1reco.E()); eta_phi_reco->Fill(p2reco.Eta() / 0.0174,
 //          p2reco.Phi() / 0.0174, p2reco.E());
 //                 Pt_reco_PhoSublead->Fill(p2reco.Pt());
 //                 Pt_reco_LeadSublead->Fill(p1reco.Pt(), p2reco.Pt());
 //                 float deta_reco = p1reco.Eta() - p2reco.Eta();
 //                 float dphi_reco = DeltaPhi(p1reco.Phi(), p2reco.Phi());
 //                 deta_dphi_reco->Fill(abs(deta_reco) / 0.0174, abs(dphi_reco) /
 //          0.0174);
 //               } */

 //       P_Pt_gen_PhoLead->Fill(p1gen.P(), p1gen.Pt());
 //       Pt_eta_gen_PhoLead->Fill(p1gen.Pt(), p1gen.Eta());
 //       P_eta_gen_PhoLead->Fill(p1gen.P(), p1gen.Eta());
 //       P_Pt_gen_PhoLead->Fill(p1gen.P(), p1gen.Pt());

 //       P_Pt_reco_PhoLead->Fill(p1reco.P(), p1reco.Pt());
 //       Pt_eta_reco_PhoLead->Fill(p1reco.Pt(), p1reco.Eta());
 //       P_eta_reco_PhoLead->Fill(p1reco.P(), p1reco.Eta());
 //       P_Pt_reco_PhoLead->Fill(p1reco.P(), p1reco.Pt());

//       P_gen_PhoLead->Fill(p1gen.P());
 //      Pt_gen_PhoLead->Fill(p1gen.Pt());
 //      Eta_gen_PhoLead->Fill(p1gen.Eta());
//
 //      P_reco_PhoLead->Fill(p1reco.P());
 //      Pt_reco_PhoLead->Fill(p1reco.Pt());
//       Eta_reco_PhoLead->Fill(p1reco.Eta());
       
 //       eta_phi_reco->Fill(p1reco.Eta(),p2reco.Phi());
 //       deta_dphi_reco->Fill(p2reco.Eta()-p1reco.Eta(),DeltaPhi(p2reco.Phi(),p1reco.Phi()));
 //     }
 //     // To quickly replace Gen_Eta->at(0) with p0.Eta() : ❯
 //     // '<,'>s/Gen_\(.[a-z]*\)->at(\([0-1]\))/p\2gen.\1()/g
 //     kgoodgen = 1;
 //   }

 //   /*     std::cout << " Run " << run << "  Event " << event
 //                 << " beamEnergy "
 //                 << lumi
 //                 //<< " tree2 ntracks " << ntracks << " impactX_HGCal_layer_1 "
 //      <<
 //                 // impactX_HGCal_layer_1
 //                 << std::endl; */
 // }
}
}
}
