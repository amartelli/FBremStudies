// g++ -Wall -o trisCheck_fBrem `root-config --cflags --glibs` Utils/setTDRStyle.cc Utils/ntupleUtils.cc Utils/stabilityUtils.cc Utils/ConvoluteTemplate.cc Utils/histoFunc.h Utils/TPileupReweighting.h  Utils/FindTrackRegion.h trisCheck_fBrem.cpp


 #include "../Utils/PUReweighting.h"

 #include "TROOT.h"
 #include "TStyle.h"
 #include "TFile.h"
 #include "TF1.h"
 #include "TH1.h"
 #include "TH2.h"
 #include "TProfile.h"
 #include "TProfile2D.h"
 #include "TCanvas.h"
 #include "TGraphErrors.h"
 #include "TGraphAsymmErrors.h"
 #include "TPaveStats.h"
 #include "TLegend.h"
 #include "TChain.h"
 #include "TRandom.h"
 #include "TVirtualFitter.h"
 #include "TLatex.h"
 #include "TLine.h"
 #include "TMath.h"

 #include <iostream>
 #include <iomanip>
 #include <string>
 #include <sstream>
 #include <ctime>
 #include <map>
 #include <algorithm>
 #include <math.h>
 #include <vector>



//void trisCheck_fBrem_DP(std::string& type, std::string& period = &("ABCD"), std::string& format = &("AOD"))
void trisCheck_fBrem_DP(std::string& type)
 {
   //-----------------
   // Input parameters

   std::cout << "\n*******************************************************************************************************************" << std::endl;

   std::string SCTK = "tk";
//   std::string SCTK = "sc";
   //std::string chargeChosen = "Pos";
   //std::string chargeChosen = "Neg";
   std::string chargeChosen = "All";

   std::string zEleSelected = "Pos";
   //std::string zEleSelected = "Neg";
   //std::string zEleSelected = "All";

   std::string period = "ABCD";
   std::string format = "AOD";

   std::cout << "type:         " << type         << std::endl;
   std::cout << "period:       " << period       << std::endl;
   std::cout << "format:       " << format       << std::endl;
   std::cout << "SCTK:         " << SCTK       << std::endl;

   // Get trees                                                          
   std::cout << std::endl;
   std::string nameFile;
   TChain* ntu;


   std::string inputFilesMC;
   inputFilesMC = "/gwterax2/users/martelli/CALIBRATION/MC/DYToEE_M20_powheg-Summer12-START53-ZSkim-runDependent-allRange.root";
   TChain* ntu_MC = new TChain("selected");
   ntu_MC->Add(inputFilesMC.c_str());

   if(type == "DA") {
     ntu = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
     ntu->Add("../DATAFolder_October2013/DoubleElectron_Run2012A-22Jan2013-v1_AOD.root");
     ntu->Add("../DATAFolder_October2013/DoubleElectron_Run2012B-22Jan2013-v1_AOD.root");
     ntu->Add("../DATAFolder_October2013/DoubleElectron_Run2012C-22Jan2013-v1_AOD.root");
     ntu->Add("../DATAFolder_October2013/DoubleElectron_Run2012D-22Jan2013-v1_AOD.root");
   }
   if(type == "F10"){
     ntu = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
     ntu->Add("../DATAFolder_October2013/DYToEE_M_20_TuneZ2star_8TeV_pythia6_Summer12_DR53X-ExtFlat10_PU_S10_START53_V15-v1_GEN-SIM-RECODEBUG.root");
   }
   if(type == "F20"){
     ntu = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
     ntu->Add("../DATAFolder_October2013/DYToEE_M_20_TuneZ2star_8TeV_pythia6_Summer12_DR53X-ExtFlat20_PU_S10_START53_V15-v1_GEN-SIM-RECODEBUG.root");
   }
   if(type == "F10S30"){
     ntu = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
     ntu->Add("../DATAFolder_October2013/DY53X-ExtFlat10S30_PU_S10_START53_V15-v1_ValidTrkHits.root");
   }
   if(type == "F20S30"){
     ntu = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
     ntu->Add("../DATAFolder_October2013/DY53X-ExtFlat20S30_PU_S10_START53_V15-v1_ValidTrkHits.root");
   }
   if(type == "F10RD"){
     ntu = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
     ntu->Add("../DATAFolder_October2013/DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-ExtFlat10_PU_RD1_START53_V7N-v1_AODSIM.root");
   }
   if(type == "F20RD"){
     ntu = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
     ntu->Add("../DATAFolder_October2013/DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-ExtFlat20_PU_RD1_START53_V7N-v1_AODSIM.root");
   }
   if(type == "MCRD"){
     ntu = new TChain("selected");
     ntu->Add("../DATAFolder_October2013/DYToEE_M20_powheg-Summer12-START53-ZSkim-runDependent-194533-194533.root");
     ntu->Add("../DATAFolder_October2013/DYToEE_M20_powheg-Summer12-START53-ZSkim-runDependent-200519-200519.root");
     ntu->Add("../DATAFolder_October2013/DYToEE_M20_powheg-Summer12-START53-ZSkim-runDependent-206859-206859.root");
   }
   if(type == "MC"){
     ntu = new TChain("simpleNtupleEoverP/SimpleNtupleEoverP");
     ntu->Add("../DATAFolder_October2013/DYToEE_M_20_TuneZ2star_8TeV_pythia6_Summer12_DR53X-NominalGeo_PU_S10_START53_V15-v1_GEN-SIM-RECODEBUG.root");
   }

   std::cout << "     type : " << type << std::setw(8) << ntu->GetEntries() << " entries" << std::endl;

   if(ntu->GetEntries() == 0)
   {
     std::cout << "Error: At least one file is empty" << std::endl; 
     return ;
   }

   std::cout << "prima di PU " << std::endl; 
   std::string PUdataFileName = "../MCWeights_MC_53X_MBstudies/pileup_69p3mb_true_Winter2013.root";

   std::map<float,float>* puReweighting = new std::map<float,float>; 
   if(type != "DA" && type != "MCRD") puReweighting = ComputePUweights(ntu, PUdataFileName, false);

   std::cout << "PU preso" << std::endl; 

   /////////////// Set branch addresses                                        
   int isW, isZ;
   int nVtx;
   float PV_z;
   float npu;
   //ele1                                                                      
   float ele1_dxy_PV;
   float ele1_dz_PV;
   float ele1_charge;
   float ele1_scEtaWidth;
   float ele1_scPhiWidth;
   float ele1_eta;
   float ele1_phi;
   float ele1_fbrem;
   float ele1_tkPt;
   float ele1_tkP;
   float ele1_scE;
   float ele1_scERaw;
   float ele1_e3x3;
   float ele1_R9;
   float ele1_scEt;
   float ele1_scEta;
   float ele1_scPhi;


   //ele2                                                                      
   float ele2_dxy_PV;
   float ele2_dz_PV;
   float  ele2_charge;
   float ele2_scEtaWidth;
   float ele2_scPhiWidth;
   float ele2_eta;
   float ele2_phi;
   float ele2_fbrem;
   float ele2_tkPt;
   float ele2_tkP;
   float ele2_scE;
   float ele2_scERaw;
   float ele2_e3x3;
   float ele2_R9;
   float ele2_scEt;
   float ele2_scEta;
   float ele2_scPhi;


   //////////////// MC_G4  || G4_X0 || Data
   //   if(type != "MCRD"){
   ntu->SetBranchStatus("*",0);
   ntu->SetBranchStatus("isW", 1);                           ntu->SetBranchAddress("isW", &isW);
   ntu->SetBranchStatus("isZ", 1);                           ntu->SetBranchAddress("isZ", &isZ);
   if(type != "DA") {
     ntu->SetBranchStatus("PUit_TrueNumInteractions", 1);      ntu->SetBranchAddress("PUit_TrueNumInteractions", &npu);
   }
   ntu->SetBranchStatus("PV_z",1);                           ntu->SetBranchAddress("PV_z",&PV_z);
   ntu->SetBranchStatus("PV_n",1);                           ntu->SetBranchAddress("PV_n",&nVtx);
   //ele1                                                                                                                 
   ntu->SetBranchStatus("ele1_dxy_PV", 1);           ntu->SetBranchAddress("ele1_dxy_PV", &ele1_dxy_PV);
   ntu->SetBranchStatus("ele1_dz_PV", 1);            ntu->SetBranchAddress("ele1_dz_PV", &ele1_dz_PV);
   ntu->SetBranchStatus("ele1_scEtaWidth", 1);       ntu->SetBranchAddress("ele1_scEtaWidth", &ele1_scEtaWidth);
   ntu->SetBranchStatus("ele1_scPhiWidth", 1);       ntu->SetBranchAddress("ele1_scPhiWidth", &ele1_scPhiWidth);

   ntu->SetBranchStatus("ele1_charge", 1);         ntu->SetBranchAddress("ele1_charge",&ele1_charge);

   ntu->SetBranchStatus("ele1_eta", 1);          ntu->SetBranchAddress("ele1_eta",&ele1_eta);
   ntu->SetBranchStatus("ele1_phi", 1);          ntu->SetBranchAddress("ele1_phi",&ele1_phi);
   ntu->SetBranchStatus("ele1_fbrem", 1);        ntu->SetBranchAddress("ele1_fbrem",&ele1_fbrem);
   ntu->SetBranchStatus("ele1_tkPt", 1);         ntu->SetBranchAddress("ele1_tkPt",&ele1_tkPt);
   ntu->SetBranchStatus("ele1_tkP", 1);         ntu->SetBranchAddress("ele1_tkP",&ele1_tkP);
   ntu->SetBranchStatus("ele1_scE", 1);          ntu->SetBranchAddress("ele1_scE",&ele1_scE);
   ntu->SetBranchStatus("ele1_scERaw", 1);          ntu->SetBranchAddress("ele1_scERaw",&ele1_scERaw);
   ntu->SetBranchStatus("ele1_e3x3", 1);          ntu->SetBranchAddress("ele1_e3x3",&ele1_e3x3);
   ntu->SetBranchStatus("ele1_scEt", 1);          ntu->SetBranchAddress("ele1_scEt",&ele1_scEt);

   ntu->SetBranchStatus("ele1_scEta", 1);          ntu->SetBranchAddress("ele1_scEta",&ele1_scEta);
   ntu->SetBranchStatus("ele1_scPhi", 1);          ntu->SetBranchAddress("ele1_scPhi",&ele1_scPhi);

   //ele2                                                                                                                  
   ntu->SetBranchStatus("ele2_dxy_PV", 1);      ntu->SetBranchAddress("ele2_dxy_PV", &ele2_dxy_PV);
   ntu->SetBranchStatus("ele2_dz_PV", 1);       ntu->SetBranchAddress("ele2_dz_PV", &ele2_dz_PV);
   ntu->SetBranchStatus("ele2_scEtaWidth", 1);       ntu->SetBranchAddress("ele2_scEtaWidth", &ele2_scEtaWidth);
   ntu->SetBranchStatus("ele2_scPhiWidth", 1);       ntu->SetBranchAddress("ele2_scPhiWidth", &ele2_scPhiWidth);

   ntu->SetBranchStatus("ele2_charge", 1);         ntu->SetBranchAddress("ele2_charge",&ele2_charge);

   ntu->SetBranchStatus("ele2_eta", 1);          ntu->SetBranchAddress("ele2_eta",&ele2_eta);
   ntu->SetBranchStatus("ele2_phi", 1);          ntu->SetBranchAddress("ele2_phi",&ele2_phi);
   ntu->SetBranchStatus("ele2_fbrem", 1);        ntu->SetBranchAddress("ele2_fbrem",&ele2_fbrem);
   ntu->SetBranchStatus("ele2_tkPt", 1);         ntu->SetBranchAddress("ele2_tkPt",&ele2_tkPt);
   ntu->SetBranchStatus("ele2_tkP", 1);         ntu->SetBranchAddress("ele2_tkP",&ele2_tkP);
   ntu->SetBranchStatus("ele2_scE", 1);          ntu->SetBranchAddress("ele2_scE",&ele2_scE);
   ntu->SetBranchStatus("ele2_scERaw", 1);          ntu->SetBranchAddress("ele2_scERaw",&ele2_scERaw);
   ntu->SetBranchStatus("ele2_e3x3", 1);          ntu->SetBranchAddress("ele2_e3x3",&ele2_e3x3);
   ntu->SetBranchStatus("ele2_scEt", 1);          ntu->SetBranchAddress("ele2_scEt",&ele2_scEt);

   ntu->SetBranchStatus("ele2_scEta", 1);          ntu->SetBranchAddress("ele2_scEta",&ele2_scEta);
   ntu->SetBranchStatus("ele2_scPhi", 1);          ntu->SetBranchAddress("ele2_scPhi",&ele2_scPhi);

   //   ntu-> SetBranchStatus("nPV",     1);          ntu->SetBranchAddress("nPV",&nVtx);
   //   }

   /*
   // electron variables
   float scEta[2];
   float scPhi[2];
   float eta[2];
   float phi[2];
   float scE[2];
   float scERaw[2];
   float scEReg[2];
   float R9[2];
   int eleID[2];
   int chargeEle[2];

   if( type == "MCRD" )
     {
       ntu -> SetBranchStatus("*",0);
       ntu -> SetBranchStatus("nPU",      1);   ntu -> SetBranchAddress("nPU",&npu);
       ntu -> SetBranchStatus("nPV",      1);   ntu_MC -> SetBranchAddress("nPV",&nVtx);
     

       ntu -> SetBranchStatus("etaSCEle",      1);   ntu -> SetBranchAddress("etaSCEle",scEta);
       ntu -> SetBranchStatus("phiSCEle",      1);   ntu -> SetBranchAddress("phiSCEle",scPhi);
       ntu -> SetBranchStatus("etaEle",        1);   ntu -> SetBranchAddress("etaEle",eta);
       ntu -> SetBranchStatus("phiEle",        1);   ntu -> SetBranchAddress("phiEle",phi);
     }
   */

   /////////// Histos                                                                       

   //   std::vector<TLine*>  tangenti; 

   //fbrem inner->middle                                                                    
   TProfile* fBremVsEta = new TProfile("fBremVsEta", "", 240, -3., 3.);
   TProfile* fBremVsPhi = new TProfile("fBremVsPhi", "", 72, -3.2, 3.2);

   TProfile* fBremVsEta_fold = new TProfile("fBremVsEta_fold", "", 120, 0., 3.);
   TProfile* fBremVsPhi_fold = new TProfile("fBremVsPhi_fold", "", 36, 0., 3.2);

   TProfile* R9VsEta = new TProfile("R9VsEta", "", 60, -3., 3.);
   TProfile* R9VsEta_fold = new TProfile("R9VsEta_fold", "", 30, 0., 3.);
   TProfile* SPoSEVsEta = new TProfile("SPoSEVsEta", "", 60, -3., 3.);
   TProfile2D* SPoSEVsFBremVsEta = new TProfile2D("SPoSEVsFBremVsEta", "", 60, -3., 3., 500, -1., 4.);

   fBremVsEta->Sumw2();
   fBremVsPhi->Sumw2();

   fBremVsEta_fold->Sumw2();
   fBremVsPhi_fold->Sumw2();

   R9VsEta->Sumw2();
   R9VsEta_fold->Sumw2();

   SPoSEVsEta->Sumw2();
   SPoSEVsFBremVsEta->Sumw2();

   TH1F* h_R9 = new TH1F("h_R9", "", 1000, 0., 1.05);
   h_R9->Sumw2();

   TH2F* BSvsPhi = new TH2F("BSvsPhi", "", 60, -4., 4., 1000, -1., 1.);
   BSvsPhi->Sumw2();

   ////////vtx                                                                       
   TH1F* h_VtxZ = new TH1F("h_VtxZ", "", 1000, -20., 20.);
   TH1F* h_Vtx = new TH1F("h_Vtx", "", 100, 0., 100.);
   TH1F* h_dxyPV = new TH1F("h_dxyPV", "", 1000, -0.02, 0.02);
   TH1F* h_dzPV = new TH1F("h_dzPV", "", 1000, -0.02, 0.02);

   h_VtxZ->Sumw2();
   h_Vtx->Sumw2();
   h_dxyPV->Sumw2();
   h_dzPV->Sumw2();
   ////////vtx                                                                       


   /*
   std::vector<TProfile*> trendPerEtaBin;                           
   for(int xBin=0; xBin<30.; ++xBin){       
       float xMin = 0. + 0.1 * xBin;
       float xMax = 0.1 + 0.1 * xBin;

       char pippo[50];
       sprintf(pippo, "etaBin_%f-%f",xMin,xMax);                          
       TProfile* TLoc = new TProfile(pippo,"", 1500, -5., 10.);              
       trendPerEtaBin.push_back(TLoc);
     }

   std::vector<TProfile*> trendPerEtaBin_SPvsPt;                           
   for(int xBin=0; xBin<30.; ++xBin){       
       float xMin = 0. + 0.1 * xBin;
       float xMax = 0.1 + 0.1 * xBin;

       char pippo[50];
       sprintf(pippo, "SP_vsPt_etaBin_%f-%f",xMin,xMax);                          
       TProfile* TLoc = new TProfile(pippo,"", 100, 0., 100.);              
       trendPerEtaBin_SPvsPt.push_back(TLoc);
     }

   std::vector<TProfile*> trendPerEtaBin_SPvsR;                           
   for(int xBin=0; xBin<30.; ++xBin){       
       float xMin = 0. + 0.1 * xBin;
       float xMax = 0.1 + 0.1 * xBin;

       char pippo[50];
       sprintf(pippo, "SP_vsR_etaBin_%f-%f",xMin,xMax);                          
       TProfile* TLoc = new TProfile(pippo,"", 130, 0., 130.);              
       trendPerEtaBin_SPvsR.push_back(TLoc);
     }

   std::vector<TProfile*> trendPerEtaBin_SEvsPt;                           
   for(int xBin=0; xBin<30.; ++xBin){       
       float xMin = 0. + 0.1 * xBin;
       float xMax = 0.1 + 0.1 * xBin;

       char pippo[50];
       sprintf(pippo, "SE_vsEta_etaBin_%f-%f",xMin,xMax);                          
       TProfile* TLoc = new TProfile(pippo,"", 100, 0., 100.);              
       trendPerEtaBin_SEvsPt.push_back(TLoc);
     }

    */



   std::cout << " N Events: " << ntu->GetEntries() << std::endl;
   float ww = 1;
   for(int entry = 0; entry < ntu->GetEntries(); ++entry)
     {
       if( entry%1000 == 0 ) std::cout << " >>> reading entry " << entry << " / " << ntu->GetEntries() << "\r" << std::endl;
       ntu->GetEntry(entry);
       if(isZ == 0) continue;

       ele1_R9 = ele1_e3x3 / ele1_scERaw;
       ele2_R9 = ele2_e3x3 / ele2_scERaw;
       float effectiveFBrem_1 = -1. * log(1. - ele1_fbrem);
       float effectiveFBrem_2 = -1. * log(1. - ele2_fbrem);

       float effectiveEta_1 = ele1_scEta;
       float effectiveEta_2 = ele2_scEta;
       float effectivePhi_1 = ele1_scPhi;
       float effectivePhi_2 = ele2_scPhi;
       if(SCTK == "tk"){
	 effectiveEta_1 = ele1_eta;
	 effectiveEta_2 = ele2_eta;
	 effectivePhi_1 = ele1_phi;
	 effectivePhi_2 = ele2_phi;
       }

       if( type != "DA" && type != "MCRD") ww = (*puReweighting)[float(int(npu))];

       /*
       scEta1 = scEta[0];
       scEta2 = scEta[1];
       scPhi1 = scPhi[0];
       scPhi2 = scPhi[1];
       eta1 = eta[0];
       eta2 = eta[1];
       phi1 = phi[0];
       phi2 = phi[1];
       scEReg1 = scEReg[0];
       scEReg2 = scEReg[1];
       eleID1 = eleID[0];
       eleID2 = eleID[1];

       float theta1 = 2*atan(exp(-eta1));
       float theta2 = 2*atan(exp(-eta2));
       float Rt1 = sin(theta1);
       float Rt2 = sin(theta2);
       float Et1 = scEReg1*Rt1;
       float Et2 = scEReg2*Rt2;
       */

       //        std::cout << " npu = " << npu << std::endl;
       //        std::cout << " ww = " << ww << std::endl;


       h_Vtx->Fill(nVtx, ww);
       h_VtxZ->Fill(PV_z, ww);

       h_dxyPV->Fill(ele1_dxy_PV, ww);
       h_dzPV->Fill(ele1_dz_PV, ww);

       h_dxyPV->Fill(ele2_dxy_PV, ww);
       h_dzPV->Fill(ele2_dz_PV, ww);

       ///////////////////////////////////////////firstEle                    

//        std::cout << " >>>> ele1_charge = " << ele1_charge << std::endl;
//        std::cout << " >>>> ele2_charge = " << ele2_charge << std::endl;

       if( (chargeChosen == "Pos" && ele1_charge == 1) || (chargeChosen == "Neg" && ele1_charge == -1) || chargeChosen == "All" ) {
	 if( (zEleSelected == "Pos" && ele1_dz_PV > 0.) || (zEleSelected == "Neg" && ele1_dz_PV < 0.) || zEleSelected == "All"){

	 h_R9->Fill(ele1_R9, ww);
	 fBremVsEta_fold->Fill(fabs(effectiveEta_1), effectiveFBrem_1, ww);
	 fBremVsPhi_fold->Fill(fabs(effectivePhi_1), effectiveFBrem_1, ww);
	 fBremVsEta->Fill(effectiveEta_1, effectiveFBrem_1, ww);
	 fBremVsPhi->Fill(effectivePhi_1, effectiveFBrem_1, ww);

	 R9VsEta->Fill(effectiveEta_1, ele1_R9, ww);
	 R9VsEta_fold->Fill(fabs(effectiveEta_1), ele1_R9, ww);

	 SPoSEVsEta->Fill(effectiveEta_1, ele1_scPhiWidth/ele1_scEtaWidth, ww);
	 SPoSEVsFBremVsEta->Fill(effectiveEta_1, effectiveFBrem_1, ele1_scPhiWidth/ele1_scEtaWidth );			

	 BSvsPhi->Fill(effectivePhi_1, ele1_dxy_PV, ww);
	 }//dz
       }// charge 

       ///////////////////////////////////////////secondEle                    

       if( (chargeChosen == "Pos" && ele2_charge == 1) || (chargeChosen == "Neg" && ele2_charge == -1) || chargeChosen == "All" ) {
	 if( (zEleSelected == "Pos" && ele2_dz_PV > 0.) || (zEleSelected == "Neg" && ele2_dz_PV < 0.) || zEleSelected == "All"){

	 h_R9->Fill(ele2_R9, ww);
	 fBremVsEta_fold->Fill(fabs(effectiveEta_2), effectiveFBrem_2, ww);
	 fBremVsPhi_fold->Fill(fabs(effectivePhi_2), effectiveFBrem_2, ww);
	 fBremVsEta->Fill(effectiveEta_2, effectiveFBrem_2, ww);
	 fBremVsPhi->Fill(effectivePhi_2, effectiveFBrem_2, ww);

	 R9VsEta->Fill(effectiveEta_2, ele2_R9, ww);
	 R9VsEta_fold->Fill(fabs(effectiveEta_2), ele2_R9, ww);

	 SPoSEVsEta->Fill(effectiveEta_2, ele2_scPhiWidth/ele2_scEtaWidth, ww);
	 SPoSEVsFBremVsEta->Fill(effectiveEta_2, effectiveFBrem_2, ele2_scPhiWidth/ele2_scEtaWidth );			

	 BSvsPhi->Fill(effectivePhi_2, ele2_dxy_PV, ww);
	 }//dz
       }// charge
      
     } // loop over Entries


   std::string plotFolderName;
   plotFolderName = "ROOT_fBrem_"+SCTK+"_charge"+chargeChosen+"_zEle"+zEleSelected;
   //   plotFolderName = "ROOT_fBrem_sc_chargeAll_dzPos30";

   //    TFile pippo((plotFolderName+"/results_"+type+"_"+period+".root").c_str(),"recreate");
   //    std::cout << " >>>>> nomeFile = " << plotFolderName+"/results_"+type+"_"+period+".root" << std::endl;

   TFile pippo((plotFolderName+"/results_"+type+"_"+period+".root").c_str(),"recreate");
   std::cout << " >>>>> nomeFile = " << plotFolderName+"/results_"+type+"_"+period+".root" << std::endl;
   fBremVsEta_fold->Write();
   fBremVsPhi_fold->Write();
   fBremVsEta->Write();
   fBremVsPhi->Write();

   R9VsEta->Write();
   R9VsEta_fold->Write();

   SPoSEVsEta->Write();
   SPoSEVsFBremVsEta->Write();

   BSvsPhi->Write();

   ////////vtx                                                                     
   h_Vtx->Write();
   h_dxyPV->Write();
   h_dzPV->Write();
   h_VtxZ->Write();
   h_R9->Write();


   /*
   for(int xBin=0; xBin<30.; ++xBin){
     trendPerEtaBin.at(xBin)->Write();
     trendPerEtaBin_SPvsPt.at(xBin)->Write();
     trendPerEtaBin_SEvsPt.at(xBin)->Write();
     trendPerEtaBin_SPvsR.at(xBin)->Write();
   }
   */

  std::cout << " TUTTO FATTO " << std::endl;
}
