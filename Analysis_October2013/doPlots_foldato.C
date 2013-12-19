#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TPaveStats.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

void doPlots_foldato(){
  gROOT->ProcessLine(".x /Users/Arabella/Public/style.C");   
  std::cout << " ci sono " << std::endl;

  int R9_rebin = 8;
  int Profile_Rebin = 1;
  int Profile_bins = 30/Profile_Rebin;
  //  int Profile_bins = 120/Profile_Rebin;

  bool doEstimateMB = true;
  //    bool doEstimateMB = false;

  bool printNew = true;

  bool PV_Z_pos = false;
  bool PV_Z_neg = false;

  bool PV_Z_pos4 = false;
  bool PV_Z_neg4 = false;

  bool PV_Z_pos8 = false;
  bool PV_Z_neg8 = false;

  bool PV_Z_pos0 = false;
  bool PV_Z_neg0 = false;

  std::string period  = "ABCD";
  //    std::string period  = "AB";
  //    std::string period  = "C";
  //      std::string period  = "C2";
  //  std::string period  = "D";
  
  //  std::string type = "sc";
  std::string type = "tk";
  
  //  std::string charge = "chargePos";
  //std::string charge = "chargeNeg";
  std::string charge = "chargeAll";

  std::string zElePos = "zEleAll";
  //  std::string zElePos = "zEleNeg";
  //  std::string zElePos = "zElePos";
  
  std::string plotDir = "ROOT_fBrem_"+type+"_"+charge+"_"+zElePos+"/";
  std::string plotDirOut = "PLOT_foldato_"+type+"_"+charge+"_"+zElePos;


  //  if(Profile_Rebin == 1) plotDirOut = "PLOT_foldato_noRebin_"+period+"_"+type+"_"+charge;

  if(PV_Z_pos == true) {
    plotDir = "ROOT_fBrem_sc_chargeAll_dzPos";
    plotDirOut = "PLOT_foldato_sc_dzPos";
  }
  if(PV_Z_neg == true) {
    plotDir = "ROOT_fBrem_sc_chargeAll_dzNeg";
    plotDirOut = "PLOT_foldato_sc_dzNeg";
  }
  if(PV_Z_pos4 == true) {
    plotDir = "ROOT_fBrem_sc_chargeAll_dzPos4";
    plotDirOut = "PLOT_foldato_sc_dzPos4";
  }
  if(PV_Z_neg4 == true) {
    plotDir = "ROOT_fBrem_sc_chargeAll_dzNeg4";
    plotDirOut = "PLOT_foldato_sc_dzNeg4";
  }
  if(PV_Z_pos8 == true) {
    plotDir = "ROOT_fBrem_sc_chargeAll_dzPos8";
    plotDirOut = "PLOT_foldato_sc_dzPos8";
  }
  if(PV_Z_neg8 == true) {
    plotDir = "ROOT_fBrem_sc_chargeAll_dzNeg8";
    plotDirOut = "PLOT_foldato_sc_dzNeg8";
  }
  if(PV_Z_pos0 == true) {
    plotDir = "ROOT_fBrem_sc_chargeAll_dzPos0";
    plotDirOut = "PLOT_foldato_sc_dzPos0";
  }
  if(PV_Z_neg0 == true) {
    plotDir = "ROOT_fBrem_sc_chargeAll_dzNeg0";
    plotDirOut = "PLOT_foldato_sc_dzNeg0";
  }

  std::cout << " plotDir = " << plotDir << std::endl;
  std::cout << " plotDirOut = " << plotDirOut << std::endl;

  int nSets = 6;
  std::vector<string> typeSet;
  typeSet.push_back("DA");
  typeSet.push_back("MC");
  typeSet.push_back("F10");
  typeSet.push_back("F20");
  typeSet.push_back("F10S30");
  typeSet.push_back("F20S30");
//   typeSet.push_back("F10RD");
//   typeSet.push_back("F20RD");

  TFile* infile;

  TProfile** fBremVsEta = new TProfile*[nSets];
  TProfile** fBremVsPhi = new TProfile*[nSets];
  TH1F** PV_z = new TH1F*[nSets];
  TH1F** h_Vtx = new TH1F*[nSets];
  TH1F** h_dxyPV = new TH1F*[nSets];
  TH1F** h_dzPV = new TH1F*[nSets];
  TH1F** h_R9 = new TH1F*[nSets];

  std::cout << " Input from  " << plotDir << std::endl;
  std::cout << " Output from " << plotDirOut << std::endl;
  for(unsigned int ItSet=0; ItSet<typeSet.size(); ++ItSet){
  // for(unsigned int ItSet=0; ItSet<1; ++ItSet){
    std::cout << " da prendere  = " << (plotDir+"/results_"+typeSet.at(ItSet)+"_"+period+".root").c_str() << std::endl;
    infile = TFile::Open((plotDir+"/results_"+typeSet.at(ItSet)+"_"+period+".root").c_str(), "read");
    std::cout << " TFile = " << infile->GetName() << std::endl;

    fBremVsEta[ItSet] = (TProfile*)(infile->Get("fBremVsEta_fold"));
    fBremVsEta[ItSet]->SetName(("fBremVsEta_"+typeSet.at(ItSet)).c_str());
    fBremVsEta[ItSet]->SetDirectory(0);
    std::cout << " TFile = " << fBremVsEta[ItSet]->GetName() << std::endl;

    fBremVsPhi[ItSet] = (TProfile*)(infile->Get("fBremVsPhi_fold"));
    fBremVsPhi[ItSet]->SetName(("fBremVsPhi_"+typeSet.at(ItSet)).c_str());
    fBremVsPhi[ItSet]->SetDirectory(0);

    PV_z[ItSet] = (TH1F*)(infile->Get("h_VtxZ"));
    PV_z[ItSet]->SetName(("PV_z_"+typeSet.at(ItSet)).c_str());
    PV_z[ItSet]->SetDirectory(0);    
    std::cout << " TFile = " << PV_z[ItSet]->GetName() << std::endl;

    h_Vtx[ItSet] = (TH1F*)(infile->Get("h_Vtx"));
    h_Vtx[ItSet]->SetName(("h_Vtx_"+typeSet.at(ItSet)).c_str());
    h_Vtx[ItSet]->SetDirectory(0);

    h_dxyPV[ItSet] = (TH1F*)(infile->Get("h_dxyPV"));
    h_dxyPV[ItSet]->SetName(("h_dxyPV_"+typeSet.at(ItSet)).c_str());
    h_dxyPV[ItSet]->SetDirectory(0);

    h_dzPV[ItSet] = (TH1F*)(infile->Get("h_dzPV"));
    h_dzPV[ItSet]->SetName(("h_dzPV_"+typeSet.at(ItSet)).c_str());
    h_dzPV[ItSet]->SetDirectory(0);

    h_R9[ItSet] = (TH1F*)(infile->Get("h_R9"));
    h_R9[ItSet]->SetName(("h_R9_"+typeSet.at(ItSet)).c_str());
    h_R9[ItSet]->SetDirectory(0);
    std::cout << " TFile = " << h_R9[ItSet]->GetName() << std::endl;
   
    infile->Close();
  }

  std::cout << " ci sono 2 " << std::endl;

  int markerStyleExt[6] = {4, 21, 25, 25, 22, 22};
  float markerSizeExt[6] = {1, 0.7, 1, 1, 1, 1};
  int lineWidthExt[6] = {2, 2, 2, 2, 2, 2};
  //  int colorsExt[6] = {kBlack, kCyan, kOrange+8, kGreen+8, kOrange-2, kGreen+2};
  int colorsExt[6] = {kBlack, kCyan+2, kRed+2, kBlue+2, kOrange+7, kAzure+7};

  std::vector<int> markerStyle;
  std::vector<float> markerSize;
  std::vector<int> lineWidth;
  std::vector<int> colors;
  for(int posVec = 0; posVec<nSets; ++posVec){
    markerStyle.push_back(markerStyleExt[posVec]);
    markerSize.push_back(markerSizeExt[posVec]);
    lineWidth.push_back(lineWidthExt[posVec]);
    colors.push_back(colorsExt[posVec]);
  }

  for(unsigned int ItSet=0; ItSet<typeSet.size(); ++ItSet){
    fBremVsEta[ItSet]->SetMarkerStyle(markerStyle.at(ItSet));
    fBremVsEta[ItSet]->SetMarkerSize(markerSize.at(ItSet));
    fBremVsEta[ItSet]->SetLineWidth(lineWidth.at(ItSet));
    fBremVsEta[ItSet]->SetLineColor(colors.at(ItSet));
    fBremVsEta[ItSet]->SetMarkerColor(colors.at(ItSet));
  }

  TLegend *legCuts1 = new TLegend(0.70,0.80,0.99,1.,NULL,"brNDC");
  legCuts1->SetTextFont(42);
  legCuts1->SetFillColor(kWhite);
  legCuts1->SetLineColor(kWhite);
  legCuts1->SetShadowColor(kWhite);
  for(unsigned int ItSet=0; ItSet<typeSet.size(); ++ItSet)
    legCuts1->AddEntry(fBremVsEta[ItSet], (typeSet.at(ItSet)).c_str(), "pl");
  

  TCanvas* cfBrem_Eta = new TCanvas();
  cSingle  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
  cSingle->SetBottomMargin(0.1);  cSingle->SetTopMargin(0.05);
  cSingle->Draw();   cSingle->cd();   gPad->SetGrid();
  fBremVsEta[0]->GetXaxis()->SetRangeUser(0, 3);
  fBremVsEta[0]->GetYaxis()->SetRangeUser(0, 1.6);
  fBremVsEta[0]->GetXaxis()->SetTitle((type+" #eta").c_str());
  fBremVsEta[0]->GetYaxis()->SetTitle("fBrem");
  fBremVsEta[0]->Draw();
  for(int posVec = 1; posVec<nSets; ++posVec)
    fBremVsEta[posVec]->Draw("same");
  legCuts1->Draw("same");
  cfBrem_Eta->Print((plotDirOut+"/fBremVsEta.png").c_str(),".png");
  cfBrem_Eta->Print((plotDirOut+"/fBremVsEta.pdf").c_str(),".pdf");
  
  //fBrem ratio
  TH1F** fBremRatio = new TH1F*[nSets];
  for(unsigned int ItSet=0; ItSet<typeSet.size(); ++ItSet){
    fBremRatio[ItSet] = new TH1F(("fBremRatio_"+typeSet.at(ItSet)).c_str(), "", Profile_bins, 0., 3.);                  
    fBremRatio[ItSet]->Sumw2();
    fBremVsEta[ItSet]->Rebin(Profile_Rebin);
  }

  // CAVEAT change Profile_bins;
  float RatioNum[6][30];
  float RatioNumE[6][30];
  for(unsigned int ItSetR=0; ItSetR<typeSet.size(); ++ItSetR){                                               
    std::cout << " ItSetR = " << ItSetR << " fBremRatio[ItSetR]->GetNbinsX() " << fBremRatio[ItSetR]->GetNbinsX() << std::endl;
    for(int bin = 1; bin <= fBremRatio[ItSetR]->GetNbinsX(); ++bin){                                         
      std::cout << " fBremVsEta[ItSet]->GetBinContent(bin) = " << fBremVsEta[ItSetR]->GetBinContent(bin) << std::endl;
      RatioNum[ItSetR][bin-1] = fBremVsEta[ItSetR]->GetBinContent(bin);
      RatioNumE[ItSetR][bin-1] = fBremVsEta[ItSetR]->GetBinError(bin);
      
      //       std::cout << " SET = " << ItSetR << " bin = " << bin-1 << std::endl;
      //       std::cout << " RatioNum[ItSetR][bin-1] = " << RatioNum[ItSetR][bin-1] << std::endl;

      if(RatioNum[0][bin-1] != 0. && ItSetR != 0){
 	double val = RatioNum[ItSetR][bin-1]/RatioNum[0][bin-1];
 	double relEr = pow(RatioNumE[ItSetR][bin-1]/RatioNum[ItSetR][bin-1],2.) + pow(RatioNumE[0][bin-1]/RatioNum[0][bin-1],2.);
	//std::cout << " val = " << val << std::endl;
	//  	std::cout << " relErr = " << relEr << std::endl;
 	fBremRatio[ItSetR]->SetBinContent(bin, val);
 	fBremRatio[ItSetR]->SetBinError(bin, val * sqrt(relEr));
      }
      else{
	fBremRatio[ItSetR]->SetBinContent(bin, 1.);
	fBremRatio[ItSetR]->SetBinError(bin, 0.);
      }

    }
  }

//   for(int bin = 0; bin < Profile_bins; ++bin){                                         
//     std::cout << " [bin] = " << bin << std::endl;
//     std::cout << " RatioNum[1][bin] = " << RatioNum[1][bin] << std::endl;
//     std::cout << " RatioNum[2][bin] = " << RatioNum[2][bin] << std::endl;
//     std::cout << " RatioNum[3][bin] = " << RatioNum[3][bin] << std::endl;
//   }

//  return;
  for(unsigned int ItSet=1; ItSet<typeSet.size(); ++ItSet){
    fBremRatio[ItSet]->SetMarkerStyle(markerStyle.at(ItSet));
    fBremRatio[ItSet]->SetMarkerSize(markerSize.at(ItSet));
    fBremRatio[ItSet]->SetLineWidth(lineWidth.at(ItSet));
    fBremRatio[ItSet]->SetLineColor(colors.at(ItSet));
    fBremRatio[ItSet]->SetMarkerColor(colors.at(ItSet));
  }

  TLegend *legCutsRatio = new TLegend(0.70,0.80,0.99,1.,NULL,"brNDC");
  legCutsRatio->SetTextFont(42);
  legCutsRatio->SetFillColor(kWhite);
  legCutsRatio->SetLineColor(kWhite);
  legCutsRatio->SetShadowColor(kWhite);
  for(unsigned int ItSet=0; ItSet<typeSet.size(); ++ItSet)
    legCutsRatio->AddEntry(fBremRatio[ItSet], (typeSet.at(ItSet)).c_str(), "pl");
  

  TCanvas* cfBrem_Ratio = new TCanvas();
  cSingle  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
  cSingle->SetBottomMargin(0.1);  cSingle->SetTopMargin(0.05);
  cSingle->Draw();   cSingle->cd();   gPad->SetGrid();
  fBremRatio[1]->GetXaxis()->SetRangeUser(0, 3);
  fBremRatio[1]->GetYaxis()->SetRangeUser(0.7, 1.3);
  fBremRatio[1]->GetXaxis()->SetTitle((type+" #eta").c_str());
  fBremRatio[1]->GetYaxis()->SetTitle("fBrem MC/DA");
  fBremRatio[1]->Draw();
  for(int posVec = 2; posVec<nSets; ++posVec)
    fBremRatio[posVec]->Draw("same");
  legCutsRatio->Draw("same");
  cfBrem_Ratio->Print((plotDirOut+"/fBremRatio.png").c_str(),".png");
  cfBrem_Ratio->Print((plotDirOut+"/fBremRatio.pdf").c_str(),".pdf");


  TCanvas* cfBremAll = new TCanvas();
  cLower  = new TPad("pad_0","pad_0",0.00,0.00,1.00,0.35);
  cUpper  = new TPad("pad_2","pad_2",0.00,0.35,1.00,1.00);
  cLower->SetBottomMargin(0.1);  cLower->SetTopMargin(0.05);
  cUpper->SetBottomMargin(0.1);  cUpper->SetTopMargin(0.05);

  cLower->Draw();   cUpper->Draw();
  cUpper->cd();   gPad->SetGrid();

  fBremVsEta[0]->GetXaxis()->SetRangeUser(0, 3);
  fBremVsEta[0]->GetYaxis()->SetRangeUser(0, 1.6);
  fBremVsEta[0]->GetXaxis()->SetTitle((type+" #eta").c_str());
  fBremVsEta[0]->GetYaxis()->SetTitle("fBrem");
  fBremVsEta[0]->Draw();
  for(int posVec = 1; posVec<nSets; ++posVec)
    fBremVsEta[posVec]->Draw("same");
  legCuts1->Draw("same");

  cLower->cd();   gPad->SetGrid();
  fBremRatio[1]->GetXaxis()->SetRangeUser(0, 3);
  fBremRatio[1]->GetYaxis()->SetRangeUser(0.7, 1.3);
  fBremRatio[1]->GetXaxis()->SetTitle((type+" #eta").c_str());
  fBremRatio[1]->GetYaxis()->SetTitle("fBrem MC/DA");
  fBremRatio[1]->Draw();
  for(int posVec = 2; posVec<nSets; ++posVec)
    fBremRatio[posVec]->Draw("same");
  legCutsRatio->Draw("same");
  cfBremAll->Print((plotDirOut+"/fBremAllVsEta.png").c_str(),".png");
  cfBremAll->Print((plotDirOut+"/fBremAllVsEta.pdf").c_str(),".pdf");

  /////////////// calibration:
  TGraphErrors** calibration_Eta = new TGraphErrors*[Profile_bins];
  TGraphErrors** calibration_Eta_F20 = new TGraphErrors*[Profile_bins];

  for(int nEtaBin=0; nEtaBin<Profile_bins; ++nEtaBin){
    calibration_Eta[nEtaBin] = new TGraphErrors();
    calibration_Eta[nEtaBin]->SetName(Form("calib_Eta_bin%d",nEtaBin));
    calibration_Eta_F20[nEtaBin] = new TGraphErrors();
    calibration_Eta_F20[nEtaBin]->SetName(Form("calib_Eta_F20_bin%d",nEtaBin));
    //     std::cout << " calibration_Eta[nEtaBin]->this = " << calibration_Eta[nEtaBin] << std::endl;   
    
    if(bin == 1) 
      //      std::cout << " bin = " << bin  << std::endl; 
    
    calibration_Eta[nEtaBin]->SetPoint(0, 0, 0);
    calibration_Eta[nEtaBin]->SetPoint(1, 1., RatioNum[1][nEtaBin]);
    calibration_Eta[nEtaBin]->SetPoint(2, 1.1, RatioNum[2][nEtaBin]);
    calibration_Eta[nEtaBin]->SetPoint(3, 1.2, RatioNum[3][nEtaBin]);
    calibration_Eta[nEtaBin]->SetPoint(4, 5., 5.);

    //     std::cout << " RatioNum[1][bin] = " << RatioNum[1][bin] << std::endl;    
    //     std::cout << " RatioNum[2][bin] = " << RatioNum[2][bin] << std::endl;    
    //     std::cout << " RatioNum[3][bin] = " << RatioNum[3][bin] << std::endl;    

    calibration_Eta_F20[nEtaBin]->SetPoint(0, 0, 0);
    calibration_Eta_F20[nEtaBin]->SetPoint(1, 1., RatioNum[1][nEtaBin]);
    calibration_Eta_F20[nEtaBin]->SetPoint(2, 1.2, RatioNum[3][nEtaBin]);
    
    calibration_Eta[nEtaBin]->SetPointError(0, 0, 0);
    calibration_Eta[nEtaBin]->SetPointError(1, 0., RatioNumE[1][nEtaBin]);
    calibration_Eta[nEtaBin]->SetPointError(2, 0., RatioNumE[2][nEtaBin]);
    calibration_Eta[nEtaBin]->SetPointError(3, 0., RatioNumE[3][nEtaBin]);
    calibration_Eta[nEtaBin]->SetPointError(4, 0., 0);
    
    calibration_Eta_F20[nEtaBin]->SetPointError(0, 0, 0);
    calibration_Eta_F20[nEtaBin]->SetPointError(1, 0., RatioNumE[1][nEtaBin]);
    calibration_Eta_F20[nEtaBin]->SetPointError(2, 0., RatioNumE[3][nEtaBin]);

    //     TCanvas* c1 = new TCanvas();
    //     calibration_Eta[nEtaBin]->Draw("AP");
    //     return;
  }

  if(doEstimateMB){
    //DA: errStat - ChargeP - ChargeN - zEleP - zEleN 
    std::vector<std::vector<float> > error_DA;
    std::vector<std::vector<float> > error_F10S30;
    std::vector<std::vector<float> > error_F20S30;

    std::vector<std::string> myfileN;
    if(Profile_Rebin == 1){
      myfileN.push_back("PLOT_foldato_tk_chargeAll_zEleAll/ErrorBands_tk_chargeAll_zEleAll.txt");
      myfileN.push_back("PLOT_foldato_tk_chargeNeg_zEleAll/ErrorBands_tk_chargeNeg_zEleAll.txt");
      myfileN.push_back("PLOT_foldato_tk_chargePos_zEleAll/ErrorBands_tk_chargePos_zEleAll.txt");
      myfileN.push_back("PLOT_foldato_tk_chargeAll_zEleNeg/ErrorBands_tk_chargeAll_zEleNeg.txt");
      myfileN.push_back("PLOT_foldato_tk_chargeAll_zElePos/ErrorBands_tk_chargeAll_zElePos.txt");
    }
    for(unsigned int countF = 0; countF<myfileN.size(); ++ countF){
      std::vector<float> errors(0);
      error_DA.push_back(errors);
      error_F10S30.push_back(errors);
      error_F20S30.push_back(errors);
    }
    
    int caso = 0;
    int count = 0;
    if(printNew == false){
      for(unsigned int countF = 0; countF<myfileN.size(); ++ countF){
	string line;
	ifstream myfile(myfileN.at(countF).c_str());
	if (myfile.is_open()){
	  while( myfile.good() ){
	    
	    //	    std::cout << " before caso = " << caso << std::endl;
	    getline(myfile,line);
	    //	    cout << line << endl;
	    if(std::string(line) == "DA")  {caso = 1; count = 0;}
	    if(std::string(line) == "MC")  {caso = 2; count = 0;}
	    if(std::string(line) == "F10")  {caso = 3;count = 0;}
	    if(std::string(line) == "F20")  {caso = 4;count = 0;}
	    if(std::string(line) == "F10S30")  {caso = 5; count = 0;}
	    if(std::string(line) == "F20S30")  {caso = 6; count = 0;}
	    //	    std::cout << " after caso = " << caso << std::endl;
	    
	    if(caso != 0 && count < Profile_bins-1){
	      ++count;
	      float a1, a2;
	      myfile >> a1 >> a2;
	      //	      std::cout << " a1 = " << a1 << " a2 = " << a2 << std::endl;
	      
	      if(countF == 0){
		if(caso == 1) (error_DA.at(countF)).push_back(float(a2));
		if(caso == 5) (error_F10S30.at(countF)).push_back(float(a2));
		if(caso == 6) (error_F20S30.at(countF)).push_back(float(a2));
		}
	      else{
		if(caso == 1) (error_DA.at(countF)).push_back(float(a1));
		if(caso == 5) (error_F10S30.at(countF)).push_back(float(a1));
		if(caso == 6) (error_F20S30.at(countF)).push_back(float(a1));
	      }
	    }
	  }
	}
	myfile.close();
      }
    }

    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////

    std::cout << " prima dei Canvas " << std::endl;
    TH1F* MBVsEta_Closure = new TH1F("MBVsEta_Closure", "", Profile_bins, 0., 3.);
    
    for(unsigned int ItSetM=2; ItSetM<3; ++ItSetM){ 
      std::cout << " ItSetM = " << ItSetM << std::endl;

      for(int nBinE = 0; nBinE<Profile_bins; ++nBinE){
	//for(int nBinE = 0; nBinE<1; ++nBinE){
	//	std::cout << " nBinE = " << nBinE << std::endl;

	TCanvas* cCalib = new TCanvas();
	cCalib->cd();

	TGraphErrors* calibration_Eta_Measured = new TGraphErrors();
	TGraphErrors* grint = new TGraphErrors(2500);
	
	TF1* fitFunc = new TF1("fitFunc", "[0] + [1] * x", 0.5, 2.);
	fitFunc->SetParameters(0.2, 1.);
	fitFunc->SetLineColor(kBlue);
	fitFunc->SetLineWidth(1);

	calibration_Eta_F20[nBinE]->Fit(fitFunc, "QRS+");
	calibration_Eta_F20[nBinE]->Draw("AP");

	for(int i=0; i<2500; i++) grint->SetPoint(i, 0.5 + 1.5/2500.*i, 0);
	//Compute the confidence intervals at the x points of the created graph
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);

	float DA_x = 0.;
	float DA_ex = 0.;

	float P0 = fitFunc->GetParameter(0);
	float P1 = fitFunc->GetParameter(1);
	float EP0 = fitFunc->GetParError(0);
	float EP1 = fitFunc->GetParError(1);

	if(P1 != 0.){
	  DA_x = (RatioNum[ItSetM][nBinE] - P0) / P1;
	  DA_ex = sqrt( pow(RatioNumE[ItSetM][nBinE]/P1, 2)  );
	  //RARA da sistemare
	  //DA_ex = sqrt( pow(RatioNumE[ItSetM][nBinE]/P1, 2) + pow(EP0/P1, 2) + pow((RatioNum[ItSetM][nBinE] - P0)/P1/P1*EP1, 2)  );
	}

	float DA_x_L = 100.;
	float DA_x_M = -100.;
	for(unsigned int iGr=0; iGr<2500; iGr++){
	  double x,y;
	  double Ex,Ey;
	  grint->GetPoint(iGr, x, y);
	  Ey = grint->GetErrorY(iGr);
	  if(RatioNum[ItSetM][nBinE] < y+Ey && DA_x_L == 100.) DA_x_L = x;
	  if(RatioNum[ItSetM][nBinE] < y-Ey && DA_x_M == -100.) DA_x_M = x;
	  if(DA_x_L != 100. && DA_x_M != -100.) break;
	}
	
	//	std::cout << " >>>> prima di SetPoint " << std::endl;
	calibration_Eta_Measured->SetPoint(0, DA_x, RatioNum[ItSetM][nBinE]);
	calibration_Eta_Measured->SetPointError(0, DA_ex, RatioNumE[ItSetM][nBinE]);
	
	if(RatioNum[ItSetM][nBinE] != 0){
	  MBVsEta_Closure->SetBinContent(nBinE+1, DA_x);
	  MBVsEta_Closure->SetBinError(nBinE+1, DA_ex);
	}

	//    	std::cout << " >>>> prima di DrawSetting " << std::endl;
	calibration_Eta_F20[nBinE]->SetMaximum(1.8);
	calibration_Eta_F20[nBinE]->SetMinimum(0.);
	calibration_Eta_F20[nBinE]->GetXaxis()->SetRangeUser(0.5, 2.);
	calibration_Eta_F20[nBinE]->GetXaxis()->SetTitle(Form("iEta bin %d",nBinE));
	calibration_Eta_F20[nBinE]->SetMarkerStyle(7);
	calibration_Eta_Measured->SetMarkerStyle(7);
	calibration_Eta_Measured->SetMarkerColor(kRed);
	calibration_Eta_Measured->SetLineColor(kRed);
    
	calibration_Eta_F20[nBinE]->Draw("AP");

	grint->SetFillColor(kMagenta);
	grint->SetFillStyle(3001);
	grint->Draw("E3, same");

	calibration_Eta_F20[nBinE]->Draw("P,same");
	calibration_Eta_Measured->Draw("P,same");

	cCalib->Print((plotDirOut+Form("/Calibration_Closure/CaliEtaBin_Closure_Set%d_bin_%d",ItSetM, nBinE)+".gif").c_str(),".gif");

	delete cCalib;
	delete fitFunc;
	delete calibration_Eta_Measured;
	delete grint;
      }
      
      TFile outFileResult(("outFileResult_Closure_"+typeSet.at(ItSetM)).c_str(), "recreate");
      MBVsEta_Closure->Write(("MBVsEta_Closure"+typeSet.at(ItSetM)).c_str());
      outFileResult.Close();
    }

    //measure from DATA
   
    std::cout << " prima dei Canvas " << std::endl;
    //rarara
    //     TGraphErrors** calibration_Eta_Measured = new TGraphErrors*[Profile_bins];
    //     TCanvas** cCalib = new TCanvas*[Profile_bins];
    TH1F** MBVsEta = new TH1F*[nSets];
    TH1F** MBVsEta_CL = new TH1F*[nSets];
    TH1F** MBVsEta_Stat = new TH1F*[nSets];
    // confidence level fit
    // charge
    TH1F** MBVsEta_Charge = new TH1F*[nSets];
    //zEle position
    TH1F** MBVsEta_zEle = new TH1F*[nSets];
    TH1F** MBVsEta_eRR = new TH1F*[nSets];

    
    for(unsigned int ItSetM=0; ItSetM<typeSet.size(); ++ItSetM){ 
      //      std::cout << " ItSetM = " << ItSetM << std::endl;

      MBVsEta[ItSetM] = new TH1F(("MBVsEta_"+typeSet.at(ItSetM)).c_str(), "", Profile_bins, 0., 3.);
      MBVsEta_CL[ItSetM] = new TH1F(("MBVsEta_CL_"+typeSet.at(ItSetM)).c_str(), "", Profile_bins, 0., 3.);
      MBVsEta_Stat[ItSetM] = new TH1F(("MBVsEta_Stat_"+typeSet.at(ItSetM)).c_str(), "", Profile_bins, 0., 3.);
      MBVsEta_Charge[ItSetM] = new TH1F(("MBVsEta_Charge_"+typeSet.at(ItSetM)).c_str(), "", Profile_bins, 0., 3.);
      MBVsEta_zEle[ItSetM] = new TH1F(("MBVsEta_zEle_"+typeSet.at(ItSetM)).c_str(), "", Profile_bins, 0., 3.);
      MBVsEta_eRR[ItSetM] = new TH1F(("MBVsEta_eRR_"+typeSet.at(ItSetM)).c_str(), "", Profile_bins, 0., 3.);
    }

    std::cout << " prima degli outout " << std::endl;                                               
    ofstream myfileOu;                                                                              
    std::string myfileNameOu = plotDirOut+"/ErrorBands_"+type+"_"+charge+"_"+zElePos+".txt";                                                         
    if(printNew) myfileOu.open(myfileNameOu.c_str(), ios::out);
    std::cout << " dopo gli outout " << std::endl;                                               

    for(unsigned int ItSetM=0; ItSetM<typeSet.size(); ++ItSetM){ 
      std::cout << " ItSetM = " << ItSetM << std::endl;
      if(printNew) myfileOu << typeSet.at(ItSetM) << std::endl;
      
      for(int nBinE = 0; nBinE<Profile_bins; ++nBinE){
	//for(int nBinE = 0; nBinE<1; ++nBinE){
	//	std::cout << " nBinE = " << nBinE << std::endl;

	TCanvas* cCalib = new TCanvas();
	cCalib->cd();

	TGraphErrors* calibration_Eta_Measured = new TGraphErrors();
	TGraphErrors* grint = new TGraphErrors(2500);
	
	TF1* fitFunc = new TF1("fitFunc", "[0] + [1] * x", 0.5, 2.);
	fitFunc->SetParameters(0.2, 1.);
	fitFunc->SetLineColor(kBlue);
	fitFunc->SetLineWidth(1);

	//	fitFunc->Draw();
	//      calibration_Eta[nBinE]->Draw("AP");
	//	return;
	//	std::cout << " >>>> prima di fit  " << std::endl;

	calibration_Eta[nBinE]->Fit(fitFunc, "QRS+");
	calibration_Eta[nBinE]->Draw("AP");

	//	std::cout << " >>>> dopo fit prima di grindt " << std::endl;
	// 	return;

	for(int i=0; i<2500; i++) grint->SetPoint(i, 0.5 + 1.5/2500.*i, 0);
	//Compute the confidence intervals at the x points of the created graph
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);

	float DA_x = 0.;
	float DA_ex = 0.;

	float P0 = fitFunc->GetParameter(0);
	float P1 = fitFunc->GetParameter(1);
	float EP0 = fitFunc->GetParError(0);
	float EP1 = fitFunc->GetParError(1);

	if(P1 != 0.){
	  DA_x = (RatioNum[ItSetM][nBinE] - P0) / P1;
	  DA_ex = sqrt( pow(RatioNumE[ItSetM][nBinE]/P1, 2)  );
	  //RARA da sistemare
	  //DA_ex = sqrt( pow(RatioNumE[ItSetM][nBinE]/P1, 2) + pow(EP0/P1, 2) + pow((RatioNum[ItSetM][nBinE] - P0)/P1/P1*EP1, 2)  );
	}

	//	std::cout << " >>>> prima di float DA_x_L = 100. " << std::endl;
	float DA_x_L = 100.;
	float DA_x_M = -100.;
	for(unsigned int iGr=0; iGr<2500; iGr++){
	  double x,y;
	  double Ex,Ey;
	  grint->GetPoint(iGr, x, y);
	  Ey = grint->GetErrorY(iGr);
	  if(RatioNum[ItSetM][nBinE] < y+Ey && DA_x_L == 100.) DA_x_L = x;
	  if(RatioNum[ItSetM][nBinE] < y-Ey && DA_x_M == -100.) DA_x_M = x;
	  if(DA_x_L != 100. && DA_x_M != -100.) break;
	}
	
	//	std::cout << " >>>> prima di SetPoint " << std::endl;
	calibration_Eta_Measured->SetPoint(0, DA_x, RatioNum[ItSetM][nBinE]);
	calibration_Eta_Measured->SetPointError(0, DA_ex, RatioNumE[ItSetM][nBinE]);

	//	if(printNew) myfileOu << DA_x[nBinE] << " " << DA_ex[nBinE] << std::endl;
	if(printNew) myfileOu << DA_x << " " << DA_ex << std::endl;
	
	if(RatioNum[ItSetM][nBinE] != 0){
	  MBVsEta[ItSetM]->SetBinContent(nBinE+1, DA_x);
	  MBVsEta[ItSetM]->SetBinError(nBinE+1, DA_ex);

	  MBVsEta_CL[ItSetM]->SetBinContent(nBinE+1, DA_x);
	  MBVsEta_CL[ItSetM]->SetBinError(nBinE+1, fabs(DA_x_M - DA_x_L)/2.);

	  if(printNew == false){
	    MBVsEta_Stat[ItSetM]->SetBinContent(nBinE+1, DA_x);
	    MBVsEta_Charge[ItSetM]->SetBinContent(nBinE+1, DA_x);
	    MBVsEta_zEle[ItSetM]->SetBinContent(nBinE+1, DA_x);
	    MBVsEta_eRR[ItSetM]->SetBinContent(nBinE+1, DA_x);

	    if(ItSetM == 0){
	      MBVsEta_Stat[ItSetM]->SetBinError(nBinE+1, (error_DA.at(0)).at(nBinE));
	      float Err_Ch_max = TMath::Max( DA_x, TMath::Max( (error_DA.at(1)).at(nBinE), (error_DA.at(2)).at(nBinE)) );
	      float Err_Ch_min = TMath::Min( DA_x, TMath::Min( (error_DA.at(1)).at(nBinE), (error_DA.at(2)).at(nBinE)) );
	      float Err_Ch = fabs(Err_Ch_max - Err_Ch_min);
// 	      std::cout << " >>>> Err_Ch_max  012 = " << Err_Ch_max << " = DA_x = " << DA_x 
// 			<< " (error_DA.at(1)).at(nBinE) = " << (error_DA.at(1)).at(nBinE) 
// 			<< " (error_DA.at(2)).at(nBinE) = " << (error_DA.at(2)).at(nBinE) << std::endl;

	      MBVsEta_Charge[ItSetM]->SetBinError(nBinE+1, Err_Ch/2.);
	      float Err_zEle_max = TMath::Max( DA_x, TMath::Max( (error_DA.at(3)).at(nBinE), (error_DA.at(4)).at(nBinE)) );
	      float Err_zEle_min = TMath::Min( DA_x, TMath::Min( (error_DA.at(3)).at(nBinE), (error_DA.at(4)).at(nBinE)) );
	      float Err_zEle = fabs(Err_zEle_max - Err_zEle_min);
// 	      std::cout << " >>>> Err_Ch_max  034 = " << Err_Ch_max << " = DA_x = " << DA_x 
// 			<< " (error_DA.at(3)).at(nBinE) = " << (error_DA.at(3)).at(nBinE) 
// 			<< " (error_DA.at(4)).at(nBinE) = " << (error_DA.at(4)).at(nBinE) << std::endl;

	      MBVsEta_zEle[ItSetM]->SetBinError(nBinE+1, Err_zEle/2.);
	      float Err_Max = TMath::Max(Err_Ch_max, Err_zEle_max);
	      float Err_Min = TMath::Min(Err_Ch_min, Err_zEle_min);
	      float Err_Sys = fabs(Err_Max - Err_Min);
	      float Err_Tot = sqrt(pow(Err_Sys/2.,2) + pow(fabs(DA_x_M - DA_x_L)/2., 2));
// 	      std::cout << " >>>> Err_Max = " << Err_Max << " Err_Min = " << Err_Min 
// 			<< " Err_Sys = " << Err_Sys 
// 			<< " Err_Tot = " << Err_Tot << std::endl;

	      MBVsEta_eRR[ItSetM]->SetBinError(nBinE+1, Err_Tot);
	    }
	    if(ItSetM == 4){
	      MBVsEta_Stat[ItSetM]->SetBinError(nBinE+1, (error_DA.at(0)).at(nBinE));
	      float Err_Ch_max = TMath::Max( DA_x, TMath::Max( (error_F10S30.at(1)).at(nBinE), (error_F10S30.at(2)).at(nBinE)) );
	      float Err_Ch_min = TMath::Min( DA_x, TMath::Min( (error_F10S30.at(1)).at(nBinE), (error_F10S30.at(2)).at(nBinE)) );
	      float Err_Ch = fabs(Err_Ch_max - Err_Ch_min);

	      MBVsEta_Charge[ItSetM]->SetBinError(nBinE+1, Err_Ch/2.);
	      float Err_zEle_max = TMath::Max( DA_x, TMath::Max( (error_F10S30.at(3)).at(nBinE), (error_F10S30.at(4)).at(nBinE)) );
	      float Err_zEle_min = TMath::Min( DA_x, TMath::Min( (error_F10S30.at(3)).at(nBinE), (error_F10S30.at(4)).at(nBinE)) );
	      float Err_zEle = fabs(Err_zEle_max - Err_zEle_min);

	      MBVsEta_zEle[ItSetM]->SetBinError(nBinE+1, Err_zEle/2.);
	      float Err_Max = TMath::Max(Err_Ch_max, Err_zEle_max);
	      float Err_Min = TMath::Min(Err_Ch_min, Err_zEle_min);
	      float Err_Sys = fabs(Err_Max - Err_Min);
	      float Err_Tot = sqrt(pow(Err_Sys/2.,2) + pow(fabs(DA_x_M - DA_x_L)/2., 2));
	      MBVsEta_eRR[ItSetM]->SetBinError(nBinE+1, Err_Tot);
	    }
	    if(ItSetM == 5){
	      MBVsEta_Stat[ItSetM]->SetBinError(nBinE+1, (error_DA.at(0)).at(nBinE));
	      float Err_Ch_max = TMath::Max( DA_x, TMath::Max( (error_F20S30.at(1)).at(nBinE), (error_F20S30.at(2)).at(nBinE)) );
	      float Err_Ch_min = TMath::Min( DA_x, TMath::Min( (error_F20S30.at(1)).at(nBinE), (error_F20S30.at(2)).at(nBinE)) );
	      float Err_Ch = fabs(Err_Ch_max - Err_Ch_min);

	      MBVsEta_Charge[ItSetM]->SetBinError(nBinE+1, Err_Ch/2.);
	      float Err_zEle_max = TMath::Max( DA_x, TMath::Max( (error_F20S30.at(3)).at(nBinE), (error_F20S30.at(4)).at(nBinE)) );
	      float Err_zEle_min = TMath::Min( DA_x, TMath::Min( (error_F20S30.at(3)).at(nBinE), (error_F20S30.at(4)).at(nBinE)) );
	      float Err_zEle = fabs(Err_zEle_max - Err_zEle_min);

	      MBVsEta_zEle[ItSetM]->SetBinError(nBinE+1, Err_zEle/2.);
	      float Err_Max = TMath::Max(Err_Ch_max, Err_zEle_max);
	      float Err_Min = TMath::Min(Err_Ch_min, Err_zEle_min);
	      float Err_Sys = fabs(Err_Max - Err_Min);
	      float Err_Tot = sqrt(pow(Err_Sys/2.,2) + pow(fabs(DA_x_M - DA_x_L)/2., 2));
	      MBVsEta_eRR[ItSetM]->SetBinError(nBinE+1, Err_Tot);
	    }
	  }//printNew == false

	}

	//    	std::cout << " >>>> prima di DrawSetting " << std::endl;
	calibration_Eta[nBinE]->SetMaximum(1.8);
	calibration_Eta[nBinE]->SetMinimum(0.);
	calibration_Eta[nBinE]->GetXaxis()->SetRangeUser(0.5, 2.);
	calibration_Eta[nBinE]->GetXaxis()->SetTitle(Form("iEta bin %d",nBinE));
	calibration_Eta[nBinE]->SetMarkerStyle(7);
	calibration_Eta_Measured->SetMarkerStyle(7);
	calibration_Eta_Measured->SetMarkerColor(kRed);
	calibration_Eta_Measured->SetLineColor(kRed);
    
	calibration_Eta[nBinE]->Draw("AP");

	grint->SetFillColor(kMagenta);
	grint->SetFillStyle(3001);
	grint->Draw("E3, same");

	calibration_Eta[nBinE]->Draw("P,same");
	calibration_Eta_Measured->Draw("P,same");

	cCalib->Print((plotDirOut+Form("/Calibration/CaliEtaBin_Set%d_bin_%d",ItSetM, nBinE)+".gif").c_str(),".gif");

	delete cCalib;
	delete fitFunc;
	delete calibration_Eta_Measured;
	delete grint;
      }
      
      TFile outFileResult(("outFileResult_"+typeSet.at(ItSetM)).c_str(), "recreate");
      MBVsEta[ItSetM]->Write(("MBVsEta_"+typeSet.at(ItSetM)).c_str());
      outFileResult.Close();
    }
    myfileOu.close();
    //    return;
  
    ////////////  MB stimato per i dati
    int TProfileColor = kBlue;
    int TProfileErr = kBlue;
    float TProfileMarkerStyle = 21;
    float TProfileErrStyle = 3001; //3002;  //3017
    MBVsEta[0]->SetMarkerColor(TProfileColor);
    MBVsEta[0]->SetMarkerStyle(TProfileMarkerStyle);
    MBVsEta_CL[0]->SetFillColor(TProfileErr);
    MBVsEta_CL[0]->SetFillStyle(TProfileErrStyle);
    MBVsEta_Stat[0]->SetFillColor(TProfileErr);
    MBVsEta_Stat[0]->SetFillStyle(TProfileErrStyle);
    MBVsEta_Charge[0]->SetFillColor(TProfileErr);
    MBVsEta_Charge[0]->SetFillStyle(TProfileErrStyle);
    MBVsEta_zEle[0]->SetFillColor(TProfileErr);
    MBVsEta_zEle[0]->SetFillStyle(TProfileErrStyle);
    MBVsEta_eRR[0]->SetFillColor(TProfileErr);
    MBVsEta_eRR[0]->SetFillStyle(TProfileErrStyle);
    
    MBVsEta[2]->SetMarkerColor(kMagenta);
    MBVsEta[2]->SetMarkerStyle(20);
    
    MBVsEta_Closure->SetMarkerColor(kMagenta);
    MBVsEta_Closure->SetMarkerStyle(20);
  
    MBVsEta[3]->SetMarkerColor(kMagenta);
    MBVsEta[3]->SetMarkerStyle(20);
    
    //   MBVsEta_F10S30_C2->SetMarkerColor(kViolet);
    //   MBVsEta_F10S30_C2->SetMarkerStyle(20);
    
    TProfileColor = kGreen;
    TProfileErr = kGreen;
    TProfileMarkerStyle = 20;
    TProfileErrStyle = 3001; 
    MBVsEta[4]->SetMarkerColor(TProfileColor);
    MBVsEta[4]->SetMarkerStyle(TProfileMarkerStyle);
    MBVsEta_CL[4]->SetFillColor(TProfileErr);
    MBVsEta_CL[4]->SetFillStyle(TProfileErrStyle);
    MBVsEta_Stat[4]->SetFillColor(TProfileErr);
    MBVsEta_Stat[4]->SetFillStyle(TProfileErrStyle);
    MBVsEta_Charge[4]->SetFillColor(TProfileErr);
    MBVsEta_Charge[4]->SetFillStyle(TProfileErrStyle);
    MBVsEta_zEle[4]->SetFillColor(TProfileErr);
    MBVsEta_zEle[4]->SetFillStyle(TProfileErrStyle);
    MBVsEta_eRR[4]->SetFillColor(TProfileErr);
    MBVsEta_eRR[4]->SetFillStyle(TProfileErrStyle);


    TProfileColor = kRed;
    TProfileErr = kRed;
    TProfileMarkerStyle = 20;
    TProfileErrStyle = 3001;

    MBVsEta[5]->SetMarkerColor(TProfileColor);
    MBVsEta[5]->SetMarkerStyle(TProfileMarkerStyle);
    MBVsEta_CL[5]->SetFillColor(TProfileErr);
    MBVsEta_CL[5]->SetFillStyle(TProfileErrStyle);
    MBVsEta_Stat[5]->SetFillColor(TProfileErr);
    MBVsEta_Stat[5]->SetFillStyle(TProfileErrStyle);
    MBVsEta_Charge[5]->SetFillColor(TProfileErr);
    MBVsEta_Charge[5]->SetFillStyle(TProfileErrStyle);
    MBVsEta_zEle[5]->SetFillColor(TProfileErr);
    MBVsEta_zEle[5]->SetFillStyle(TProfileErrStyle);
    MBVsEta_eRR[5]->SetFillColor(TProfileErr);
    MBVsEta_eRR[5]->SetFillStyle(TProfileErrStyle);

    
    //   MBVsEta_F20S30_CS->SetMarkerColor(kRed+2);
    //   MBVsEta_F20S30_CS->SetMarkerStyle(20);
  
    
    
    TLegend* legCutsMB = new TLegend(0.20,0.65,0.49,0.85,NULL,"brNDC");
    legCutsMB->SetTextFont(42);
    legCutsMB->SetFillColor(kWhite);
    legCutsMB->SetLineColor(kWhite);
    legCutsMB->SetShadowColor(kWhite);
    
    legCutsMB->AddEntry(MBVsEta[0], (typeSet.at(0)).c_str(),"pl");
    legCutsMB->AddEntry(MBVsEta[4], (typeSet.at(4)).c_str(),"pl");
    legCutsMB->AddEntry(MBVsEta[5], (typeSet.at(5)).c_str(),"pl");
    //  legCutsMBB->AddEntry(MBVsEta_F10, "MC +F10 (calib. with Std. +F20%)","pl");
    //  legCutsMBB->AddEntry(MBVsEta_F10S30_C2, "MC +F10S30 (calib. with Std. +F20%S30%)","pl");
    //  legCutsMBB->AddEntry(MBVsEta_F20S30_CS, "MC +F20%S30% (calib. with Std. F10S30 F20S30)","pl");
    
    
    TLatex* latexLabel = new TLatex();
    latexLabel->SetTextSize(0.04);
    latexLabel->SetTextFont(42);
    latexLabel->SetLineWidth(2);
    latexLabel->SetNDC();
    
    TCanvas* cMaterialBudget = new TCanvas();
    cSingle  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
    cSingle->SetBottomMargin(0.1);  cSingle->SetTopMargin(0.05);
    cSingle->Draw();   cSingle->cd();   gPad->SetGrid();
    MBVsEta[0]->GetXaxis()->SetTitle((type+" #eta").c_str());
    MBVsEta[0]->GetYaxis()->SetTitle("MB_{estim.}/MB_{MC}");
    MBVsEta[0]->GetYaxis()->SetRangeUser(0.9, 1.6);
    MBVsEta[0]->GetXaxis()->SetRangeUser(0., 2.5);

    MBVsEta[0]->DrawCopy();
    //   MBVsEta_D->SetFillColor(kBlue);
    //   MBVsEta_D->SetFillStyle(3001);
    //   MBVsEta_D->Draw("e2same");
  
    /*
    MBVsEta_DM->SetFillColor(kAzure);
    MBVsEta_DM->SetFillStyle(3001); 
    MBVsEta_DM->DrawCopy("e2same");

    MBVsEta_Ch->SetFillColor(kBlue-7);
    MBVsEta_Ch->SetFillStyle(3018);
    MBVsEta_Ch->Draw("e2same");
    
    MBVsEta_Per->SetFillColor(kBlue-7);
    MBVsEta_Per->SetFillStyle(3017);
    MBVsEta_Per->Draw("e2same");
    MBVsEta_D->DrawCopy("same");
    */

    MBVsEta[4]->DrawCopy("same");
    MBVsEta[5]->DrawCopy("same");
    
    //  MBVsEta_F10->DrawCopy("same");
    //  MBVsEta_F10S30_C2->DrawCopy("same");
    //  MBVsEta_F10S30_C2->DrawCopy();
    
    
    MBVsEta_CL[0]->DrawCopy("e2same");
//      MBVsEta_Charge[0]->DrawCopy("same");
//      MBVsEta_zEle[0]->DrawCopy("same");
//      MBVsEta_Stat[0]->DrawCopy("e2same");
//      MBVsEta_eRR[0]->DrawCopy("e2same");

    MBVsEta_CL[4]->DrawCopy("e2same");
//      MBVsEta_Charge[4]->DrawCopy("same");
//      MBVsEta_zEle[4]->DrawCopy("same");
//      MBVsEta_Stat[4]->DrawCopy("e2same");
//      MBVsEta_eRR[4]->DrawCopy("e2same");

    MBVsEta_CL[5]->DrawCopy("e2same");
//      MBVsEta_Charge[5]->DrawCopy("same");
//      MBVsEta_zEle[5]->DrawCopy("same");
//      MBVsEta_Stat[5]->DrawCopy("e2same");
//      MBVsEta_eRR[5]->DrawCopy("e2same");


    MBVsEta[0]->DrawCopy("same");
    MBVsEta[4]->DrawCopy("same");
    MBVsEta[5]->DrawCopy("same");
    
    latexLabel->DrawLatex(0.20, 0.90, "Calibration from MC Std +F10% +F20%");
    legCutsMB->Draw("same");
    cMaterialBudget->Print((plotDirOut+"/MB_Estimated_Bande.png").c_str(),".png");
    cMaterialBudget->Print((plotDirOut+"/MB_Estimated_Bande.pdf").c_str(),".pdf");
    cMaterialBudget->Print((plotDirOut+"/MB_Estimated_Bande.C").c_str(),".C");
  
    TFile perGiacomo("perGiacomo.root", "recreate");
    perGiacomo.cd();
    MBVsEta[0]->Write();
    MBVsEta[4]->Write();
    MBVsEta[5]->Write();
    MBVsEta_CL[0]->Write();
    MBVsEta_CL[4]->Write();
    MBVsEta_CL[5]->Write();
    legCutsMB->Write();
    perGiacomo.Close();


    //MC closure
    TLegend* legCutsC = new TLegend(0.20,0.65,0.49,0.85,NULL,"brNDC");
    legCutsC->SetTextFont(42);
    legCutsC->SetFillColor(kWhite);
    legCutsC->SetLineColor(kWhite);
    legCutsC->SetShadowColor(kWhite);
    
    legCutsC->AddEntry(MBVsEta_Closure, (typeSet.at(2)).c_str(),"pl");
    
    TLatex* latexLabelC = new TLatex();
    latexLabelC->SetTextSize(0.04);
    latexLabelC->SetTextFont(42);
    latexLabelC->SetLineWidth(2);
    latexLabelC->SetNDC();
    
    TCanvas* cClosure = new TCanvas();
    cSingle  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
    cSingle->SetBottomMargin(0.1);  cSingle->SetTopMargin(0.05);
    cSingle->Draw();   cSingle->cd();   gPad->SetGrid();
    MBVsEta_Closure->GetXaxis()->SetTitle((type+" #eta").c_str());
    MBVsEta_Closure->GetYaxis()->SetTitle("MB_{estim.}/MB_{MC}");
    MBVsEta_Closure->GetYaxis()->SetRangeUser(0.9, 1.6);
    MBVsEta_Closure->GetXaxis()->SetRangeUser(0., 2.5);
    MBVsEta_Closure->DrawCopy();
    
    latexLabelC->DrawLatex(0.20, 0.90, "Calibration from MC Std +F20%");
    legCutsC->Draw("same");
    cClosure->Print((plotDirOut+"/MB_Closure.png").c_str(),".png");
    cClosure->Print((plotDirOut+"/MB_Closure.pdf").c_str(),".pdf");
  }
  




  return;

  if(pippo)  {


  TLegend *legCutsMB = new TLegend(0.70,0.80,0.99,1.,NULL,"brNDC");
  legCutsMB->SetTextFont(42);
  legCutsMB->SetFillColor(kWhite);
  legCutsMB->SetLineColor(kWhite);
  legCutsMB->SetShadowColor(kWhite);
  legCutsMB->AddEntry(MBVsEta, "DA 2012ABCD","pl");
  legCutsMB->AddEntry(MBVsEta_F10S30, "MC +F10%S30%","pl");
  legCutsMB->AddEntry(MBVsEta_F20S30, "MC +F20%S30%","pl");
  //   legCutsMB->AddEntry(MBVsEta_F10S30, "MB(F10S30) calib. from Std,F10,F20","pl");
  //   legCutsMB->AddEntry(MBVsEta_F20S30, "MB(F20S30) from Std,F10,F20","pl");
  //   legCutsMB->AddEntry(MBVsEta_F20S30_CS, "MB(F20S30) calib. from Std,F10S30,F20S30","pl");
  
  //  legCutsMBB->AddEntry(MBVsEta_D, "DA 2012ABCD","pl");
  //  legCutsMBB->AddEntry(MBVsEta_F10, "MC +F10 (calib. with Std. +F20%)","pl");                    
  //  legCutsMBB->AddEntry(MBVsEta_F10S30_C2, "MC +F10S30 (calib. with Std. +F20%S30%)","pl");       
  //  legCutsMBB->AddEntry(MBVsEta_F10S30, "MC +F10%S30%","pl");
  //  legCutsMBB->AddEntry(MBVsEta_F20S30, "MC +F20%S30%","pl");
  //  legCutsMBB->AddEntry(MBVsEta_F20S30_CS, "MC +F20%S30% (calib. with Std. F10S30 F20S30)","pl"); 


  TCanvas* cMB_Est = new TCanvas();
  cSingle  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
  cSingle->SetBottomMargin(0.1);  cSingle->SetTopMargin(0.05);
  cSingle->Draw();   cSingle->cd();   gPad->SetGrid();
  MBVsEta->GetXaxis()->SetTitle("sc #eta");
  if(type == "tkEle")  MBVsEta->GetXaxis()->SetTitle("ele #eta");
  MBVsEta->GetYaxis()->SetTitle("MB estimated (%)");
  MBVsEta->GetYaxis()->SetRangeUser(0.9, 1.6);
  MBVsEta->GetXaxis()->SetRangeUser(0., 2.5);
  MBVsEta->Draw();
  MBVsEta_F10S30->Draw("same");
  MBVsEta_F20S30->Draw("same");
  legCutsMB->Draw("same");
  cMB_Est->Print((plotDirOut+"/MB_Est.png").c_str(),".png");
  

  //////////////////////////////////////////////////////////////////////////

  //1)
  TLegend *legCutsMB_1 = new TLegend(0.70,0.80,0.99,1.,NULL,"brNDC");
  legCutsMB_1->SetTextFont(42);
  legCutsMB_1->SetFillColor(kWhite);
  legCutsMB_1->SetLineColor(kWhite);
  legCutsMB_1->SetShadowColor(kWhite);
  legCutsMB_1->AddEntry(MBVsEta_F10S30, "MB(F10S30) calib. from Std,F10,F20","pl");
  legCutsMB_1->AddEntry(MBVsEta_F20S30, "MB(F20S30) from Std,F10,F20","pl");
  legCutsMB_1->AddEntry(MBVsEta_F20S30_CS, "MB(F20S30) calib. from Std,F10S30,F20S30","pl");

  TCanvas* cMB_Est_1 = new TCanvas();
  cSingle  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
  cSingle->SetBottomMargin(0.1);  cSingle->SetTopMargin(0.05);
  cSingle->Draw();   cSingle->cd();   gPad->SetGrid();

  MBVsEta_F10S30->GetXaxis()->SetTitle("sc #eta");
  if(type == "tkEle")  MBVsEta_F10S30->GetXaxis()->SetTitle("ele #eta");
  MBVsEta_F10S30->GetYaxis()->SetTitle("MB estimated (%)");
  MBVsEta_F10S30->GetYaxis()->SetRangeUser(0.9, 1.5);
  MBVsEta_F10S30->GetXaxis()->SetRangeUser(0., 2.5);
  MBVsEta_F10S30->Draw();
  MBVsEta_F20S30->Draw("same");
  MBVsEta_F20S30_CS->Draw("same");
  legCutsMB_1->Draw("same");
  cMB_Est_1->Print((plotDirOut+"/MB_Estimated.pdf").c_str(),".pdf");
  cMB_Est_1->Print((plotDirOut+"/MB_Estimated.png").c_str(),".png");
 

  //2)
  TLegend *legCutsMB_2 = new TLegend(0.70,0.80,0.99,1.,NULL,"brNDC");
  legCutsMB_2->SetTextFont(42);
  legCutsMB_2->SetFillColor(kWhite);
  legCutsMB_2->SetLineColor(kWhite);
  legCutsMB_2->SetShadowColor(kWhite);
  //  legCutsMB->AddEntry(MBVsEta, "DA_MB from F10 F20","pl");
   legCutsMB_2->AddEntry(MBVsEta_F10, "MB(F10) calib. from Std,F20","pl");
   legCutsMB_2->AddEntry(MBVsEta_F10S30_C2, "MB(F10S30) calib. from F20,F20S30","pl");


  TCanvas* cMB_Est_2 = new TCanvas();
  cSingle  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
  cSingle->SetBottomMargin(0.1);  cSingle->SetTopMargin(0.05);
  cSingle->Draw();   cSingle->cd();   gPad->SetGrid();
  MBVsEta_F10->GetXaxis()->SetTitle("sc #eta");
  if(type == "tkEle")  MBVsEta_F10S30->GetXaxis()->SetTitle("ele #eta");
  MBVsEta_F10->GetYaxis()->SetTitle("MB estimated (%)");
  MBVsEta_F10->GetYaxis()->SetRangeUser(0.2, 1.4);
  MBVsEta_F10->GetXaxis()->SetRangeUser(0., 2.5);
  MBVsEta_F10->Draw();
    MBVsEta_F10S30_C2->Draw("same");
  legCutsMB_2->Draw("same");
  cMB_Est_2->Print((plotDirOut+"/MB_Estimated_2.pdf").c_str(),".png");
  cMB_Est_2->Print((plotDirOut+"/MB_Estimated_2.png").c_str(),".png");
  //////////////////////////////////////////////////////////////////////////

  //3)                                                                                        
  TLegend *legCutsMB_3 = new TLegend(0.70,0.80,0.99,1.,NULL,"brNDC");
  legCutsMB_3->SetTextFont(42);
  legCutsMB_3->SetFillColor(kWhite);
  legCutsMB_3->SetLineColor(kWhite);
  legCutsMB_3->SetShadowColor(kWhite);
  legCutsMB_3->AddEntry(MBVsEta_F10S30, "MB(F10S30) calib. from Std,F10,F20","pl");
  legCutsMB_3->AddEntry(MBVsEta_F20S30, "MB(F20S30) from Std,F10,F20","pl");
  legCutsMB_3->AddEntry(MBVsEta_F20S30_CS, "MB(F20S30) calib. from Std,F10S30,F20S30","pl");
  legCutsMB_3->AddEntry(MBVsEta_F10, "MB(F10) calib. from Std,F20","pl");

  TCanvas* cMB_Est_3 = new TCanvas();
  cSingle  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
  cSingle->SetBottomMargin(0.1);  cSingle->SetTopMargin(0.05);
  cSingle->Draw();   cSingle->cd();   gPad->SetGrid();

  MBVsEta_F10S30->GetXaxis()->SetTitle("sc #eta");
  if(type == "tkEle")  MBVsEta_F10S30->GetXaxis()->SetTitle("ele #eta");
  MBVsEta_F10S30->GetYaxis()->SetTitle("MB estimated (%)");
  MBVsEta_F10S30->GetYaxis()->SetRangeUser(0.9, 1.5);
  MBVsEta_F10S30->GetXaxis()->SetRangeUser(0., 2.5);
  MBVsEta_F10S30->Draw();
  MBVsEta_F20S30->Draw("same");
  MBVsEta_F20S30_CS->Draw("same");
  MBVsEta_F10->Draw("same");
  legCutsMB_3->Draw("same");
  cMB_Est_3->Print((plotDirOut+"/MB_Estimated_3.pdf").c_str(),".pdf");
  cMB_Est_1->Print((plotDirOut+"/MB_Estimated_3.png").c_str(),".png");

  //////////////////////////////////////////////////////////////////////////


  TLegend *legCutsMBB = new TLegend(0.20,0.65,0.49,0.85,NULL,"brNDC");
  legCutsMBB->SetTextFont(42);
  legCutsMBB->SetFillColor(kWhite);
  legCutsMBB->SetLineColor(kWhite);
  legCutsMBB->SetShadowColor(kWhite);
  legCutsMBB->AddEntry(MBVsEta_D, "DA 2012ABCD","pl");
  //  legCutsMBB->AddEntry(MBVsEta_F10, "MC +F10 (calib. with Std. +F20%)","pl");
  //  legCutsMBB->AddEntry(MBVsEta_F10S30_C2, "MC +F10S30 (calib. with Std. +F20%S30%)","pl");
  legCutsMBB->AddEntry(MBVsEta_F10S30, "MC +F10%S30%","pl");
  legCutsMBB->AddEntry(MBVsEta_F20S30, "MC +F20%S30%","pl");
  //  legCutsMBB->AddEntry(MBVsEta_F20S30_CS, "MC +F20%S30% (calib. with Std. F10S30 F20S30)","pl");


  TLatex *latexLabel = new TLatex();
  latexLabel->SetTextSize(0.04);
  latexLabel->SetTextFont(42);
  latexLabel->SetLineWidth(2);
  latexLabel->SetNDC();

  TCanvas* cMB_Est_Bande = new TCanvas();
  cSingle  = new TPad("pad_0","pad_0",0.00,0.00,1.00,1.00);
  cSingle->SetBottomMargin(0.1);  cSingle->SetTopMargin(0.05);
  cSingle->Draw();   cSingle->cd();   gPad->SetGrid();
  MBVsEta_D->GetXaxis()->SetTitle("sc #eta");
  if(type == "tkEle")  MBVsEta_D->GetXaxis()->SetTitle("ele #eta");
  MBVsEta_D->GetYaxis()->SetTitle("MB_{estim.}/MB_{MC Std.}");
  MBVsEta_D->GetYaxis()->SetRangeUser(0.9, 1.6);
  MBVsEta_D->GetXaxis()->SetRangeUser(0., 2.5);

  
  MBVsEta_D->DrawCopy();
//   MBVsEta_D->SetFillColor(kBlue);
//   MBVsEta_D->SetFillStyle(3001);
//   MBVsEta_D->Draw("e2same");

  MBVsEta_DM->SetFillColor(kAzure);
  MBVsEta_DM->SetFillStyle(3001); 
  MBVsEta_DM->DrawCopy("e2same");

  MBVsEta_Ch->SetFillColor(kBlue-7);
  MBVsEta_Ch->SetFillStyle(3018);
  MBVsEta_Ch->Draw("e2same");

  MBVsEta_Per->SetFillColor(kBlue-7);
  MBVsEta_Per->SetFillStyle(3017);
  MBVsEta_Per->Draw("e2same");
  MBVsEta_D->DrawCopy("same");

  //  MBVsEta_F10->DrawCopy("same");

  //  MBVsEta_F10S30_C2->DrawCopy("same");
  //  MBVsEta_F10S30_C2->DrawCopy();


  MBVsEta_F10S30->DrawCopy("same");
  MBVsEta_F10S30_Ch->SetFillColor(kGreen-7);  
  MBVsEta_F10S30_Ch->SetFillStyle(3018);
  MBVsEta_F10S30_Ch->Draw("e2same");
  MBVsEta_F10S30_Per->SetFillColor(kGreen-7);  
  MBVsEta_F10S30_Per->SetFillStyle(3017);
  MBVsEta_F10S30_Per->Draw("e2same");

  MBVsEta_F10S30_DM->SetFillColor(kSpring);
  MBVsEta_F10S30_DM->SetFillStyle(3001);
  MBVsEta_F10S30_DM->DrawCopy("e2same");

  MBVsEta_F20S30->DrawCopy("same");
  MBVsEta_F20S30_Ch->SetFillColor(kRed-9);  
  MBVsEta_F20S30_Ch->SetFillStyle(3018);
  MBVsEta_F20S30_Ch->Draw("e2same");
  MBVsEta_F20S30_Per->SetFillColor(kRed+2);  
  MBVsEta_F20S30_Per->SetFillStyle(3017);
  MBVsEta_F20S30_Per->Draw("e2same");
  
  MBVsEta_F20S30_DM->SetFillColor(kPink);
  MBVsEta_F20S30_DM->SetFillStyle(3001);
  MBVsEta_F20S30_DM->DrawCopy("e2same");
  MBVsEta_F20S30->DrawCopy("same");

  //  MBVsEta_F20S30_CS->DrawCopy("same");

  latexLabel->DrawLatex(0.20, 0.90, "Calibration from MC Std +F10% +F20%");
  legCutsMBB->Draw("same");
  cMB_Est_Bande->Print((plotDirOut+"/MB_Estimated_Bande.png").c_str(),".png");
  cMB_Est_Bande->Print((plotDirOut+"/MB_Estimated_Bande.pdf").c_str(),".pdf");
  }

  return;


  PV_z_DA->Rebin(2);
  PV_z_MC_New->Rebin(2);
  PV_z_X0_F10S30->Rebin(2);
  PV_z_X0_F20S30->Rebin(2);
  PV_z_X0_F10->Rebin(2);
  PV_z_X0_F20->Rebin(2);

  PV_z_DA->SetLineColor(kBlack);
  PV_z_DA->SetLineWidth(2);

  PV_z_MC_New->SetLineColor(kCyan);
  PV_z_MC_New->SetLineWidth(2);
  PV_z_MC_New->Scale(1./PV_z_MC_New->Integral()*PV_z_DA->Integral());

  PV_z_X0_F10S30->SetLineColor(kOrange-2);
  PV_z_X0_F10S30->SetLineWidth(2);
  PV_z_X0_F10S30->Scale(1./PV_z_X0_F10S30->Integral()*PV_z_DA->Integral());

  PV_z_X0_F20S30->SetLineColor(kGreen+2);
  PV_z_X0_F20S30->SetLineWidth(2);
  PV_z_X0_F20S30->Scale(1./PV_z_X0_F20S30->Integral()*PV_z_DA->Integral());

  PV_z_X0_F10->SetLineColor(kOrange+8);
  PV_z_X0_F10->SetLineWidth(2);
  PV_z_X0_F10->Scale(1./PV_z_X0_F10->Integral()*PV_z_DA->Integral());

  PV_z_X0_F20->SetLineColor(kGreen+8);
  PV_z_X0_F20->SetLineWidth(2);
  PV_z_X0_F20->Scale(1./PV_z_X0_F20->Integral()*PV_z_DA->Integral());

  TLegend *EasyLegDE = new TLegend(0.70,0.80,0.8,1.,NULL,"brNDC");
  EasyLegDE->SetTextFont(42);
  EasyLegDE->SetFillColor(kWhite);
  EasyLegDE->SetLineColor(kWhite);
  EasyLegDE->SetShadowColor(kWhite);

  EasyLegDE->AddEntry(PV_z_DA, "DA", "l");
  EasyLegDE->AddEntry(PV_z_MC_New, "MC V15", "l");
  EasyLegDE->AddEntry(PV_z_X0_F10S30, "F10S30 V15", "l");
  EasyLegDE->AddEntry(PV_z_X0_F20S30, "F20S30 V15", "l");
  EasyLegDE->AddEntry(PV_z_X0_F10, "F10 V15", "l");
  EasyLegDE->AddEntry(PV_z_X0_F20, "F20 V15", "l");

  TCanvas* cPV_z_DA = new TCanvas();
  gPad->SetGrid();
  PV_z_DA->GetXaxis()->SetTitle("PV_z");
  PV_z_DA->Draw("hist");
  PV_z_MC_New->Draw("hist,same");
  PV_z_X0_F10S30->Draw("hist,same");
  PV_z_X0_F20S30->Draw("hist,same");
  PV_z_X0_F10->Draw("hist,same");
  PV_z_X0_F20->Draw("hist,same");
  EasyLegDE->Draw("same");
  cPV_z_DA->Print((plotDirOut+"/cPV_z_DA.png").c_str(),".png");





  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  h_DEtaIn_vs_bcN_DA->SetLineColor(kBlack);
  h_DEtaIn_vs_bcN_DA->SetLineWidth(2);

  h_DEtaIn_vs_bcN_MC_New->SetLineColor(kRed);
  h_DEtaIn_vs_bcN_MC_New->SetLineWidth(2);
  h_DEtaIn_vs_bcN_MC_New->Scale(1./h_DEtaIn_vs_bcN_MC_New->Integral()*h_DEtaIn_vs_bcN_DA->Integral());

//   h_DEtaIn_vs_bcN_X0_F10S30->SetLineColor(kOrange-2);
//   h_DEtaIn_vs_bcN_X0_F10S30->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_X0_F10S30->Scale(1./h_DEtaIn_vs_bcN_X0_F10S30->Integral()*h_DEtaIn_vs_bcN_DA->Integral());

//   h_DEtaIn_vs_bcN_X0_F20S30->SetLineColor(kGreen+2);
//   h_DEtaIn_vs_bcN_X0_F20S30->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_X0_F20S30->Scale(1./h_DEtaIn_vs_bcN_X0_F20S30->Integral()*h_DEtaIn_vs_bcN_DA->Integral());

//   h_DEtaIn_vs_bcN_X0_F10->SetLineColor(kOrange+8);
//   h_DEtaIn_vs_bcN_X0_F10->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_X0_F10->Scale(1./h_DEtaIn_vs_bcN_X0_F10->Integral()*h_DEtaIn_vs_bcN_DA->Integral());

//   h_DEtaIn_vs_bcN_X0_F20->SetLineColor(kGreen+8);
//   h_DEtaIn_vs_bcN_X0_F20->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_X0_F20->Scale(1./h_DEtaIn_vs_bcN_X0_F20->Integral()*h_DEtaIn_vs_bcN_DA->Integral());

  TLegend *EasyLegDE = new TLegend(0.70,0.80,0.8,1.,NULL,"brNDC");
  EasyLegDE->SetTextFont(42);
  EasyLegDE->SetFillColor(kWhite);
  EasyLegDE->SetLineColor(kWhite);
  EasyLegDE->SetShadowColor(kWhite);

  EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_DA, "DA", "l");
  EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_MC_New, "MC V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_X0_F10S30, "F10S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_X0_F20S30, "F20S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_X0_F10, "F10 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_X0_F20, "F20 V15", "l");

  TCanvas* cDEtaIn_vs_bcN = new TCanvas();
  gPad->SetGrid();
  h_DEtaIn_vs_bcN_DA->GetXaxis()->SetTitle("nDEtaIn_vs_bcN");
  h_DEtaIn_vs_bcN_DA->Draw("hist");
  h_DEtaIn_vs_bcN_MC_New->Draw("hist,same");
//   h_DEtaIn_vs_bcN_X0_F10S30->Draw("hist,same");
//   h_DEtaIn_vs_bcN_X0_F20S30->Draw("hist,same");
//   h_DEtaIn_vs_bcN_X0_F10->Draw("hist,same");
//   h_DEtaIn_vs_bcN_X0_F20->Draw("hist,same");
  EasyLegDE->Draw("same");
  cDEtaIn_vs_bcN->Print((plotDirOut+"/DEtaIn_vs_bcN.png").c_str(),".png");


  ////////////////////
  h_DEtaIn_vs_bcN_EN_DA->SetLineColor(kBlack);
  h_DEtaIn_vs_bcN_EN_DA->SetLineWidth(2);

  h_DEtaIn_vs_bcN_EN_MC_New->SetLineColor(kRed);
  h_DEtaIn_vs_bcN_EN_MC_New->SetLineWidth(2);
  h_DEtaIn_vs_bcN_EN_MC_New->Scale(1./h_DEtaIn_vs_bcN_EN_MC_New->Integral()*h_DEtaIn_vs_bcN_EN_DA->Integral());

//   h_DEtaIn_vs_bcN_EN_X0_F10S30->SetLineColor(kOrange-2);
//   h_DEtaIn_vs_bcN_EN_X0_F10S30->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EN_X0_F10S30->Scale(1./h_DEtaIn_vs_bcN_EN_X0_F10S30->Integral()*h_DEtaIn_vs_bcN_EN_DA->Integral());

//   h_DEtaIn_vs_bcN_EN_X0_F20S30->SetLineColor(kGreen+2);
//   h_DEtaIn_vs_bcN_EN_X0_F20S30->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EN_X0_F20S30->Scale(1./h_DEtaIn_vs_bcN_EN_X0_F20S30->Integral()*h_DEtaIn_vs_bcN_EN_DA->Integral());

//   h_DEtaIn_vs_bcN_EN_X0_F10->SetLineColor(kOrange+8);
//   h_DEtaIn_vs_bcN_EN_X0_F10->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EN_X0_F10->Scale(1./h_DEtaIn_vs_bcN_EN_X0_F10->Integral()*h_DEtaIn_vs_bcN_EN_DA->Integral());

//   h_DEtaIn_vs_bcN_EN_X0_F20->SetLineColor(kGreen+8);
//   h_DEtaIn_vs_bcN_EN_X0_F20->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EN_X0_F20->Scale(1./h_DEtaIn_vs_bcN_EN_X0_F20->Integral()*h_DEtaIn_vs_bcN_EN_DA->Integral());

  TLegend *EasyLegDE = new TLegend(0.70,0.80,0.8,1.,NULL,"brNDC");
  EasyLegDE->SetTextFont(42);
  EasyLegDE->SetFillColor(kWhite);
  EasyLegDE->SetLineColor(kWhite);
  EasyLegDE->SetShadowColor(kWhite);

  EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_DA, "DA", "l");
  EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_MC_New, "MC V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_X0_F10S30, "F10S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_X0_F20S30, "F20S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_X0_F10, "F10 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_X0_F20, "F20 V15", "l");

  TCanvas* cDEtaIn_vs_bcN_EN = new TCanvas();
  gPad->SetGrid();
  h_DEtaIn_vs_bcN_EN_DA->GetXaxis()->SetTitle("nDEtaIn_vs_bcN_EN");
  h_DEtaIn_vs_bcN_EN_DA->Draw("hist");
  h_DEtaIn_vs_bcN_EN_MC_New->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EN_X0_F10S30->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EN_X0_F20S30->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EN_X0_F10->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EN_X0_F20->Draw("hist,same");
  EasyLegDE->Draw("same");
  cDEtaIn_vs_bcN_EN->Print((plotDirOut+"/DEtaIn_vs_bcN_EN.png").c_str(),".png");


  /////////////////////////////////////////////////////////////////////////
  h_DEtaIn_vs_bcN_EP_DA->SetLineColor(kBlack);
  h_DEtaIn_vs_bcN_EP_DA->SetLineWidth(2);

  h_DEtaIn_vs_bcN_EP_MC_New->SetLineColor(kRed);
  h_DEtaIn_vs_bcN_EP_MC_New->SetLineWidth(2);
  h_DEtaIn_vs_bcN_EP_MC_New->Scale(1./h_DEtaIn_vs_bcN_EP_MC_New->Integral()*h_DEtaIn_vs_bcN_EP_DA->Integral());

//   h_DEtaIn_vs_bcN_EP_X0_F10S30->SetLineColor(kOrange-2);
//   h_DEtaIn_vs_bcN_EP_X0_F10S30->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EP_X0_F10S30->Scale(1./h_DEtaIn_vs_bcN_EP_X0_F10S30->Integral()*h_DEtaIn_vs_bcN_EP_DA->Integral());

//   h_DEtaIn_vs_bcN_EP_X0_F20S30->SetLineColor(kGreen+2);
//   h_DEtaIn_vs_bcN_EP_X0_F20S30->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EP_X0_F20S30->Scale(1./h_DEtaIn_vs_bcN_EP_X0_F20S30->Integral()*h_DEtaIn_vs_bcN_EP_DA->Integral());

//   h_DEtaIn_vs_bcN_EP_X0_F10->SetLineColor(kOrange+8);
//   h_DEtaIn_vs_bcN_EP_X0_F10->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EP_X0_F10->Scale(1./h_DEtaIn_vs_bcN_EP_X0_F10->Integral()*h_DEtaIn_vs_bcN_EP_DA->Integral());

//   h_DEtaIn_vs_bcN_EP_X0_F20->SetLineColor(kGreen+8);
//   h_DEtaIn_vs_bcN_EP_X0_F20->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EP_X0_F20->Scale(1./h_DEtaIn_vs_bcN_EP_X0_F20->Integral()*h_DEtaIn_vs_bcN_EP_DA->Integral());

  TLegend *EasyLegDE = new TLegend(0.70,0.80,0.8,1.,NULL,"brNDC");
  EasyLegDE->SetTextFont(42);
  EasyLegDE->SetFillColor(kWhite);
  EasyLegDE->SetLineColor(kWhite);
  EasyLegDE->SetShadowColor(kWhite);

  EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_DA, "DA", "l");
  EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_MC_New, "MC V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_X0_F10S30, "F10S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_X0_F20S30, "F20S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_X0_F10, "F10 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_X0_F20, "F20 V15", "l");

  TCanvas* cDEtaIn_vs_bcN_EP = new TCanvas();
  gPad->SetGrid();
  h_DEtaIn_vs_bcN_EP_DA->GetXaxis()->SetTitle("nDEtaIn_vs_bcN_EP");
  h_DEtaIn_vs_bcN_EP_DA->Draw("hist");
  h_DEtaIn_vs_bcN_EP_MC_New->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EP_X0_F10S30->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EP_X0_F20S30->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EP_X0_F10->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EP_X0_F20->Draw("hist,same");
  EasyLegDE->Draw("same");
  cDEtaIn_vs_bcN_EP->Print((plotDirOut+"/DEtaIn_vs_bcN_EP.png").c_str(),".png");
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  h_DEtaIn_vs_bcN_EN_Vtx_DA->SetLineColor(kBlack);
  h_DEtaIn_vs_bcN_EN_Vtx_DA->SetLineWidth(2);

  h_DEtaIn_vs_bcN_EN_Vtx_MC_New->SetLineColor(kRed);
  h_DEtaIn_vs_bcN_EN_Vtx_MC_New->SetLineWidth(2);
  h_DEtaIn_vs_bcN_EN_Vtx_MC_New->Scale(1./h_DEtaIn_vs_bcN_EN_Vtx_MC_New->Integral()*h_DEtaIn_vs_bcN_EN_Vtx_DA->Integral());

//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F10S30->SetLineColor(kOrange-2);
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F10S30->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F10S30->Scale(1./h_DEtaIn_vs_bcN_EN_Vtx_X0_F10S30->Integral()*h_DEtaIn_vs_bcN_EN_Vtx_DA->Integral());

//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F20S30->SetLineColor(kGreen+2);
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F20S30->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F20S30->Scale(1./h_DEtaIn_vs_bcN_EN_Vtx_X0_F20S30->Integral()*h_DEtaIn_vs_bcN_EN_Vtx_DA->Integral());

//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F10->SetLineColor(kOrange+8);
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F10->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F10->Scale(1./h_DEtaIn_vs_bcN_EN_Vtx_X0_F10->Integral()*h_DEtaIn_vs_bcN_EN_Vtx_DA->Integral());

//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F20->SetLineColor(kGreen+8);
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F20->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F20->Scale(1./h_DEtaIn_vs_bcN_EN_Vtx_X0_F20->Integral()*h_DEtaIn_vs_bcN_EN_Vtx_DA->Integral());

  TLegend *EasyLegDE = new TLegend(0.70,0.80,0.8,1.,NULL,"brNDC");
  EasyLegDE->SetTextFont(42);
  EasyLegDE->SetFillColor(kWhite);
  EasyLegDE->SetLineColor(kWhite);
  EasyLegDE->SetShadowColor(kWhite);

  EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_Vtx_DA, "DA", "l");
  EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_Vtx_MC_New, "MC V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_Vtx_X0_F10S30, "F10S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_Vtx_X0_F20S30, "F20S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_Vtx_X0_F10, "F10 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EN_Vtx_X0_F20, "F20 V15", "l");

  TCanvas* cDEtaIn_vs_bcN_EN_Vtx = new TCanvas();
  gPad->SetGrid();
  h_DEtaIn_vs_bcN_EN_Vtx_DA->GetXaxis()->SetTitle("nDEtaIn_vs_bcN_EN_Vtx");
  h_DEtaIn_vs_bcN_EN_Vtx_DA->Draw("hist");
  h_DEtaIn_vs_bcN_EN_Vtx_MC_New->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F10S30->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F20S30->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F10->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EN_Vtx_X0_F20->Draw("hist,same");
  EasyLegDE->Draw("same");
  cDEtaIn_vs_bcN_EN_Vtx->Print((plotDirOut+"/DEtaIn_vs_bcN_EN_Vtx.png").c_str(),".png");


  /////////////////////////////////////////////////////////////////////////
  h_DEtaIn_vs_bcN_EP_Vtx_DA->SetLineColor(kBlack);
  h_DEtaIn_vs_bcN_EP_Vtx_DA->SetLineWidth(2);

  h_DEtaIn_vs_bcN_EP_Vtx_MC_New->SetLineColor(kRed);
  h_DEtaIn_vs_bcN_EP_Vtx_MC_New->SetLineWidth(2);
  h_DEtaIn_vs_bcN_EP_Vtx_MC_New->Scale(1./h_DEtaIn_vs_bcN_EP_Vtx_MC_New->Integral()*h_DEtaIn_vs_bcN_EP_Vtx_DA->Integral());

//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F10S30->SetLineColor(kOrange-2);
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F10S30->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F10S30->Scale(1./h_DEtaIn_vs_bcN_EP_Vtx_X0_F10S30->Integral()*h_DEtaIn_vs_bcN_EP_Vtx_DA->Integral());

//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F20S30->SetLineColor(kGreen+2);
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F20S30->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F20S30->Scale(1./h_DEtaIn_vs_bcN_EP_Vtx_X0_F20S30->Integral()*h_DEtaIn_vs_bcN_EP_Vtx_DA->Integral());

//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F10->SetLineColor(kOrange+8);
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F10->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F10->Scale(1./h_DEtaIn_vs_bcN_EP_Vtx_X0_F10->Integral()*h_DEtaIn_vs_bcN_EP_Vtx_DA->Integral());

//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F20->SetLineColor(kGreen+8);
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F20->SetLineWidth(2);
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F20->Scale(1./h_DEtaIn_vs_bcN_EP_Vtx_X0_F20->Integral()*h_DEtaIn_vs_bcN_EP_Vtx_DA->Integral());

  TLegend *EasyLegDE = new TLegend(0.70,0.80,0.8,1.,NULL,"brNDC");
  EasyLegDE->SetTextFont(42);
  EasyLegDE->SetFillColor(kWhite);
  EasyLegDE->SetLineColor(kWhite);
  EasyLegDE->SetShadowColor(kWhite);

  EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_Vtx_DA, "DA", "l");
  EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_Vtx_MC_New, "MC V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_Vtx_X0_F10S30, "F10S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_Vtx_X0_F20S30, "F20S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_Vtx_X0_F10, "F10 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_bcN_EP_Vtx_X0_F20, "F20 V15", "l");

  TCanvas* cDEtaIn_vs_bcN_EP_Vtx = new TCanvas();
  gPad->SetGrid();
  h_DEtaIn_vs_bcN_EP_Vtx_DA->GetXaxis()->SetTitle("nDEtaIn_vs_bcN_EP_Vtx");
  h_DEtaIn_vs_bcN_EP_Vtx_DA->Draw("hist");
  h_DEtaIn_vs_bcN_EP_Vtx_MC_New->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F10S30->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F20S30->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F10->Draw("hist,same");
//   h_DEtaIn_vs_bcN_EP_Vtx_X0_F20->Draw("hist,same");
  EasyLegDE->Draw("same");
  cDEtaIn_vs_bcN_EP_Vtx->Print((plotDirOut+"/DEtaIn_vs_bcN_EP_Vtx.png").c_str(),".png");


  /////////////////

  h_DEtaIn_vs_Eta_DA->SetLineColor(kBlack);
  h_DEtaIn_vs_Eta_DA->SetLineWidth(2);

  h_DEtaIn_vs_Eta_MC_New->SetLineColor(kRed);
  h_DEtaIn_vs_Eta_MC_New->SetLineWidth(2);
  h_DEtaIn_vs_Eta_MC_New->Scale(1./h_DEtaIn_vs_Eta_MC_New->Integral()*h_DEtaIn_vs_Eta_DA->Integral());

//   h_DEtaIn_vs_Eta_X0_F10S30->SetLineColor(kOrange-2);
//   h_DEtaIn_vs_Eta_X0_F10S30->SetLineWidth(2);
//   h_DEtaIn_vs_Eta_X0_F10S30->Scale(1./h_DEtaIn_vs_Eta_X0_F10S30->Integral()*h_DEtaIn_vs_Eta_DA->Integral());

//   h_DEtaIn_vs_Eta_X0_F20S30->SetLineColor(kGreen+2);
//   h_DEtaIn_vs_Eta_X0_F20S30->SetLineWidth(2);
//   h_DEtaIn_vs_Eta_X0_F20S30->Scale(1./h_DEtaIn_vs_Eta_X0_F20S30->Integral()*h_DEtaIn_vs_Eta_DA->Integral());

//   h_DEtaIn_vs_Eta_X0_F10->SetLineColor(kOrange+8);
//   h_DEtaIn_vs_Eta_X0_F10->SetLineWidth(2);
//   h_DEtaIn_vs_Eta_X0_F10->Scale(1./h_DEtaIn_vs_Eta_X0_F10->Integral()*h_DEtaIn_vs_Eta_DA->Integral());

//   h_DEtaIn_vs_Eta_X0_F20->SetLineColor(kGreen+8);
//   h_DEtaIn_vs_Eta_X0_F20->SetLineWidth(2);
//   h_DEtaIn_vs_Eta_X0_F20->Scale(1./h_DEtaIn_vs_Eta_X0_F20->Integral()*h_DEtaIn_vs_Eta_DA->Integral());

  TLegend *EasyLegDE = new TLegend(0.70,0.80,0.8,1.,NULL,"brNDC");
  EasyLegDE->SetTextFont(42);
  EasyLegDE->SetFillColor(kWhite);
  EasyLegDE->SetLineColor(kWhite);
  EasyLegDE->SetShadowColor(kWhite);

  EasyLegDE->AddEntry(h_DEtaIn_vs_Eta_DA, "DA", "l");
  EasyLegDE->AddEntry(h_DEtaIn_vs_Eta_MC_New, "MC V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_Eta_X0_F10S30, "F10S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_Eta_X0_F20S30, "F20S30 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_Eta_X0_F10, "F10 V15", "l");
//   EasyLegDE->AddEntry(h_DEtaIn_vs_Eta_X0_F20, "F20 V15", "l");

  TCanvas* cDEtaIn_vs_Eta = new TCanvas();
  gPad->SetGrid();
  h_DEtaIn_vs_Eta_DA->GetXaxis()->SetTitle("nDEtaIn_vs_Eta");
  h_DEtaIn_vs_Eta_DA->Draw("hist");
  h_DEtaIn_vs_Eta_MC_New->Draw("hist,same");
//   h_DEtaIn_vs_Eta_X0_F10S30->Draw("hist,same");
//   h_DEtaIn_vs_Eta_X0_F20S30->Draw("hist,same");
//   h_DEtaIn_vs_Eta_X0_F10->Draw("hist,same");
//   h_DEtaIn_vs_Eta_X0_F20->Draw("hist,same");
  EasyLegDE->Draw("same");
  cDEtaIn_vs_Eta->Print((plotDirOut+"/DEtaIn_vs_Eta.png").c_str(),".png");




  /////////////////////////////////////////////////////////////////////////////////////////
  h_Vtx_DA->SetLineColor(kBlack);
  h_Vtx_DA->SetLineWidth(2);

  h_Vtx_MC_New->SetLineColor(kCyan);
  h_Vtx_MC_New->SetLineWidth(2);
  h_Vtx_MC_New->Scale(1./h_Vtx_MC_New->Integral()*h_Vtx_DA->Integral());

  h_Vtx_X0_F10S30->SetLineColor(kOrange-2);
  h_Vtx_X0_F10S30->SetLineWidth(2);
  h_Vtx_X0_F10S30->Scale(1./h_Vtx_X0_F10S30->Integral()*h_Vtx_DA->Integral());

  h_Vtx_X0_F20S30->SetLineColor(kGreen+2);
  h_Vtx_X0_F20S30->SetLineWidth(2);
  h_Vtx_X0_F20S30->Scale(1./h_Vtx_X0_F20S30->Integral()*h_Vtx_DA->Integral());

  h_Vtx_X0_F10->SetLineColor(kOrange+8);
  h_Vtx_X0_F10->SetLineWidth(2);
  h_Vtx_X0_F10->Scale(1./h_Vtx_X0_F10->Integral()*h_Vtx_DA->Integral());

  h_Vtx_X0_F20->SetLineColor(kGreen+8);
  h_Vtx_X0_F20->SetLineWidth(2);
  h_Vtx_X0_F20->Scale(1./h_Vtx_X0_F20->Integral()*h_Vtx_DA->Integral());

  if(period == "doNewD"){
  h_Vtx_MC_runD->SetLineColor(kBlue);
  h_Vtx_MC_runD->SetLineWidth(2);
  h_Vtx_MC_runD->Scale(1./h_Vtx_MC_runD->Integral()*h_Vtx_DA->Integral());

  h_Vtx_MC_runD_OOTm200->SetLineColor(kGreen);
  h_Vtx_MC_runD_OOTm200->SetLineWidth(2);
  h_Vtx_MC_runD_OOTm200->Scale(1./h_Vtx_MC_runD_OOTm200->Integral()*h_Vtx_DA->Integral());

  h_Vtx_MC_OOTm200->SetLineColor(kRed);
  h_Vtx_MC_OOTm200->SetLineWidth(2);
  h_Vtx_MC_OOTm200->Scale(1./h_Vtx_MC_OOTm200->Integral()*h_Vtx_DA->Integral());
  }

//   h_Vtx_DA->Rebin(R9_rebin);
//   h_Vtx_MC_New->Rebin(R9_rebin);
//   h_Vtx_X0_F10S30->Rebin(R9_rebin);
//   h_Vtx_X0_F20S30->Rebin(R9_rebin);
//   h_Vtx_X0_F10->Rebin(R9_rebin);
//   h_Vtx_X0_F20->Rebin(R9_rebin);
//   if(period == "doNewD"){
//   h_Vtx_MC_runD->Rebin(R9_rebin);
//   h_Vtx_MC_runD_OOTm200->Rebin(R9_rebin);
//   h_Vtx_MC_OOTm200->Rebin(R9_rebin);
//   }

  TLegend *EasyLeg = new TLegend(0.70,0.80,0.8,1.,NULL,"brNDC");
  EasyLeg->SetTextFont(42);
  EasyLeg->SetFillColor(kWhite);
  EasyLeg->SetLineColor(kWhite);
  EasyLeg->SetShadowColor(kWhite);

  EasyLeg->AddEntry(h_Vtx_DA, "DA", "l");
  EasyLeg->AddEntry(h_Vtx_MC_New, "MC V15", "l");
  EasyLeg->AddEntry(h_Vtx_X0_F10S30, "F10S30 V15", "l");
  EasyLeg->AddEntry(h_Vtx_X0_F20S30, "F20S30 V15", "l");
  EasyLeg->AddEntry(h_Vtx_X0_F10, "F10 V15", "l");
  EasyLeg->AddEntry(h_Vtx_X0_F20, "F20 V15", "l");
  if(period == "doNewD"){
  EasyLeg->AddEntry(h_Vtx_MC_runD, "MC V15-runD", "l");
  EasyLeg->AddEntry(h_Vtx_MC_runD_OOTm200, "V15-runD-OOTm200", "l");
  EasyLeg->AddEntry(h_Vtx_MC_OOTm200, "MC V15-OOTm200", "l");
  }

  TCanvas* cVtx = new TCanvas();
  gPad->SetGrid();
  h_Vtx_DA->GetXaxis()->SetTitle("nVtx");
  h_Vtx_DA->Draw("hist");
  h_Vtx_MC_New->Draw("hist,same");
  h_Vtx_X0_F10S30->Draw("hist,same");
  h_Vtx_X0_F20S30->Draw("hist,same");
  h_Vtx_X0_F10->Draw("hist,same");
  h_Vtx_X0_F20->Draw("hist,same");
  if(period == "doNewD"){
    h_Vtx_MC_runD->Draw("hist,same"); 
    h_Vtx_MC_runD_OOTm200->Draw("hist,same");
    h_Vtx_MC_OOTm200->Draw("hist,same");
  }
  EasyLeg->Draw("same");
  cVtx->Print((plotDirOut+"/Vtx.png").c_str(),".png");

//   nTangHits->SetLineColor(kBlue);
//   nTangHits->SetMarkerColor(kBlue);
//   nTangHits->Scale(1./nTangHits->Integral()*nTangHits_DA->Integral());
//   nTangHits_X0->SetLineColor(kRed);
//   nTangHits_X0->SetMarkerColor(kRed);
//   nTangHits_X0->Scale(1./nTangHits_X0->Integral()*nTangHits_DA->Integral());
//   TCanvas* cnTangHits = new TCanvas();
//   nTangHits->GetXaxis()->SetTitle("nTangHits");
//   nTangHits->Draw();
//   nTangHits_X0->Draw("same");
//   nTangHits_DA->Draw("same");
//   EasyLeg->Draw("same");

  h_dxyPV_DA->SetLineWidth(2);
  h_dzPV_DA->SetLineWidth(2);

  h_dxyPV_MC_New->SetLineColor(kCyan);
  h_dxyPV_MC_New->SetLineWidth(2);
  h_dzPV_MC_New->SetLineColor(kCyan);
  h_dzPV_MC_New->SetLineWidth(2);
  h_dxyPV_MC_New->Scale(1./h_dxyPV_MC_New->Integral()*h_dxyPV_DA->Integral());
  h_dzPV_MC_New->Scale(1./h_dzPV_MC_New->Integral()*h_dzPV_DA->Integral());

  h_dxyPV_X0_F10S30->SetLineColor(kOrange-2);
  h_dxyPV_X0_F10S30->SetLineWidth(2);
  h_dzPV_X0_F10S30->SetLineColor(kOrange-2);
  h_dzPV_X0_F10S30->SetLineWidth(2);
  h_dxyPV_X0_F10S30->Scale(1./h_dxyPV_X0_F10S30->Integral()*h_dxyPV_DA->Integral());
  h_dzPV_X0_F10S30->Scale(1./h_dzPV_X0_F10S30->Integral()*h_dzPV_DA->Integral());

  h_dxyPV_X0_F20S30->SetLineColor(kGreen+2);
  h_dxyPV_X0_F20S30->SetLineWidth(2);
  h_dzPV_X0_F20S30->SetLineColor(kGreen+2);
  h_dzPV_X0_F20S30->SetLineWidth(2);
  h_dxyPV_X0_F20S30->Scale(1./h_dxyPV_X0_F20S30->Integral()*h_dxyPV_DA->Integral());
  h_dzPV_X0_F20S30->Scale(1./h_dzPV_X0_F20S30->Integral()*h_dzPV_DA->Integral());

  h_dxyPV_X0_F10->SetLineColor(kOrange+8);
  h_dxyPV_X0_F10->SetLineWidth(2);
  h_dzPV_X0_F10->SetLineColor(kOrange+8);
  h_dzPV_X0_F10->SetLineWidth(2);
  h_dxyPV_X0_F10->Scale(1./h_dxyPV_X0_F10->Integral()*h_dxyPV_DA->Integral());
  h_dzPV_X0_F10->Scale(1./h_dzPV_X0_F10->Integral()*h_dzPV_DA->Integral());

  h_dxyPV_X0_F20->SetLineColor(kGreen+8);
  h_dxyPV_X0_F20->SetLineWidth(2);
  h_dzPV_X0_F20->SetLineColor(kGreen+8);
  h_dzPV_X0_F20->SetLineWidth(2);
  h_dxyPV_X0_F20->Scale(1./h_dxyPV_X0_F20->Integral()*h_dxyPV_DA->Integral());
  h_dzPV_X0_F20->Scale(1./h_dzPV_X0_F20->Integral()*h_dzPV_DA->Integral());

  if(period == "doNewD"){
  h_dxyPV_MC_runD->SetLineColor(kBlue);
  h_dxyPV_MC_runD->SetLineWidth(2);
  h_dxyPV_MC_runD->Scale(1./h_dxyPV_MC_runD->Integral()*h_dxyPV_DA->Integral());

  h_dxyPV_MC_runD_OOTm200->SetLineColor(kGreen);
  h_dxyPV_MC_runD_OOTm200->SetLineWidth(2);
  h_dxyPV_MC_runD_OOTm200->Scale(1./h_dxyPV_MC_runD_OOTm200->Integral()*h_dxyPV_DA->Integral());

  h_dxyPV_MC_OOTm200->SetLineColor(kRed);
  h_dxyPV_MC_OOTm200->SetLineWidth(2);
  h_dxyPV_MC_OOTm200->Scale(1./h_dxyPV_MC_OOTm200->Integral()*h_dxyPV_DA->Integral());

  h_dzPV_MC_runD->SetLineColor(kBlue);
  h_dzPV_MC_runD->SetLineWidth(2);
  h_dzPV_MC_runD->Scale(1./h_dzPV_MC_runD->Integral()*h_dzPV_DA->Integral());

  h_dzPV_MC_runD_OOTm200->SetLineColor(kGreen);
  h_dzPV_MC_runD_OOTm200->SetLineWidth(2);
  h_dzPV_MC_runD_OOTm200->Scale(1./h_dzPV_MC_runD_OOTm200->Integral()*h_dzPV_DA->Integral());

  h_dzPV_MC_OOTm200->SetLineColor(kRed);
  h_dzPV_MC_OOTm200->SetLineWidth(2);
  h_dzPV_MC_OOTm200->Scale(1./h_dzPV_MC_OOTm200->Integral()*h_dzPV_DA->Integral());
  }


  h_dzPV_DA->Rebin(R9_rebin);
  h_dzPV_MC_New->Rebin(R9_rebin);
  h_dzPV_X0_F10S30->Rebin(R9_rebin);
  h_dzPV_X0_F20S30->Rebin(R9_rebin);
  h_dzPV_X0_F10->Rebin(R9_rebin);
  h_dzPV_X0_F20->Rebin(R9_rebin);
  if(period == "doNewD"){
  h_dzPV_MC_runD->Rebin(R9_rebin);
  h_dzPV_MC_runD_OOTm200->Rebin(R9_rebin);
  h_dzPV_MC_OOTm200->Rebin(R9_rebin);
  }

  h_dxyPV_DA->Rebin(R9_rebin);
  h_dxyPV_MC_New->Rebin(R9_rebin);
  h_dxyPV_X0_F10S30->Rebin(R9_rebin);
  h_dxyPV_X0_F20S30->Rebin(R9_rebin);
  h_dxyPV_X0_F10->Rebin(R9_rebin);
  h_dxyPV_X0_F20->Rebin(R9_rebin);
  if(period == "doNewD"){
  h_dxyPV_MC_runD->Rebin(R9_rebin);
  h_dxyPV_MC_runD_OOTm200->Rebin(R9_rebin);
  h_dxyPV_MC_OOTm200->Rebin(R9_rebin);
  }



  TCanvas* cdxyPV = new TCanvas();
  gPad->SetGrid();
  gPad->SetLogy();
  h_dxyPV_DA->GetXaxis()->SetTitle("dxyPV");
  h_dxyPV_DA->Draw("hist");
  h_dxyPV_MC_New->Draw("hist,same");
  h_dxyPV_X0_F10S30->Draw("hist,same");
  h_dxyPV_X0_F20S30->Draw("hist,same");
  h_dxyPV_X0_F10->Draw("hist,same");
  h_dxyPV_X0_F20->Draw("hist,same");
  if(period == "doNewD"){
    h_dxyPV_MC_runD->Draw("hist,same"); 
    h_dxyPV_MC_runD_OOTm200->Draw("hist,same");
    h_dxyPV_MC_OOTm200->Draw("hist,same");
  }
  EasyLeg->Draw("same");
  cdxyPV->Print((plotDirOut+"/dxyPV.png").c_str(),".png");

  TCanvas* cdzPV = new TCanvas();
  gPad->SetGrid();
  gPad->SetLogy();
  h_dzPV_DA->GetXaxis()->SetTitle("dzPV");
  h_dzPV_DA->Draw("hist");
  h_dzPV_MC_New->Draw("hist,same");
  h_dzPV_X0_F10S30->Draw("hist,same");
  h_dzPV_X0_F20S30->Draw("hist,same");
  h_dzPV_X0_F10->Draw("hist,same");
  h_dzPV_X0_F20->Draw("hist,same");
  if(period == "doNewD"){
    h_dzPV_MC_runD->Draw("hist,same"); 
    h_dzPV_MC_runD_OOTm200->Draw("hist,same");
    h_dzPV_MC_OOTm200->Draw("hist,same");
  }
  EasyLeg->Draw("same");
  cdzPV->Print((plotDirOut+"/dzPV.png").c_str(),".png");

  ///////////////////////////////////////////////////////////////////////
  DP_DA->SetLineWidth(2);

  DP_MC_New->SetLineColor(kCyan);
  DP_MC_New->SetLineWidth(2);
  DP_MC_New->Scale(1./DP_MC_New->Integral()*DP_DA->Integral());

  DP_X0_F10S30->SetLineColor(kOrange-2);
  DP_X0_F10S30->SetLineWidth(2);
  DP_X0_F10S30->Scale(1./DP_X0_F10S30->Integral()*DP_DA->Integral());

  DP_X0_F20S30->SetLineColor(kGreen+2);
  DP_X0_F20S30->SetLineWidth(2);
  DP_X0_F20S30->Scale(1./DP_X0_F20S30->Integral()*DP_DA->Integral());

  DP_X0_F10->SetLineColor(kOrange+8);
  DP_X0_F10->SetLineWidth(2);
  DP_X0_F10->Scale(1./DP_X0_F10->Integral()*DP_DA->Integral());

  DP_X0_F20->SetLineColor(kGreen+8);
  DP_X0_F20->SetLineWidth(2);
  DP_X0_F20->Scale(1./DP_X0_F20->Integral()*DP_DA->Integral());

  TCanvas* cDP = new TCanvas();
  gPad->SetGrid();
  DP_DA->GetXaxis()->SetTitle("DP_{tg}");
  DP_DA->Draw("hist");
  DP_MC_New->Draw("hist,same");
  DP_X0_F10S30->Draw("hist,same");
  DP_X0_F20S30->Draw("hist,same");
  DP_X0_F10->Draw("hist,same");
  DP_X0_F20->Draw("hist,same");
  EasyLeg->Draw("same");
  cDP->Print((plotDirOut+"/DP.png").c_str(),".png");
  ///////////////////////////////////////////////////////////////////////

  DPVsEta_DA->SetLineColor(kBlack);
  
  DPVsEta_MC_New->SetLineColor(kRed);
  DPVsEta_MC_New->Scale(1./DPVsEta_MC_New->Integral()*DPVsEta_DA->Integral());

  DPVsEta_X0_F10S30->SetLineColor(kOrange-2);
  DPVsEta_X0_F10S30->Scale(1./DPVsEta_X0_F10S30->Integral()*DPVsEta_DA->Integral());

  DPVsEta_X0_F20S30->SetLineColor(kGreen+2);
  DPVsEta_X0_F20S30->Scale(1./DPVsEta_X0_F20S30->Integral()*DPVsEta_DA->Integral());

  DPVsEta_X0_F10->SetLineColor(kOrange+8);
  DPVsEta_X0_F10->Scale(1./DPVsEta_X0_F10->Integral()*DPVsEta_DA->Integral());

  DPVsEta_X0_F20->SetLineColor(kGreen+8);
  DPVsEta_X0_F20->Scale(1./DPVsEta_X0_F20->Integral()*DPVsEta_DA->Integral());

  TCanvas* cDPVsEta = new TCanvas();
  gPad->SetGrid();
  DPVsEta_DA->GetXaxis()->SetTitle("DPVsEta_{tg}");
  DPVsEta_DA->Draw("hist");
  DPVsEta_MC_New->Draw("hist, same");
  DPVsEta_X0_F10S30->Draw("hist, same");
  DPVsEta_X0_F20S30->Draw("hist, same");
  DPVsEta_X0_F10->Draw("hist, same");
  DPVsEta_X0_F20->Draw("hist, same");
  EasyLeg->Draw("same");
  cDPVsEta->Print((plotDirOut+"/DPVsEta.png").c_str(),".png");

  ///////////////////////////////////////////////////////////////////////
  h_R9_DA->SetLineWidth(2);

  h_R9_MC_New->SetLineColor(kRed);
  h_R9_MC_New->SetLineWidth(2);
  h_R9_MC_New->Scale(1./h_R9_MC_New->Integral()*h_R9_DA->Integral());

  h_R9_X0_F10S30->SetLineColor(kOrange-2);
  h_R9_X0_F10S30->SetLineWidth(2);
  h_R9_X0_F10S30->Scale(1./h_R9_X0_F10S30->Integral()*h_R9_DA->Integral());

  h_R9_X0_F20S30->SetLineColor(kGreen+2);
  h_R9_X0_F20S30->SetLineWidth(2);
  h_R9_X0_F20S30->Scale(1./h_R9_X0_F20S30->Integral()*h_R9_DA->Integral());

  h_R9_X0_F10->SetLineColor(kOrange+8);
  h_R9_X0_F10->SetLineWidth(2);
  h_R9_X0_F10->Scale(1./h_R9_X0_F10->Integral()*h_R9_DA->Integral());

  h_R9_X0_F20->SetLineColor(kGreen+8);
  h_R9_X0_F20->SetLineWidth(2);
  h_R9_X0_F20->Scale(1./h_R9_X0_F20->Integral()*h_R9_DA->Integral());


  h_R9_DA->Rebin(R9_rebin);
  h_R9_MC_New->Rebin(R9_rebin);
  h_R9_X0_F10S30->Rebin(R9_rebin);
  h_R9_X0_F20S30->Rebin(R9_rebin);
  h_R9_X0_F10->Rebin(R9_rebin);
  h_R9_X0_F20->Rebin(R9_rebin);

  TLegend *EasyLeg2 = new TLegend(0.70,0.80,0.8,1.,NULL,"brNDC");
  EasyLeg2->SetTextFont(42);
  EasyLeg2->SetFillColor(kWhite);
  EasyLeg2->SetLineColor(kWhite);
  EasyLeg2->SetShadowColor(kWhite);

  EasyLeg2->AddEntry(h_Vtx_DA, "DA", "l");
  EasyLeg2->AddEntry(h_Vtx_MC_New, "MC V15", "l");
  EasyLeg2->AddEntry(h_Vtx_X0_F10S30, "F10S30 V15", "l");
  EasyLeg2->AddEntry(h_Vtx_X0_F20S30, "F20S30 V15", "l");
  EasyLeg2->AddEntry(h_Vtx_X0_F10, "F10 V15", "l");
  EasyLeg2->AddEntry(h_Vtx_X0_F20, "F20 V15", "l");

  if(period == "doNewD"){
  h_R9_MC_runD->SetLineColor(kBlue);
  h_R9_MC_runD->SetLineWidth(2);
  h_R9_MC_runD->Scale(1./h_R9_MC_runD->Integral()*h_R9_DA->Integral());

  h_R9_MC_runD_OOTm200->SetLineColor(kGreen);
  h_R9_MC_runD_OOTm200->SetLineWidth(2);
  h_R9_MC_runD_OOTm200->Scale(1./h_R9_MC_runD_OOTm200->Integral()*h_R9_DA->Integral());

  h_R9_MC_OOTm200->SetLineColor(kRed);
  h_R9_MC_OOTm200->SetLineWidth(2);
  h_R9_MC_OOTm200->Scale(1./h_R9_MC_OOTm200->Integral()*h_R9_DA->Integral());

  h_R9_MC_runD->Rebin(R9_rebin);
  h_R9_MC_runD_OOTm200->Rebin(R9_rebin);
  h_R9_MC_OOTm200->Rebin(R9_rebin);

  EasyLeg2->AddEntry(h_R9_MC_runD, "MC V15-runD", "l");
  EasyLeg2->AddEntry(h_R9_MC_runD_OOTm200, "V15-runD-OOTm200", "l");
  EasyLeg2->AddEntry(h_R9_MC_OOTm200, "MC V15-OOTm200", "l");
  }

  TCanvas* ch_R9 = new TCanvas();
  gPad->SetGrid();
  h_R9_DA->GetXaxis()->SetTitle("h_R9_{tg}");
  h_R9_DA->GetXaxis()->SetRangeUser(0.4, 1.01);
  h_R9_DA->Draw("hist");
  h_R9_MC_New->Draw("hist,same");
  h_R9_X0_F10S30->Draw("hist,same");
  h_R9_X0_F20S30->Draw("hist,same");
  h_R9_X0_F10->Draw("hist,same");
  h_R9_X0_F20->Draw("hist,same");
  if(period == "doNewD"){
  h_R9_MC_runD->Draw("hist,same");
  h_R9_MC_runD_OOTm200->Draw("hist,same");
  h_R9_MC_OOTm200->Draw("hist,same");
  }
  EasyLeg2->Draw("same");
  ch_R9->Print((plotDirOut+"/R9.png").c_str(),".png");
  ///////////////////////////////////////////////////////////////////////

//   h_pT_DA->SetLineWidth(2);
//   //  return;
//   h_pT_MC_New->SetLineColor(kRed);
//   h_pT_MC_New->SetLineWidth(2);
//   h_pT_MC_New->Scale(1./h_pT_MC_New->Integral()*h_pT_DA->Integral());

//   h_pT_X0_F10S30->SetLineColor(kOrange-2);
//   h_pT_X0_F10S30->SetLineWidth(2);
//   h_pT_X0_F10S30->Scale(1./h_pT_X0_F10S30->Integral()*h_pT_DA->Integral());

//   h_pT_X0_F20S30->SetLineColor(kGreen+2);
//   h_pT_X0_F20S30->SetLineWidth(2);
//   h_pT_X0_F20S30->Scale(1./h_pT_X0_F20S30->Integral()*h_pT_DA->Integral());

//   h_pT_X0_F10->SetLineColor(kOrange+8);
//   h_pT_X0_F10->SetLineWidth(2);
//   h_pT_X0_F10->Scale(1./h_pT_X0_F10->Integral()*h_pT_DA->Integral());

//   h_pT_X0_F20->SetLineColor(kGreen+8);
//   h_pT_X0_F20->SetLineWidth(2);
//   h_pT_X0_F20->Scale(1./h_pT_X0_F20->Integral()*h_pT_DA->Integral());


// //   h_pT_DA->Rebin(pT_rebin);
// //   h_pT_MC_New->Rebin(pT_rebin);
// //   h_pT_X0_F10S30->Rebin(pT_rebin);
// //   h_pT_X0_F20S30->Rebin(pT_rebin);
// //   h_pT_X0_F10->Rebin(pT_rebin);
// //   h_pT_X0_F20->Rebin(pT_rebin);

//   TLegend *EasyLeg2 = new TLegend(0.70,0.80,0.8,1.,NULL,"brNDC");
//   EasyLeg2->SetTextFont(42);
//   EasyLeg2->SetFillColor(kWhite);
//   EasyLeg2->SetLineColor(kWhite);
//   EasyLeg2->SetShadowColor(kWhite);

//   EasyLeg2->AddEntry(h_Vtx_DA, "DA", "l");
//   EasyLeg2->AddEntry(h_Vtx_MC_New, "MC V15", "l");
//   EasyLeg2->AddEntry(h_Vtx_X0_F10S30, "F10S30 V15", "l");
//   EasyLeg2->AddEntry(h_Vtx_X0_F20S30, "F20S30 V15", "l");
//   EasyLeg2->AddEntry(h_Vtx_X0_F10, "F10 V15", "l");
//   EasyLeg2->AddEntry(h_Vtx_X0_F20, "F20 V15", "l");


//   TCanvas* ch_pT = new TCanvas();
//   gPad->SetGrid();
//   h_pT_DA->GetXaxis()->SetTitle("h_pT_{tg}");
//   //  h_pT_DA->GetXaxis()->SetRangeUser(0.4, 1.01);
//   h_pT_DA->Draw("hist");
//   h_pT_MC_New->Draw("hist,same");
//   h_pT_X0_F10S30->Draw("hist,same");
//   h_pT_X0_F20S30->Draw("hist,same");
//   h_pT_X0_F10->Draw("hist,same");
//   h_pT_X0_F20->Draw("hist,same");
//   EasyLeg2->Draw("same");
//   ch_pT->Print((plotDirOut+"/pT.png").c_str(),".png");


 return;

  for(int fBremSubl_It=0; fBremSubl_It<27; ++fBremSubl_It){
    fBremSubl_DA[fBremSubl_It]->SetLineColor(kBlack);
    //    fBremSubl_DA[fBremSubl_It]->SetMarkerColor(kBlack);
    fBremSubl_DA[fBremSubl_It]->SetLineWidth(2);
    //    fBremSubl_DA[fBremSubl_It]->SetMarkerStyle(7);

    fBremSubl_MC_New[fBremSubl_It]->SetLineColor(kCyan);
    //    fBremSubl_MC_New[fBremSubl_It]->SetMarkerColor(kCyan);
    fBremSubl_MC_New[fBremSubl_It]->SetLineWidth(2);
    //    fBremSubl_MC_New[fBremSubl_It]->SetMarkerStyle(7);

    fBremSubl_X0_F10S30[fBremSubl_It]->SetLineColor(kOrange-2);
    //    fBremSubl_X0_F10S30[fBremSubl_It]->SetMarkerColor(kOrange-2);
    fBremSubl_X0_F10S30[fBremSubl_It]->SetLineWidth(2);
    //    fBremSubl_X0_F10S30[fBremSubl_It]->SetMarkerStyle(7);

    fBremSubl_X0_F20S30[fBremSubl_It]->SetLineColor(kGreen+2);
    //    fBremSubl_X0_F20S30[fBremSubl_It]->SetMarkerColor(kGreen+2);
    fBremSubl_X0_F20S30[fBremSubl_It]->SetLineWidth(2);
    //    fBremSubl_X0_F20S30[fBremSubl_It]->SetMarkerStyle(7);

    fBremSubl_X0_F10[fBremSubl_It]->SetLineColor(kOrange+8);
    //    fBremSubl_X0_F10[fBremSubl_It]->SetMarkerColor(kOrange+8);
    fBremSubl_X0_F10[fBremSubl_It]->SetLineWidth(2);
    //    fBremSubl_X0_F10[fBremSubl_It]->SetMarkerStyle(7);

    fBremSubl_X0_F20[fBremSubl_It]->SetLineColor(kGreen+8);
    //    fBremSubl_X0_F20[fBremSubl_It]->SetMarkerColor(kGreen+8);
    fBremSubl_X0_F20[fBremSubl_It]->SetLineWidth(2);
    //    fBremSubl_X0_F20[fBremSubl_It]->SetMarkerStyle(7);
  }

  TCanvas** cfBrem_vsEta_layers_BP = new TCanvas*[3];
  for(int TC_it=0; TC_it<3; ++TC_it){
    cfBrem_vsEta_layers_BP[TC_it] = new TCanvas();
    cfBrem_vsEta_layers_BP[TC_it]->cd();
    gPad->SetGrid();

    char nome[50];
    sprintf (nome, "BP%d", TC_it+1);
    std::string etaBin = std::string(nome);
    
    fBremSubl_DA[TC_it]->GetXaxis()->SetTitle((etaBin).c_str());
    //    fBremSubl_DA[TC_it]->GetYaxis()->SetRangeUser(-15., 2.);
    //  fBremSubl_isBP2->GetXaxis()->SetRangeUser(-2.5, 2.5);
    fBremSubl_DA[TC_it]->Draw("hist");
    fBremSubl_MC_New[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20[TC_it]->Draw("hist, same");
    cfBrem_vsEta_layers_BP[TC_it]->Print((plotDirOut+"/cfBrem_vsEta_layers_"+etaBin+".png").c_str(), "png");
  }

  TCanvas** cfBrem_vsEta_layers_FP = new TCanvas*[2];
  for(int TC_it=3; TC_it<5; ++TC_it){
    cfBrem_vsEta_layers_FP[TC_it] = new TCanvas();
    cfBrem_vsEta_layers_FP[TC_it]->cd();
    gPad->SetGrid();

    char nome[50];
    sprintf (nome, "FP%d", TC_it-2);
    std::string etaBin = std::string(nome);
    
    fBremSubl_DA[TC_it]->GetXaxis()->SetTitle((etaBin).c_str());
    //    fBremSubl_DA[TC_it]->GetYaxis()->SetRangeUser(-15., 2.);
    //  fBremSubl_isBP2->GetXaxis()->SetRangeUser(-2.5, 2.5);
    fBremSubl_DA[TC_it]->Draw("hist");
    fBremSubl_MC_New[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20[TC_it]->Draw("hist, same");
    cfBrem_vsEta_layers_FP[TC_it]->Print((plotDirOut+"/cfBrem_vsEta_layers_"+etaBin+".png").c_str(), "png");
  }

  TCanvas** cfBrem_vsEta_layers_TIB = new TCanvas*[4];
  for(int TC_it=5; TC_it<9; ++TC_it){
    cfBrem_vsEta_layers_TIB[TC_it] = new TCanvas();
    cfBrem_vsEta_layers_TIB[TC_it]->cd();
    gPad->SetGrid();

    char nome[50];
    sprintf (nome, "TIB%d", TC_it-4);
    std::string etaBin = std::string(nome);
    
    fBremSubl_DA[TC_it]->GetXaxis()->SetTitle((etaBin).c_str());
    //    fBremSubl_DA[TC_it]->GetYaxis()->SetRangeUser(-15., 2.);
    //  fBremSubl_isBP2->GetXaxis()->SetRangeUser(-2.5, 2.5);
    fBremSubl_DA[TC_it]->Draw("hist");
    fBremSubl_MC_New[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20[TC_it]->Draw("hist, same");
    cfBrem_vsEta_layers_TIB[TC_it]->Print((plotDirOut+"/cfBrem_vsEta_layers_"+etaBin+".png").c_str(), "png");
  }

  TCanvas** cfBrem_vsEta_layers_TID = new TCanvas*[3];
  for(int TC_it=9; TC_it<12; ++TC_it){
    cfBrem_vsEta_layers_TID[TC_it] = new TCanvas();
    cfBrem_vsEta_layers_TID[TC_it]->cd();
    gPad->SetGrid();

    char nome[50];
    sprintf (nome, "TID%d", TC_it-8);
    std::string etaBin = std::string(nome);
    
    fBremSubl_DA[TC_it]->GetXaxis()->SetTitle((etaBin).c_str());
    //    fBremSubl_DA[TC_it]->GetYaxis()->SetRangeUser(-15., 2.);
    //  fBremSubl_isBP2->GetXaxis()->SetRangeUser(-2.5, 2.5);
    fBremSubl_DA[TC_it]->Draw("hist");
    fBremSubl_MC_New[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20[TC_it]->Draw("hist, same");
    cfBrem_vsEta_layers_TID[TC_it]->Print((plotDirOut+"/cfBrem_vsEta_layers_"+etaBin+".png").c_str(), "png");
  }

  TCanvas** cfBrem_vsEta_layers_TOB = new TCanvas*[6];
  for(int TC_it=12; TC_it<18; ++TC_it){
    cfBrem_vsEta_layers_TOB[TC_it] = new TCanvas();
    cfBrem_vsEta_layers_TOB[TC_it]->cd();
    gPad->SetGrid();

    char nome[50];
    sprintf (nome, "TOB%d", TC_it-11);
    std::string etaBin = std::string(nome);
    
    fBremSubl_DA[TC_it]->GetXaxis()->SetTitle((etaBin).c_str());
    //    fBremSubl_DA[TC_it]->GetYaxis()->SetRangeUser(-15., 2.);
    //  fBremSubl_isBP2->GetXaxis()->SetRangeUser(-2.5, 2.5);
    fBremSubl_DA[TC_it]->Draw("hist");
    fBremSubl_MC_New[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20[TC_it]->Draw("hist, same");
    cfBrem_vsEta_layers_TOB[TC_it]->Print((plotDirOut+"/cfBrem_vsEta_layers_"+etaBin+".png").c_str(), "png");
  }

  TCanvas** cfBrem_vsEta_layers_TEC = new TCanvas*[9];
  for(int TC_it=18; TC_it<27; ++TC_it){
    cfBrem_vsEta_layers_TEC[TC_it] = new TCanvas();
    cfBrem_vsEta_layers_TEC[TC_it]->cd();
    gPad->SetGrid();

    char nome[50];
    sprintf (nome, "TEC%d", TC_it-17);
    std::string etaBin = std::string(nome);
    
    fBremSubl_DA[TC_it]->GetXaxis()->SetTitle((etaBin).c_str());
    //    fBremSubl_DA[TC_it]->GetYaxis()->SetRangeUser(-15., 2.);
    //  fBremSubl_isBP2->GetXaxis()->SetRangeUser(-2.5, 2.5);
    fBremSubl_DA[TC_it]->Draw("hist");
    fBremSubl_MC_New[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20S30[TC_it]->Draw("hist, same");
    fBremSubl_X0_F10[TC_it]->Draw("hist, same");
    fBremSubl_X0_F20[TC_it]->Draw("hist, same");
    cfBrem_vsEta_layers_TEC[TC_it]->Print((plotDirOut+"/cfBrem_vsEta_layers_"+etaBin+".png").c_str(), "png");
  }



  return;




  std::cout << " CIAO " << std::endl;
}
