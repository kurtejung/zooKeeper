
// drawZoo.C
// Kurt Jung, UIC
// v1, Sept. 2016

// List-based method of drawing the zoo plots
// Simply make a list of plots you want, like:
// *************************
// HadronRAA  <-tab->  InclJetRAA
// BJpsiRAA   <-tab->  BJetRAA
// *************************
// and it will be drawn (either one or two panel)
// If you have more objects to plot on the right side than on the left, you can organize your input file like:
// *************************
// HadronRAA  <-tab->  InclJetRAA
// BJpsiRAA   <-tab->  BJetRAA <-tab-> BMesonRpA
// *************************
// Only the objects without a tab or space in front will be drawn on the left side
// Objects are also FIFO - the earlier you ask for them in your input file, the earlier they are drawn

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "TROOT.h"
#include "TStyle.h"

#include "TFile.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TH1.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"

using namespace std;
const int nObjects = 13;
//keywords for input list
const string objectsToPlot[nObjects] = {"HadronRAA", "HadronRpA", "InclJetRpA", "InclJetRAA", "BJetRAA", "PhotonRAA", "ZRpA", "WRAA", "BJpsiRAA", "BMesonRpA", "BJetRpA", "CJetRpA", "ZRAA"};


//**********************  USER PARAMETERS ******************************//

//number of panels are determined automatically by the input list

const string inputList = "plot.txt";

//set x-axis for left plot
const double xMax = 100;

//set x-axis maximum for right plot;
const double xMax2 = 400;

//set maximum for y-axis value
const double yMax = 3;

//split the legend into two pieces (only for 2-panel plots!)
const bool splitLegend = true;

//draw the global uncertainties?
const bool drawLumiUnc = false;
const bool drawGlauberUnc = false;

//redraw the points on top of all the systematic uncertainties?
const bool redrawPoints = false;

//**********************  END USER PARAMETERS ******************************//


void drawZoo(){

  bool debug = false;

  ifstream instr(inputList);
  if(!instr){ cout << "ERROR! Cannot find " << inputList << endl; exit(0); }
  vector<string> leftPlot;
  vector<string> rightPlot;
  string obj;
  while(instr >> obj){
    leftPlot.push_back(obj);
    char ch;
    while((ch = instr.peek()) != std::char_traits<char>::eof()){
      if(ch == '\t' || ch == ' '){
         ch = instr.get();
         instr >> obj;
         rightPlot.push_back(obj);
      }
      if(ch == '\n') break;
    }
  }
  
  for(unsigned int i=0; i<leftPlot.size(); i++){
    if(debug) cout << "left side: "<< leftPlot.at(i) << endl;
  }
  for(unsigned int i=0; i<rightPlot.size(); i++){
    if(debug) cout << "right side: "<< rightPlot.at(i) << endl;
  }

  TGraphAsymmErrors *systUnc[nObjects];
  TGraphErrors *plots[nObjects];
  int runPeriod[nObjects]; // 1 = 2.76 TeV PbPb, 2 = 5.02 TeV pPb, 3 = 5.02 TeV PbPb
  bool isPrelim[nObjects];
  for(int i=0; i<nObjects; i++){
    isPrelim[i] = false;
  }


  double mymarkersize = 1.;
  
  //----------------------------------------- charged hadrons (RAA)
  double ptBins_h[27]={0};
  double ptError_h[27]={0};
  
  double raa_h[27];
  double raaStat_h[27];
  double raaSyst_h[27];

  double ptSystXlow_h[27];
  double ptSystXhigh_h[27];

  // reading numbers
  ifstream inData_h; 
  TString inFile("raa_h05.dat");
  inData_h.open(inFile);
  if(inData_h.fail()) {
    cerr << "unable to open file raa_h05.dat for reading" << endl;
    exit(1);
  }
  Int_t j=0;
  Double_t xx_low, xx_high,raa,raa_syst, raa_stat, raa_stat2;
  while(!inData_h.eof())
  { 
    inData_h >> xx_low >> xx_high >> raa >> raa_syst >> raa_stat;
    ptBins_h[j]   = xx_low+ (xx_high-xx_low)/2;
    raa_h[j]      = raa;
    raaStat_h[j]  = raa_syst;
    raaSyst_h[j]  = raa_stat;
    
    ptSystXlow_h[j] = (xx_high-xx_low)/2;
    ptSystXhigh_h[j]= (xx_high-xx_low)/2;
    // cout<<"pT"<<ptBins_h[j]<<"\t raa= "<<raa_h[j]<<"\t syst= "<<raaStat_h[j]<<"\t stat ="<<raaSyst_h[j]<<"\t x_low= "<<ptSystXlow_h[j]<<"\t high= "<<ptSystXhigh_h[j]<<endl;
    j++;
  }     
  inData_h.close();
  // done reading numebrs
  
  TGraphErrors *pgRaa_h          = new TGraphErrors(27, ptBins_h, raa_h, ptError_h, raaStat_h);
  TGraphAsymmErrors *pgRaaSyst_h = new TGraphAsymmErrors(27, ptBins_h, raa_h, ptSystXlow_h,ptSystXhigh_h,raaSyst_h,raaSyst_h);
  pgRaa_h->SetName("pgRaa_h");
  pgRaa_h->SetMarkerStyle(20);
  pgRaa_h->SetMarkerSize(1.);

  plots[0] = new TGraphErrors(*pgRaa_h);
  
  //systm error
  pgRaaSyst_h->SetName("pgRaaSyst_h");
  pgRaaSyst_h->SetFillColor(TColor::GetColor("#33ccff"));
  pgRaaSyst_h->SetFillStyle(1001);

  systUnc[0] = new TGraphAsymmErrors(*pgRaaSyst_h);
  systUnc[0]->SetTitle("Charged particle R_{AA}  (0-5%)  |#eta| < 1");
  runPeriod[0] = 1;


//----------------------------------------- charged hadrons (RpA)
  double ptBins_hp[33]={0};
  double ptError_hp[33]={0};

  double raa_hp[33];
  double raaStat_hp[33];
  double raaSyst_hp[33];

  double ptSystXlow_hp[33];
  double ptSystXhigh_hp[33];

  // reading numbers
  ifstream inData_hp;
  TString inFilep("rpa_h.dat");
  inData_hp.open(inFilep);
  if(inData_hp.fail()) {
    cerr << "unable to open file rpa_h.dat for reading" << endl;
    exit(1);
  }
  j=0;
  while(!inData_hp.eof())
  {
    inData_hp >> xx_low >> xx_high >> raa >> raa_syst >> raa_stat >> raa_stat2;
    ptBins_hp[j]   = xx_low+ (xx_high-xx_low)/2;
    raa_hp[j]      = raa;
    raaStat_hp[j]  = raa_syst;
    raaSyst_hp[j]  = sqrt(raa_stat*raa_stat + raa_stat2*raa_stat2);

    ptSystXlow_hp[j] = (xx_high-xx_low)/2;
    ptSystXhigh_hp[j]= (xx_high-xx_low)/2;
    // cout<<"pT"<<ptBins_h[j]<<"\t raa= "<<raa_h[j]<<"\t syst= "<<raaStat_h[j]<<"\t stat ="<<raaSyst_h[j]<<"\t x_low= "<<ptSystXlow_h[j]<<"\t high= "<<ptSystXhigh_h[j]<<endl;
    j++;
  }
  inData_hp.close();
  // done reading numebrs

  TGraphErrors *pgRaa_hp          = new TGraphErrors(33, ptBins_hp, raa_hp, ptError_hp, raaStat_hp);
  TGraphAsymmErrors *pgRaaSyst_hp = new TGraphAsymmErrors(33, ptBins_hp, raa_hp, ptSystXlow_hp,ptSystXhigh_hp,raaSyst_hp,raaSyst_hp);
  pgRaa_hp->SetName("pgRaa_hp");
  pgRaa_hp->SetMarkerStyle(21);
  pgRaa_hp->SetMarkerSize(1.);

  plots[1] = new TGraphErrors(*pgRaa_hp);

  //systm error
  pgRaaSyst_hp->SetName("pgRaaSyst_hp");
  pgRaaSyst_hp->SetFillColor(kViolet+6);
  pgRaaSyst_hp->SetFillStyle(1001);

  systUnc[1] = new TGraphAsymmErrors(*pgRaaSyst_hp);
  systUnc[1]->SetTitle("Charged particle R_{pA}  |#eta_{CM}| < 1");
  runPeriod[1] = 2;


//----------------------------------------- jets (RpA)
  double ptBins_jp[19]={0};
  double ptError_jp[19]={0};

  double raa_jp[19];
  double raaStat_jp[19];
  double raaSyst_jp[19];
  double raaSyst2_jp[19];

  double ptSystXlow_jp[19];
  double ptSystXhigh_jp[19];

  // reading numbers
  ifstream inData_jp;
  TString inFilejp("rpa_j.dat");
  inData_jp.open(inFilejp);
  if(inData_jp.fail()) {
    cerr << "unable to open file rpa_j.dat for reading" << endl;
    exit(1);
  }
  j=0;
  while(!inData_jp.eof())
  {
    inData_jp >> xx_low >> xx_high >> raa >> raa_syst >> raa_stat >> raa_stat2;
    if( xx_low < 30. ) continue;
    ptBins_jp[j]   = xx_low+ (xx_high-xx_low)/2;
    raa_jp[j]      = raa;
    raaStat_jp[j]  = raa_syst;
    raaSyst_jp[j]  = raa_stat/100.*raa;
    raaSyst2_jp[j]  =  raa_stat2/100.*raa; //sqrt( raa_stat*raa_stat + raa_stat2*raa_stat2)/100.* raa;

    ptSystXlow_jp[j] = (xx_high-xx_low)/2;
    ptSystXhigh_jp[j]= (xx_high-xx_low)/2;
    // cout<<"pT"<<ptBins_h[j]<<"\t raa= "<<raa_h[j]<<"\t syst= "<<raaStat_h[j]<<"\t stat ="<<raaSyst_h[j]<<"\t x_low= "<<ptSystXlow_h[j]<<"\t high= "<<ptSystXhigh_h[j]<<endl;
    j++;
  }
  inData_jp.close();
  // done reading numebrs

  TGraphErrors *pgRaa_jp          = new TGraphErrors(19, ptBins_jp, raa_jp, ptError_jp, raaStat_jp);
  TGraphAsymmErrors *pgRaaSyst_jp = new TGraphAsymmErrors(19, ptBins_jp, raa_jp, ptSystXlow_jp,ptSystXhigh_jp,raaSyst2_jp,raaSyst2_jp);
  TGraphAsymmErrors *pgRaaSyst2_jp = new TGraphAsymmErrors(19, ptBins_jp, raa_jp, ptSystXlow_jp,ptSystXhigh_jp,raaSyst_jp,raaSyst_jp);
  pgRaa_jp->SetName("pgRaa_jp");
  pgRaa_jp->SetMarkerStyle(34);
  pgRaa_jp->SetMarkerSize(1.2);

  plots[2] = new TGraphErrors(*pgRaa_jp);

  //systm error
  pgRaaSyst_jp->SetName("pgRaaSyst_jp");
  pgRaaSyst_jp->SetFillColor(TColor::GetColor("#FFBF00"));
  pgRaaSyst2_jp->SetName("pgRaaSyst2_jp");
  pgRaaSyst2_jp->SetFillColor(TColor::GetColor("#FFBF00"));
  //pgRaaSyst2_jp->SetFillStyle(1);
  pgRaaSyst2_jp->SetLineWidth(1);
  //pgRaaSyst2_jp->SetFillStyle(3001);

  systUnc[2] = new TGraphAsymmErrors(*pgRaaSyst_jp);
  systUnc[2]->SetTitle("Inclusive jet  R_{pA} |#eta_{CM}| < 0.5");
  isPrelim[2] = true;
  runPeriod[2] = 2;

  //----------------------------------------- jet RAA

  TFile *jetRAAFile = new TFile("JetRAA_datapoints.root");
  TGraphAsymmErrors *pgRaa_jet = (TGraphAsymmErrors*)jetRAAFile->Get("RAA_R3_staterr_cent0");
  TGraphAsymmErrors *pgRaaSyst_jet = (TGraphAsymmErrors*)jetRAAFile->Get("RAA_R3_syserr_cent0");
  pgRaa_jet->SetName("pgRaa_jet");
  pgRaa_jet->SetMarkerStyle(33);
  pgRaa_jet->SetMarkerSize(1.4);

  plots[3] = new TGraphErrors(pgRaa_jet->GetN(), pgRaa_jet->GetX(), pgRaa_jet->GetY(), pgRaa_jet->GetEX(), pgRaa_jet->GetEY());
  
  //systm error
  pgRaaSyst_jet->SetName("pgRaaSyst_jet");
  pgRaaSyst_jet->SetFillStyle(1001);
  pgRaaSyst_jet->SetFillColor(TColor::GetColor("#00FF60"));

  systUnc[3] = new TGraphAsymmErrors(*pgRaaSyst_jet);
  systUnc[3]->SetTitle("Inclusive jet  R_{AA}  (0-5%)  |#eta| < 2");
  //isPrelim[3] = true;
  runPeriod[3] = 1;


  //----------------------------------------- b-jet RAA (UPDATED with Erratum included)
  double ptBins_bjet[5]   = {85,100,120,150,210}; 
  double ptError_bjet[5]  = {0,0,0,0,0};
  
  double raa_bjet[5]      = {0.3927,0.3277,0.4598,0.4034,0.4214}; 
  double raaStat_bjet[5]  = {0.0672, 0.0567, 0.0849, 0.1038, 0.201};
  double raaSyst_bjet[5]  = {0.0917,0.0923,0.1272,0.1529,0.2771};

  double ptSystXlow_bjet[5]      = {5,10,10,20,40};
  double ptSystXhigh_bjet[5]     = {5,10,10,20,40};
  
  TGraphErrors *pgRaa_bjet          = new TGraphErrors(5, ptBins_bjet, raa_bjet, ptError_bjet, raaStat_bjet);
  TGraphAsymmErrors *pgRaaSyst_bjet = new TGraphAsymmErrors(5, ptBins_bjet, raa_bjet, ptSystXlow_bjet,ptSystXhigh_bjet,raaSyst_bjet,raaSyst_bjet);
  pgRaa_bjet->SetName("pgRaa_bjet");
  pgRaa_bjet->SetMarkerStyle(21);
  pgRaa_bjet->SetMarkerColor(kRed);
  pgRaa_bjet->SetMarkerSize(1.);

  plots[4] = new TGraphErrors(*pgRaa_bjet);
  
  //systm error
  pgRaaSyst_bjet->SetName("pgRaaSyst_bjet");
  pgRaaSyst_bjet->SetFillColor(kGray);//TColor::GetColor("#FFBF00"));

  systUnc[4] = new TGraphAsymmErrors(*pgRaaSyst_bjet);
  systUnc[4]->SetTitle("b-jet R_{AA} (0-10%)  |#eta| < 2");
  runPeriod[4] = 1;

  //----------------------------------------- photon
  double ptBins_photon[]   = {22.26, 27.30, 34.35, 44.45, 61.72}; 
  double ptError_photon[]  = {0.0, 0.0, 0.0, 0.0, 0.0};
  
  double raa_photon[]      = {1.03, 0.84, 1.37, 0.98, 0.99}; 
  double raaStat_photon[]  = {0.12 ,0.14, 0.22, 0.24, 0.31 };
  double raaSyst_photon[]  = {0.29, 0.27,  0.25, 0.22, 0.21};

  double ptSystXlow_photon[]      = {2.5, 2.5, 5, 5, 15};
  double ptSystXhigh_photon[]     = {2.5, 2.5, 5, 5, 15};
  
  TGraphErrors *pgRaa_photon          = new TGraphErrors(5, ptBins_photon, raa_photon, ptError_photon, raaStat_photon);
  TGraphAsymmErrors *pgRaaSyst_photon = new TGraphAsymmErrors(5, ptBins_photon, raa_photon, ptSystXlow_photon,ptSystXhigh_photon,raaSyst_photon,raaSyst_photon);
  pgRaa_photon->SetName("pgRaa_photon");
  pgRaa_photon->SetMarkerStyle(23);
  pgRaa_photon->SetMarkerSize(1.);

  plots[5] = new TGraphErrors(*pgRaa_photon);

  //systm error
  pgRaaSyst_photon->SetName("pgRaaSyst_photon");
  pgRaaSyst_photon->SetFillColor(TColor::GetColor("#ffff00"));

  systUnc[5] = new TGraphAsymmErrors(*pgRaaSyst_photon);
  systUnc[5]->SetTitle("Isolated photon  (0-10%)  |#eta| < 1.44");

  runPeriod[5] = 1;


  //----------------------------------------- Z RpA (arXiv: 1512.06461v2)
  double ptBins_zpPb[13]   = {1.25,3.75, 6.25, 8.75, 11.25, 13.75, 17.5, 25., 35., 45., 60., 85., 125.};
  double ptError_zpPb[13]  = {0.0};
  
  double raa_zpPb[13]        = {0}; //PLACEHOLDERS - RATIO/POWHEG IS NOT ON HEPDATA!!
  double raaStat_zpPb[13]    = {0};
  double raaSyst_zpPb[13]    = {0};
  double ptSystXlow_zpPb[13] = {1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 2.5, 5, 5, 5, 10, 15, 25};
  double ptSystXhigh_zpPb[13] = {1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 2.5, 5, 5, 5, 10, 15, 25};

  TGraphErrors *pgRpa_z          = new TGraphErrors(13, ptBins_zpPb, raa_zpPb, ptError_zpPb, raaStat_zpPb);
  TGraphAsymmErrors *pgRpaSyst_z = new TGraphAsymmErrors(13, ptBins_zpPb, raa_zpPb, ptSystXlow_zpPb,ptSystXhigh_zpPb,raaSyst_zpPb,raaSyst_zpPb);
  pgRpa_z->SetName("pgRpa_z");   pgRpa_z->SetMarkerStyle(kFullCircle);
  pgRpa_z->SetMarkerSize(1.);
  pgRpa_z->SetMarkerStyle(22);

  plots[6] = new TGraphErrors(*pgRpa_z);

  //systm error
  pgRpaSyst_z->SetName("pgRpaSyst_z");
  pgRpaSyst_z->SetFillColor(TColor::GetColor("#009999"));

  systUnc[6] = new TGraphAsymmErrors(*pgRpaSyst_z);
  systUnc[6]->SetTitle("Z Boson R_{pA} (0-100%) |y| < 2");
  runPeriod[6] = 1;


  // ----------------------------------------- W raa
  double ptBins_w[]   = {80.38}; 
  double ptError_w[]  = {0.0};
  
  double raa_w[]      = {1.04}; 
  double raaStat_w[]  = {0.07};
  double raaSyst_w[]  = {0.12};
  double ptSystXlow_w[]      = {4};
  double ptSystXhigh_w[]     = {4};

  TGraphErrors *pgRaa_w          = new TGraphErrors(1, ptBins_w, raa_w, ptError_w, raaStat_w);
  TGraphAsymmErrors *pgRaaSyst_w = new TGraphAsymmErrors(1, ptBins_w, raa_w, ptSystXlow_w,ptSystXhigh_w,raaSyst_w,raaSyst_w);
  pgRaa_w->SetName("pgRaa_w");
  pgRaa_w->SetMarkerStyle(kFullCircle);
  pgRaa_w->SetMarkerSize(1.);
  pgRaa_w->SetMarkerStyle(21);

  plots[7] = new TGraphErrors(*pgRaa_w);
  
  //systm error
  pgRaaSyst_w->SetName("pgRaaSyst_w");
  pgRaaSyst_w->SetFillColor(TColor::GetColor("#ff88ff"));

  systUnc[7] = new TGraphAsymmErrors(*pgRaaSyst_w);
  systUnc[7]->SetTitle("W Boson R_{AA} (0-100%) p_{T}^{#mu} > 25 GeV/c, |#eta^{#mu}| < 2.1");
  runPeriod[7] = 1;


  // ----------------------------------------- 2012 non-prompt Jpsi
  double ptBins_npjpsi[]   = {7.31,8.97,11.32,16.52}; 
  double ptError_npjpsi[]  = {0.0,0.0,0.0,0.0};
  
  double raa_npjpsi[]      = {0.52,0.43,0.43,0.34};  
  double raaStat_npjpsi[]  = {0.12,0.08,0.09,0.07};
  double raaSyst_npjpsi[]  = {0.06,0.05,0.05,0.04};

  double ptSystXlow_npjpsi[]      = {0.81, 0.97, 1.32, 3.52};
  double ptSystXhigh_npjpsi[]     = {0.69, 1.03, 1.68, 13.48};

  TGraphErrors *pgRaa_npjpsi          = new TGraphErrors(4, ptBins_npjpsi, raa_npjpsi, ptError_npjpsi, raaStat_npjpsi);
  TGraphAsymmErrors *pgRaaSyst_npjpsi = new TGraphAsymmErrors(4, ptBins_npjpsi, raa_npjpsi, ptSystXlow_npjpsi,ptSystXhigh_npjpsi,raaSyst_npjpsi,raaSyst_npjpsi);
  pgRaa_npjpsi->SetName("pgRaa_npjpsi");
  pgRaa_npjpsi->SetMarkerStyle(29);
  pgRaa_npjpsi->SetMarkerSize(1.2);
  pgRaa_npjpsi->SetMarkerColor(kRed);

  plots[8] = (TGraphErrors*)(pgRaa_npjpsi->Clone(pgRaa_npjpsi->GetName()));
  
  //systm error
  pgRaaSyst_npjpsi->SetName("pgRaaSyst_npjpsi");
  // pgRaaSyst_npjpsi->SetFillColor(TColor::GetColor("#ee7711"));
  pgRaaSyst_npjpsi->SetFillColor(TColor::GetColor("#ba8a98"));

  systUnc[8] = (TGraphAsymmErrors*)(pgRaaSyst_npjpsi->Clone(pgRaaSyst_npjpsi->GetName()));
  systUnc[8]->SetTitle("B #rightarrow J/#psi (0-100%) |y| < 2.4");
  runPeriod[8] = 1;

  //--------------------------------------------- B Mesons (2015)
  double ptBins_bMeson[] = {12.5,17.5,22.5,27.5,45};
  double ptError_bMeson[] = {2.5,2.5,2.5,2.5,15};

  double raa_bMeson[] = {1.11, 1.19, 0.91, 0.86,1.14};
  double raaStat_bMeson[] = {0.08,0.1,0.12,0.18,0.19};
  double raaSyst_bMeson[] = {0.17,0.17,0.13,0.12,0.17};

  TGraphErrors *pgRaa_bMeson = new TGraphErrors(5, ptBins_bMeson, 
    raa_bMeson, ptError_bMeson, raaStat_bMeson);
  TGraphAsymmErrors *pgRaaSyst_bMeson = new TGraphAsymmErrors(5, ptBins_bMeson,
    raa_bMeson, ptError_bMeson, ptError_bMeson, raaSyst_bMeson, raaSyst_bMeson);
  pgRaa_bMeson->SetMarkerStyle(20);
  pgRaaSyst_bMeson->SetFillColor(kGreen+2);

  plots[9] = new TGraphErrors(*pgRaa_bMeson);
  systUnc[9] = new TGraphAsymmErrors(*pgRaaSyst_bMeson);
  systUnc[9]->SetTitle("B^{+} Meson R_{pA} |#eta| < 2");
  isPrelim[9] = true;
  runPeriod[9] = 2;
  
  //---------------------------------------------- B jet RpA (Final, 2015)

  TFile *fin = new TFile("RpA_BJet_Output.root");
  TH1D *bjetRpA = (TH1D*)fin->Get("RpA");
  TGraphAsymmErrors *pgRpASyst_bjet = new TGraphAsymmErrors(7);
  for(int i=0; i<7; i++){
  	TGraphErrors *temp = (TGraphErrors*)fin->Get(Form("RpA_SystErr_bin%d",i));
  	double xpoint, ypoint;
  	temp->GetPoint(0,xpoint,ypoint);
  	pgRpASyst_bjet->SetPoint(i,xpoint,ypoint);
  	pgRpASyst_bjet->SetPointError(i,bjetRpA->GetBinWidth(i+10)/2.,bjetRpA->GetBinWidth(i+10)/2., temp->GetErrorY(0), temp->GetErrorY(0));
  	temp->Delete();
  }
  pgRpASyst_bjet->SetMarkerStyle(20);
  pgRpASyst_bjet->SetFillColor(kCyan-7);

  plots[10] = new TGraphErrors(bjetRpA);
  systUnc[10] = new TGraphAsymmErrors(*pgRpASyst_bjet);
  systUnc[10]->SetTitle("b jet R_{pA}^{Pythia} |#eta_{CM}| < 2");
  runPeriod[10] = 2;

   //------------------------------------------- C jet RpA (Prelim, 2015)

  TFile *finc = new TFile("dataOverPythia_CJet_PRELIM.root");
  TH1D *cjetRpA = (TH1D*)finc->Get("hpARatio");
  TGraphErrors *pgRpASyst_cjet = (TGraphErrors*)finc->Get("systErrRatiopA");
  TGraphAsymmErrors *temp = new TGraphAsymmErrors(pgRpASyst_cjet->GetN(), pgRpASyst_cjet->GetX(), pgRpASyst_cjet->GetY(), pgRpASyst_cjet->GetEX(), pgRpASyst_cjet->GetEX(), pgRpASyst_cjet->GetEY(), pgRpASyst_cjet->GetEY());
 
  cjetRpA->SetMarkerStyle(25);
  temp->SetMarkerStyle(25);
  temp->SetFillColor(kMagenta-4);

  plots[11] = new TGraphErrors(cjetRpA);
  systUnc[11] = new TGraphAsymmErrors(*temp);
  systUnc[11]->SetTitle("c jet R_{pA}^{Pythia} |#eta_{CM}| < 2");
  isPrelim[11] = true;
  runPeriod[11] = 2;

 //------------------------------------------- Z RAA (arXiv: 1410.4825v2, Z->ll)

  //double ptBins_z[]   = {91.19};
  double ptBins_z[7]   = {2.5,7.5,15,25,35,45,75};
  //double ptError_z[]  = {0.0};
  double ptError_z[7]  = {0.0};
  
  //double raa_z[]      = {1.06}; 
  //double raaStat_z[]  = {0.05};
  //double raaSyst_z[]  = {0.08};
  //double ptSystXlow_z[]      = {4};
  //double ptSystXhigh_z[]     = {4};
  double raa_z[7]        = {0.99,1.29,0.93,1.27,1.18,1.28,0.89}; 
  double raaStat_z[7]    = {0.09,0.14,0.10,0.20,0.31,0.40,0.28};
  double raaSyst_z[7]    = {0.08,0.11,0.08,0.11,0.10,0.11,0.07};
  double ptSystXlow_z[7] = {2.5,2.5,5,5,5,5,25};
  double ptSystXhigh_z[7] = {2.5,2.5,5,5,5,5,25};

  TGraphErrors *pgRaa_z          = new TGraphErrors(7, ptBins_z, raa_z, ptError_z, raaStat_z);
  TGraphAsymmErrors *pgRaaSyst_z = new TGraphAsymmErrors(7, ptBins_z, raa_z, ptSystXlow_z,ptSystXhigh_z,raaSyst_z,raaSyst_z);
  pgRaa_z->SetName("pgRaa_z");   pgRaa_z->SetMarkerStyle(kFullCircle);
  pgRaa_z->SetMarkerSize(1.);
  pgRaa_z->SetMarkerStyle(22);

  plots[12] = new TGraphErrors(*pgRaa_z);

  //systm error
  pgRaaSyst_z->SetName("pgRpaSyst_z");
  pgRaaSyst_z->SetFillColor(TColor::GetColor("#ff8888"));

  systUnc[12] = new TGraphAsymmErrors(*pgRaaSyst_z);
  systUnc[12]->SetTitle("Z Boson R_{AA} (0-100%) |y| < 2");
  runPeriod[12] = 1;

  //--------------------------------------- Global uncertainties 

  // b-jet RpA Pythia uncert.
  TGraphErrors *pythiaErr = new TGraphErrors(1);
  pythiaErr->SetPoint(0,4,1);
  pythiaErr->SetPointError(0,4,0.22);
  pythiaErr->SetFillColor(kRed-9);

  // pPb lumi uncert.
  TGraphErrors *pPblumiErr = new TGraphErrors(1);
  pPblumiErr->SetPoint(0,8,1);
  pPblumiErr->SetPointError(0,2,0.039);
  pPblumiErr->SetFillColor(kGreen+2);
  
  // PbPb lumi uncert.
  TBox *PbPbLumi = new TBox(0,0.9568966,3,1.043103);
  PbPbLumi->SetFillColor(kBlue+2);
  PbPbLumi->SetFillStyle(1001);

  //TAA uncertainty
  TBox *taaUncert = new TBox(3,1-0.0726,6,1.0726);
  taaUncert->SetFillColor(kBlue-7);
  taaUncert->SetLineColor(kBlue-7);

      //TpA uncertainty
  TBox *tpaUncert = new TBox(10,1-0.035,14,1.035);
  tpaUncert->SetFillColor(kGreen-7);
  tpaUncert->SetLineColor(kGreen-7);


  for(int i=0; i<nObjects; i++){
    systUnc[i]->SetMarkerColor(plots[i]->GetMarkerColor());
    systUnc[i]->SetMarkerStyle(plots[i]->GetMarkerStyle());
    systUnc[i]->SetMarkerSize(plots[i]->GetMarkerSize());
  }

  //*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
  //*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*

  //   HadronRAA, HadronRpA, InclJetRpA, InclJetRAA, BJetRAA, PhotonRAA, ZRAA, WRAA, BJpsiRAA, BMesonRAA, BJetRpA, CJetRpA

  bool leftFlags[nObjects];
  bool rightFlags[nObjects];
  for(int i=0; i<nObjects; i++){
    leftFlags[i] = false;
    rightFlags[i] = false;
  }

  vector<int> plottingOrderLeft;
  vector<int> plottingOrderRight;

  for(unsigned int i=0; i<leftPlot.size(); i++){
    bool tempCheck = false;
    for(int j=0; j<nObjects; j++){
      if(leftPlot.at(i) == objectsToPlot[j]){ leftFlags[j] = true; plottingOrderLeft.push_back(j); tempCheck = true; break; }
    }
    if(leftPlot.at(i) != "placeholder" && !tempCheck) cout << "Warning! " << leftPlot.at(i) << " is not a valid plot option" << endl;
  }
  for(unsigned int i=0; i<rightPlot.size(); i++){
    bool tempCheck = false;
    for(int j=0; j<nObjects; j++){
      if(rightPlot.at(i) == objectsToPlot[j]){ rightFlags[j] = true; plottingOrderRight.push_back(j); tempCheck = true; break; }
    }
    if(rightPlot.at(i) != "placeholder" && !tempCheck) cout << "Warning! " << rightPlot.at(i) << " is not a valid plot option!" << endl;
  }

  TLatex *cmsP = new TLatex(0.03*xMax,0.915*yMax,"CMS ");
  cmsP->SetTextFont(42);
  cmsP->SetTextSize(0.0558);

  TLatex *cmsPrelim = new TLatex(0.2*xMax,0.915*yMax,"Preliminary");
  cmsPrelim->SetTextFont(42);
  cmsPrelim->SetTextSize(0.055);
  cmsPrelim->SetTextColor(kRed+2);

  bool run1flag[2] = {false, false};
  bool run2flag[2] = {false, false};

  TCanvas *cplot;
  if(rightPlot.size()>0){ cplot = new TCanvas("cplot","",1200,600); cplot->Divide(2,1); }
  else cplot = new TCanvas("cplot","",600,600);

  TLegend *leg1 = new TLegend(0.2,0.62,0.88,0.92);
  if(rightPlot.size()==0 || splitLegend) leg1->SetY2(0.851);
  if(drawGlauberUnc || drawLumiUnc) leg1->SetY1(0.679);
  TLegend *leg2 = new TLegend(0.2,0.62,0.88,0.92);
  for(int i=0; i<nObjects; i++){
    if(leftFlags[i]){
      TLegendEntry *le = leg1->AddEntry(systUnc[i],systUnc[i]->GetTitle(), "fpl");
      if(isPrelim[i]) le->SetTextColor(kRed+2);
      if(runPeriod[i] == 1) run1flag[0] = true;
      if(runPeriod[i] == 2) run2flag[0] = true;
    }
  }
  for(int i=0; i<nObjects; i++){
    if(rightFlags[i]){
      TLegendEntry *le;
      if(splitLegend) le =  leg2->AddEntry(systUnc[i],systUnc[i]->GetTitle(), "fpl");
      else le = leg1->AddEntry(systUnc[i],systUnc[i]->GetTitle(), "fpl");
      if(isPrelim[i]) le->SetTextColor(kRed+2);
      if(runPeriod[i] == 1) run1flag[1] = true;
      if(runPeriod[i] == 2) run2flag[1] = true;
    }
  }
  TLegend *glauberLeg = new TLegend(0.2,0.588,0.5,0.665);
  if(run1flag[1] && run2flag[1]){
    glauberLeg->AddEntry(taaUncert, "T_{AA} Unc.","fl");
    glauberLeg->AddEntry(tpaUncert, "T_{pA} Unc.","fl");
  }
  else if(run1flag[1]){ 
    glauberLeg->AddEntry(taaUncert, "T_{AA} Unc.","fl");
  }
  else if(run2flag[1]){
    glauberLeg->AddEntry(tpaUncert, "T_{pA} Unc.","fl");
  }

  TLegend *lumiLeg = new TLegend(0.6,0.588,0.9,0.665);
  if(run1flag[1] && run2flag[1]){
    lumiLeg->AddEntry(PbPbLumi, "PbPb Lumi. Unc.","fl");
    lumiLeg->AddEntry(pPblumiErr, "pPb Lumi. Unc.","fl");
    if(leftFlags[10] || leftFlags[11] || rightFlags[10] || rightFlags[11]) lumiLeg->AddEntry(pythiaErr, "Pythia Unc.","fl"); 
  }
  else if(run1flag[1]){ 
    lumiLeg->AddEntry(PbPbLumi, "PbPb Lumi. Unc.","fl");
  }
  else if(run2flag[1]){
    lumiLeg->AddEntry(pPblumiErr, "pPb Lumi. Unc.","fl");
    if(leftFlags[10] || leftFlags[11] || rightFlags[10] || rightFlags[11]) lumiLeg->AddEntry(pythiaErr, "Pythia Unc.","fl"); 
  }

  TLine *lineLeft = new TLine(0,1,xMax,1);
  lineLeft->SetLineStyle(7);
  TLine *lineRight = new TLine(0,1,xMax2,1);
  lineRight->SetLineStyle(7);

  TLatex *cmsLumi1 = new TLatex(49.24,yMax*1.01,"166 #mub^{-1} (PbPb 2.76 TeV)");
  cmsLumi1->SetTextFont(43);
  cmsLumi1->SetTextSize(22);
  TLatex *cmsLumi2 = new TLatex(212.5,yMax*1.01,"35 nb^{-1} (pPb, 5.02 TeV)");
  cmsLumi2->SetTextFont(43);
  cmsLumi2->SetTextSize(25);
  TLatex *cmsLumi3 = new TLatex(0.01,yMax*1.01,"166 #mub^{-1} (PbPb 2.76 TeV), 35 nb^{-1} (pPb, 5.02 TeV)");
  cmsLumi3->SetTextFont(43);
  cmsLumi3->SetTextSize(22);
  TLatex *cmsLumi1_2 = (TLatex*)cmsLumi1->Clone();
  TLatex *cmsLumi2_2 = (TLatex*)cmsLumi2->Clone();
  TLatex *cmsLumi3_2 = (TLatex*)cmsLumi3->Clone();


  cplot->cd(1);
  TH1D *thefirstAxis = new TH1D("thefirstAxis","",1,0,xMax);
  thefirstAxis->SetMaximum(yMax);
  thefirstAxis->SetMinimum(0);
  thefirstAxis->SetXTitle("p_{T} (GeV/c)");
  thefirstAxis->SetYTitle("Nuclear Modification Factor");
  thefirstAxis->Draw();
  lineLeft->Draw("same");
  cmsP->Draw("same");
  for(int i=0; i<nObjects; i++){
    if((leftFlags[i] || rightFlags[i]) && isPrelim[i]){ cmsPrelim->Draw("same"); break; }
  }
  for(unsigned int i=0; i<plottingOrderLeft.size(); i++){
    int iter = plottingOrderLeft.at(i);
    if(debug) cout << "left plotting " << objectsToPlot[iter] << endl;
    systUnc[iter]->Draw("2,same");
    plots[iter]->Draw("P z same");
    if(drawLumiUnc) if(iter==10 || iter==11) pythiaErr->Draw("2,5,Same");
  }
  if(redrawPoints){
    for(unsigned int i=0; i<plottingOrderLeft.size(); i++){
      int iter = plottingOrderLeft.at(i);
      plots[iter]->Draw("P z same");
    }
  }
  if(splitLegend) leg1->Draw("same");
  if(run1flag[1] && run2flag[1]){ 
    cmsLumi3->SetX(0.0*xMax); 
    cmsLumi3->Draw("same"); 
    if(drawLumiUnc) PbPbLumi->Draw("same"); 
    if(drawLumiUnc) pPblumiErr->Draw("2,5,same");
    if(drawGlauberUnc) taaUncert->Draw("same");
    if(drawGlauberUnc) tpaUncert->Draw("same");
  }
  else if(run1flag[1]){ 
    cmsLumi1->SetX(0.48*xMax); 
    cmsLumi1->Draw("same"); 
    if(drawLumiUnc) PbPbLumi->Draw("2,5,same");
    if(drawGlauberUnc) taaUncert->Draw("same");
  }
  else if(run2flag[1]){
    cmsLumi2->SetX(0.45*xMax); 
    cmsLumi2->Draw("same"); 
    if(drawLumiUnc) pPblumiErr->Draw("2,5,same"); 
    if(drawGlauberUnc) tpaUncert->Draw("same");
  }
  if(drawLumiUnc) lumiLeg->Draw("same");
  if(drawGlauberUnc) glauberLeg->Draw("same");


  if(rightPlot.size()>0){
    cplot->cd(2);
    TH1D *theSecondAxis = new TH1D("theSecondAxis","",1,0,xMax2);
    theSecondAxis->SetMaximum(yMax);
    theSecondAxis->SetMinimum(0);
    theSecondAxis->SetXTitle("p_{T} (GeV/c)");
    theSecondAxis->SetYTitle("Nuclear Modification Factor");
    theSecondAxis->Draw();
    lineRight->Draw("same");
    for(unsigned int i=0; i<plottingOrderRight.size(); i++){
      int iter = plottingOrderRight.at(i);
      if(debug) cout << "right plotting " << objectsToPlot[iter] << endl;
      systUnc[iter]->Draw("2,same");
      plots[iter]->Draw("P z same");
      if(drawLumiUnc) if(iter==10 || iter==11) pythiaErr->Draw("2,5,Same");
    }
    if(redrawPoints){
      for(unsigned int i=0; i<plottingOrderRight.size(); i++){
        int iter = plottingOrderRight.at(i);
        plots[iter]->Draw("P z same");
      }
    }
    if(run1flag[1] && run2flag[1]){ 
      cmsLumi3_2->SetX(0.0*xMax2); 
      cmsLumi3_2->Draw("same"); 
      if(drawLumiUnc) PbPbLumi->Draw("same"); 
      if(drawLumiUnc) pPblumiErr->Draw("2,5,same");
      if(drawGlauberUnc) taaUncert->Draw("same");
      if(drawGlauberUnc) tpaUncert->Draw("same");
    }
    else if(run1flag[1]){ 
      cmsLumi1_2->SetX(0.48*xMax2); 
      cmsLumi1_2->Draw("same"); 
      if(drawLumiUnc) PbPbLumi->Draw("2,5,same");
      if(drawGlauberUnc) taaUncert->Draw("same");
    }
    else if(run2flag[1]){
      cmsLumi2_2->SetX(0.45*xMax2); 
      cmsLumi2_2->Draw("same"); 
      if(drawLumiUnc) pPblumiErr->Draw("2,5,same"); 
      if(drawGlauberUnc) tpaUncert->Draw("same");
    }
  }
  if(!splitLegend) leg1->Draw("same");
  else if(rightPlot.size()>0) leg2->Draw("same");


}