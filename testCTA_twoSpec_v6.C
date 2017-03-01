#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <stdio.h>
#include <TStyle.h>

/*
 * =====================================================================================
 *
 *       Filename:  testCTA.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/07/09 15:10:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Daniel Mazin (dm), mazin@ifae.es
 *        Company:  IFAE, Barcelona, Spain
 *
 * =====================================================================================
 */

#include "makeCTAspec_v6.C"

double OnTimeDefault = 10.0*60.*60.; //10 hours //The number is in seconds
//TString ctaPerformDefault = "SubarrayEtoytest.root";
//TString ctaPerformDefault = "SubarrayEkb1.root";
//TString ctaPerformDefault = "kb_J_v3.root";
//TString ctaPerformDefault = "SubarrayE_IFAE_20101021.root";
//TString ctaPerformDefault = "SubarrayX_IFAE_DATE.root";
//TString ctaPerformDefault = "/Users/mazin/physics/CTA/taskB/IFAE_Nov2010/SubarrayE_IFAE_50hours_20101102.root";
TString ctaPerformDefault = "SubarrayE_IFAE_50hours_20101102.root";
//TString ctaPerformDefault = "SubarrayE_IFAE_20101030.root";
//TString specstrDef = "[0]*(2.0e-10)*pow((x/1.),-1.6)*exp(-x/12.)"; //Mkn501 stecker 
//TString specstrDef = "[0]*(1.3e-10)*pow((x/1.),-2.0)*exp(-x/55.)"; //Mkn501 franceschini 

//TString specstrDef = "[0]*(6.0e-10)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))";  //Crab Nebula (MAGIC paper) with [0]=1. 
//TString specstrDef2 = "[0]*(6.0e-10)*pow((x/.3),(-2.36+(-0.26*TMath::Log10(x/.3))))";  //Crab Nebula (MAGIC paper) with [0]=1. 

//TString specstrDef = "([0]*(5.83e-12)*pow((x/1.0),-2.62)+(2e-10)*exp(-(log10(x*3)*log10(x*3))/(0.002)))*exp(-x/1.0)";  //Crab Nebula (HEGRA paper) with [0]=1 + feature at 0.1TeV + cutoff at 1 TeV
TString specstrDef = "([0]*(2.83e-11)*pow((x/1.0),-2.62))";  //Crab Nebula (HEGRA paper) with [0]=1
//TString specstrDef2 = "([0]*(2.83e-11)*pow((x/1.0),-2.62))";  //Crab Nebula (HEGRA paper) with [0]=1 + feature at 0.1TeV + cutoff at 1 TeV
TString specstrDef2 = "([0]*(2.83e-11)*pow((x/1.0),-2.62))*exp(-x/15.0)";  //Crab Nebula (HEGRA paper) with [0]=1 + feature at 0.1TeV + cutoff at 1 TeV
//
//
//


//TString specstrDef = "[0]*(2.83e-11)*pow((x/1.0),-2.62)";  //Crab Nebula (HEGRA paper) with [0]=1
//TString specstrDef2 = "[0]*(2.83e-11)*pow((x/1.0),-2.12)";  //Crab Nebula (HEGRA paper) with [0]=1

//TString specstrDef = "[0]*(2.83e-11)*pow((x/1.0),-2.62)+(3e-11)*exp(-(log10(x)*log10(x))/(0.001))";  //Crab Nebula (HEGRA paper) with [0]=1 + feature 
//TString specstrDef = "([0]*(2.83e-11)*pow((x/1.0),-2.62)+(3e-11)*exp(-(log10(x)*log10(x))/(0.001)))*exp(-x/1.0)";  //Crab Nebula (HEGRA paper) with [0]=1 + cutoff + feature 
//TString specstrDef = "[0]*(2.83e-11)*pow((x/1.0),-2.62)+(3e-9)*exp(-(log10(x*10)*log10(x*10))/(0.001))";  //Crab Nebula (HEGRA paper) with [0]=1 + feature 
//TString specstrDef = "[0]*(2.83e-11)*pow((x/1.0),-2.62)+(3e-9)*exp(-(log10(x*10)*log10(x*10))/(0.01))";  //Crab Nebula (HEGRA paper) with [0]=1 + feature 
//TString specstrDef = "([0]*(2.83e-11)*pow((x/1.0),-2.62)+(3e-9)*exp(-(log10(x*10)*log10(x*10))/(0.01)))*exp(-x/1.0)";  //Crab Nebula (HEGRA paper) with [0]=1 + feature at 0.1TeV + cutoff at 1 TeV

//TString specstrDef = "([0]*(2.83e-11)*pow((x/1.0),-2.62)+(5e-10)*exp(-(log10(x*3)*log10(x*3))/(0.002)))*exp(-x/1.0)";  //Crab Nebula (HEGRA paper) with [0]=1 + feature at 0.1TeV + cutoff at 1 TeV
//
//TString specstrDef = "([0]*(5.83e-12)*pow((x/1.0),-2.62)+(5e-10)*exp(-(log10(x*3)*log10(x*3))/(0.002)))*exp(-x/1.0)";  //Crab Nebula (HEGRA paper) with [0]=1 + feature at 0.1TeV + cutoff at 1 TeV
//
//TString specstrDef = "([0]*(2.83e-11)*pow((x/1.0),-2.62)+(2e-10)*exp(-(log10(x*3)*log10(x*3))/(0.0005)))*exp(-x/1.0)";  //Crab Nebula (HEGRA paper) with [0]=1 + feature at 0.1TeV + cutoff at 1 TeV


//TString specstrDef = "([0]*(2.83e-12)*pow((x/1.0),-2.62)+(2e-10)*exp(-(log10(x*3)*log10(x*3))/(0.002)))*exp(-x/1.0)";  //Crab Nebula (HEGRA paper) with [0]=1 + feature at 0.1TeV + cutoff at 1 TeV
//
//
//
//TString specstrDef = "pow(10,(0.0693527*pow(log10(x),4)-0.100257*pow(log10(x),3)-0.488969*pow(log10(x),2)-1.78426*log10(x)-10.952))"; //d=50pc, t=2000yr

Double_t threshold = 0.0; //TeV (dummy, will be overwritten);
//define range where the chi2 test will be performed
Double_t EminForTest = 0.3; //TeV 
Double_t EmaxForTest = 100.0; //TeV 

bool testCTA_twoSpec_v6(const char* filename=ctaPerformDefault,
                      double crabunits=1.00, 
                      Double_t obsTime = OnTimeDefault,
                      TString specstr1 = specstrDef,
                      TString specstr2 = specstrDef2,
                      Double_t size_source = 0.01,
                      Bool_t kApplyRebinning = kTRUE,
                      Bool_t kApplyEnRes = kTRUE,
                      Bool_t kApplyAtten = kFALSE)
{

  FILE *fp;
//  static const Int_t nH1=20;
  Bool_t fE2=kTRUE;  //if you wish to plot E2 dN/dE set to kTRUE

  const char *specstrS1 = specstr1;
  TF1 *origSpec1 = new TF1("origSpec1",specstrS1,0.01,200.);  
  origSpec1->SetNpx(10000);
  origSpec1->FixParameter(0,crabunits);
  origSpec1->SetLineWidth(2);
  origSpec1->SetLineColor(4);

  TF1 *origSpec1E2 = new TF1("origSpec1E2","origSpec1*x*x",0.01,200.);  
  origSpec1E2->SetNpx(10000);
  origSpec1E2->SetLineWidth(2);
  origSpec1E2->SetLineColor(4);

  TCanvas *c1 = new TCanvas("spectrum","my CTA spectrum",800,600);
  c1->SetFillColor(10);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetLeftMargin(0.128571);
  c1->SetRightMargin(0.0777778);
  c1->SetTopMargin(0.118217);
  c1->SetBottomMargin(0.128915);
  c1->SetFrameFillColor(0);
  c1->SetFrameLineWidth(2);
  c1->SetFrameBorderMode(0);
  c1->SetFrameLineWidth(3);
  c1->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);

  //make an array of low edges (in log) between emin and emax
  const int nbins = 36; //18; //30;
  const double emin = 0.03; //0.05; //1.122e-02; //0.01122; // 0.05; //TeV
  const double emax = 120; //50.; //1.122e+02; //112.2; // 50.; //TeV
  const double eminLog = log10(emin);
  const double emaxLog = log10(emax);
  double steplog = (emaxLog-eminLog)/(nbins);
  Double_t lowedges[151];
  for (int i=0;i<=nbins;i++)
  {
    double e = eminLog + i*steplog;
    lowedges[i] = pow(10,e);
  }

  TH1D *spObserved1 = new TH1D("spObserved1","",nbins,lowedges);
  spObserved1->Sumw2();
  spObserved1->GetXaxis()->SetTitle("energy E(TeV)");
  spObserved1->GetXaxis()->SetTitleOffset(1.2);
  spObserved1->SetMaximum(1e-9);
  spObserved1->SetMinimum(1e-14);
  spObserved1->GetYaxis()->SetTitleOffset(1.4);
  spObserved1->GetXaxis()->SetMoreLogLabels();
  spObserved1->GetXaxis()->SetNoExponent();
  spObserved1->GetXaxis()->SetLabelOffset(0.015);
  spObserved1->GetYaxis()->SetTitle("dN/dE (ph TeV^{-1} cm^{-2} s^{-1})");
  if (fE2)
    spObserved1->GetYaxis()->SetTitle("E^{2} dN/dE (TeV cm^{-2} s^{-1})");

  spObserved1->SetLineWidth(2);
  spObserved1->SetLineColor(1);
  spObserved1->SetMarkerStyle(20);

  //new: use TGraphAsymmErrors to store the result
  TGraphAsymmErrors *specgr1 = new TGraphAsymmErrors(1);
  specgr1->SetLineWidth(2);
  specgr1->SetLineColor(4);
  specgr1->SetMarkerColor(4);
  specgr1->SetMarkerStyle(20);

  //new: use TGraphAsymmErrors to store the integral flux  
  TGraphAsymmErrors *ifluxgr1 = new TGraphAsymmErrors(1);

  //new: use TGraphAsymmErrors to store the number and error on excess events  in Erec
  TGraphAsymmErrors *excessgraph1 = new TGraphAsymmErrors(1);

  //CALL MAKECTASPEC
  bool rc = makeCTAspec(spObserved1,specgr1,ifluxgr1,excessgraph1,filename,kApplyRebinning,kApplyEnRes,kApplyAtten,origSpec1,obsTime,size_source);
  if (!rc) {
    cout << " errors detected in makeCRAspec routine, aborting ... " << endl;
    return 0;
  }

  if(fE2)
  {
    for (int i=1;i<=nbins;i++)
    {
      double flux   = spObserved1->GetBinContent(i);
      double fluxer = spObserved1->GetBinError(i);
      double xLedge = spObserved1->GetBinLowEdge(i);  //OK
      double xHedge = spObserved1->GetBinLowEdge(i+1);  //OK
      double energy = exp(log(xLedge*xHedge)/2.);  // mean log 
      spObserved1->SetBinContent(i,energy*energy*flux);  
      spObserved1->SetBinError(i,energy*energy*fluxer);  

    }
      //graph
    double x,y;
    for (int i=0;i<specgr1->GetN();i++)
    {
      specgr1->GetPoint(i,x,y);
      y *= x*x;
      specgr1->SetPoint(i,x,y);
      double eyl = specgr1->GetErrorYlow(i);
      double eyh = specgr1->GetErrorYhigh(i);
      eyl *= x*x;
      eyh *= x*x;
      specgr1->SetPointEYlow(i,eyl);
      specgr1->SetPointEYhigh(i,eyh);
    }
  }

  TH2D *sp;
  if(fE2) sp = new TH2D("sp","",1000,0.01,100,1000,1e-14,1e-9);
  else  sp = new TH2D("sp","",1000,0.01,100,1000,1e-16,1e-6);
  sp->Sumw2();
  sp->GetXaxis()->SetTitle("energy E(TeV)");
  sp->GetXaxis()->SetTitleOffset(1.2);
  sp->GetYaxis()->SetTitleOffset(1.4);
  sp->GetXaxis()->SetMoreLogLabels();
  sp->GetXaxis()->SetNoExponent();
  sp->GetXaxis()->SetLabelOffset(0.015);
  sp->GetYaxis()->SetTitle("dN/dE (ph TeV^{-1} cm^{-2} s^{-1})");
  if (fE2)
    sp->GetYaxis()->SetTitle("E^{2} dN/dE (TeV cm^{-2} s^{-1})");
  sp->Draw();

  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();

//  spObserved1->Draw("same");
  specgr1->Draw("P");

  double chi2;
  double ndf;
  if (fE2)
  {
    specgr1->Fit(origSpec1E2,"","",emin,emax);
    chi2=origSpec1E2->GetChisquare();
    ndf=origSpec1E2->GetNDF();
    origSpec1E2->Draw("same");
  }
  else
  {
    specgr1->Fit(origSpec1,"","",emin,emax);
    chi2=origSpec1->GetChisquare();
    ndf=origSpec1->GetNDF();
    origSpec1->Draw("same");
  }



  TF1 *origSpecCrab;
  if(fE2) origSpecCrab = new TF1("origSpecCrab","(6.0e-10)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))*x*x",0.01,200.);  //Crab Nebula (MAGIC paper) with [0]=1. 
  else origSpecCrab = new TF1("origSpecCrab","(6.0e-10)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))",0.01,200.);  //Crab Nebula (MAGIC paper) with [0]=1. 
  origSpecCrab->SetLineStyle(7);
  origSpecCrab->SetLineColor(36);
  origSpecCrab->Draw("same");

  TLegend *leg = new TLegend(0.50,0.82,0.99,0.99);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  if (fE2) leg->AddEntry(origSpec1E2,Form("intrinsic spectrum 1"),"l");
  else leg->AddEntry(origSpec1,Form("intrinsic spectrum 1"),"l");
//  leg->AddEntry(spObserved1,Form("expected spectrum (CTA), %.2lfh",obsTime/3600),"pl");
  leg->AddEntry(specgr1,Form("expected spectrum 1 (CTA), %.2lfh",obsTime/3600),"pl");

  TF1 * pl = new TF1("pl","[0]*pow(x,[1])",0.01,200);
  pl->SetParameters(1e-12,-2);
  pl->SetLineColor(8);
  spObserved1->Fit(pl,"REN","same",0.02,50.);
  cout << " (HIST) Fit probability " << pl->GetProb() << " chi2/ndf " << pl->GetChisquare() << "," << pl->GetNDF() << endl;
  specgr1->Fit(pl,"REN","same",0.02,5.);
  cout << " (GRAPH) Fit probability " << pl->GetProb() << " chi2/ndf " << pl->GetChisquare() << "," << pl->GetNDF() << endl;
  cout << " (GRAPH) slope " << pl->GetParameter(1)-2 << " +/- " << pl->GetParError(1) << endl;

  cout << endl;
  cout << "####################################################################################################"<<endl;
  cout << "      Differential energy spectrum "<<endl;
  cout << "####################################################################################################"<<endl;
  cout << Form ("Bin\tEnergy\t\tE_low\t\tE_high\t\tFlux\t\tFlux_neg_error\tFlux_pos_error") << endl;
  if (fE2) cout << Form ("#\tTeV\t\tTeV\t\tTeV\t\tTeV cm^-2 s^-1\t") << endl;
  else cout << Form ("#\tTeV\t\tTeV\t\tTeV\t\tcm^-2 s^-1 TeV^-1\t") << endl;
  double ener, enl, enh, iflux, ifluxerp, ifluxern;
  for (int j=0; j<specgr1->GetN(); j++)
  {
    specgr1->GetPoint(j,ener,iflux);
    ifluxerp = specgr1->GetErrorYhigh(j);
    ifluxern = specgr1->GetErrorYlow(j);
    enh = specgr1->GetErrorXhigh(j);
    enl = specgr1->GetErrorXlow(j);
    cout << Form ("%d\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e",j,ener,ener-enl,ener+enh,iflux,ifluxern,ifluxerp)<< endl;
  }
  cout << "####################################################################################################"<<endl;

  ifluxgr1->GetPoint(0,ener,iflux);
  cout << endl;
  cout << "####################################################################################################"<<endl;
  cout << "      Integral flux "<<endl;
  cout << "####################################################################################################"<<endl;
  cout << " Integral flux above "<< ener <<" TeV is " << endl;
  cout << Form("\t F_meas = %.2e (+%.2e) (-%0.2e) cm^{-2} s^{-1}",iflux,ifluxgr1->GetErrorYhigh(0),ifluxgr1->GetErrorYlow(0)) << endl;
  if (ifluxgr1->GetErrorYlow(0)>0.) cout << Form("\t Significance of F_meas is %.2lf sigma",iflux/ifluxgr1->GetErrorYlow(0)) << endl;
  cout << " The real integrated flux above "<< emin<<" TeV is " << endl;
  double F_real = origSpec1->Integral(emin,emax);
  cout << Form("\t F_real = %0.2e cm^{-2} s^{-1}",F_real) << endl;
  if (iflux>F_real) cout << Form("INFO: F_meas is %.2lf sigma above F_real",(iflux-F_real)/ifluxgr1->GetErrorYhigh(0)) << endl;
  else
    cout << Form("INFO: F_meas is %.2lf sigma below F_real",(F_real-iflux)/ifluxgr1->GetErrorYlow(0)) << endl;
  cout << "####################################################################################################"<<endl;

  cout << " chi2/NDF " << chi2/ndf << endl;

  /*****************************************************/
  /*****************************************************/
  /*****************************************************/
  /*****************************************************/
// prepare for the call of the second spectrum
  /*****************************************************/
  /*****************************************************/
  /*****************************************************/
  TH1D *spObserved2 = new TH1D("spObserved2","",nbins,lowedges);
  spObserved2->Sumw2();
  spObserved2->GetXaxis()->SetTitle("energy E(TeV)");
  spObserved2->GetXaxis()->SetTitleOffset(1.2);
  spObserved2->SetMaximum(1e-9);
  spObserved2->SetMinimum(1e-14);
  spObserved2->GetYaxis()->SetTitleOffset(1.4);
  spObserved2->GetXaxis()->SetMoreLogLabels();
  spObserved2->GetXaxis()->SetNoExponent();
  spObserved2->GetXaxis()->SetLabelOffset(0.015);
  spObserved2->GetYaxis()->SetTitle("dN/dE (ph TeV^{-1} cm^{-2} s^{-1})");
  if (fE2)
    spObserved2->GetYaxis()->SetTitle("E^{2} dN/dE (TeV cm^{-2} s^{-1})");

  spObserved2->SetLineWidth(2);
  spObserved2->SetLineColor(1);
  spObserved2->SetMarkerStyle(20);

  TGraphAsymmErrors *specgr2 = new TGraphAsymmErrors(1);
  specgr2->SetLineWidth(2);
  specgr2->SetLineColor(kGreen+3);
  specgr2->SetMarkerColor(kGreen+3);
  specgr2->SetMarkerStyle(20);

  TGraphAsymmErrors *ifluxgr2 = new TGraphAsymmErrors(1);
  TGraphAsymmErrors *excessgraph2 = new TGraphAsymmErrors(1);

  const char *specstrS2 = specstr2;
  TF1 *origSpec2 = new TF1("origSpec2",specstrS2,0.01,200.);  
  origSpec2->SetNpx(10000);
  origSpec2->FixParameter(0,crabunits);
  origSpec2->SetLineWidth(2);
  origSpec2->SetLineColor(kGreen+3);

  TF1 *origSpec2E2 = new TF1("origSpec2E2","origSpec2*x*x",0.01,200.);  
  origSpec2E2->SetNpx(10000);
  origSpec2E2->SetLineWidth(2);
  origSpec2E2->SetLineColor(kGreen+3);

  bool rc2 = makeCTAspec(spObserved2,specgr2,ifluxgr2,excessgraph2,filename,kApplyRebinning,kApplyEnRes,kApplyAtten,origSpec2,obsTime,size_source);
  if (!rc2) {
    cout << " errors detected in makeCRAspec routine (second), aborting ... " << endl;
    return 0;
  }

  if(fE2)
  {
    for (int i=1;i<=nbins;i++)
    {
      double flux   = spObserved2->GetBinContent(i);
      double fluxer = spObserved2->GetBinError(i);
      double xLedge = spObserved2->GetBinLowEdge(i);  //OK
      double xHedge = spObserved2->GetBinLowEdge(i+1);  //OK
      double energy = exp(log(xLedge*xHedge)/2.);  // mean log 
      spObserved2->SetBinContent(i,energy*energy*flux);  
      spObserved2->SetBinError(i,energy*energy*fluxer);  

    }
      //graph
    double x,y;
    for (int i=0;i<specgr2->GetN();i++)
    {
      specgr2->GetPoint(i,x,y);
      y *= x*x;
      specgr2->SetPoint(i,x,y);
      double eyl = specgr2->GetErrorYlow(i);
      double eyh = specgr2->GetErrorYhigh(i);
      eyl *= x*x;
      eyh *= x*x;
      specgr2->SetPointEYlow(i,eyl);
      specgr2->SetPointEYhigh(i,eyh);
    }
  }

  if (fE2)
  {
    origSpec2E2->Draw("same");
  }
  else
  {
    origSpec2->Draw("same");
  }
  specgr2->Draw("P");

  if (fE2) leg->AddEntry(origSpec2E2,Form("intrinsic spectrum 2"),"l");
  else leg->AddEntry(origSpec2,Form("intrinsic spectrum 2"),"l");
  leg->AddEntry(specgr2,Form("expected spectrum 2 (CTA), %.2lfh",obsTime/3600),"pl");
  leg->AddEntry(origSpecCrab,"Crab (MAGIC)","pl");
  leg->Draw();

  // compare the two excess event distributions
  TCanvas *c2 = new TCanvas("excessdd","differences in the excess",800,600);
  c2->SetFillColor(10);
  c2->SetFillColor(10);
  c2->SetBorderMode(0);
  c2->SetBorderSize(0);
  c2->SetLeftMargin(0.128571);
  c2->SetRightMargin(0.0777778);
  c2->SetTopMargin(0.118217);
  c2->SetBottomMargin(0.128915);
  c2->SetFrameFillColor(0);
  c2->SetFrameLineWidth(2);
  c2->SetFrameBorderMode(0);
  c2->SetFrameLineWidth(3);
  c2->SetFrameBorderMode(0);
  c2->Divide(1,2);

  c2->cd(1);
  double max1=0.; 
  double max2=0.;
  double minr=1.;
  double maxr=0.;
  double xe1;
  double ye1;
  double xe2;
  double ye2;
  double ratio;
  double r2;
  double su1;
  double su2;
  double sd1;
  double sd2;
  double siU;
  double siD;
  int jk=0;
  double  chi2r=0.;
  int  ndf2;
  TGraphAsymmErrors *excessRatio = new TGraphAsymmErrors(1);
  for (int k=0;k<excessgraph1->GetN();k++)
  {
    excessgraph1->GetPoint(k,xe1,ye1);
    if (ye1>max1)
      max1 = ye1;
    excessgraph2->GetPoint(k,xe2,ye2);
    if (ye2>max2)
      max2 = ye2;

    if (ye1>0. && ye2 >0. && xe1 == xe2 && xe1 >= EminForTest && xe1 <= EmaxForTest)
    {
//      cout << " UHU " << k << endl;
      sd1 = excessgraph1->GetErrorYlow(k);
      sd1 *= sd1;
      sd2 = excessgraph2->GetErrorYlow(k);
      sd2 *= sd2;
      su1 = excessgraph1->GetErrorYhigh(k);
      su1 *= su1;
      su2 = excessgraph2->GetErrorYhigh(k);
      su2 *= su2;
      ratio = ye1/ye2;
      if (minr > ratio)
        minr = ratio;
      if (maxr < ratio)
        maxr = ratio;

      r2 = ratio*ratio;
      siU = r2/ye1/ye1*su1 + r2/ye2/ye2*su2;
      siU = sqrt(siU);
      siD = r2/ye1/ye1*sd1 + r2/ye2/ye2*sd2;
      siD = sqrt(siD);
      excessRatio->SetPoint(jk,xe1,ratio);
      excessRatio->SetPointError(jk,0.,0.,siD,siU);

      //get chi2
      chi2r += (1.- ratio)*(1.-ratio)/siD/siU;
      jk++;
    }
  }

  ndf2 = jk;

  double maxm = max2>max1 ? max2 : max1;

  gPad->SetBottomMargin(0.15);

  TH2D *exh = new TH2D("sp","",100,0.01,200,100,0.5,maxm*2); 
  exh->GetXaxis()->SetTitle("Rec. energy E(TeV)");
  exh->GetXaxis()->SetLabelSize(0.05);
  exh->GetXaxis()->SetTitleSize(0.05);
  exh->GetYaxis()->SetLabelSize(0.05);
  exh->GetYaxis()->SetTitleSize(0.05);
  exh->GetXaxis()->SetTitleOffset(1.2);
  exh->GetYaxis()->SetTitleOffset(1.0);
  exh->GetXaxis()->SetMoreLogLabels();
  exh->GetXaxis()->SetNoExponent();
  exh->GetXaxis()->SetLabelOffset(0.015);
  exh->GetYaxis()->SetTitle("excess events");
  exh->Draw();

  gPad->SetLogy();
  gPad->SetLogx();

  excessgraph1->SetLineWidth(2);
  excessgraph1->SetLineColor(kBlue);
  excessgraph1->SetMarkerStyle(20);
  excessgraph1->SetMarkerColor(kBlue);
  excessgraph1->Draw("P");

  excessgraph2->SetLineWidth(2);
  excessgraph2->SetLineColor(kGreen+3);
  excessgraph2->SetMarkerStyle(20);
  excessgraph2->SetMarkerColor(kGreen+3);
  excessgraph2->Draw("P");

  TLegend *leg1a = new TLegend(0.50,0.82,0.99,0.99);
  leg1a->SetFillColor(10);
  leg1a->SetBorderSize(1);
  leg1a->AddEntry(excessgraph1,"Excess events, spectrum 1","pl");
  leg1a->AddEntry(excessgraph2,"Excess events, spectrum 2","pl");
  leg1a->Draw();

  c2->cd(2);
  gPad->SetLogx();
  gPad->SetBottomMargin(0.15);

  TH2D *exhr = new TH2D("sp","",100,0.01,200,100,minr-0.5,maxr+0.5); 
  exhr->GetXaxis()->SetTitle("Rec. energy E(TeV)");
  exhr->GetXaxis()->SetLabelSize(0.05);
  exhr->GetXaxis()->SetTitleSize(0.05);
  exhr->GetYaxis()->SetLabelSize(0.05);
  exhr->GetYaxis()->SetTitleSize(0.05);
  exhr->GetXaxis()->SetTitleOffset(1.3);
  exhr->GetYaxis()->SetTitleOffset(1.0);
  exhr->GetXaxis()->SetMoreLogLabels();
  exhr->GetXaxis()->SetNoExponent();
  exhr->GetXaxis()->SetLabelOffset(0.015);
  exhr->GetYaxis()->SetTitle("Ratio in excess events");
  exhr->Draw();

  excessRatio->SetLineWidth(2);
  excessRatio->SetLineColor(kOrange+10);
  excessRatio->SetMarkerStyle(20);
  excessRatio->SetMarkerColor(kOrange+10);
  excessRatio->Draw("P");

  TF1 * fr = new TF1 ("fr","[0]",emin,emax);
  fr->FixParameter(0,1.);
  fr->Draw("same");
//  excessRatio->Fit(fr);

  TLegend *leg1b = new TLegend(0.60,0.82,0.99,0.99);
  leg1b->SetFillColor(10);
  leg1b->SetBorderSize(1);
  leg1b->AddEntry(excessRatio,"Excess events ratio","pl");
  leg1b->Draw();

  cout << " chi2/NDF " << chi2r/(double)ndf2 << endl;
  cout << " Probability to be consistent: " << TMath::Prob(chi2r,ndf2) << endl;

  //KS test
  //stupid trick
  TH1D* exc1 = new TH1D(*spObserved1);
  exc1->Reset();
  TH1D* exc2 = new TH1D(*spObserved1);
  exc2->Reset();
  for (int k=0;k<excessgraph1->GetN();k++)
  {
    excessgraph1->GetPoint(k,xe1,ye1);
    if (xe1 >= EminForTest && xe1 <= EmaxForTest)
    {
      sd1 = excessgraph1->GetErrorYlow(k);
      su1 = excessgraph1->GetErrorYhigh(k);
      exc1->SetBinContent(k+1,ye1);
      exc1->SetBinError(k+1,(sd1+su1)/2.);

      excessgraph2->GetPoint(k,xe2,ye2);
      sd2 = excessgraph2->GetErrorYlow(k);
      su2 = excessgraph2->GetErrorYhigh(k);
      exc2->SetBinContent(k+1,ye2);
      exc2->SetBinError(k+1,(sd2+su2)/2.);
    }

  }
  double KS1 = exc1->KolmogorovTest(exc2);
  cout << "####################################################################################################"<<endl;
  cout << " Kolmogorov Smirnov  agreement (no N option) is " << KS1 << endl;
  double KS2 = exc1->KolmogorovTest(exc2,"N");
  cout << " Kolmogorov Smirnov  agreement (with N option) is " << KS2 << endl;
  cout << "####################################################################################################"<<endl;

  TLatex *lat = new TLatex();
  lat->SetTextSize(0.06);
  lat->SetNDC();
  lat->DrawLatex(0.2,0.8,Form("#chi^{2}/NDF = %.2lf/%d, Prob = %.2e",chi2r,ndf2,TMath::Prob(chi2r,ndf2)));
  lat->DrawLatex(0.2,0.6,Form("Chosen energy range is from %.1lf TeV to %.1lf TeV",EminForTest, EmaxForTest));

  c2->cd(1);
  lat->DrawLatex(0.15,0.3,Form("KolmogorovSmirnovTest Probability = %.2e",KS1));

  c1->cd();

  double Chi2T = exc1->Chi2Test(exc2);
  cout << "####################################################################################################"<<endl;
  cout << " Chi2 Test agreement between the two histograms is " << Chi2T << endl;
  cout << "####################################################################################################"<<endl;

  return 1;

}
