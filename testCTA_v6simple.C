#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
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

double OnTimeDefault = 50.*60.*60.; //10 hours //The number is in seconds
TString ctaPerformDefault = "SubarrayE_IFAE_50hours_20101102.root";

TString specstrDef = "[0]*(6.0e-10)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))";  //Crab Nebula (MAGIC paper) with [0]=1. 


//TString attenFileNameDefaultTest = "exptau_z0.859_modelFranceschini.dat";
Double_t threshold = 0.1; //TeV (dummy, will be overwritten);
TString attenFileNameDefaultTest = "exptau_z0.030_modelFranceschini.dat";

bool testCTA_v6simple(const char* filename=ctaPerformDefault,
                      double crabunits=0.1, //1.00, 
                      Double_t obsTime = OnTimeDefault,
                      TString specstr = specstrDef,
                      Double_t size_source = 0.,
                      Bool_t kApplyRebinning = kFALSE,
                      Bool_t kApplyEnRes = kTRUE,
                      Bool_t kApplyAtten = kFALSE)
{

  FILE *fp;
  Bool_t fE2=kTRUE;  //if you wish to plot E2 dN/dE set to kTRUE

  const char *specstrS = specstr;
  TF1 *origSpec = new TF1("origSpec",specstrS,0.01,200.);  
  origSpec->SetNpx(10000);
  origSpec->FixParameter(0,crabunits);
  origSpec->SetLineWidth(2);
  origSpec->SetLineColor(1);

  TF1 *origSpecE2 = new TF1("origSpecE2","origSpec*x*x",0.01,200.);  
  origSpecE2->SetNpx(10000);
  origSpecE2->SetLineWidth(2);
  origSpecE2->SetLineColor(1);

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
  const int nbins = 18; //32; //18; //18; //30;
  const double emin = 0.03; //0.05; //1.122e-02; //0.01122; // 0.05; //TeV
  const double emax = 50; //50.; //1.122e+02; //112.2; // 50.; //TeV
  const double eminLog = log10(emin);
  const double emaxLog = log10(emax);
  double steplog = (emaxLog-eminLog)/(nbins);
  Double_t lowedges[151];
  for (int i=0;i<=nbins;i++)
  {
    double e = eminLog + i*steplog;
    lowedges[i] = pow(10,e);
  }

  TH1D *spObserved = new TH1D("spObserved","",nbins,lowedges);
  spObserved->Sumw2();
  spObserved->GetXaxis()->SetTitle("energy E(TeV)");
  spObserved->GetXaxis()->SetTitleOffset(1.2);
  spObserved->SetMaximum(1e-9);
  spObserved->SetMinimum(1e-14);
  spObserved->GetYaxis()->SetTitleOffset(1.4);
  spObserved->GetXaxis()->SetMoreLogLabels();
  spObserved->GetXaxis()->SetNoExponent();
  spObserved->GetXaxis()->SetLabelOffset(0.015);
  spObserved->GetYaxis()->SetTitle("dN/dE (ph TeV^{-1} cm^{-2} s^{-1})");
  if (fE2)
    spObserved->GetYaxis()->SetTitle("E^{2} dN/dE (TeV cm^{-2} s^{-1})");

  spObserved->SetLineWidth(2);
  spObserved->SetLineColor(1);
  spObserved->SetMarkerStyle(20);

  //get attenuation if needed
  TSpline3 *SplineEnVsAtt2 =NULL;
  double energyMin = 1e16;  //dummy very big
  double energyMax = 1e-16;  //dummy very small
  if (kApplyAtten)
  {
    //read attenuation factors and fill them into a spline
    cout << "reading attenuation factors from " << attenFileNameDefaultTest << endl;
    FILE *fp = NULL;
    const char *cname = attenFileNameDefaultTest;
    if (!(fp = fopen(cname, "r")))
    {
      cout << "ERROR: could not open attenuation file " << cname << endl;
      return 0;
    }
    char line[99];
    TString strline;
    Int_t k=0;

    TArrayD AttenCoef(NUMPOINT);
    TArrayD EnergyPoi(NUMPOINT);

    while(1==1)
    {
      if (fgets(line,199,fp)==NULL)
        break;
      strline = line;

      if (!strline.Contains("#"))
      {
        sscanf(line,"%lf%lf", &EnergyPoi[k], &AttenCoef[k]);
        if (EnergyPoi.At(k) < energyMin)
          energyMin = EnergyPoi.At(k);
        if (EnergyPoi.At(k) > energyMax)
          energyMax = EnergyPoi.At(k);

        cout << EnergyPoi[k] << " " << AttenCoef[k] << endl;
        k++;

        if (k>=NUMPOINT)
          break;
      }
    }

    fclose(fp);

    SplineEnVsAtt2 = new TSpline3("attenuation2",
        EnergyPoi.GetArray(), AttenCoef.GetArray(), k);

  }

  //new: use TGraphAsymmErrors to store the result
  TGraphAsymmErrors *specgr = new TGraphAsymmErrors(1);
  specgr->SetLineWidth(2);
  specgr->SetLineColor(4);
  specgr->SetMarkerColor(4);
  specgr->SetMarkerStyle(20);

  //new: use TGraphAsymmErrors to store the integral flux  
  TGraphAsymmErrors *ifluxgr = new TGraphAsymmErrors(1);

  //new: use TGraphAsymmErrors to store the number and error on excess events  in Erec
  TGraphAsymmErrors *excessgraph = new TGraphAsymmErrors(1);


  bool rc = makeCTAspec(spObserved,specgr,ifluxgr,excessgraph,filename,kApplyRebinning,kApplyEnRes,kApplyAtten,origSpec,obsTime,size_source,attenFileNameDefaultTest, &threshold);
  if (!rc) {
    cout << " errors detected in makeCRAspec routine, aborting ... " << endl;
    return 0;
  }

  //Graph to contain absorbed expected spectrum 
  TGraph *specAbsorbed = new TGraph(1);
  specAbsorbed->SetLineWidth(2);
  specAbsorbed->SetLineColor(kBlue+2);

  if(fE2)
  {
    for (int i=1;i<=nbins;i++)
    {
      double flux   = spObserved->GetBinContent(i);
      double fluxer = spObserved->GetBinError(i);
      double xLedge = spObserved->GetBinLowEdge(i);  //OK
      double xHedge = spObserved->GetBinLowEdge(i+1);  //OK
      double energy = exp(log(xLedge*xHedge)/2.);  // mean log 
      spObserved->SetBinContent(i,energy*energy*flux);  
      spObserved->SetBinError(i,energy*energy*fluxer);  

    }
      //graph
    double x,y;
    for (int i=0;i<specgr->GetN();i++)
    {
      specgr->GetPoint(i,x,y);
      y *= x*x;
      specgr->SetPoint(i,x,y);
      double eyl = specgr->GetErrorYlow(i);
      double eyh = specgr->GetErrorYhigh(i);
      eyl *= x*x;
      eyh *= x*x;
      specgr->SetPointEYlow(i,eyl);
      specgr->SetPointEYhigh(i,eyh);
    }
  }

  TH2D *sp;
  if(fE2) sp = new TH2D("sp","",1000,0.01,100,1000,1e-14,1e-9);
  else  sp = new TH2D("sp","",1000,0.01,100,1000,1e-17,1e-6);
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

//  spObserved->Draw("same");
  specgr->Draw("P");

  double chi2=0.;
  double ndf=0;
  double atten=1.;
  if (fE2)
  {
    cout << endl << " FITTING " << endl << endl;
    for (int i=0; i<specgr->GetN();i++)
    {
      double xsp;
      double ysp;
      double ytrue;
      specgr->GetPoint(i,xsp,ysp);
      if (kApplyAtten) 
        atten = SplineEnVsAtt2->Eval(xsp);
//      cout << endl << " E = " << xsp << ", " << atten << " flux " <<  origSpecE2->Eval(xsp)*atten <<  endl;
      ytrue = origSpecE2->Eval(xsp)*atten;
      chi2 += (ysp-ytrue)*(ysp-ytrue)/specgr->GetErrorYlow(i)/specgr->GetErrorYhigh(i);
      cout << i << " chi2 " << chi2 << endl;
      ndf++;

      specAbsorbed->SetPoint(i,xsp,origSpecE2->Eval(xsp)*atten);
    }
//    specgr->Fit(origSpecE2,"","",emin,emax);
//    chi2=origSpecE2->GetChisquare();
//    ndf=origSpecE2->GetNDF();
    origSpecE2->Draw("same");
    cout << endl << " END OF FITTING, chi2 = " << chi2 << ", NDF = " << ndf << endl << endl;
  }
  else
  {
    cout << endl << " FITTING " << endl << endl;
    for (int i=0; i<specgr->GetN();i++)
    {
      double xsp;
      double ysp;
      double ytrue;
      specgr->GetPoint(i,xsp,ysp);
      if (kApplyAtten) 
        atten = SplineEnVsAtt2->Eval(xsp);
      ytrue = origSpec->Eval(xsp)*atten;
      chi2 += (ysp-ytrue)*(ysp-ytrue)/specgr->GetErrorYlow(i)/specgr->GetErrorYhigh(i);
      cout << i << " chi2 " << chi2 << endl;
      ndf++;

      specAbsorbed->SetPoint(i,xsp,ytrue);
    }
    origSpec->Draw("same");
    cout << endl << " END OF FITTING, chi2 = " << chi2 << ", NDF = " << ndf << endl << endl;
  }

//  specAbsorbed->SetMarkerStyle(20);
  if (kApplyAtten)
    specAbsorbed->Draw("L"); 

    cout << endl << " END 1 " << endl; 

  TF1 *origSpecCrab=NULL;
  if(fE2) origSpecCrab = new TF1("origSpecCraab","(6.0e-10)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))*x*x",0.01,200.);  //Crab Nebula (MAGIC paper) with [0]=1. 
  else origSpecCrab = new TF1("origSpecCraab","(6.0e-10)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))",0.01,200.);  //Crab Nebula (MAGIC paper) with [0]=1. 
  origSpecCrab->SetLineStyle(7);
  origSpecCrab->SetLineColor(36);
  origSpecCrab->Draw("same");

  TLegend *leg = new TLegend(0.50,0.82,0.99,0.99);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  if (fE2) leg->AddEntry(origSpecE2,Form("intrinsic spectrum, scaled by %.3f",crabunits),"l");
  else leg->AddEntry(origSpec,Form("intrinsic spectrum, scaled by %.3f",crabunits),"l");
//  leg->AddEntry(spObserved,Form("expected spectrum (CTA), %.2lfh",obsTime/3600),"pl");
  leg->AddEntry(specgr,Form("expected spectrum (CTA), %.2lfh",obsTime/3600),"pl");
  leg->AddEntry(origSpecCrab,"Crab (MAGIC)","pl");
  leg->Draw();

  TF1 * pl = new TF1("pl","[0]*pow(x,[1])",0.01,200);
  pl->SetParameters(1e-12,-2);
  pl->SetLineColor(8);
  spObserved->Fit(pl,"REN","same",0.02,50.);
  cout << " (HIST) Fit probability " << pl->GetProb() << " chi2/ndf " << pl->GetChisquare() << "," << pl->GetNDF() << endl;
  specgr->Fit(pl,"REN","same",0.02,5.);
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
  for (int j=0; j<specgr->GetN(); j++)
  {
    specgr->GetPoint(j,ener,iflux);
    ifluxerp = specgr->GetErrorYhigh(j);
    ifluxern = specgr->GetErrorYlow(j);
    enh = specgr->GetErrorXhigh(j);
    enl = specgr->GetErrorXlow(j);
    if (ener>threshold)
    cout << Form ("%d\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e",j,ener,ener-enl,ener+enh,iflux,ifluxern,ifluxerp)<< endl;
    else
    cout << Form ("%d\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t(point below threshold)",j,ener,ener-enl,ener+enh,iflux,ifluxern,ifluxerp)<< endl;
  }
  cout << "####################################################################################################"<<endl;

  ifluxgr->GetPoint(0,ener,iflux);
  cout << endl;
  cout << "####################################################################################################"<<endl;
  cout << "      Integral flux "<<endl;
  cout << "####################################################################################################"<<endl;
  cout << " Integral flux above "<< ener <<" TeV is " << endl;
  cout << Form("\t F_meas = %.2e (+%.2e) (-%0.2e) cm^{-2} s^{-1}",iflux,ifluxgr->GetErrorYhigh(0),ifluxgr->GetErrorYlow(0)) << endl;
  if (ifluxgr->GetErrorYlow(0)>0.) cout << Form("\t Significance of F_meas is %.2lf sigma",iflux/ifluxgr->GetErrorYlow(0)) << endl;
  cout << " The real integrated flux above "<< emin<<" TeV is " << endl;
  double F_real = origSpec->Integral(emin,emax);
  cout << Form("\t F_real = %0.2e cm^{-2} s^{-1}",F_real) << endl;
  if (iflux>F_real) cout << Form("INFO: F_meas is %.2lf sigma above F_real",(iflux-F_real)/ifluxgr->GetErrorYhigh(0)) << endl;
  else
    cout << Form("INFO: F_meas is %.2lf sigma below F_real",(F_real-iflux)/ifluxgr->GetErrorYlow(0)) << endl;

  cout << "####################################################################################################"<<endl;
  cout << " Comparison between the input spectrum and the simulated spectrum leads to " << endl;
  cout << " chi2/NDF = " << chi2/ndf << endl;
  cout << "####################################################################################################"<<endl;
//  cout << " REAL " << origSpec->Integral(emin,emax) << endl;

//  for (int i=0; i<specAbsorbed->GetN();i++)
//  {
//    double xsp;
//    double ysp;
//    double ytrue;
//    specAbsorbed->GetPoint(i,xsp,ysp);
//    cout << " E = " << xsp << ", flux " <<  ysp <<  endl;
//  }


  cout << " THRESHOLD of the analysis is "  << threshold  << " TeV "<< endl;

  return 1;

}
