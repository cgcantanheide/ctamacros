#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TF3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <stdio.h>
#include <TStyle.h>
#include <TMarker.h>
#include <TEllipse.h>

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

#include "makeCTAspec_extended_v6.C"

Double_t getflux (double energy)
{

  double res = ((6.0e-10)*pow((energy/.3),(-2.31+(-0.26*TMath::Log10(energy/.3)))));
//  Double_t res = (6.8e-12*pow(energy/1.0,-2.60));
//  cout <<  "gf= " << res << endl;
  return res;

}


Double_t getflux_cooling (double energy)
{

//  return ((6.0e-10)*pow((energy/.3),(-2.31+(-0.26*TMath::Log10(energy/.3)))));
  return (6.0e-10*pow(energy/0.3,-2.0));

}

double OnTimeDefault = 50.*60.*60.; //50 hours //The number is in seconds
//TString ctaPerformDefault = "SubarrayE_IFAE_50hours_11032011_offaxis.root";
TString ctaPerformDefault = "SubarrayB_IFAE_50hours_09052011_offaxis.root";

#define SIGMAMIN 0.035
#define DISTMAX  4.5
#define INTMAX  5.0

//dimension of map (in degrees):
double mapXmin = -5.;
double mapXmax = 5.;
double mapYmin = -5.;
double mapYmax = 5.;
//energy limits (in TeV):
double mapEmin = 0.015; 
double mapEmax = 200.; 

Double_t threshold = 0.1; //TeV (dummy, will be overwritten);

Double_t gausfunctionPSF(Double_t *x, Double_t *par);
Double_t gausfunctionPSF_cool(Double_t *x, Double_t *par);
Double_t gausfunctionPSF_XY(Double_t *x, Double_t *par);
Double_t ringfunction(Double_t *x, Double_t *par);
Double_t fadding(Double_t *x, Double_t *par);
Double_t fadding2gauss(Double_t *x, Double_t *par);
Double_t fadding2Gaussian(Double_t *x, Double_t *par);
TF3 *morphospec1;
TF3 *morphospec2;

bool testCTA_v6extended(const char* filename=ctaPerformDefault,
                        double crabunits=0.1, //1.00, 
                        Double_t obsTime = OnTimeDefault,
                        Int_t whichshape=3,
                        Double_t size_source = 1.0,  //integration radius in degrees + 0.14 deg to catch the tails
                        Double_t Xpos = 0.0, //center in X
                        Double_t X_sigma = 0.0, // sigma in X
                        Double_t Ypos = 2.0, //center in Y
                        Double_t Y_sigma = 0.55, //sigma in Y
                        Double_t fluxInner = 1.0, //used for some of the morphologies to determine the flux fraction in the inner part, in units of parameter 'crabunits' 
                        Bool_t kApplyRebinning = kTRUE)
{

  Bool_t fE2=kTRUE;  //if you wish to plot E2 dN/dE set to kTRUE

  if (whichshape < 0 || whichshape > 4)
  {
    std::cout << "ERROR: no morphology specified. Select 'whichshape' from the list below:" << std::endl;
//    std::cout << " Avaliable shapes: 0 = Gaussian, 1 = Cooling with brightness - TBD, 2 = Ring, 3 = Ring+Gaussian PL Source, 4=2 Gaussians - Cooling" << std::endl;
    std::cout << "******************************************************************************************" << std::endl;
    std::cout << " PLEASE SELECT ONE OF THE FOLLOWING SHAPES:               "                                 << std::endl;
    std::cout << " shape 0: GAUSSIAN SHAPE                                  "                                 << std::endl;
    std::cout << " shape 1: SHELL SHAPE                                     "                                 << std::endl;
    std::cout << " shape 2: SHELL+Gauss SHAPE                               "                                 << std::endl;
    std::cout << " shape 3: Cooling effect                                  "                                 << std::endl;
    std::cout << " shape 4: GAUSS WITH ASYMMETRIC X and Y EXTENSIONS        "                                 << std::endl;
    std::cout << "******************************************************************************************" << std::endl;
    return 0;
  }

  if (size_source > INTMAX)
  {
    cout << "ERROR: integration radius is chosen to be " << size_source << " deg." << endl;
    cout << "This is too large. Currently values below " << INTMAX << " deg are supported only." << endl;
    return 0;
  }

  if ((Ypos*Ypos+Xpos*Xpos)>DISTMAX*DISTMAX)
  {
    cout << "ERROR: source is too far away from camera center, please use less than " << DISTMAX << " deg" << endl;
    return 0;
  }

  if (X_sigma < SIGMAMIN)
  {
    cout << "WARNING: sigma X of the source is chosen to be too small: " << X_sigma << " deg. Will assume point source using " << SIGMAMIN << " deg instead " << endl;
    X_sigma = SIGMAMIN;
  }
  else if (whichshape==4) //increase intrinsic sigma by VERY ROUGH telescope resolution
  {
      X_sigma = sqrt(X_sigma*X_sigma+SIGMAMIN*SIGMAMIN);
  }

  if (Y_sigma < SIGMAMIN)
  {
    cout << "WARNING: sigma Y of the source is chosen to be too small: " << Y_sigma << " deg. Will iassume a point source using " << SIGMAMIN << " deg instead " << endl;
    Y_sigma = SIGMAMIN;
  }
  else if (whichshape==4) //increase intrinsic sigma by VERY ROUGH telescope resolution
  {
      Y_sigma = sqrt(Y_sigma*Y_sigma+SIGMAMIN*SIGMAMIN);
  }


  //switch
  //TF3 stuff
  TF3 *morphospec = NULL;
  switch (whichshape)
  {
    case 0:                                 // ---------------------------------------------------- > Gaussian
      {
        std::cout << "******************************************************************************************" << std::endl;
        std::cout << " YOU SELECTED A GAUSSIAN SHAPE "                                                            << std::endl;
        std::cout << " parameter X_sigma " << X_sigma << " will be used for the sigma, Y_sigma unused "           << std::endl;
        std::cout << " parameter crabunits " << crabunits << " specifies the integral flux in Crab Units"         << std::endl;
        std::cout << " parameter fluxInner is unused                                                  "           << std::endl;
        std::cout << " Used spectrum: see function 'getflux'                                           "           << std::endl;
        std::cout << "******************************************************************************************" << std::endl;
        morphospec = new TF3("morphospec",gausfunctionPSF,mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax,5);  // from -5 to 5 deg in X and Y
//        morphospec->SetNpx(100000);
        morphospec->FixParameter(2,SIGMAMIN);   // 
        morphospec->FixParameter(1,X_sigma);
        morphospec->FixParameter(3,Xpos);   // 
        morphospec->FixParameter(4,Ypos);   // 
        morphospec->SetParameter(0,crabunits);
      }
        break;
    case 1:                                  // ---------------------------------------------------- > Shell
      {
        std::cout << "******************************************************************************************" << std::endl;
        std::cout << " YOU SELECTED A SHELL SHAPE    "                                                            << std::endl;
        std::cout << " parameter X_sigma " << X_sigma << " used as sigma of ExternalGauss "                       << std::endl;
        std::cout << " parameter Y_sigma " << Y_sigma << " used as sigma of InternalGauss "                       << std::endl;
        std::cout << " WARNING: Please make sure Y_sigma is smaller than X_sigma          "                       << std::endl;
        std::cout << " parameter crabunits " << crabunits << " specifies the integral flux in Crab Units"         << std::endl;
        std::cout << " parameter fluxInner is unused                                                  "           << std::endl;
        std::cout << " Used spectrum: see function 'getflux'                                           "           << std::endl;
        std::cout << "******************************************************************************************" << std::endl;
        morphospec = new TF3("morphospec",ringfunction,mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax,6);  // from -5 to 5 deg in X and Y
//        morphospec->SetNpx(100000);
        morphospec->FixParameter(2,SIGMAMIN);   // 
        morphospec->FixParameter(1,X_sigma); //ExternalGaus
        morphospec->FixParameter(3,Y_sigma); //InternalGaus
        morphospec->FixParameter(4,Xpos);   // 
        morphospec->FixParameter(5,Ypos);   // 
        morphospec->SetParameter(0,crabunits);
      }
      break;
    case 2:                              // ---------------------------------------------------- > Shell+Gauss
      {
        std::cout << "******************************************************************************************" << std::endl;
        std::cout << " YOU SELECTED A SHELL+Gauss SHAPE  "                                                        << std::endl;
        std::cout << " parameter X_sigma " << X_sigma << " used as sigma of ExternalGauss "                       << std::endl;
        std::cout << " parameter Y_sigma " << Y_sigma << " used as sigma of InternalGauss "                       << std::endl;
        std::cout << " WARNING: Please make sure Y_sigma is smaller than X_sigma           "                       << std::endl;
        std::cout << " parameter fluxInner " << fluxInner << " is unused for the strength of the inner Gauss "    << std::endl;
        std::cout << " parameter crabunits " << crabunits << " specifies the integral flux in Crab Units"         << std::endl;
        std::cout << " Used spectrum: see function 'getflux'                                           "           << std::endl;
        std::cout << "******************************************************************************************" << std::endl;
        morphospec = new TF3("morphospec",fadding,mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax,7);  // from -5 to 5 deg in X and Y
//        morphospec->SetNpx(100000);
        morphospec->FixParameter(1,SIGMAMIN);   // 
        morphospec->FixParameter(2,X_sigma); //ExternalGaus
        morphospec->FixParameter(3,Y_sigma); //InternalGaus
        morphospec->FixParameter(4,fluxInner); 
        morphospec->FixParameter(5,Xpos);   // 
        morphospec->FixParameter(6,Ypos);   // 
        morphospec->SetParameter(0,crabunits);
     }
      break;
    case 3:                          // ---------------------------------------------------- > Cooling
      {
        std::cout << "******************************************************************************************" << std::endl;
        std::cout << " YOU SELECTED A COOLING EFFECT     "                                                        << std::endl;
        std::cout << " parameter X_sigma " << X_sigma << " used as sigma of Gaus-like emission "                  << std::endl;
        std::cout << " parameter Y_sigma " << Y_sigma << " used as sigma of Colling              "                << std::endl;
        std::cout << " WARNING: Please make sure X_sigma is sufficiently smaller than Y_sigma    "                << std::endl;
        std::cout << " parameter fluxInner " << fluxInner << " is unused for the strength of the inner Gauss "    << std::endl;
        std::cout << " parameter crabunits " << crabunits << " is unused for the strength of the cooling emission "    << std::endl;
        std::cout << " Used spectrum for the inner part: see function 'getflux'                      "           << std::endl;
        std::cout << " Used spectrum for the inner part: see function 'getflux_cooling'              "           << std::endl;
        std::cout << "******************************************************************************************" << std::endl;
        morphospec1 = new TF3("morphospec1",gausfunctionPSF,mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax,5);  // from -5 to 5 deg in X and Y
        morphospec1->FixParameter(2,SIGMAMIN);   // 
        morphospec1->FixParameter(1,X_sigma);
        morphospec1->FixParameter(3,Xpos);   // 
        morphospec1->FixParameter(4,Ypos);   // 
        morphospec1->SetParameter(0,fluxInner);

        morphospec2 = new TF3("morphospec2",gausfunctionPSF_cool,mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax,5);  // from -5 to 5 deg in X and Y
        morphospec2->FixParameter(2,SIGMAMIN);   // 
        morphospec2->FixParameter(1,Y_sigma);
        morphospec2->FixParameter(3,Xpos);   // 
        morphospec2->FixParameter(4,Ypos);   // 
        morphospec2->SetParameter(0,crabunits);

        morphospec = new TF3("morphospec",fadding2Gaussian,mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax);  // from -5 to 5 deg in X and Y
//        morphospec->FixParameter(2,SIGMAMIN);   // 
//        morphospec->FixParameter(1,X_sigma);
//        morphospec->FixParameter(3,Xpos);   // 
//        morphospec->FixParameter(4,Ypos);   // 
//        morphospec->SetParameter(0,crabunits);
//        morphospec->FixParameter(2+5,SIGMAMIN);   // 
//        morphospec->FixParameter(1+5,Y_sigma);
//        morphospec->FixParameter(3+5,Xpos);   // 
//        morphospec->FixParameter(4+5,Ypos);   // 
//        morphospec->SetParameter(0+5,fluxInner);
//        morphospec->SetNpx(100000);
//        morphospec->FixParameter(1,X_sigma);   // 
//        morphospec->FixParameter(2,Y_sigma); //
//        morphospec->FixParameter(3,fluxInner); // in untins of par[0]
//        morphospec->FixParameter(4,Xpos);   // 
//        morphospec->FixParameter(5,Ypos);   // 
//        morphospec->SetParameter(0,crabunits);
      }
      break;

    case 4:                          // ---------------------------------------------------- > Gaus with asymmetric X and Y
      {
        std::cout << "******************************************************************************************" << std::endl;
        std::cout << " YOU SELECTED A GAUSS WITH ASYMMETRIC X and Y EXTENSIONS  "                                 << std::endl;
        std::cout << " parameter X_sigma " << X_sigma << " used as sigma in X direction  "                      << std::endl;
        std::cout << " parameter Y_sigma " << Y_sigma << " used as sigma in Y direction  "                      << std::endl;
        std::cout << " parameter crabunits " << crabunits << " specifies the integral flux in Crab Units"         << std::endl;
        std::cout << " Used spectrum: see function 'getflux'                                           "           << std::endl;
        std::cout << " parameter fluxInner is unused                                     "                      << std::endl;
        std::cout << "******************************************************************************************" << std::endl;
        morphospec = new TF3("morphospec",gausfunctionPSF_XY,mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax,5);  
//        morphospec->SetNpx(100000);
        // par 0 : normalization
        // par 1: X0
        // par 2: sigma X
        // par 3: Y0
        // par 4: sigma Y
        morphospec->SetParameter(0,crabunits);
        morphospec->FixParameter(1,Xpos);
        morphospec->FixParameter(2,X_sigma);
        morphospec->FixParameter(3,Ypos);
        morphospec->FixParameter(4,Y_sigma);
      }

      break;

    default:   // -1 :-)
        std::cout << "******************************************************************************************" << std::endl;
        std::cout << " PLEASE SELECT ONE OF THE FOLLOWING SHAPES:               "                                 << std::endl;
        std::cout << " shape 0: GAUSSIAN SHAPE                                  "                                 << std::endl;
        std::cout << " shape 1: SHELL SHAPE                                     "                                 << std::endl;
        std::cout << " shape 2: SHELL+Gauss SHAPE                               "                                 << std::endl;
        std::cout << " shape 3: Cooling effect                                  "                                 << std::endl;
        std::cout << " shape 4: GAUSS WITH ASYMMETRIC X and Y EXTENSIONS        "                                 << std::endl;
        std::cout << "******************************************************************************************" << std::endl;
      break;
  }

  //make an array of low edges (in log) between emin and emax
  const int nbins = 15; //18; //32; //18; //30;
  const double emin = 0.05; //0.05; //1.122e-02; //0.01122; // 0.05; //TeV
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

//  cout << " chosen function is : " << endl;
//  cout << morphospec->GetExpFormula() << endl;

//  TF1 *origSpec = new TF1("origSpec",specstrS,mapEmin,mapEmax);  
  TF1 *origSpec;
  TF1 *origSpec1=NULL;
  TF1 *origSpec2=NULL; 
  if (whichshape!=3)
    origSpec = new TF1("origSpec","getflux(x)*[0]",mapEmin,mapEmax); 
  else
  {
    origSpec1 = new TF1("origSpec1","[0]*getflux(x)",mapEmin,mapEmax); 
    origSpec2 = new TF1("origSpec2","[0]*getflux_cooling(x)",mapEmin,mapEmax); 
    origSpec = new TF1("origSpec","[0]*getflux(x)+[1]*getflux_cooling(x)",mapEmin,mapEmax); 
  }
  origSpec->SetNpx(10000);
  origSpec->SetParameter(0,crabunits*BOOST);
  if (whichshape==3)
  {
    origSpec1->SetParameter(0,fluxInner*BOOST);
    origSpec2->SetParameter(0,crabunits*BOOST);
  }
  origSpec->SetLineWidth(2);
  origSpec->SetLineColor(1);

  double fluxIntegral = morphospec->Integral(mapXmin,mapXmax,mapYmin,mapYmax,emin, emax);
  double flb = origSpec->Integral(emin,emax);
  double fluxIntegral1, fluxIntegral2, flb1, flb2;
  if (whichshape==3) {
    fluxIntegral1 = morphospec1->Integral(mapXmin,mapXmax,mapYmin,mapYmax,emin, emax);
    fluxIntegral2 = morphospec2->Integral(mapXmin,mapXmax,mapYmin,mapYmax,emin, emax);
    flb1 = origSpec1->Integral(emin,emax);
    flb2 = origSpec2->Integral(emin,emax);
  }
//  cout << " check prior " << fluxIntegral << " UUUUUUUU  " << flb << endl;
//  cout << morphospec1->Integral(mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax) << endl;
//  cout << morphospec2->Integral(mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax) << endl; 
//  cout << morphospec1->Integral(mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax)+morphospec2->Integral(mapXmin,mapXmax,mapYmin,mapYmax,mapEmin,mapEmax) << endl; 
  if (whichshape!=3)
    morphospec->FixParameter(0,crabunits*flb/fluxIntegral/BOOST);
  else
  {
    morphospec1->FixParameter(0,fluxInner*flb1/fluxIntegral1/BOOST);
    morphospec2->FixParameter(0,crabunits*flb2/fluxIntegral2/BOOST);
//    cout << Form(" check prior 2: flb1 %.2e, flb2 %.2e, p1 %.2e, p2 %.2e\n",flb1,flb2,fluxInner*flb1/fluxIntegral1/BOOST,crabunits*flb2/fluxIntegral2/BOOST) ;
  }

  if (whichshape==2)
    morphospec->FixParameter(4,fluxInner*flb/fluxIntegral/BOOST);
//  else if (whichshape==3)
//    morphospec->FixParameter(3,fluxInner*flb/fluxIntegral/BOOST);

//  morphospec2->FixParameter(0,10.);

  if (whichshape==3)
  {
    origSpec->SetParameter(0,fluxInner);
    origSpec->SetParameter(1,crabunits);
  }
  else
    origSpec->SetParameter(0,crabunits);

//  double fluxIntegral_check = morphospec->Integral(mapXmin,mapXmax,mapYmin,mapYmax,emin,emax);
//  double flb_check = origSpec->Integral(emin,emax);
//  cout << " check integral: " << fluxIntegral_check << " " << flb_check << endl;

  TF1 *origSpecE2 = new TF1("origSpecE2","origSpec*x*x",mapEmin,mapEmax);  
  origSpecE2->SetNpx(10000);
  origSpecE2->SetLineWidth(2);
  origSpecE2->SetLineColor(1);

  TCanvas *c1 = new TCanvas("spectrum","my CTA spectrum",1200,600);
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
  c1->Divide(2,1);
  gStyle->SetOptStat(0);


//  double fluxIntegral_check2 = morphospec->Integral(mapXmin,mapXmax,mapYmin,mapYmax,emin,emax);
//  double flb_check2 = origSpec->Integral(emin,emax);
//  cout << " check integral: " << fluxIntegral_check2 << " " << flb_check2 << endl;

  TH1D *spObserved = new TH1D("spObserved","",nbins,lowedges);
  spObserved->Sumw2();
  spObserved->GetXaxis()->SetTitle("energy E(TeV)");
  spObserved->GetXaxis()->SetTitleOffset(1.2);
  spObserved->SetMaximum(1e-9);
  spObserved->SetMinimum(1e-14);
  spObserved->GetYaxis()->SetTitleOffset(1.4);
//  spObserved->GetXaxis()->SetMoreLogLabels();
  spObserved->GetXaxis()->SetNoExponent();
  spObserved->GetXaxis()->SetLabelOffset(0.015);
  spObserved->GetYaxis()->SetTitle("dN/dE (ph TeV^{-1} cm^{-2} s^{-1})");
  if (fE2)
    spObserved->GetYaxis()->SetTitle("E^{2} dN/dE (TeV cm^{-2} s^{-1})");

  spObserved->SetLineWidth(2);
  spObserved->SetLineColor(1);
  spObserved->SetMarkerStyle(20);

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

  bool rc = makeCTAspec_extended(spObserved,specgr,ifluxgr,excessgraph,filename,kApplyRebinning,morphospec,obsTime,size_source,Xpos,Ypos, &threshold);
  if (!rc) {
    cout << " errors detected in makeCRAspec routine, aborting ... " << endl;
    return 0;
  }

  //fill 2D histogram with fluxes
  TH2D *myflux = new TH2D("myflux",Form("ideal fluxmap from %.3lf TeV to %.3lf TeV",emin,emax),NBINMORPHO,mapXmin,mapXmax,NBINMORPHO,mapYmin,mapYmax);
  myflux->SetXTitle("x(deg)");
  myflux->SetYTitle("y(deg)");
  myflux->SetZTitle("Integral flux (cm^{-2} s^{-1})");
  double x1,x2,x3,x4;
  morphospec->SetNpx(10);
  morphospec->SetNpy(10);
  morphospec->SetNpz(100);
  for (int ik = 1; ik<= myflux->GetNbinsX(); ik++)
  {
      x1 = myflux->GetXaxis()->GetBinLowEdge(ik);
      x2 = myflux->GetXaxis()->GetBinLowEdge(ik+1);
      for (int jk = 1; jk<= myflux->GetNbinsY(); jk++)
      {
          x3 = myflux->GetYaxis()->GetBinLowEdge(jk);
          x4 = myflux->GetYaxis()->GetBinLowEdge(jk+1);
          myflux->SetBinContent(ik,jk,morphospec->Integral(x1,x2,x3,x4,emin,emax));
      }
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

  c1->cd(1);
  gPad->SetLeftMargin(0.12);
  TH2D *sp;
  if(fE2) sp = new TH2D("sp","",1000,0.01,100,1000,1e-14,1e-9);
  else  sp = new TH2D("sp","",1000,0.01,100,1000,1e-17,1e-6);
  sp->Sumw2();
  sp->GetXaxis()->SetTitle("energy E(TeV)");
  sp->GetXaxis()->SetTitleOffset(1.2);
  sp->GetYaxis()->SetTitleOffset(1.5);
//  sp->GetXaxis()->SetMoreLogLabels();
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

    cout << endl << " END 1 " << endl; 

  TF1 *origSpecCrab=NULL;
  if(fE2) origSpecCrab = new TF1("origSpecCraab","(6.0e-10)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))*x*x",mapEmin,mapEmax);  //Crab Nebula (MAGIC paper) with [0]=1. 
  else origSpecCrab = new TF1("origSpecCraab","(6.0e-10)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))",mapEmin,mapEmax);  //Crab Nebula (MAGIC paper) with [0]=1. 
  origSpecCrab->SetLineStyle(7);
  origSpecCrab->SetLineColor(36);
  origSpecCrab->Draw("same");

  TLegend *leg = new TLegend(0.50,0.82,0.99,0.99);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  if (fE2) leg->AddEntry(origSpecE2,"my source, intrinsic spectrum","l");
  else leg->AddEntry(origSpec,"my source, intrinsic spectrum","l");
//  leg->AddEntry(spObserved,Form("expected spectrum (CTA), %.2lfh",obsTime/3600),"pl");
  leg->AddEntry(specgr,Form("expected spectrum (CTA), %.2lfh",obsTime/3600),"pl");
  leg->AddEntry(origSpecCrab,"Crab (MAGIC)","pl");
  leg->Draw();

  TF1 * pl = new TF1("pl","[0]*pow(x,[1])",mapEmin,mapEmax);
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

  cout << "INFO: F_real_2 is " << morphospec->Integral(mapXmin,mapXmax,mapYmin,mapYmax,emin,emax) << endl; 
  cout << "INFO: F_real_3 is " << morphospec->Integral(-2,2,mapYmin,mapYmax,emin,emax) << endl; 
  cout << "INFO: F_real_4 is " << morphospec->Integral(-2,2,0,4,emin,emax) << endl; 

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

  gStyle->SetPalette(60);
  c1->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetRightMargin(0.17);
  gPad->SetFillColor(10);
  gPad->SetBorderMode(0);
  gPad->SetBorderSize(0);
  gPad->SetLogz();
  double maxf = myflux->GetMaximum();
  myflux->SetMinimum(maxf*1e-3);
  myflux->GetZaxis()->SetTitleOffset(1.5);
  myflux->Draw("colz");
//  morphospec->Draw("lego2");
  double r = size_source + sqrt(2.)*BINSIZEMORPHO/2;
  TEllipse *elli = new TEllipse(Xpos,Ypos,r);
  elli->SetFillStyle(0);
  elli->SetLineColor(kOrange+2);
  elli->SetLineWidth(3);
  elli->Draw();

  TMarker *ma = new TMarker(Xpos,Ypos,5);
  ma->Draw();

  return 1;


}

// **************************************************************************************************************
// Power Law function
// **************************************************************************************************************
Double_t pwlfunction(Double_t* x, Double_t* par)
{
  return (Double_t) par[0]*TMath::Power(x[0],-par[1]);
}
// **************************************************************************************************************
// Power Law + exp function
// **************************************************************************************************************
Double_t expfunction(Double_t* x, Double_t* par)
{
  return (Double_t) par[0]*TMath::Power(x[0],-par[1])*exp(-x[0]/par[2]);
}
// **************************************************************************************************************
// Composite 
// **************************************************************************************************************
Double_t fadding(Double_t *x, Double_t *par)
{
  Double_t X=x[0];
  Double_t Y=x[1];
  Double_t E=x[2];
//  Double_t Xpos=0; // position in camera coordinates
//  Double_t Ypos=0; 
  Double_t spl=par[1]; // point like source
  Double_t sinner=par[3];
  Double_t souter=par[2];
  Double_t enterm = getflux(E);
  Double_t r2 = pow(X-par[5],2)+pow(Y-par[6],2);
  Double_t z  = ((par[0]*(exp(-r2/(2*souter*souter))-exp(-r2/(2*sinner*sinner))))+par[4]*exp(-r2/(2*spl*spl)))*enterm;
  return z;
}
// **************************************************************************************************************
// Adding two Gaussian 
// **************************************************************************************************************
Double_t fadding2gauss(Double_t *x, Double_t *par)
{
  Double_t X=x[0];
  Double_t Y=x[1];
  Double_t E=x[2];
//  Double_t Xpos=0; // position in camera coordinates
//  Double_t Ypos=0; 
  Double_t s1=par[1];
  Double_t s2=par[2];
  Double_t r2 = pow(X-par[4],2)+pow(Y-par[5],2);
  Double_t enterm = getflux(E);
  Double_t z  = (par[0]*(exp(-r2/(2*s1*s1))) + par[3]*exp(-r2/(2*s2*s2)))*enterm;
  return z;
}
// **************************************************************************************************************
// Shell-type function 
// **************************************************************************************************************
Double_t ringfunction(Double_t *x, Double_t *par)
{
  
  Double_t X=x[0];
  Double_t Y=x[1];
  Double_t E=x[2];
//  Double_t Xpos=0; // position in camera coordinates
//  Double_t Ypos=0; 
  Double_t s=par[1];
  Double_t sinner=par[3];
  Double_t spsf=par[2];
  Double_t r2 = pow(X-par[4],2)+pow(Y-par[5],2);

  Double_t enterm = getflux(E);

  Double_t z  = par[0]*(exp(-r2/(2*((s*s)+(spsf*spsf))))-exp(-r2/(2*((sinner*sinner)+(spsf*spsf)))))*enterm;
    return z;
}
// **************************************************************************************************************
// Gaus Function and PSF 
// **************************************************************************************************************
Double_t gausfunctionPSF(Double_t *x, Double_t *par)
{
  
  Double_t X=x[0];
  Double_t Y=x[1];
  Double_t E=x[2];
//  Double_t Xpos=0; // position in camera coordinates
//  Double_t Ypos=0; 
  Double_t s=par[1];
  Double_t spsf=par[2];
  Double_t r2 = pow(X-par[3],2)+pow(Y-par[4],2);

  Double_t enterm = getflux(E);

  Double_t nenner = 2*((s*s)+(spsf*spsf));

  Double_t z  = par[0]*exp(-r2/nenner)*enterm;

//  cout << "gausfunctionPSF " << r2 << " " << E << " " << enterm << " " << z << endl;
//  if (r2<0.01)
//    cout << Form("gausfunctionPSF [0]=%lf [1]=%.2lf [2]=%.2lf [3]=%.1lf [4]=%.1lf r2=%e nen=%e e=%lf enterm=%e z=%e",par[0],par[1],par[2],par[3],par[4],r2,nenner,x[2],enterm,z) << endl;
  return z;
}

// **************************************************************************************************************
// Gaus Function and PSF  and cooling spectrum
// **************************************************************************************************************
Double_t gausfunctionPSF_cool(Double_t *x, Double_t *par)
{
  
  Double_t X=x[0];
  Double_t Y=x[1];
  Double_t E=x[2];
//  Double_t Xpos=0; // position in camera coordinates
//  Double_t Ypos=0; 
  Double_t s=par[1];
  Double_t spsf=par[2];
  Double_t r2 = pow(X-par[3],2)+pow(Y-par[4],2);

  Double_t enterm = getflux_cooling(E);

  Double_t z  = par[0]*exp(-r2/(2*((s*s)+(spsf*spsf))))*enterm;
//  cout << "gausfunctionPSF_cool " << z << endl;
    return z;
}

Double_t fadding2Gaussian(Double_t *x, Double_t *par)
{
  
  return (morphospec1->EvalPar(x,par) + morphospec2->EvalPar(x,par));
//  return (morphospec2->EvalPar(x,par) + morphospec1->EvalPar(x,par));
//  double a = morphospec1->Eval(x[0],x[1],x[2]);
//  double b = morphospec2->Eval(x[0],x[1],x[2]);
////  cout << x[2] << " " <<  a << " " << b << endl;
//  return (a+b);
}

// **************************************************************************************************************
// Gauss Function with different sigmas in X and Y
// **************************************************************************************************************
Double_t gausfunctionPSF_XY(Double_t *x, Double_t *par)
{
  
  Double_t X=x[0];
  Double_t Y=x[1];
  Double_t E=x[2];
//  Double_t Xpos=0; // position in camera coordinates
//  Double_t Ypos=0; 
  Double_t r2 = 0.5*pow((X-par[1])/par[2],2)+0.5*pow((Y-par[3])/par[4],2);

  Double_t enterm = getflux(E);

  Double_t z  = par[0]*exp(-r2)*enterm;
    return z;
}
