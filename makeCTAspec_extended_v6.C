#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPad.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TF3.h>
#include <TString.h>
#include <TSpline.h>
#include <TArrayD.h>
#include <TRandom3.h>

/*
 * =====================================================================================
 *
 *       Filename:  makeCTAspec_extended.C
 *
 *    Description:  
 * The program  gets as _input_:
 *   (1) Intrinsic spectral shape and flux
 *   (2) Distance //obsolete, the distance for extragal. sources is provided through (5)
 *   (3) CTA effective area
 *   (4) Observation time
 *   (5) EBL attenuation file, i.e. taus must be calculated by a different software.
 *   (6) Binning of the final spectrum // provided through a TH1D histogram

 * _Output_ is:
 *   (i)  Binned spectrum (events, flux) as TH1D (obsolete)
 *   (ii) Binned spectrum (events, flux) as TGraphAsymmErrors
 *   (iii) Integral flux and corresponding error in a given energy range (the range is the same as for the energy spectrum), aslo as TGraphAsymmErrors
 *
 *        Version:  6.0
 *        Created:  04/04/11 
 *
 *         Author:  Daniel Mazin (dm), mazin@ifae.es
 *        Company:  IFAE, Spain                
 *
 * =====================================================================================
 */

#define GAUSL 20   //limit f the gauss statistics. if lower, asymmetric will be used DO NOT USE LARGER THAN 25
#define NBINS 151  //maximum number of energy bins
#define MINEVT 7    //min events
#define MINAREA 1e4  //min area in cm2 (used to decide if rebinning should be done)
#define MINSIG 3.0  //min significance 
#define MINBKG 0.03  //min above background level (in %)
#define SAVEEMIN 0.029 //min energy in TeV above which spectra can be made
#define NEVT 10000    //number of events per energy bin (Etrue) of migration matrix
#define SMALL 1e-3   //used for minimum background
#define MINETAU  0.01 //min energy for which tau is available
#define MAXETAU  50.  //max energy for which tau is available
#define BOOST 1 // 1000 //boost factor to make TF1::Integral() method work better. May need to change it

#define NUMPOINT 200

#define NBINMORPHO 51  
#define BINSIZEMORPHO 0.2  

#define kNODEBUG

//defaults
TF3 *spIntrDefault = new TF3("spIntrDefault","(6.0e-11)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))*exp(-(z*z/2.+y*y/2))",0.01,100.,-3.,3.,-3.,3.);  //Crab Nebula (MAGIC paper)
//TF1 *spIntrDef = new TF1("spIntrDef","(6.0e-11)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))",0.01,200);
Double_t threshDefault = 0.1;

double effOnTimeDefault = 20.*60*60; //20 hours //The number is in seconds

Double_t xMin=0.01; //TeV
Double_t xMax=100.; //TeV

Double_t XSourceDefault = 0.0; //source position, deg
Double_t YSourceDefault = 0.0; //source position, deg
Double_t RSourceDefault = 1.0; //integraton radius, deg

const int minPointNotRebinned = 3; //min number of spectral points which must be valid before trying to rebin

  //background function:
//differential background in 50h of observation (x in TeV)
TF1 *bg = new TF1("bg","300*pow(x,-3.2)",xMin,xMax); //in events
TF1 *bgSpec = new TF1("bgSpec","1.5e-13*pow(x,-3.2)",xMin,xMax); //in ph/cm2/s/TeV 

//flag between BG from spectrum and fixed BG in events
//Bool_t kUseBGFromSpec = kTRUE;

//flag between point-like source and extended source 
Bool_t kUseExtended = kTRUE;

//flag between BG from spectrum and fixed BG in events
Bool_t kUseRandom = kTRUE;

//flag if you want to apply spill-over correction (default kFALSE)
Bool_t kUseSpillOver = kFALSE; //obsolete

double efficgampsf = 1.2; //scale factor to account for efficiency of the theta2 cut in case of an extended source
double alpha = 0.2;  // 0.2 = 5 times better bg estimation; in Li&Ma alpha = t_on / t_off
const double fT50 = 50.*60*60; //sec in 50 hours
const int fIntSteps = 3;
Double_t fFractionArea = 0.00; //0.1; //relative uncertainty of the collection area
Double_t fFractionSystBkg = 0.01; // NOT USED

TRandom3 rnd(2749);
TRandom3 rnd2(5427);

Double_t SignificanceLiMa(Double_t s, Double_t b, Double_t alpha) //taking from the standard MAGIC Mars code
{
    const Double_t sum = s+b;

    if (s==0 || b==0 || sum==0)
        return 0;

    if (sum<0 || alpha<=0)
        return -1;

    const Double_t l = s*TMath::Log(s/sum*(alpha+1)/alpha);
    const Double_t m = b*TMath::Log(b/sum*(alpha+1)      );

    return l+m<0 ? -1 : TMath::Sqrt((l+m)*2);
}

Double_t IntegrateFromSpline(TSpline3 *spline, double Xlow, double Xhigh)
{
  //simpson rule:
  // Int_a^b f(x) dx = (b-a)/6*(f(a) + 4*f((a+b)/2) + f(b))

  if (Xlow >= Xhigh) 
    return (1e-16);

  double x, h, h2, s;
  double F1_x, F1_x_h, F1_a, F1_b;

  int n2 = 2 * fIntSteps;
  double integral = 0.;

  h = (Xhigh - Xlow) / (double)n2;
  h2 = h + h;
  s = 0.0;
  x = Xlow + h;

  for(int j = 0; j < fIntSteps; j++) {

    F1_x = spline->Eval(x);
    F1_x_h = spline->Eval(x+h);

    s = 4. * F1_x + 2. * F1_x_h + s;
    x = x + h2;
  }


  F1_a = spline->Eval(Xlow);
  F1_b = spline->Eval(Xhigh);

  s = F1_a + s - F1_b;
//  cout << " F1_a, F1_b, s  " << F1_a << "  " << F1_b << "  " << s ;

  if(s < 1e-4) {
//    std::cout << "WARNING <MREBLIntegrator2::SIMPSRedshift> " << std::endl
//      << "    s = " << s << " Endsumme bei Simpson3 ist kleiner als 1e-4" << std::endl;
    s=0.;
  }

  integral = s * h/3.;

  if (integral<0)
    integral=0.;

  return integral;

}

Double_t IntegrateFromSplineLog(TSpline3 *spline, double Xlow, double Xhigh)
{
  //simpson rule:
  // Int_a^b f(x) dx = (b-a)/6*(f(a) + 4*f((a+b)/2) + f(b))

  if (Xlow >= Xhigh) 
    return (1e-16);

  double x, h, h2, s;
  double F1_x, F1_x_h, F1_a, F1_b;

  int n2 = 2 * fIntSteps;
  double integral = 0.;

  h = (Xhigh - Xlow) / (double)n2;
  h2 = h + h;
  s = 0.0;
  x = Xlow + h;

  for(int j = 0; j < fIntSteps; j++) {

    F1_x = spline->Eval(x) * log(10) * pow(10,x);
    F1_x_h = spline->Eval(x+h) * log(10) * pow(10,x+h);

    s = 4. * F1_x + 2. * F1_x_h + s;
    x = x + h2;
  }


  F1_a = spline->Eval(Xlow) * log(10.) * pow(10,Xlow);
  F1_b = spline->Eval(Xhigh) * log(10.) * pow(10,Xhigh);

  s += F1_a - F1_b;

  if(s < 1e-4) {
    s=0.;
  }

  integral = s * h/3.;

  if (integral<0)
    integral=0.;

  return integral;

}

Double_t IntegrateFromGraphLog(TGraph *graph, double Xlow, double Xhigh)
{
  //simpson rule:
  // Int_a^b f(x) dx = (b-a)/6*(f(a) + 4*f((a+b)/2) + f(b))

//  cout << " in the IntegrateFromGraphLog ..." << endl;
//  cout << " %d points " << graph->GetN() << endl;
  for (int g=0; g<graph->GetN();g++)
  {
    double x,y;
    graph->GetPoint(g,x,y);
//    cout << Form("Graph, point %d, X=%.2e, Y=%.2e",g,x,y) << endl;
  }

  if (Xlow >= Xhigh) 
    return (1e-16);

  double x, h, h2, s;
  double F1_x, F1_x_h, F1_a, F1_b;

  int n2 = 2 * fIntSteps;
  double integral = 0.;

  h = (Xhigh - Xlow) / (double)n2;
  h2 = h + h;
  s = 0.0;
  x = Xlow + h;

  for(int j = 0; j < fIntSteps; j++) {

    F1_x = graph->Eval(x) * log(10) * pow(10,x);
    F1_x_h = graph->Eval(x+h) * log(10) * pow(10,x+h);

    s = 4. * F1_x + 2. * F1_x_h + s;
    x = x + h2;
  }


  F1_a = graph->Eval(Xlow) * log(10.) * pow(10,Xlow);
  F1_b = graph->Eval(Xhigh) * log(10.) * pow(10,Xhigh);

  s += F1_a - F1_b;

  if(s < 1e-4) {
    s=0.;
  }

  integral = s * h/3.;

  if (integral<0)
    integral=0.;

  return integral;

}

Double_t IntegrateForTau(TSpline3 *spline, TF1 *flux, double Xlow, double Xhigh)
{
  //simpson rule:
  // Int_a^b f(x) dx = (b-a)/6*(f(a) + 4*f((a+b)/2) + f(b))

  if (Xlow < MINETAU && Xhigh < MINETAU)
    return 1.;
  if (Xlow > MAXETAU && Xhigh > MAXETAU)
    return 0.;

  if (Xlow >= Xhigh) 
    return (1e-16);

  double x, h, h2, s, s2;
  double F1_x, F1_x_h, F1_a, F1_b;
  double F2_x, F2_x_h, F2_a, F2_b;
  double atten;

  int n2 = 2 * fIntSteps;

  h = (Xhigh - Xlow) / (double)n2;
  h2 = h + h;
  s = 0.0;
  s2 = 0.0;
  x = Xlow + h;

  for(int j = 0; j < fIntSteps; j++) {

    //zaehler:
    //WORKAROUND:
    if (x < MINETAU)
      atten = 1.;
    else if (x > MAXETAU)
      atten = 0.;
    else
      atten = spline->Eval(x);

    if (atten<0.)
      return 0.;

    F1_x = atten * flux->Eval(x);

    if (F1_x<=0.)
      return 0.;

    if ((x+h) < MINETAU)
      atten = 1.;
    else if ((x+h) > MAXETAU)
      atten = 0.;
    else
      atten = spline->Eval(x+h);

    F1_x_h = atten * flux->Eval(x+h);

    if (F1_x_h<=0.)
      return 0.;

    s = 4. * F1_x + 2. * F1_x_h + s;

    //nenner:
    F2_x = flux->Eval(x);
    F2_x_h = flux->Eval(x+h);

    s2 = 4. * F2_x + 2. * F2_x_h + s2;

    //increase x
    x = x + h2;
  }

  if (Xlow < MINETAU)
    atten = 1.;
  else if (Xlow > MAXETAU)
    atten = 0.;
  else
    atten = spline->Eval(Xlow);
  F1_a = atten * flux->Eval(Xlow);

  if (Xhigh < MINETAU)
    atten = 1.;
  else if (Xhigh > MAXETAU)
    atten = 0.;
  else
    atten = spline->Eval(Xhigh);
  F1_b = atten * flux->Eval(Xhigh);

  s = F1_a + s - F1_b;

  F2_a = flux->Eval(Xlow);
  F2_b = flux->Eval(Xhigh);
  s2 = F2_a + s2 - F2_b;

//  cout << endl << " F1_a, F1_b, s  " << F1_a << "  " << F1_b << "  " << s ;
//  cout << " F2_a, F2_b, s2  " << F2_a << "  " << F2_b << "  " << s2 << endl << endl;

  if(s < 0.) {
    s=0.;
  }
  if (s2<0)
    s2=1e-7;  //small number

  double result = s/s2;

  if (result<0)
    result=0.;

  return result;
}


//have to provide:
//1. center of the source (x,y in degree): XSource, YSource
//2. integration radius (in degree): RSource

//main function
bool makeCTAspec_extended(TH1D *spObserved, //expected spectrum to be observed with CTA
                          TGraphAsymmErrors *specgraph, //the same as above
                          TGraphAsymmErrors *ifluxgraph,  //store here the integral flux and its error
                          TGraphAsymmErrors *excessgraph, //distribution of excess events
                          const char* filename,          //root format response file
                          Bool_t kRebin=kFALSE,         //should rebinning be applied if the number of excess events too low?
                          TF3 *spIntrMain=spIntrDefault,   //function of the source shape
                          Double_t effOnTime=effOnTimeDefault, // observation time
                          double RSource = RSourceDefault, //integration radius, deg
                          double XSource = XSourceDefault,
                          double YSource = YSourceDefault,
                          Double_t *threshold = &threshDefault)
{
//  double par0 = spIntrMain->GetParameter(0);
//  spIntrMain->SetParameter(0,par0);

  //poisson errors for small numbers (up to 25)
  double poisUp[] = {1.0, 1.4314, 1.793, 2.0894, 2.3452, 
                     2.5732, 2.7808, 2.9726, 3.1518,  3.3204,
                     3.4802, 3.6324, 3.778, 3.9178, 4.0526,
                     4.1826,  4.3084, 4.4306, 4.5492, 4.6646, 
                     4.777,  4.8866,  4.9938, 5.0986, 5.201,  5.3016}; 
  double poisLo[] = {0.0, 0.8486, 1.2268, 1.5288, 1.7872,
                     2.017, 2.2256, 2.4182, 2.5978, 2.767,    
                     2.9272, 3.0798, 3.2256, 3.3656, 3.5006,
                     3.6308, 3.7568, 3.879, 3.9976, 4.1132,
                     4.2256, 4.3354, 4.4426, 4.5474, 4.65, 4.7506};  

//  rnd.SetSeed(0);
//  rnd.SetSeed(0);
//  double H72 = 72./72.;  //change second number if you want to change the Hubble constant!!


  cout << "Entering the makeCTAspec (version for extended sources) macro ... " << endl;

  if (XSource*XSource+YSource*YSource+RSource*RSource > pow(NBINMORPHO*BINSIZEMORPHO/2.,2) )
  {
    cout << "ERROR: source integration radius R = " << RSource;
    cout << " combine with the source position X = " << XSource << " and Y = " << YSource;
    cout << " is too large " << endl;
    cout << "Tipp: choose X*X+Y*Y+R*R < " << pow(NBINMORPHO*BINSIZEMORPHO/2.,2)  << endl;
    return 0;
  }

//  double fEfficgampsf = 1.2; //efficiency of gamma-rays for an extended source

  //create fMorpho: X,Y histogram with dimensions of skymap (currently -2, 2 degrees)
  TH2F *fMorpho = new TH2F ("fMorpho","Morpho", NBINMORPHO, -1*NBINMORPHO*BINSIZEMORPHO/2.,NBINMORPHO*BINSIZEMORPHO/2.,NBINMORPHO, -1*NBINMORPHO*BINSIZEMORPHO/2.,NBINMORPHO*BINSIZEMORPHO/2.);

  //open file
  TFile* f = NULL;
  if (!(f = TFile::Open(filename)))
  {
    cout << "ERROR: could not open root performance file " << filename << endl;
    return 0;
  }
  
  //get histograms
  TH2* effColArea = NULL;
  if (!(effColArea = (TH2*)f->Get("EffectiveArea_offaxis")))
  {
    cout << "ERROR: did not find histogram EffectiveArea in the provided root performance file " << filename << endl;
    return 0;
  }

  TH2* effColAreaTrue = NULL;
  if (!(effColAreaTrue = (TH2*)f->Get("EffectiveAreaEtrue")))
  {
    cout << "ERROR: did not find histogram EffectiveAreaEtrue in the provided root performance file " << filename << endl;
    return 0;
  }

  TH2* effColArea80 = NULL;
  if (!(effColArea80 = (TH2*)f->Get("EffectiveArea80_offaxis")))
  {
    cout << "ERROR: did not find histogram EffectiveArea80 in the provided root performance file " << filename << endl;
    return 0;
  }

//  TH1* res = (TH1*)f->Get("AngRes");
  TH2* res = NULL;
  if (!(res = (TH2*)f->Get("AngRes80_offaxis")))
  {
    cout << "ERROR: did not find histogram AngRes80 in the provided root performance file " << filename << endl;
    return 0;
  }
  TH2* bgdeg = NULL;
  if (!( bgdeg = (TH2*)f->Get("BGRatePerSqDeg_offaxis")))
  {
    cout << "ERROR: did not find histogram BGRatePerSqDeg in the provided root performance file " << filename << endl;
    return 0;
  }
//  TH1* bgrate = (TH1*)f->Get("BGRatePerSqDeg");
  TH2* bgrate = NULL;
  if (!(bgrate = (TH2*)f->Get("BGRate_offaxis")))
  {
    cout << "ERROR: did not find histogram BGRate in the provided root performance file " << filename << endl;
    return 0;
  }
  TH2* enres = NULL;
  if (!(enres = (TH2*)f->Get("ERes")))
  {
    cout << "ERROR: did not find histogram ERes in the provided root performance file " << filename << endl;
    return 0;
  }
  TH3F* migmat = NULL;
  if (!(migmat = (TH3F*)(f->Get("MigMatrix")->Clone("mymigmat"))))
  {
    cout << "ERROR: did not find migration matrix in the provided root performance file " << filename << endl;
    return 0;
  }

  //get number of bins and bin width 
  //for offset in degrees from the camera center
  static const int nbdeg = enres->GetNbinsY();
  static const float degbin = enres->GetYaxis()->GetBinWidth(1);

  //check min and max energies
  double emin = pow(10,bgrate->GetXaxis()->GetBinLowEdge(1));
  double emax = pow(10,bgrate->GetXaxis()->GetBinLowEdge(bgrate->GetNbinsX()+1));
  double useremin = spObserved->GetBinLowEdge(1);
  double useremax = spObserved->GetBinLowEdge(spObserved->GetNbinsX()+1);
  if (emin > useremin)
  {
    cout << "ERROR: min energy provided by user " << useremin << " TeV is below min energy in performance file: " << emin << " TeV" << endl;
    cout << "Tip: choose min energy above " << emin << " TeV" << endl;
    return 0;
  }
  else if (SAVEEMIN > useremin)
  {
    cout << "ERROR: min energy provided by user " << useremin << " TeV is below SAVEEMIN for spectrum and lightcurve (" << SAVEEMIN << " TeV)" << endl;
    cout << "Tip: either choose min energy above " << SAVEEMIN << " TeV or switch off the energy resolution option" << endl;
    return 0;
  }
  if (emax < useremax)
  {
    cout << "ERROR: max energy provided by user " << useremax << " TeV is above max energy in performance file: " << emax << " TeV" << endl;
    cout << "Tip: choose max energy below " << emax << " TeV" << endl;
    return 0;
  }

  //energy of the intrinsic spectrum _spIntr_ is in the observer frame!!!!

  //prepare loop over the bins of the observed spectrum
  Int_t nbins = spObserved->GetNbinsX();

  cout << "filling collection areas ... " << endl;
  //fill the effective collection area into a tgraph
  TGraph *GraphEnVsArea[nbdeg];
#ifdef kDEBUG
  cout << "  effColArea->GetNbinsX() " << effColArea->GetNbinsX() << endl;
#endif
  for (int k=0; k<effColArea->GetNbinsY(); k++)
  {
    GraphEnVsArea[k] = new TGraph(1);
    GraphEnVsArea[k]->SetName(Form("GraphEnVsArea_%d",k));
    for (int i=0;i<effColArea->GetNbinsX();i++)
    {
      //in log energy (TeV)
      GraphEnVsArea[k]->SetPoint(i,effColArea->GetXaxis()->GetBinCenter(i+1),effColArea->GetBinContent(i+1,k+1));
#ifdef kDEBUG
      cout << " bin GraphEnVsArea  " << k << ": " << pow(10,effColArea->GetXaxis()->GetBinCenter(i+1)) << "   " << effColArea->GetBinContent(i+1, k+1) << endl;
#endif
    }
  }

  TGraph *GraphEnVsAreaTrue[nbdeg];
  cout << " Collection area in Etrue: " <<  endl;
  for (int k=0; k<effColArea->GetNbinsY(); k++)
  {
    GraphEnVsAreaTrue[k] = new TGraph(1);
    GraphEnVsAreaTrue[k]->SetName(Form("GraphEnVsAreaTrue_%d",k));
    for (int i=0;i<effColAreaTrue->GetNbinsX();i++)
    {
      //in log energy (TeV)
      GraphEnVsAreaTrue[k]->SetPoint(i,effColAreaTrue->GetXaxis()->GetBinCenter(i+1),effColAreaTrue->GetBinContent(i+1,k+1));
#ifdef kDEBUG
      cout << " bin GraphEnVsAreaTrue " << k << ": " << pow(10,effColAreaTrue->GetXaxis()->GetBinCenter(i+1)) << "   " << effColAreaTrue->GetBinContent(i+1, k+1) << endl;
#endif
    }
  }

  TGraph *GraphEnVsArea80[nbdeg];
  cout << " Collection area, 80per cent containment: " <<  endl;
  for (int k=0; k<effColArea->GetNbinsY(); k++)
  {
    GraphEnVsArea80[k] = new TGraph(1);
    GraphEnVsArea80[k]->SetName(Form("GraphEnVsArea80_%d",k));
    for (int i=0;i<effColArea80->GetNbinsX();i++)
    {
      //in log energy (TeV)
      GraphEnVsArea80[k]->SetPoint(i,effColArea80->GetXaxis()->GetBinCenter(i+1),effColArea80->GetBinContent(i+1,k+1));
#ifdef kDEBUG
      cout << " bin GraphEnVsArea80 " << k << ": " << pow(10,effColArea80->GetXaxis()->GetBinCenter(i+1)) << "   " << effColArea80->GetBinContent(i+1,k+1) << endl;
#endif
    }
  }

  cout << "filling background rates ... " << endl;
  //fill the background rates into a graph
  TGraph *GraphEnVsBg[nbdeg];
#ifdef kDEBUG
  cout << "  bgrate->GetNbinsX() " << bgrate->GetNbinsX() << endl;
#endif

  for (int k=0; k<effColArea->GetNbinsY(); k++)
  {
    GraphEnVsBg[k] = new TGraph(1);
    GraphEnVsBg[k]->SetName(Form("GraphEnVsBg_%d",k));
    for (int i=0;i<bgrate->GetNbinsX();i++)
    {
      double le = bgrate->GetXaxis()->GetBinLowEdge(i+1);
      double he = bgrate->GetXaxis()->GetBinLowEdge(i+2);
      double enr = pow(10,he) - pow(10,le); //bin width in TeV
      //X: in log10 energy (TeV)
      //Y: in Hz per TeV
      GraphEnVsBg[k]->SetPoint(i,bgrate->GetXaxis()->GetBinCenter(i+1),bgrate->GetBinContent(i+1,k+1)/enr);

      //    cout << pow(10,bgrate->GetBinCenter(i+1)) << "   TeV, " << bgrate->GetBinContent(i+1) << " Hz " << endl;
#ifdef kDEBUG
      cout << " bin GraphEnVsBg " << k << ": " << i+1 << ", " << enr << ", " << bgrate->GetBinContent(i+1,k+1)/enr << ", " << bgrate->GetXaxis()->GetBinCenter(i+1) << endl;
#endif
    }
  }

  cout << "filling background rates per sqdegree ... " << endl;
  //fill the background rates per deg2 into a graph
  TGraph *GraphEnVsBgDeg[nbdeg];
#ifdef kDEBUG
  cout << "  bgdeg->GetNbinsX() " << bgdeg->GetNbinsX() << endl;
#endif

  for (int k=0; k<effColArea->GetNbinsY(); k++)
  {
    GraphEnVsBgDeg[k] = new TGraph(1);
    GraphEnVsBgDeg[k]->SetName(Form("GraphEnVsBgDeg_%d",k));
    for (int i=0;i<bgdeg->GetNbinsX();i++)
    {
      double le = bgdeg->GetXaxis()->GetBinLowEdge(i+1);
      double he = bgdeg->GetXaxis()->GetBinLowEdge(i+2);
      double enr = pow(10,he) - pow(10,le);  //bin width in TeV
#ifdef kDEBUG
      cout << " le " << le  << endl;
      cout << " he " << he  << endl;
      cout << " bin " << k << ": " << i+1 << ", " << enr << endl;
#endif
      //X: in log10 energy (TeV)
      //Y: in Hz per TeV
      GraphEnVsBgDeg[k]->SetPoint(i,bgdeg->GetXaxis()->GetBinCenter(i+1),bgdeg->GetBinContent(i+1,k+1)/enr);

      //    cout << i+1 << "  bgperdeg : " << enr << ", " << bgdeg->GetBinContent(i+1)/enr << ", " << bgdeg->GetBinCenter(i+1) << endl;
    }
  }

  cout << "filling angular resolution ... " << endl;
  //fill the angular resolution into a graph 
  TGraph *GraphEnVsAngRes[nbdeg];
#ifdef kDEBUG
  cout << "  res->GetNbinsX() " << res->GetNbinsX() << endl;
#endif

  for (int k=0; k<effColArea->GetNbinsY(); k++)
  {
    GraphEnVsAngRes[k] = new TGraph(1);
    GraphEnVsAngRes[k]->SetName(Form("GraphEnVsAngRes_%d",k));
    for (int i=0;i<res->GetNbinsX();i++)
    {
      GraphEnVsAngRes[k]->SetPoint(i,res->GetXaxis()->GetBinCenter(i+1),res->GetBinContent(i+1,k+1));
      //    cout << i+1 << " angres : " << res->GetBinContent(i+1) << ", " << res->GetBinCenter(i+1) << endl;
    }
  }

  cout << "filling energy resolution ... " << endl;
  //fill the angular resolution into a graph 
  TGraph *GraphEnVsEnRes[nbdeg];
#ifdef kDEBUG
  cout << "  enres->GetNbinsX() " << enres->GetNbinsX() << endl;
#endif

  for (int k=0; k<effColArea->GetNbinsY(); k++)
  {
    GraphEnVsEnRes[k] = new TGraph(1);
    GraphEnVsEnRes[k]->SetName(Form("GraphEnVsEnRes_%d",k));
    for (int i=0;i<enres->GetNbinsX();i++)
    {
      GraphEnVsEnRes[k]->SetPoint(i,enres->GetXaxis()->GetBinCenter(i+1),enres->GetBinContent(i+1,k+1));
#ifdef kDEBUG
      cout << " bin " << k << ": " << i+1 << " enres : " << enres->GetBinContent(i+1,k+1) << ", " << enres->GetXaxis()->GetBinCenter(i+1) << endl;
#endif
    }
  }

  //needed for integral flux (integral over the sky map)
  TH1D *excessrate = new TH1D(*spObserved);
  excessrate->Reset();
  excessrate->SetName("excessrate");
  TH1D *bkgrate = new TH1D(*spObserved);
  bkgrate->Reset();
  bkgrate->SetName("bkgrate");

  //////////////////////////////////////////////////////////
  // find bins in fMorpho you want to integrate over 
  // (square around XSource, YSource with RSource "radius")
  //////////////////////////////////////////////////////////
//  int bincntr = fMorpho->FindBin(XSource,YSource);
  TArrayI binsmorp;
  binsmorp.Set(NBINMORPHO*NBINMORPHO);
  TArrayI binsrad;
  binsrad.Set(NBINMORPHO*NBINMORPHO);
  int bm = 0;
  for (int ii=1; ii<= NBINMORPHO; ii++) //X
  {
    double x = fMorpho->GetXaxis()->GetBinCenter(ii);
    for (int kk=1; kk<= NBINMORPHO; kk++)  //Y
    {
      double y = fMorpho->GetYaxis()->GetBinCenter(kk);
      if (( (XSource-x)*(XSource-x) + (YSource-y)*(YSource-y) ) <= (RSource*RSource + BINSIZEMORPHO*BINSIZEMORPHO/2) ) 
      {
        //this bin I want!
        int nbi = fMorpho->FindBin(x,y);
        double rm = sqrt(x*x+y*y);  //distance to camera center
        int rb = (int)(rm/degbin) + 1;  //bin in offset dependent histograms
        binsmorp.AddAt(nbi,bm);
        binsrad.AddAt(rb,bm);
        bm++;
      }
    }
  }
  binsmorp.Set(bm);
  binsrad.Set(bm);

#ifdef kDEBUG
  cout << " Number of usefull bins is :" << bm << endl;
#endif
 
  //////////////////////////////////////////////////////////
  ////////// FIRST LOOP OVER THE ENERGY BINS ///////////////
  //////////////////////////////////////////////////////////
  TH1D *excessSmeared = new TH1D(*spObserved); //with energy smearing
  excessSmeared->Reset();
  excessSmeared->SetDirectory(NULL);
  excessSmeared->SetName("excessSmeared");

  //rebin the histogram allowing for events with true energy outside of the specified bounds
  Double_t ledges[NBINS]; 
  Double_t ledgesUser[NBINS]; 
  int kk=0;
  ledges[0] = emin;
  double x1,x2,x3,x4;  //needed later
  int xb1,yb1,zb1;  //needed later

#ifdef kDEBUG
  cout << " energy edges " << endl;
  cout <<  " true bin 0, E=" << ledges[0] << "TeV"<< endl;
#endif
  for (kk=1;kk<=nbins;kk++)
  {
    double xLedge = spObserved->GetBinLowEdge(kk);  //OK
    ledges[kk] = xLedge;
    ledgesUser[kk-1] = xLedge;
#ifdef kDEBUG
    cout <<  " true  bin " << kk << ", E= " << ledges[kk] << "TeV, user bin" << kk-1 << ", E= " << ledgesUser[kk-1] << "TeV"<<endl;
#endif
  }
  ledges[kk] = spObserved->GetBinLowEdge(kk);
  ledgesUser[nbins] = useremax;
#ifdef kDEBUG
  cout <<  " true  bin " << kk << ", E= " << ledges[kk] << "TeV, user bin " << nbins << ", E= " << ledgesUser[nbins] << "TeV"<< endl;
#endif
  if (useremax < emax) 
  {
    kk++;
    ledges[kk] = emax;
  } 
#ifdef kDEBUG
  cout <<  " true bin " << kk << ", E= " << ledges[kk] << "TeV"<< endl;
#endif
  const int nbins2 = kk;

  TH1D *excessIdeal = new TH1D();
  excessIdeal->SetBins(nbins2,ledges);
  excessIdeal->Reset();
  excessIdeal->SetDirectory(NULL);
  excessIdeal->SetName("excessIdeal");

  // for the assumed TF1 spectrum dN/dE and observation time T,
  // calculate the weights in every Etrue bin as
  // Integral(dN/dE dE) * Aeff * T

  //first normalize the Etrue lines to 1
  cout << " normalizing migration matrix along Etrue " << endl;
  double sume = 0.;
  for (int k=1; k<=migmat->GetNbinsZ(); k++) //over offset bins
  {
    for (int j=1;j<=migmat->GetNbinsY();j++) //rows (Etrue)
    {
      sume=0.;
      for (int i=1;i<=migmat->GetNbinsX();i++) //columns (Eest)
        sume += migmat->GetBinContent(i,j,k);
      if (sume>0.)
      {
        for (int i=1;i<=migmat->GetNbinsX();i++) //columns (Eest)
        {
          if (migmat->GetBinContent(i,j,k)>0.) 
          {
#ifdef kDEBUG
            cout << Form("  i %d, j %d, k %d, sum = %.2e, v = %.2e",i,j,k,sume,migmat->GetBinContent(i,j,k)/sume);
#endif
            migmat->SetBinContent(i,j,k,migmat->GetBinContent(i,j,k)/sume);
          }
        }
      }
    }
    cout << "." ; fflush(stdout);
  }
  cout << "  DONE " << endl;

  //now, integrate the spectrum, multiply by  Aeff and  T
  //FOR EACH BIN IN X,Y PLANE!!
  // create migration matrix for each bin in MORPHO
  TH2D *migmorp[NBINMORPHO*NBINMORPHO];
  int mnx = migmat->GetNbinsX();
  int mny = migmat->GetNbinsY();
//  TArrayF *migmorpAr[NBINMORPHO*NBINMORPHO][mny]; //first index: space bin in Morpho, 
//  float *migmorpF[NBINMORPHO*NBINMORPHO][mny][mnx];   
                                                  //second index: y bin migmat  
  double mstartx = migmat->GetXaxis()->GetBinLowEdge(1);
  double mstopx = migmat->GetXaxis()->GetBinLowEdge(mnx+1);
  double mstarty = migmat->GetYaxis()->GetBinLowEdge(1);
  double mstopy = migmat->GetYaxis()->GetBinLowEdge(mny+1);
//  for (int i=0; i< NBINMORPHO*NBINMORPHO; i++)
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    migmorp[bb] = new TH2D (Form("migmorp%d",bb),Form("migmorp%d",bb),mnx,mstartx,mstopx,mny,mstarty,mstopy); 
//    for (int ik=0; ik < mny ; ik++)
//    {
//        //        migmorpAr[bb][ik] = new TArrayF(mnx); //mnx bins in X in migmat
//        for (int ix=0; ix < mnx ; ix++)
//        {
//            migmorpF[bb][ik][ix] = new float (0.);
//        }
//    }
#ifdef kDEBUG
    cout << " bin " << ii << " " << bb << endl;
#endif
  }

  double areaef =0.;
  double areaefRec80 =0.;
  double areaefRec =0.;
  double fintegral = 0.;
  double lowe=0.;
  double highe=0.;
  double meane=0;
  double atten3 = 1.;

  cout << " folding the migration matrices with Teff, spectrum and effective area " << endl; 
  // scan over all space bins MORPHO
//  for (int ii=0; ii< NBINMORPHO; ii++)  //in X
//  {
//    double xm = fMorpho->GetXaxis()->GetBinCenter(ii+1);
//    x1 = fMorpho->GetXaxis()->GetBinLowEdge(ii+1);
//    x2 = fMorpho->GetXaxis()->GetBinLowEdge(ii+2);
//    for (int kk=0; kk< NBINMORPHO; kk++)  // in Y
//    {
  //scam over fMorpho bins which are inside the ROI
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    //    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+2)  - (int) (binsmorp.GetAt(ii)/(NBINMORPHO+2));
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    int glbb = binsmorp.GetAt(ii);
    fMorpho->GetBinXYZ(glbb,xb1,yb1,zb1);
    double xm = fMorpho->GetXaxis()->GetBinCenter(xb1);
    x1 = fMorpho->GetXaxis()->GetBinLowEdge(xb1);
    x2 = fMorpho->GetXaxis()->GetBinLowEdge(xb1+1);
    double ym = fMorpho->GetYaxis()->GetBinCenter(yb1);
    x3 = fMorpho->GetYaxis()->GetBinLowEdge(yb1);
    x4 = fMorpho->GetYaxis()->GetBinLowEdge(yb1+1);
    double rm = sqrt(xm*xm+ym*ym);  //distance to camera center
    int rb = (int)(rm/degbin) + 1;  //bin in offset dependent histograms
#ifdef kDEBUG
    cout << Form("bb=%d, x1=%lf, x2=%lf, x3=%lf, x4=%lf, rm=%lf, rb=%d, ", bb,x1,x2,x3,x4,rm,rb); fflush(stdout);
#endif
    //find out which distance from camera center we have
    //      cout << endl << Form("j, energy, Flux, atten, time, effarea, expected excess") << endl;
    for (int j=1;j<=migmat->GetNbinsY();j++) //rows (Etrue)
    {
      lowe = migmat->GetYaxis()->GetBinLowEdge(j);
      highe = migmat->GetYaxis()->GetBinLowEdge(j+1);
      meane = (highe+lowe)/2.;
      lowe = pow(10,lowe);
      highe = pow(10,highe);

#ifdef kDEBUG
      cout << " lowe " << lowe << ", highe " << highe << ": "; 
      cout << Form("x1=%lf, x2=%lf, x3=%lf, x4=%lf ", x1,x2,x3,x4); 
#endif

      areaef = GraphEnVsAreaTrue[rb]->Eval(meane); //in true energy

     //correct for cut efficiency
      areaefRec80 = 0.;
      areaefRec80 = GraphEnVsArea80[rb]->Eval(meane); //in true energy
      areaefRec = 0.;
      areaefRec = GraphEnVsArea[rb]->Eval(meane); //in true energy
      if (areaefRec80>0. && areaefRec>0.)
      {
        areaef *= 1.25*areaefRec80/areaefRec;
#ifdef kDEBUG
        cout << Form(" areaf %.2e, corr = %.2e ", areaef, 1.25*areaefRec80/areaefRec); 
#endif
      }


      if (areaef<1e-4)
        areaef = 1e-4; //in m2

#ifdef kDEBUG
      cout << " areaef " << areaef; fflush(stdout); 
#endif
      //        cout << "  INTEGRAL " << spIntrMain->Integral(lowe,highe,x3,x4,x1,x2) << Form(", coord %lf-%lf;%lf-%lf",x3,x4,x1,x2) << ", atten=" << atten3 << " effOnTime=" << effOnTime << ", areaef=" << areaef;
      //      fintegral = spIntrMain->Integral(lowe,highe,-5,5,-5,5) * atten3 * effOnTime * areaef; //ideal number of gammas in this Etrue bin
      fintegral = 0.;
      if(lowe > spIntrMain->GetZmin() && highe < spIntrMain->GetZmax()) //inside the energy limits
        fintegral = spIntrMain->Integral(x1,x2,x3,x4,lowe,highe) * atten3 * effOnTime * areaef; //ideal number of gammas in this Etrue bin; 

      fintegral *=1e4; //area is in m2 (m2 -> cm2)
#ifdef kDEBUG
      cout << " fintegral: " << fintegral << endl;
#endif
      for (int i=1;i<=migmat->GetNbinsX();i++) //columns (Eest)
      {
        double val = migmat->GetBinContent(i,j,rb);
        val *= fintegral;
        migmorp[bb]->SetBinContent(i,j,val); // or should it be shuffled?
//        migmorpAr[bb][j-1]->AddAt(migmat->GetBinContent(i,j,rb)*fintegral,i-1); //mnx bins in X in mimat
      }
    }
    cout << "." ; fflush(stdout);
  }
  cout << "  DONE " << endl;


  // we now rebin the migration matrix in Etrue and in Eest
  // the user specified binning is:
#ifdef kDEBUG
  cout << endl << "User specified binning is: " << endl;
  cout << Form("\tLowEdge(TeV)\tUpperEdge(TeV)") << endl; 
  for (int i=1; i<=spObserved->GetNbinsX(); i++)
    cout << Form("%d\t%.2e\t%.2e",i, spObserved->GetBinLowEdge(i), spObserved->GetBinLowEdge(i+1)) << endl; 
#endif

  Double_t ledgesCoarse[NBINS]; 
  for  (int i=1; i<=nbins; i++)
  {
    double lelim = log10(spObserved->GetBinLowEdge(i));
    double uelim = log10(spObserved->GetBinLowEdge(i+1));
#ifdef kDEBUG
    cout << "i=" << i  << " " << lelim << " :" ; 
#endif
    for (int j=1; j<=migmat->GetNbinsX();j++)
    {
      //       cout << "  " << migmat->GetXaxis()->GetBinLowEdge(j);
      if (lelim >= migmat->GetXaxis()->GetBinLowEdge(j) &&
          lelim < migmat->GetXaxis()->GetBinLowEdge(j+1)) //found low edge
      {
        ledgesCoarse[i-1] = pow(10,migmat->GetXaxis()->GetBinLowEdge(j));
#ifdef kDEBUG
        cout << " new: " << ledgesCoarse[i-1] << " :" << endl;; 
#endif
        break;
      }
    }
    if (i==nbins) //last bin
    {
#ifdef kDEBUG
      cout << endl << "ulim  " << uelim << endl;;
#endif
      for (int jk=1; jk<=migmat->GetNbinsX();jk++)
      {
        if (uelim <= migmat->GetXaxis()->GetBinLowEdge(jk)) //upperedge found
        {
          ledgesCoarse[i] = pow(10,migmat->GetXaxis()->GetBinLowEdge(jk));
#ifdef kDEBUG
          cout <<  " last  " <<  ledgesCoarse[i] << endl;;
#endif
          break;
        }
      }
    }
  }
  
  Double_t ledgesCoarseTRUE[NBINS]; 
  Int_t nbinsCoarseTrue=nbins+2;
  ledgesCoarseTRUE[0] = pow(10,mstarty);
  for  (int i=1; i<=(nbins+1); i++)
  {
    ledgesCoarseTRUE[i] = ledgesCoarse[i-1];
  }
  ledgesCoarseTRUE[nbinsCoarseTrue]=pow(10,mstopy);
   

#ifdef kDEBUG
  cout << " Coarse True: 0 " << ledgesCoarseTRUE[0] << endl;
  for  (int i=0; i<=nbins; i++)
  {
    cout << "Coarse True  " << i+1 << " " << ledgesCoarseTRUE[i+1] << " " <<  ledgesCoarse[i] <<  " " << spObserved->GetBinLowEdge(i+1) << endl;
  }
  cout << " Coarse True: last " << ledgesCoarseTRUE[nbinsCoarseTrue] << endl;
#endif

  cout << " Filling in Coarse migration matrices " << endl;
  //define coarse response matrix PER X,Y BIN iN MORPHO
  TH2F *migmorpCoarse[NBINMORPHO*NBINMORPHO];
//  TArrayF *migmorpCoarseAr[NBINMORPHO*NBINMORPHO][nbinsCoarseTrue]; //first index: space bin in Morpho,
                                                                    //second index: nbinsCoarseTrue 
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    migmorpCoarse[bb] = new TH2F (Form("migmorpCoarse%d",bb),Form("migmorpCoarse%d",bb),nbins,ledgesCoarse, nbinsCoarseTrue, ledgesCoarseTRUE); 
  }
//  TH2F* migmatCoarse = new TH2F("migmatCoarse","",nbins,ledgesCoarse,nbinsCoarseTrue,ledgesCoarseTRUE);
  //fill the coarse response matrix
  //only bins which are inside ROI
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    int glbb = binsmorp.GetAt(ii);
    fMorpho->GetBinXYZ(glbb,xb1,yb1,zb1);
    //      int nbm = ii*NBINMORPHO+kk;// num bin in morpho (pay attention here!)
#ifdef kDEBUG
    cout << endl << " NBM (BB) and bins : " << Form("%d, %d, %d, %d, %d ", bb, glbb, xb1, yb1, zb1) ;
#endif
    //get distance from the center!
//    double xm = fMorpho->GetXaxis()->GetBinCenter(xb1);
//    double ym = fMorpho->GetYaxis()->GetBinCenter(yb1);
//    double rm = sqrt(xm*xm+ym*ym);  //distance to camera center
//    int rb = (int)(rm/degbin) + 1;  //bin in offset dependent histograms
    for (int j=1; j<=migmat->GetNbinsY();j++) //over etrue
    {
      double bcenty = pow(10,migmat->GetYaxis()->GetBinCenter(j));
#ifdef kDEBUG
      cout << endl << " Y " << bcenty ;
#endif
      for (int i=1; i<=migmat->GetNbinsX();i++) //over erec
      {
        double bcentx = pow(10,migmat->GetXaxis()->GetBinCenter(i));
        //          double value = migmat->GetBinContent(i,j,rb);
        double value = migmorp[bb]->GetBinContent(i,j);
//        double value = migmorpAr[bb][j-1]->GetAt(i-1);
#ifdef kDEBUG
        if (value>0.) 
          cout << ", X " << bcentx << ", V=" <<value ;
#endif
        migmorpCoarse[bb]->Fill(bcentx,bcenty,value);
      }
    }
    cout << "." ; fflush(stdout);
  }
  cout << "  DONE " << endl;

  //print the coarse matrix
//  for (int j=1; j<=migmatCoarse->GetNbinsY();j++) //over etrue
//  {
//    cout << endl << " line " << j ;
//    for (int i=1; i<=migmatCoarse->GetNbinsX();i++) //over erec
//    {
//      cout << "  " << migmatCoarse->GetBinContent(i,j);
//    }
//  }
//  cout << endl;

  // create TH1 with ideal number of gamma events in Etrue
  // FOR EACH BIN X,Y IN MORPHO
  TH1D *gammaIdeal[NBINMORPHO*NBINMORPHO];
  TH1D *gammaIdealRec[NBINMORPHO*NBINMORPHO];
  TH1D *spillOver[NBINMORPHO*NBINMORPHO];
//  for (int i=0; i< NBINMORPHO*NBINMORPHO; i++)
//  {
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    gammaIdeal[bb] = migmorpCoarse[bb]->ProjectionY(Form("gammaIdeal_%d",bb)); //ideal in Etrue
    gammaIdealRec[bb] = migmorpCoarse[bb]->ProjectionX(Form("gammaIdealec_%d",bb)); //ideal in Erec
    spillOver[bb]  = new TH1D(*(gammaIdealRec[bb]));
    spillOver[bb]->SetName(Form("spillOver_%d",bb));
    spillOver[bb]->Reset();
  }
  //get spillover factors
  cout << " ************************************************" << endl;
  cout << " SPILL OVER (NOT USED) " << endl;
  cout << " ************************************************" << endl;
//  for (int ii=0; ii< NBINMORPHO*NBINMORPHO; ii++)
//  {
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    for (int i=1; i<=spillOver[bb]->GetNbinsX(); i++)
    {
      double v1 = gammaIdeal[bb]->GetBinContent(i+1);
      double v2 = gammaIdealRec[bb]->GetBinContent(i);
      if (v2>0. && v1>0.)
      {
#ifdef kDEBUG
        cout << " SPILL OVER " << i << " " << v1/v2 << endl;
#endif
        spillOver[bb]->SetBinContent(i,v1/v2);  // true / reconstructed
      }
    }
  }


  TH1D *gammaIdealShuffled[NBINMORPHO*NBINMORPHO];
//  for (int ii=0; ii< NBINMORPHO*NBINMORPHO; ii++)
//  {
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    gammaIdealShuffled[bb] = new TH1D(*(gammaIdeal[bb]));
    gammaIdealShuffled[bb]->SetName(Form("gammaIdealShuffled_%d",bb));
    gammaIdealShuffled[bb]->Reset();
    for (int i=1; i<=nbinsCoarseTrue; i++)  //nbinsCoarseTrue=nbins+2;
    {
      double value =gammaIdeal[bb]->GetBinContent(i);
      if (kUseRandom) value = rnd.PoissonD(value); 
      gammaIdealShuffled[bb]->SetBinContent(i,value);
    }
  }

  //adjust bins of the output histogram
  spObserved->Reset();
  spObserved->SetBins(nbins,ledgesCoarse);

  //shuffle now the 2d migmat
  //shuffle the number of gamma events in Etrue according to Poission distribition
#ifdef kDEBUG
  cout << " Shuffeling magration matrix:" << endl;
#endif
  TH2F *migmatCoarseShuffled[NBINMORPHO*NBINMORPHO];
//  for (int ii=0; ii< NBINMORPHO*NBINMORPHO; ii++)
//  {
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    migmatCoarseShuffled[bb] = new TH2F(*(migmorpCoarse[bb]));
    migmatCoarseShuffled[bb]->SetName(Form("migmatCoarseShuffled_%d",bb));
    migmatCoarseShuffled[bb]->Reset();

    for (int j=1; j<=nbinsCoarseTrue; j++)  //nbinsCoarseTrue=nbins+2;
    {
      double value = gammaIdealShuffled[bb]->GetBinContent(j);
      double valueId = gammaIdeal[bb]->GetBinContent(j);
      for (int i=1; i<nbins; i++)
      {
        double v2 = migmorpCoarse[bb]->GetBinContent(i,j);
        if (valueId>0.)
        {
//          cout << " BINGO " << valueId << " " << v2 << " " << value << endl;
          v2 = v2/valueId*value;
        }
        migmatCoarseShuffled[bb]->SetBinContent(i,j,v2);
      }
    }

    //print the coarse shuffled matrix
#ifdef kDEBUG
    for (int j=1; j<=migmatCoarseShuffled[bb]->GetNbinsY();j++) //over etrue
    {
      cout << endl << " shuffled line in True Energy " << migmatCoarseShuffled[bb]->GetYaxis()->GetBinCenter(j) << " TeV:";
      for (int i=1; i<=migmatCoarseShuffled[bb]->GetNbinsX();i++) //over erec
      {
        cout << " " << migmatCoarseShuffled[bb]->GetBinContent(i,j);
      }
    }
    cout << endl;
#endif
  }

  //create distribution of expected excess events
  TH1D *gammaExpected[NBINMORPHO*NBINMORPHO];
//  for (int ii=0; ii< NBINMORPHO*NBINMORPHO; ii++)
//  {
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    gammaExpected[bb] = migmatCoarseShuffled[bb]->ProjectionX(Form("gammaExpected_%d",bb));
  }

  cout << " ************************************************" << endl;
  cout << " ************************************************" << endl;
  cout << " ************************************************" << endl;
  cout << "         LOOP OVER THE ENERGY BINS               " << endl;
  cout << " ************************************************" << endl;
  cout << " ************************************************" << endl;
  cout << " ************************************************" << endl;

  //second loop
  cout << endl << "entering the second loop over " << nbins << " energy bins of the spectrum ... " << endl;
  cout << "i, energy, sol_ang, idealExcess, excess, excess/idealExcess, BG, Flux, sigma " << endl;;
  int ik = 0;
  double fFlux=0.;
  double fFluxReb=0.;
  double fFluxRebEr =0.;
  double fExcess=0.;
  double fON=0;
  double fExcessReb=0.;
  double fFluxRebErP=0.;
  double fFluxRebErN=0.;
  bool toBeRebinned = kFALSE;
  double EnergyLowBound=0.;
  bool kPointSet = kFALSE;
  for (int i=0;i<nbins;i++)
  {
    double xLedge = spObserved->GetBinLowEdge(i+1);  //OK
    double xHedge = spObserved->GetBinLowEdge(i+2);  //OK
    double de = spObserved->GetBinWidth(i+1);
    double energy = exp(log(xLedge*xHedge)/2.);  // mean log 
    cout << i << "  " << energy << "  ";

    fExcess = 0.;
    //integrate only where you are told to:
    //i.e. only bins which are up to RSource degrees away 
    //from XSource, YSource
    for (int ii=0; ii< binsmorp.GetSize(); ii++)
    {
      int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
//      fExcess += TMath::Nint(gammaExpected[bb]->GetBinContent(i+1)); //convert to integer
      fExcess += gammaExpected[bb]->GetBinContent(i+1); //convert to integer
    }
    fExcess = TMath::Nint(fExcess); //to be sure

//    double psf = GraphEnVsAngRes->Eval(log10(energy)); 
//    cout << " " << psf;  //we do not care about PSF
//    double solid_angle = TMath::Pi()*(pow(1.6*psf,2)+pow(RSource,2)); //for Jim's toy
//    double solid_angle = TMath::Pi()*(pow(psf,2)+pow(RSource,2));  //because we use r80
    double solid_angle = BINSIZEMORPHO*BINSIZEMORPHO;  
    // we take bin size of fMorpho
    cout << "  " << solid_angle*binsrad.GetSize();
    
    //get BG for the energy bin, point-like
    double fBGpoint = 0.; 
    double fBGext = 0.;
    for (int ii=0; ii< binsrad.GetSize(); ii++)
    {
      int bb = binsrad.GetAt(ii);
      fBGpoint += effOnTime*IntegrateFromGraphLog(GraphEnVsBg[bb],log10(xLedge),log10(xHedge)); 
      fBGext += solid_angle*effOnTime*IntegrateFromGraphLog(GraphEnVsBgDeg[bb],log10(xLedge),log10(xHedge)); 
    }
   
    double fBG = kUseExtended ? fBGext : fBGpoint;
    //calculate expected BG 
    double fBGeff = fBG / alpha;
    //shuffle fBG in background region
    if (kUseRandom) fBGeff = rnd.PoissonD(fBGeff); //integer
    fBGeff = TMath::Nint(fBGeff); //to be sure
    //shuffle fBG in ON region
    if (kUseRandom) fBG = rnd2.PoissonD(fBG); //integer
    fBG = TMath::Nint(fBG); //to be sure

    fON = fExcess + fBG;   //integer
    fExcess = fON - fBGeff*alpha;  //this is the real excess

    excessrate->SetBinContent(i+1,fExcess);  //for the light curve
    bkgrate->SetBinContent(i+1,fBGeff); //for the light curve

    //calculate error on the excess:
    double fExcessError = sqrt(fBGeff*alpha*alpha+fON); //formula 5 in Li&Ma
    //calculate asymmetric errors, too
    double fOnErPos = (fON<GAUSL) ? poisUp[int(fON)] : sqrt(fON);//0.5 + sqrt(fON + 0.25); 
    fOnErPos *= fOnErPos;
    double fOnErNeg = (fON<GAUSL) ? poisLo[int(fON)] : sqrt(fON); //-0.5 + sqrt(fON + 0.25); 
    fOnErNeg *= fOnErNeg;
    double fBGerPos = (fBGeff<GAUSL) ? poisUp[int(fBGeff)] : sqrt(fBGeff); //0.5 + sqrt(fBGeff + 0.25);
    fBGerPos *= fBGerPos;
    double fBGerNeg = (fBGeff<GAUSL) ? poisLo[int(fBGeff)] : sqrt(fBGeff); //-0.5 + sqrt(fBGeff + 0.25);
    fBGerNeg *= fBGerNeg;
    double fExcessErrorP = sqrt(fBGerPos*alpha*alpha+fOnErPos); //formula 5 in Li&Ma
    double fExcessErrorN = sqrt(fBGerNeg*alpha*alpha+fOnErNeg); //formula 5 in Li&Ma

    //fill the number and error on excess events  in Erec
    excessgraph->SetPoint(i,energy,fExcess);
    excessgraph->SetPointEYlow(i,fExcessErrorN);
    excessgraph->SetPointEYhigh(i,fExcessErrorP);

    //get the relative difference between the idel expected excess in Etrue and measured fExcess and thus get measured flux
    fFlux = 0.;
    double excessideal = 0.; 
    double excessidealRec = 0;
    double sourceflux = 0.;
    for (int ii=0; ii< binsmorp.GetSize(); ii++)
    {
      int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
      excessideal += gammaIdeal[bb]->GetBinContent(i+2); // in Etrue
      excessidealRec += gammaIdealRec[bb]->GetBinContent(i+1); //in Erec 
      int glbb = binsmorp.GetAt(ii);
      fMorpho->GetBinXYZ(glbb,xb1,yb1,zb1);
      x1 = fMorpho->GetXaxis()->GetBinLowEdge(xb1);
      x2 = fMorpho->GetXaxis()->GetBinLowEdge(xb1+1);
      x3 = fMorpho->GetYaxis()->GetBinLowEdge(yb1);
      x4 = fMorpho->GetYaxis()->GetBinLowEdge(yb1+1);
      sourceflux += spIntrMain->Integral(x1,x2,x3,x4,xLedge,xHedge); 
    }
#ifdef kDEBUG
    cout << " source flux = " << sourceflux << " " ;
#endif
    if (excessideal>0.)
       fFlux = sourceflux / de * fExcess / excessidealRec;
//       fFlux = spIntr->Integral(xLedge,xHedge) / de * fExcess / excessideal * spillOver->GetBinContent(i);
//       fFlux = spIntr->Integral(xLedge,xHedge) / de * fExcess / excessideal;
//       fFlux = spIntr->Integral(xLedge,xHedge) / de;
//       fFlux = spIntr->Eval(energy);

    double AreaEff = 0.;
    for (int ii=0; ii< binsmorp.GetSize(); ii++)
    {
      int bb = binsrad.GetAt(ii);
      AreaEff += GraphEnVsArea[bb]->Eval(log10(energy));  // needed to judge if the excess can be trusted only (if Aeff too small then dont), do not care that summing up collection areas does not make sense
    }
    if (fExcess <= 0.)
    {
      cout << ", fExcess is only " << fExcess << " skipping the bin ... " << endl;
      if (AreaEff<MINAREA || !kRebin) //eff area is still too low or no rebinning wished
      {
        toBeRebinned = kFALSE;
        continue;
      }
      else //eff area is good, try to rebin
      {
        if (!toBeRebinned) EnergyLowBound = xLedge;
        toBeRebinned = kTRUE;
        fFluxReb += fFlux*de;
        fExcessReb += fExcess;
        fFluxRebEr += fFlux*fFlux*fExcessError*fExcessError*de*de; //assuming Excess of 1
        fFluxRebErP += fFlux*fFlux*fExcessErrorP*fExcessErrorP*de*de; //assuming Excess of 1
        fFluxRebErN += fFlux*fFlux*fExcessErrorN*fExcessErrorN*de*de; //assuming Excess of 1
        continue;
      }
    }

//    double fFluxEr = (fFlux/fExcess)*fExcessError;
//    double fFluxErP = (fFlux/fExcess)*fExcessErrorP;
//    double fFluxErN = (fFlux/fExcess)*fExcessErrorN;

    double relareaer = fFractionArea;
    relareaer *= fFlux*fFlux*relareaer;
    double fFluxEr = sqrt((fFlux*fFlux/fExcess/fExcess)*fExcessError*fExcessError + relareaer);
    double fFluxErP = sqrt((fFlux*fFlux/fExcess/fExcess)*fExcessErrorP*fExcessErrorP + relareaer);
    double fFluxErN = sqrt((fFlux*fFlux/fExcess/fExcess)*fExcessErrorN*fExcessErrorN + relareaer);
    if (fFluxEr==0.)
    {
      cout << ", fFluxEr is only " << fFluxEr << " (excessideal = " << excessideal << "), skipping the bin ... " << endl;
      if (kRebin && kPointSet) 
      {
        if (!toBeRebinned) EnergyLowBound = xLedge;
        toBeRebinned = kTRUE;
        fExcessReb += fExcess;
        fFluxReb += fFlux*(xHedge-xLedge);
      }
//      fFluxRebEr += fFluxEr*fFluxEr;
      continue;
    }

    cout << "  " << excessideal << " " << fExcess << "  " << fExcess/excessideal << " " << fBGeff*alpha << "  " <<  fFlux << " " << fFlux/fFluxEr << " ";

    //make bgeff positive in order to avoid problems with divisions by 0
    fBGeff = fBGeff > 0. ? fBGeff : SMALL; 

    //check:
    // if it is more than 3 sigma 
    // and number of excess events are greater than 7
    // and excess is bigger than 3% of the BG
    if (AreaEff<MINAREA) //the area is still too low, no rebinning
    {
      if (fExcess/(fBGeff*alpha) < MINBKG || 
          (fFlux/fFluxEr<MINSIG) ||
          fExcess < MINEVT) 
      {
        cout << " the point did not survive the cuts, and area too small ... " << endl;
        toBeRebinned = kFALSE;
        continue;
      }
    }
    else //area is big enough
    {
      if (fExcess < MINEVT ||             //TRY TO REBIN
            (fFlux/fFluxEr<MINSIG) || 
            fExcess/(fBGeff*alpha) < MINBKG )
      {
        if (kRebin && kPointSet)
        {
          if (!toBeRebinned) EnergyLowBound = xLedge;
          toBeRebinned = kTRUE;
          fFluxReb += fFlux*de;
          fFluxRebEr += fFluxEr*fFluxEr*de*de;
          fFluxRebErP += fFluxErP*fFluxErP*de*de;
          fFluxRebErN += fFluxErN*fFluxErN*de*de;
          //      fFlux=0.;
          cout << " the point did not survive the cuts ... ";
        }
        else
        {
          cout << " the point did not survive the cuts ... " << endl;
          continue;
        }
      }
    }

    //check if the rebinned point can survive the cuts
    if(toBeRebinned)
    {
      fFlux = fFluxReb / (xHedge-EnergyLowBound);
      fFluxEr = sqrt(fFluxRebEr) / (xHedge-EnergyLowBound);
      fFluxErP = sqrt(fFluxRebErP) / (xHedge-EnergyLowBound);
      fFluxErN = sqrt(fFluxRebErN) / (xHedge-EnergyLowBound);
      fExcessReb += fExcess;
      if (fFlux/fFluxEr<MINSIG || 
          fExcessReb < MINEVT || 
          fExcessReb/(fBGeff*alpha) < MINBKG)
      {
        cout << endl << "\t the rebinned point did not survive the cuts:  ";
        cout << Form ("%lf--%lf: %e, %e, %lf, %lf",EnergyLowBound,xHedge,fFlux,fFluxEr,fFlux/fFluxEr,fExcessReb);
      }
      else
      {
        cout << endl << "\t REBINNING SUCCESSFUL! ... : ";
        cout << Form ("in E = %lf--%lf: %e, %e, %lf, %lf",EnergyLowBound,xHedge,fFlux,fFluxEr,fFlux/fFluxEr,fExcessReb);
        energy = exp(log(EnergyLowBound*xHedge)/2.);  // mean log 
        xLedge=EnergyLowBound;
        toBeRebinned = kFALSE;
        fFluxReb = 0.; 
        fExcessReb = 0.;
        fFluxRebEr = 0.; 
        fFluxRebErP = 0.; 
        fFluxRebErN = 0.; 
      }
    }


    if (toBeRebinned)
    {
      cout << endl;
      continue;
    }

    specgraph->SetPoint(ik,energy,fFlux);
    if (fON < GAUSL || fBGeff < GAUSL) {
      specgraph->SetPointError(ik,energy-xLedge,xHedge-energy,fFluxErN,fFluxErP);
      cout << " USING ASYMMETRIC ERRORBARS ";
    }
    else
      specgraph->SetPointError(ik,energy-xLedge,xHedge-energy,fFluxEr,fFluxEr);

    ik++;
    cout << endl;

    if (ik>=minPointNotRebinned) 
      kPointSet = kTRUE;

  }

//  exit(1);

  //calculate the integral flux
  double ExcessTot = 0.;
  double ExcessIdeal = 0.;
  double BGtot = 0.;
  double FluxIdeal = 0.; 
  ExcessTot = excessrate->Integral();
  BGtot = bkgrate->Integral();
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    ExcessIdeal += gammaIdeal[bb]->Integral(2,nbins+1);
    int glbb = binsmorp.GetAt(ii);
    fMorpho->GetBinXYZ(glbb,xb1,yb1,zb1);
#ifdef kDEBUG
    cout << Form("glbb %d, bin %d, xlb = %d, ylb = %d", glbb, ii, xb1, yb1) << endl;
#endif
    x1 = fMorpho->GetXaxis()->GetBinLowEdge(xb1);
    x2 = fMorpho->GetXaxis()->GetBinLowEdge(xb1+1);
    x3 = fMorpho->GetYaxis()->GetBinLowEdge(yb1);
    x4 = fMorpho->GetYaxis()->GetBinLowEdge(yb1+1);
//    double sum = spIntrMain->Integral(useremin,useremax,x3,x4,x1,x2);
    double sum = spIntrMain->Integral(x1,x2,x3,x4,useremin,useremax);
#ifdef kDEBUG
    cout << Form("bin %d, x1 = %.2e, x2 = %0.2e, x3 = %.2e, x4 = %0.2e, sum = %.2e", ii, x1, x2, x3, x4, sum) << endl;
#endif
    FluxIdeal += sum;
#ifdef kDEBUG
    cout << Form("bin %d, ExcessTot %.2e, ExcessIdeal %.2e", ii, ExcessTot, ExcessIdeal) << endl;
#endif
  }
  double ONtot = ExcessTot + BGtot*alpha;

//  double ExcessTotEr = SignificanceLiMa(ONtot,BGtot,alpha);
  double bkgerrpos; 
  double bkgerrneg; 
  double onerrpos; 
  double onerrneg; 
  if (BGtot < GAUSL) //small numbers
  {
    bkgerrpos = poisUp[int(BGtot)]; // (0.5+sqrt(BGtot + 0.25));
    bkgerrneg = poisLo[int(BGtot)]; //(-0.5+sqrt(BGtot + 0.25));
    bkgerrpos *= bkgerrpos;
    bkgerrneg *= bkgerrneg;
  }
  else  //gaussian regime
  {
    bkgerrpos = BGtot;
    bkgerrneg = BGtot;
  }

  if (ONtot < GAUSL) //small numbers
  {
    onerrpos = poisUp[int(ONtot)]; //(0.5+sqrt(ONtot + 0.25));
    onerrneg = poisLo[int(ONtot)]; //(-0.5+sqrt(ONtot + 0.25));
    onerrpos *= onerrpos;
    onerrneg *= onerrneg;
  }
  else //gaussian regime
  {
    onerrpos = ONtot;
    onerrneg = ONtot;
  }

  double ExcessTotErPos = sqrt(alpha*alpha*bkgerrpos+onerrpos);
  double ExcessTotErNeg = sqrt(alpha*alpha*bkgerrneg+onerrneg);
  double FluxIreal = (ExcessIdeal!=0.) ? FluxIdeal*ExcessTot/ExcessIdeal : 0.;
  double FluxIrealEP = ((ONtot-BGtot*alpha)!=0.) ? FluxIreal*ExcessTotErPos/(ONtot-BGtot*alpha) : 0.;
  double FluxIrealEN = ((ONtot-BGtot*alpha)!=0.) ? FluxIreal*ExcessTotErNeg/(ONtot-BGtot*alpha) : 0.;
  ifluxgraph->SetPoint(0,useremin,FluxIreal);
  ifluxgraph->SetPointEYlow(0,FluxIrealEN);
  ifluxgraph->SetPointEYhigh(0,FluxIrealEP);

  cout << Form("Excess real = %.2e, Excess ideal = %.2e", ExcessTot, ExcessIdeal) << endl << endl;
  cout << endl << "    Integral flux "<< endl; 
  cout << Form("FluxReal = %.2e, FluxIdeal = %.2e ", FluxIreal, FluxIdeal) << endl << endl;

  //rebin the histogram
  double x,y, exl, eyl, eyh;
  Double_t lowedges[NBINS]; 
  int ikk;
  for (ikk=0;ikk<specgraph->GetN();ikk++)
  {
    specgraph->GetPoint(ikk,x,y);
    exl = specgraph->GetErrorXlow(ikk);
    lowedges[ikk] = x - exl;
  }
  double exh = specgraph->GetErrorXhigh(ikk-1);
  lowedges[ikk] = x + exh;
  spObserved->Reset();
  spObserved->SetBins(ikk,lowedges);
  for (ikk=0;ikk<specgraph->GetN();ikk++)
  {
    specgraph->GetPoint(ikk,x,y);
    eyl = specgraph->GetErrorYlow(ikk);
    eyh = specgraph->GetErrorYhigh(ikk);
    spObserved->SetBinContent(ikk+1,y);
    spObserved->SetBinError(ikk+1,(eyh+eyl)/2);
  }
//  cout << " highest energy reached: "  << endl; 
//  cout << lowedges[ikk] << endl; 


  // calculate energy threshold from gammaIdeal
  TH1D *gammaIdeal2 = migmat->ProjectionY("gammaIdeal2",0,-1,1,1); //ideal in Etrue
  gammaIdeal2->Reset();
  TH1D *gammaIdealTemp = new TH1D(*gammaIdeal2);
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    gammaIdealTemp = migmorp[bb]->ProjectionY("gammaIdealTemp"); //ideal in Etrue
    gammaIdeal2->Add(gammaIdealTemp);
    gammaIdealTemp->Reset();
//
//    for (int j=1;j<=migmat->GetNbinsY();j++) //rows (Etrue)
//    {
//      for (int i=1;i<=migmat->GetNbinsX();i++) //columns (Eest)
//      {
//          double v = migmorpAr[bb][j-1]->GetAt(i-1);
//          gammaIdeal2->AddBinContent(j,v);
//      }
//    }
  }
  (*threshold) = pow(10,gammaIdeal2->GetBinCenter(gammaIdeal2->GetMaximumBin()));
  cout << " THRESHOLD is "  << *threshold  << " TeV "<< endl;
  



  delete excessrate;
  delete bkgrate;
  for (int k=0; k<effColArea->GetNbinsY(); k++)
  {
    delete GraphEnVsArea[k];
    delete GraphEnVsAreaTrue[k];
    delete GraphEnVsArea80[k];
    delete GraphEnVsBg[k];
    delete GraphEnVsAngRes[k];
    delete GraphEnVsEnRes[k];
    delete GraphEnVsBgDeg[k];
  }
  delete excessIdeal;
  delete excessSmeared;

//  delete gammaIdealTemp;
  for (int ii=0; ii< binsmorp.GetSize(); ii++)
  {
    int bb = binsmorp.GetAt(ii) - (NBINMORPHO+1)  - ((int) (binsmorp.GetAt(ii)/(NBINMORPHO+2))) * 2;
    delete migmorp[bb];
    delete migmorpCoarse[bb];
    delete spillOver[bb];
    delete gammaIdealShuffled[bb];
    delete migmatCoarseShuffled[bb];
    delete gammaIdeal[bb];
    delete gammaIdealRec[bb];
    delete gammaExpected[bb];
//    for (int ik=0; ik < mny ; ik++) 
//        delete migmorpAr[bb][ik];
  }

  delete fMorpho;
  delete gammaIdeal2;

  f->Close();
  delete f;

  return 1;
}

