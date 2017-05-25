#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPad.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TVector.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TString.h>
#include <TSpline.h>
#include <TArrayD.h>
#include <TRandom3.h>

/*
 * =====================================================================================
 *
 *       Filename:  makeCTAspec.C
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
 *        Version:  5.0
 *        Created:  10/06/10 
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
//#define PINDEX -2.62   //photon index to use in the migration matrix to wieght the entries
#define SMALL 1e-3   //used for minimum background
#define MINETAU  0.01 //min energy for which tau is available
#define MAXETAU  50.  //max energy for which tau is available

#define NUMPOINT 200

//defaults
TF1 *spIntrDefault = new TF1("spIntr","(6.0e-11)*pow((x/.3),(-2.31+(-0.26*TMath::Log10(x/.3))))",0.01,100.);  //Crab Nebula (MAGIC paper)
Double_t zSourceDefault = 0.2;
TString attenFileNameDefault = "exptau_z0.030_modelFranceschini.dat";
//Int_t specBinningOutDefault = 100;

double effOnTimeDefault = 20.*60*60; //20 hours //The number is in seconds
double sizeDegDefault  = 0.01; // default size of the gamma-ray emission in deg

Double_t xMin=0.01; //TeV
Double_t xMax=100.; //TeV

const int minPointNotRebinned = 3; //min number of spectral points which must be valid before trying to rebin

  //background function:
//differential background in 50h of observation (x in TeV)
TF1 *bg = new TF1("bg","300*pow(x,-3.2)",xMin,xMax); //in events
TF1 *bgSpec = new TF1("bgSpec","1.5e-13*pow(x,-3.2)",xMin,xMax); //in ph/cm2/s/TeV 

//flag between BG from spectrum and fixed BG in events
//Bool_t kUseBGFromSpec = kTRUE;

//flag between point-like source and extended source 
//Bool_t kUseExtended = kTRUE;

//flag between BG from spectrum and fixed BG in events
//Bool_t kUseRandom = kTRUE;

//flag if you want to apply spill-over correction (default kFALSE)
Bool_t kUseSpillOver = kFALSE; //obsolete

double efficgampsf = 1.2; //scale factor to account for efficiency of the theta2 cut in case of an extended source
double alphaDefault = 0.2;  // 0.2 = 5 times better bg estimation; in Li&Ma alpha = t_on / t_off
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


//main function
bool makeCTAspec(TH1D *spObserved, //expected spectrum to be observed with CTA
                 TGraphAsymmErrors *specgraph, //the same as above
                 TGraphAsymmErrors *ifluxgraph,  //store here the integral flux and its error
                 TGraphAsymmErrors *excessgraph, //distribution of excess events
		 TH1D *gammaExp,	// expected gamma counts 
		 TH1D *bkgExp,		// expected bkg counts
                 const char* filename,          //root format response file
                 TSpline3 *SplineEnVsAtt, //Spline with attenuation //**new
		 TVector *Stats,		// Vector to save some additional numbers
		 Double_t energyMin=1e16,	//min energy of Spline //**new
		 Double_t energyMax=1e-16,	//max energy of Spline //**new
                 Bool_t kRebin=kFALSE,         //should rebinning be applied if the number of excess events too low?
                 Bool_t kEnRes=kFALSE,         //should the energy resolution be applied to the data?
                 Bool_t kAtten=kFALSE,         //should EBL absorption be taken into account?
                 Bool_t kVerbose=kFALSE,         //print output? //**new
                 TF1 *spIntr=spIntrDefault,   //function of the source shape
                 Double_t effOnTime=effOnTimeDefault, // observation time
                 Double_t size_deg = sizeDegDefault, //radius of the emission in degrees
                 Bool_t kUseExtended=kTRUE,         //flag between point-like source and extended source
                 Bool_t kUseRandom=kTRUE,         //flag between BG from spectrum and fixed BG in events 
		 Double_t alpha = alphaDefault, //alpha value 
		 Double_t minEvt= MINEVT,
		 Double_t minBkg= MINBKG, 
		 Double_t minAeff= MINAREA, 
		 Double_t minSig= MINSIG,
		 Int_t seedOn=0,
		 Int_t seedOff=0
		 )
{

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

  if (kUseRandom)
      rnd.SetSeed(seedOn);
      rnd2.SetSeed(seedOff);
  double H72 = 72./72.;  //change second number if you want to change the Hubble constant!!

  if (kVerbose)
      cout << "entering the makeCTAspec macro ... " << endl;

//  const float photindex = PINDEX;

//  double fEfficgampsf = 1.2; //efficiency of gamma-rays for an extended source

//  const Bool_t kIntegrateTau=kTRUE;   //used to choose between tau integration methods
  const Bool_t kIntegrateTau2=kTRUE;  //used to choose between tau integration methods
//  const Bool_t kLogEnRes=kTRUE;       //if true, energy resolution is assumed to follow a Gaussian distribution in log10(Eest/Etrue)
                                      //if false, energy resolution is assumed to follow a Gaussian distribution in (Eest/Etrue)-1

  //open file
  TFile* f = NULL;
  if (!(f = TFile::Open(filename)))
  {
    cout << "ERROR: could not open root performance file " << filename << endl;
    return 0;
  }
  
  //get histograms
  TH1* effColArea = NULL;
  if (!(effColArea = (TH1*)f->Get("EffectiveArea")))
  {
    cout << "ERROR: did not find histogram EffectiveArea in the provided root performance file " << filename << endl;
    return 0;
  }

  TH1* effColAreaTrue = NULL;
  if (!(effColAreaTrue = (TH1*)f->Get("EffectiveAreaEtrue")))
  {
    cout << "ERROR: did not find histogram EffectiveAreaEtrue in the provided root performance file " << filename << endl;
    return 0;
  }

  TH1* effColArea80 = NULL;
  if (!(effColArea80 = (TH1*)f->Get("EffectiveArea80")))
  {
    cout << "ERROR: did not find histogram EffectiveArea80 in the provided root performance file " << filename << endl;
    return 0;
  }

  TH1* res = NULL;
  if (!(res = (TH1*)f->Get("AngRes80")))
  {
    cout << "ERROR: did not find histogram AngRes80 in the provided root performance file " << filename << endl;
    cout << "No Off-axis performance available? Trying to use only On-Axis..." << endl;
    if (!(res = (TH1*)f->Get("AngRes")))
    {
	cout << "ERROR: did not find histogram AngRes either in the provided root performance file " << filename << endl;
	return 0;
    }
    else
    {
	kUseExtended = kFALSE;
	cout << "... worked" << endl;
    }
  }
  TH1* bgdeg = NULL;
  if (!( bgdeg = (TH1*)f->Get("BGRatePerSqDeg")))
  {
    cout << "ERROR: did not find histogram BGRatePerSqDeg in the provided root performance file " << filename << endl;
    return 0;
  }
//  TH1* bgrate = (TH1*)f->Get("BGRatePerSqDeg");
  TH1* bgrate = NULL;
  if (!(bgrate = (TH1*)f->Get("BGRate")))
  {
    cout << "ERROR: did not find histogram BGRate in the provided root performance file " << filename << endl;
    return 0;
  }
  TH1* enres = NULL;
  if (!(enres = (TH1*)f->Get("ERes")))
  {
    cout << "ERROR: did not find histogram ERes in the provided root performance file " << filename << endl;
    return 0;
  }
  TH2F* migmat = NULL;
  if (!(migmat = (TH2F*)(f->Get("MigMatrix")->Clone("mymigmat"))))
  {
    cout << "ERROR: did not find migration matrix in the provided root performance file " << filename << endl;
    return 0;
  }


  //check min and max energies
  double emin = pow(10,bgrate->GetBinLowEdge(1));
  double emax = pow(10,bgrate->GetBinLowEdge(bgrate->GetNbinsX()+1));
  double useremin = spObserved->GetBinLowEdge(1);
  double useremax = spObserved->GetBinLowEdge(spObserved->GetNbinsX()+1);
  if (emin > useremin)
  {
    cout << "ERROR: min energy provided by user " << useremin << " TeV is below min energy in performance file: " << emin << " TeV" << endl;
    cout << "Tip: choose min energy above " << emin << " TeV" << endl;
    return 0;
  }
  else if (kEnRes && SAVEEMIN > useremin)
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

  //prepare for EBL attenuation
  // Not needed here since it is included in funtion call!
  //energy of the intrinsic spectrum _spIntr_ is in the observer frame!!!!


  //prepare loop over the bins of the observed spectrum
  Int_t nbins = spObserved->GetNbinsX();

  //////////////////////////////////////////////////////
  // fill the effective collection area into a tgraph //
  // fills GraphEnVsArea
  // fills GraphEnVsAreaTrue	- Aeff for true energy
  // fills GraphEnVsArea80	- Aeff with 80% containment
  //////////////////////////////////////////////////////
  if (kVerbose)
      cout << "filling collection areas ... " << endl;
  //fill the effective collection area into a tgraph
  TGraph *GraphEnVsArea = new TGraph(1);
  GraphEnVsArea->SetName("GraphEnVsArea");
  if (kVerbose)
      cout << "  effColArea->GetNbinsX() " << effColArea->GetNbinsX() << endl;
  for (int i=0;i<effColArea->GetNbinsX();i++)
  {
    //in log energy (TeV)
    GraphEnVsArea->SetPoint(i,effColArea->GetBinCenter(i+1),effColArea->GetBinContent(i+1));
    if (kVerbose)
	cout << pow(10,effColArea->GetBinCenter(i+1)) << "   " << effColArea->GetBinContent(i+1) << endl;
  }

  TGraph *GraphEnVsAreaTrue = new TGraph(1);
  GraphEnVsAreaTrue->SetName("GraphEnVsAreaTrue");
  if (kVerbose)
      cout << " Collection area in Etrue: " <<  endl;
  for (int i=0;i<effColAreaTrue->GetNbinsX();i++)
  {
    //in log energy (TeV)
    GraphEnVsAreaTrue->SetPoint(i,effColAreaTrue->GetBinCenter(i+1),effColAreaTrue->GetBinContent(i+1));
    if (kVerbose)
	cout << pow(10,effColAreaTrue->GetBinCenter(i+1)) << "   " << effColAreaTrue->GetBinContent(i+1) << endl;
  }

  TGraph *GraphEnVsArea80 = new TGraph(1);
  GraphEnVsArea80->SetName("GraphEnVsArea80");
  if (kVerbose)
      cout << " Collection area, 80per cent containment: " <<  endl;
  for (int i=0;i<effColArea80->GetNbinsX();i++)
  {
    //in log energy (TeV)
    GraphEnVsArea80->SetPoint(i,effColArea80->GetBinCenter(i+1),effColArea80->GetBinContent(i+1));
    if (kVerbose)
	cout << pow(10,effColArea80->GetBinCenter(i+1)) << "   " << effColArea80->GetBinContent(i+1) << endl;
  }

  ////////////////////////////////////////////
  // fill the background rates into a graph //
  // fills GraphEnVsBg		- background rate vs energy
  // fills GraphEnVsBgDeg	- background rate vs energy per sr
  ////////////////////////////////////////////
  if (kVerbose)
      cout << "filling background rates ... " << endl;
  //fill the background rates into a graph
  TGraph *GraphEnVsBg = new TGraph(1);
  GraphEnVsBg->SetName("GraphEnVsBg");
  if (kVerbose)
      cout << "  bgrate->GetNbinsX() " << bgrate->GetNbinsX() << endl;

  for (int i=0;i<bgrate->GetNbinsX();i++)
  {
    double le = bgrate->GetBinLowEdge(i+1);
    double he = bgrate->GetBinLowEdge(i+2);
    double enr = pow(10,he) - pow(10,le); //bin width in TeV
    //X: in log10 energy (TeV)
    //Y: in Hz per TeV
    GraphEnVsBg->SetPoint(i,bgrate->GetBinCenter(i+1),bgrate->GetBinContent(i+1)/enr);

//    cout << pow(10,bgrate->GetBinCenter(i+1)) << "   TeV, " << bgrate->GetBinContent(i+1) << " Hz " << endl;
    if (kVerbose)
	cout << i+1 << ", " << enr << ", " << bgrate->GetBinContent(i+1)/enr << ", " << bgrate->GetBinCenter(i+1) << endl;
  }

  if (kVerbose)
      cout << "filling background rates per sqdegree ... " << endl;
  //fill the background rates per deg2 into a graph
  TGraph *GraphEnVsBgDeg = new TGraph(1);
  GraphEnVsBgDeg->SetName("GraphEnVsBgDeg");
  if (kVerbose)
      cout << "  bgdeg->GetNbinsX() " << bgdeg->GetNbinsX() << endl;

  for (int i=0;i<bgdeg->GetNbinsX();i++)
  {
    double le = bgdeg->GetBinLowEdge(i+1);
    double he = bgdeg->GetBinLowEdge(i+2);
    double enr = pow(10,he) - pow(10,le);  //bin width in TeV
    //X: in log10 energy (TeV)
    //Y: in Hz per TeV
    GraphEnVsBgDeg->SetPoint(i,bgdeg->GetBinCenter(i+1),bgdeg->GetBinContent(i+1)/enr);

//    cout << i+1 << "  bgperdeg : " << enr << ", " << bgdeg->GetBinContent(i+1)/enr << ", " << bgdeg->GetBinCenter(i+1) << endl;
  }

  if (kVerbose)
      cout << "filling angular resolution ... " << endl;
  //fill the angular resolution into a graph 
  TGraph *GraphEnVsAngRes = new TGraph(1);
  GraphEnVsAngRes->SetName("GraphEnVsAngRes");
  if (kVerbose)
      cout << "  res->GetNbinsX() " << res->GetNbinsX() << endl;

  for (int i=0;i<res->GetNbinsX();i++)
  {
    GraphEnVsAngRes->SetPoint(i,res->GetBinCenter(i+1),res->GetBinContent(i+1));
//    cout << i+1 << " angres : " << res->GetBinContent(i+1) << ", " << res->GetBinCenter(i+1) << endl;
  }

  if (kVerbose)
      cout << "filling energy resolution ... " << endl;
  //fill the angular resolution into a graph 
  TGraph *GraphEnVsEnRes = new TGraph(1);
  GraphEnVsEnRes->SetName("GraphEnVsEnRes");
  if (kVerbose)
      cout << "  enres->GetNbinsX() " << enres->GetNbinsX() << endl;

  for (int i=0;i<enres->GetNbinsX();i++)
  {
    GraphEnVsEnRes->SetPoint(i,enres->GetBinCenter(i+1),enres->GetBinContent(i+1));
    if (kVerbose)
	cout << i+1 << " enres : " << enres->GetBinContent(i+1) << ", " << enres->GetBinCenter(i+1) << endl;
  }

  //needed for integral flux
  TH1D *excessrate = new TH1D(*spObserved);
  excessrate->Reset();
  excessrate->SetName("excessrate");
  TH1D *bkgrate = new TH1D(*spObserved);
  bkgrate->Reset();
  bkgrate->SetName("bkgrate");

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
  if (kVerbose)
  {
      cout << " energy edges " << endl;
      cout <<  " true bin 0, E=" << ledges[0] << "TeV"<< endl;
  }
  for (kk=1;kk<=nbins;kk++)
  {
    double xLedge = spObserved->GetBinLowEdge(kk);  //OK
    ledges[kk] = xLedge;
    ledgesUser[kk-1] = xLedge;
    if (kVerbose)
	cout <<  " true  bin " << kk << ", E= " << ledges[kk] << "TeV, user bin" << kk-1 << ", E= " << ledgesUser[kk-1] << "TeV"<<endl;
  }
  ledges[kk] = spObserved->GetBinLowEdge(kk);
  ledgesUser[nbins] = useremax;
  if (kVerbose)
      cout <<  " true  bin " << kk << ", E= " << ledges[kk] << "TeV, user bin " << nbins << ", E= " << ledgesUser[nbins] << "TeV"<< endl;
  if (useremax < emax) 
  {
    kk++;
    ledges[kk] = emax;
  } 
  if (kVerbose)
      cout <<  " true bin " << kk << ", E= " << ledges[kk] << "TeV"<< endl;
  const int nbins2 = kk;

  TH1D *excessIdeal = new TH1D();
  excessIdeal->SetBins(nbins2,ledges);
  excessIdeal->Reset();
  excessIdeal->SetDirectory(NULL);
  excessIdeal->SetName("excessIdeal");

  //first normalize the Etrue lines to 1
  double sume = 0.;
  for (int j=1;j<=migmat->GetNbinsY();j++) //rows (Etrue)
  {
    sume=0.;
    for (int i=1;i<=migmat->GetNbinsX();i++) //columns (Eest)
      sume += migmat->GetBinContent(i,j);
    if (sume>0.)
    {
      for (int i=1;i<=migmat->GetNbinsX();i++) //columns (Eest)
        migmat->SetBinContent(i,j,migmat->GetBinContent(i,j)/sume);
    }
  }

  //now, integrate the spectrum, multiply by  Aeff and  T
  double areaef =0.;
  double fintegral = 0.;
  double lowe=0.;
  double highe=0.;
  double meane=0;
  double atten3 = 1.;
//  for (int j=1;j<=migmat->GetNbinsY();j++) //rows (Etrue)
  if (kVerbose)
      cout << endl << Form("j, energy, Flux, atten, time, effarea, expected excess") << endl;
  for (int j=1;j<=migmat->GetYaxis()->GetNbins();j++) //rows (Etrue)
  {
    lowe = migmat->GetYaxis()->GetBinLowEdge(j);
    highe = migmat->GetYaxis()->GetBinLowEdge(j+1);
    meane = (highe+lowe)/2.;
    lowe = pow(10,lowe);
    highe = pow(10,highe);
    if (kAtten)
    {
      if (kIntegrateTau2) //USE THIS!
      {
        double xl = lowe < energyMin ? energyMin : lowe;
        double xh = highe > energyMax ? energyMax : highe;
        atten3 = IntegrateForTau(SplineEnVsAtt,spIntr,xl,xh);
      }
      if (H72!=1.)
        atten3 = pow(atten3,H72);
    }

    areaef = GraphEnVsAreaTrue->Eval(meane); //in true energy
    if (areaef<1e-3)
      areaef = 1e-3; //cm2
    fintegral = spIntr->Integral(lowe,highe) * atten3 * effOnTime * areaef; //ideal number of gammas in this Etrue bin
    fintegral *=1e4; //area is in m2 (m2 -> cm2)
//    cout << " Ntrue: " << fintegral << endl;
    if (kVerbose)
    {
	if (kAtten)
	  cout << Form("%d %.2e %.2e %.2e %.2e %.0lf %lf %lf",j,pow(10,meane),spIntr->Integral(lowe,highe), atten3, SplineEnVsAtt->Eval(pow(10,meane)), effOnTime, areaef, fintegral) << endl;
	else
	  cout << Form("%d %.2e %.2e %.2e %.0lf %lf %lf",j,pow(10,meane),spIntr->Integral(lowe,highe), atten3, effOnTime, areaef, fintegral) << endl;
    }
    for (int i=1;i<=migmat->GetNbinsX();i++) //columns (Eest)
      {
//        cout << "  " << i << " VM=" << migmat->GetBinContent(i,j) << " INT=" << fintegral;
        migmat->SetBinContent(i,j,migmat->GetBinContent(i,j)*fintegral); // or should it be shuffled?
      }
  }
//  cout << endl;
//  cout << endl;

  // we now rebin the migration matrix in Etrue and in Eest
  // the user specified binning is:
  if (kVerbose)
  {
      cout << endl << "User specified binning is: " << endl;
      cout << Form("\tLowEdge(TeV)\tUpperEdge(TeV)") << endl; 
  }
  for (int i=1; i<=spObserved->GetNbinsX(); i++)
      if (kVerbose)
	  cout << Form("%d\t%.2e\t%.2e",i, spObserved->GetBinLowEdge(i), spObserved->GetBinLowEdge(i+1)) << endl; 

  Double_t ledgesCoarse[NBINS]; 
  for  (int i=1; i<=nbins; i++)
  {
    double lelim = log10(spObserved->GetBinLowEdge(i));
    double uelim = log10(spObserved->GetBinLowEdge(i+1));
    if (kVerbose)
	cout << "i=" << i  << " " << lelim << " :" ; 
    for (int j=1; j<=migmat->GetNbinsX();j++)
    {
//       cout << "  " << migmat->GetXaxis()->GetBinLowEdge(j);
      if (lelim >= migmat->GetXaxis()->GetBinLowEdge(j) &&
          lelim < migmat->GetXaxis()->GetBinLowEdge(j+1)) //found low edge
      {
        ledgesCoarse[i-1] = pow(10,migmat->GetXaxis()->GetBinLowEdge(j));
        if (kVerbose)
	    cout << " new: " << ledgesCoarse[i-1] << " :" << endl;; 
        break;
      }
    }
    if (i==nbins) //last bin
    {
      if (kVerbose)
	   cout << endl << "ulim  " << uelim << endl;;
      for (int jk=1; jk<=migmat->GetNbinsX();jk++)
      {
        if (uelim <= migmat->GetXaxis()->GetBinLowEdge(jk)) //upperedge found
        {
          ledgesCoarse[i] = pow(10,migmat->GetXaxis()->GetBinLowEdge(jk));
	  if (kVerbose)
	      cout <<  " last  " <<  ledgesCoarse[i] << endl;;
          break;
        }
      }
    }
  }
  
  Double_t ledgesCoarseTRUE[NBINS]; 
  Int_t nbinsCoarseTrue=nbins+2;
  ledgesCoarseTRUE[0] = pow(10,migmat->GetYaxis()->GetBinLowEdge(1));
  for  (int i=1; i<=(nbins+1); i++)
  {
    ledgesCoarseTRUE[i] = ledgesCoarse[i-1];
  }
  ledgesCoarseTRUE[nbinsCoarseTrue]=pow(10,migmat->GetYaxis()->GetBinLowEdge(migmat->GetNbinsY()+1));
   

  if (kVerbose)
      cout << " Coarse True: 0 " << ledgesCoarseTRUE[0] << endl;
  for  (int i=0; i<=nbins; i++)
  {
    if (kVerbose)
	cout << "Coarse True  " << i+1 << " " << ledgesCoarseTRUE[i+1] << " " <<  ledgesCoarse[i] <<  " " << spObserved->GetBinLowEdge(i+1) << endl;
  }
  if (kVerbose)
      cout << " Coarse True: last " << ledgesCoarseTRUE[nbinsCoarseTrue] << endl;

  //define coarse response matrix
  TH2F* migmatCoarse = new TH2F("migmatCoarse","",nbins,ledgesCoarse,nbinsCoarseTrue,ledgesCoarseTRUE);
  //fill the coarse response matrix
  for (int j=1; j<=migmat->GetNbinsY();j++) //over etrue
  {
    double bcenty = pow(10,migmat->GetYaxis()->GetBinCenter(j));
//    cout << endl << " Y " << bcenty ;
    for (int i=1; i<=migmat->GetNbinsX();i++) //over erec
    {
      double bcentx = pow(10,migmat->GetXaxis()->GetBinCenter(i));
      double value = migmat->GetBinContent(i,j);
//      cout << ", X " << bcentx << ", V=" <<value ;
      migmatCoarse->Fill(bcentx,bcenty,value);
    }
  }
//  cout << endl;
//  cout << endl;

  //print the coarse matrix
  for (int j=1; j<=migmatCoarse->GetNbinsY();j++) //over etrue
  {
    if (kVerbose)
	cout << endl << " line " << j ;
    for (int i=1; i<=migmatCoarse->GetNbinsX();i++) //over erec
    {
      if (kVerbose)
	  cout << "  " << migmatCoarse->GetBinContent(i,j);
    }
  }
  if (kVerbose)
      cout << endl;

  // create TH1 with ideal number of gamma events in Etrue
  TH1D *gammaIdeal = migmatCoarse->ProjectionY("gammaIdeal"); //ideal in Etrue
  TH1D *gammaIdealRec = migmatCoarse->ProjectionX("gammaIdealec"); //ideal in Erec
  //get spillover factors
  TH1D *spillOver = new TH1D(*gammaIdealRec);
  spillOver->SetName("spillOver");
  spillOver->Reset();
  if (kVerbose)
  {
      cout << " ************************************************" << endl;
      cout << " ************************************************" << endl;
      cout << " ************************************************" << endl;
      cout << " ************************************************" << endl;
      cout << " SPILL OVER (NOT USED) " << endl;
  }
  for (int i=1; i<=spillOver->GetNbinsX(); i++)
  {
    double v1 = gammaIdeal->GetBinContent(i+1);
    double v2 = gammaIdealRec->GetBinContent(i);
    if (v2>0. && v1>0.)
    {
      if (kVerbose)
	  cout << " SPILL OVER " << i << " " << v1/v2 << endl;
      spillOver->SetBinContent(i,v1/v2);  // true / reconstructed
    }
  }


  TH1D *gammaIdealShuffled = new TH1D(*gammaIdeal);
  gammaIdealShuffled->SetName("gammaIdealShuffled");
  gammaIdealShuffled->Reset();
  for (int i=1; i<nbinsCoarseTrue; i++)
  {
    double value =gammaIdeal->GetBinContent(i);
    if (kUseRandom) value = rnd.PoissonD(value); 
    gammaIdealShuffled->SetBinContent(i,value);
  }

  //adjust bins of the output histogram
  spObserved->Reset();
  spObserved->SetBins(nbins,ledgesCoarse);

  //shuffle now the 2d migmat
  //shuffle the number of gamma events in Etrue according to Poission distribition
  TH2F* migmatCoarseShuffled = new TH2F(*migmatCoarse);
  migmatCoarseShuffled->SetName("migmatCoarseShuffled");
  migmatCoarseShuffled->Reset();

  for (int j=1; j<nbinsCoarseTrue; j++)
  {
    double value = gammaIdealShuffled->GetBinContent(j);
    double valueId = gammaIdeal->GetBinContent(j);
    for (int i=1; i<nbins; i++)
    {
      double v2 = migmatCoarse->GetBinContent(i,j);
      if (valueId>0.)
      {
        v2 = v2/valueId*value;
      }
      migmatCoarseShuffled->SetBinContent(i,j,v2);
    }
  }

  //print the coarse shuffled matrix
  if (kVerbose)
  {
      for (int j=1; j<=migmatCoarseShuffled->GetNbinsY();j++) //over etrue
      {
	cout << endl << " shuffled line " << j ;
	for (int i=1; i<=migmatCoarseShuffled->GetNbinsX();i++) //over erec
	{
	  cout << "  " << migmatCoarseShuffled->GetBinContent(i,j);
	}
      }
      cout << endl;
  }

  //create distribution of expected excess events
  TH1D *gammaExpected = migmatCoarseShuffled->ProjectionX("gammaExpected");

  //second loop
  if (kVerbose)
  {
      cout << endl << "entering the second loop over " << nbins << " energy bins of the spectrum ... " << endl;
      cout << "i, energy, psf, sol_ang, idealExcess, excess, excess/idealExcess, BG, Flux, sigma " << endl;;
  }
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
  bkgExp->Reset();
  bkgExp->SetBins(nbins,ledgesUser);
  gammaExp->Reset();
  gammaExp->SetBins(nbins,ledgesUser);
  for (int i=0;i<nbins;i++)
  {
    gammaExp->SetBinContent(i+1,gammaExpected->GetBinContent(i+1));
  }

  for (int i=0;i<nbins;i++)
  {
    double xLedge = spObserved->GetBinLowEdge(i+1);  //OK
    double xHedge = spObserved->GetBinLowEdge(i+2);  //OK
    double de = spObserved->GetBinWidth(i+1);
    double energy = exp(log(xLedge*xHedge)/2.);  // mean log 
    if (kVerbose)
	cout << i << "  " << energy << "  ";

    fExcess = TMath::Nint(gammaExpected->GetBinContent(i+1)); //convert to integer

    double psf = GraphEnVsAngRes->Eval(log10(energy)); 
    if (kVerbose)
	cout << " " << psf;
//    double solid_angle = TMath::Pi()*(pow(1.6*psf,2)+pow(size_deg,2)); //for Jim's toy
    //Refine the formula??
    double solid_angle = TMath::Pi()*(pow(psf,2)+pow(size_deg,2));  //because we use r80
    if (kVerbose)
	cout << "  " << solid_angle;
    
    //get BG for the energy bin, point-like
    double fBGpoint = effOnTime*IntegrateFromGraphLog(GraphEnVsBg,log10(xLedge),log10(xHedge)); 
    double fBGext;
   
    fBGext = solid_angle*effOnTime*IntegrateFromGraphLog(GraphEnVsBgDeg,log10(xLedge),log10(xHedge)); 

    double fBG = kUseExtended ? fBGext : fBGpoint;
    //calculate expected BG 
    double fBGeff = fBG / alpha;
    //shuffle fBG in background region
    bkgExp->SetBinContent(i+1,fBG);
    if (kUseRandom)
    {
	fBGeff = rnd.PoissonD(fBGeff); //integer
	fBGeff = TMath::Nint(fBGeff); //to be sure // only needed if random used
    }
    //shuffle fBG in ON region
    if (kUseRandom) 
    {
	fBG = rnd2.PoissonD(fBG); //integer
	fBG = TMath::Nint(fBG); //to be sure // only needed if random used
    }

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
    double excessideal = gammaIdeal->GetBinContent(i+2); //in Etrue
    double excessidealRec = gammaIdealRec->GetBinContent(i+1); //in Erec
    if (excessideal>0.)
//       fFlux = spIntr->Integral(xLedge,xHedge) / de * fExcess / excessideal * spillOver->GetBinContent(i);
       fFlux = spIntr->Integral(xLedge,xHedge) / de * fExcess / excessidealRec;
//       fFlux = spIntr->Integral(xLedge,xHedge) / de * fExcess / excessideal;
//       fFlux = spIntr->Integral(xLedge,xHedge) / de;
//       fFlux = spIntr->Eval(energy);

    double atten3 = 1.;
    if (kAtten && kIntegrateTau2)
    {
      atten3 = IntegrateForTau(SplineEnVsAtt,spIntr,xLedge,xHedge);
      if (H72!=1.)
        atten3 = pow(atten3,H72);

      fFlux *= atten3;
//      cout << endl << " E = " << xLedge << " - " << xHedge << ", " << atten3 << " flux " <<  atten3* spIntr->Integral(xLedge,xHedge) / de * energy * energy << " " << atten3*spIntr->Eval(energy)*energy*energy << endl;
    }

    double AreaEff = GraphEnVsArea->Eval(log10(energy));  // needed to judge if the excess can be trusted only (if Aeff too small then dont)
    if (fExcess <= 0.)
    {
      if (kVerbose)
	  cout << ", fExcess is only " << fExcess << " skipping the bin ... " << endl;
      if (AreaEff<minAeff|| !kRebin) //eff area is still too low or no rebinning wished
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
      if (kVerbose)
	  cout << ", fFluxEr is only " << fFluxEr << " skipping the bin ... " << endl;
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

    if (kVerbose)
	cout << "  " << excessideal << " " << fExcess << "  " << fExcess/excessideal << " " << fBGeff*alpha << "  " <<  fFlux << " " << fFlux/fFluxEr << " ";

    //make bgeff positive in order to avoid problems with divisions by 0
    fBGeff = fBGeff > 0. ? fBGeff : SMALL; 

    //check:
    // if it is more than 3 sigma 
    // and number of excess events are greater than 7
    // and excess is bigger than 3% of the BG
    if (AreaEff<minAeff) //the area is still too low, no rebinning
    {
      if (fExcess/(fBGeff*alpha) < minBkg || 
          (fFlux/fFluxEr<minSig) ||
          fExcess < minEvt) 
      {
        if (kVerbose)
	    cout << " the point did not survive the cuts, and area too small ... " << endl;
        toBeRebinned = kFALSE;
        continue;
      }
    }
    else //area is big enough
    {
      if (fExcess < minEvt ||             //TRY TO REBIN
            (fFlux/fFluxEr<minSig) || 
            fExcess/(fBGeff*alpha) < minBkg)
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
	  if (kVerbose)
	      cout << " the point did not survive the cuts ... ";
        }
        else
        {
	  if (kVerbose)
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
      if (fFlux/fFluxEr<minSig|| 
          fExcessReb < minEvt|| 
          fExcessReb/(fBGeff*alpha) < minBkg)
      {
        if (kVerbose)
	{
	    cout << endl << "\t the rebinned point did not survive the cuts:  ";
	    cout << Form ("%lf--%lf: %e, %e, %lf, %lf",EnergyLowBound,xHedge,fFlux,fFluxEr,fFlux/fFluxEr,fExcessReb);
	}
      }
      else
      {
        if (kVerbose)
	{
	    cout << endl << "\t REBINNING SUCCESSFUL! ... : ";
	    cout << Form ("in E = %lf--%lf: %e, %e, %lf, %lf",EnergyLowBound,xHedge,fFlux,fFluxEr,fFlux/fFluxEr,fExcessReb);
	}
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
      if (kVerbose)
	  cout << endl;
      continue;
    }

    specgraph->SetPoint(ik,energy,fFlux);
    if (fON < GAUSL || fBGeff < GAUSL) {
      specgraph->SetPointError(ik,energy-xLedge,xHedge-energy,fFluxErN,fFluxErP);
      if (kVerbose)
	  cout << " USING ASYMMETRIC ERRORBARS ";
    }
    else
      specgraph->SetPointError(ik,energy-xLedge,xHedge-energy,fFluxEr,fFluxEr);

    ik++;
    if (kVerbose)
	cout << endl;

    if (ik>=minPointNotRebinned) 
      kPointSet = kTRUE;

  }

//  exit(1);

  //calculate the integral flux
  double ExcessTot = excessrate->Integral();
  double ExcessIdeal = gammaIdeal->Integral(2,nbins+1);
  double BGtot = bkgrate->Integral();
  double FluxIdeal = spIntr->Integral(useremin,useremax);
  double ONtot = ExcessTot + BGtot*alpha;

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

  if (kVerbose)
  {
      cout << endl << "    Integral flux "<< endl; 
      cout << Form("FluxReal = %.2e, FluxIdeal = %.2e ", FluxIreal, FluxIdeal) << endl << endl;
  }

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
  ((*Stats))[0] = SignificanceLiMa(ONtot,BGtot,alpha);

  if (kVerbose)
      cout << Form("Total Li&Ma significance = %.3e", ((*Stats))[0]);

  TH1D *gammaIdeal2 = migmat->ProjectionY("gammaIdeal2"); //ideal in Etrue
  ((*Stats))[1] = pow(10,gammaIdeal2->GetBinCenter(gammaIdeal2->GetMaximumBin()));
  //Stats->Print();

  if (kVerbose)
      cout << " THRESHOLD is "  << ((*Stats))[1]  << " TeV "<< endl;
  



//  delete mig;
//  delete weightfunc;
//  delete gaum;
  delete excessrate;
//  cout << "free 1 " << endl;
  delete bkgrate;
//  cout << "free 2 " << endl;
//  cout << "free 3 " << endl;
  f->Close();
//  cout << "free 4 " << endl;
  delete f;
//  cout << "free 5 " << endl;
  if (SplineEnVsAtt) delete SplineEnVsAtt;
//  cout << "free 6 " << endl;
  delete GraphEnVsArea;
//  cout << "free 7 " << endl;
  delete GraphEnVsBg;
//  cout << "free 8 " << endl;
  delete GraphEnVsAngRes;
//  cout << "free 9 " << endl;
  delete GraphEnVsEnRes;
//  cout << "free 9a " << endl;
  delete GraphEnVsBgDeg;
//  cout << "free 10" << endl;
  delete excessIdeal;
//  cout << "free 11" << endl;
  delete excessSmeared;
//  cout << "free 12" << endl;

//  delete gammaIdeal;
//  delete gammaIdealRec;
//  delete gammaIdeal2;

  return 1;
}

