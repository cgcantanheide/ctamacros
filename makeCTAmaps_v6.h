#include <TRandom.h>
#include <iostream>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPad.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TF2.h>
#include <TString.h>
#include <TSpline.h>
#include <TArrayD.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TEllipse.h>


#define noDEBUG
#define noDEBUG_SHAPE
std::vector<TH2D*> CreateShape(double percentofCrab,double Gaussian_sigma, int selected_shape, int selected_spectrum, double EExpcut);
std::vector<TH2D*>  ConstructPhotonMap(std::vector<TH2D*> fluxMap, TString subarray,  bool IsEx);
std::vector<TH2D*> ReadFromFile(double percentofCrab,TString fshapename,TString fhistoname,int selected_spectrum);


// Under Testing (GP)
TString sigmar(double gamma_zero, double radius, double index_change);

Double_t pwlfunction(Double_t* x, Double_t* par);
Double_t gausfunctionPSF(Double_t *x, Double_t *par);
Double_t expfunction(Double_t* x, Double_t* par);
Double_t ringfunction(Double_t *x, Double_t *par);
Double_t fadding(Double_t *x, Double_t *par);
Double_t fadding2gauss(Double_t *x, Double_t *par);
double GetEffectiveArea( double logEMin, double logEMax, TString cnfname );
double GetBgRate( double logEMin, double logEMax, TString cnfname );
double GetPSF( double logEMin, double logEMax , TString cnfname);
Double_t gausfunctionPSF(Double_t *x, Double_t *par);
float GetBinArea(TH2D *histo,int xbin, int ybin);
TH2D* GausSmooth(TH2D *histo, float sigma);
void call_fit(std::vector<TH2D*> photonMap, TString cnfname, TString psf_outfile, double Gaussian_sigma );
void call_radialprofile(std::vector<TH2D*> photonMap );
void call_radialprofile2(std::vector<TH2D*> photonMap, bool xslice,double halfwidth,double range1,double range2,bool average);
int SetOffsetBins(float offset);
double LiMa(float Non, float Noff, double Alpha);
void help();



