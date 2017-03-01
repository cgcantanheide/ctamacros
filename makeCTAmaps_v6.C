/*
 * =====================================================================================
 *
 *       Filename:  makeCTAmaps.C
 *
 *   Description:  
 *   Simulation of skymaps (excess & flux maps) for CTA
 *   Version:  1.0
 *   Created:  24/01/11 
 *
 *   Author:  Emma de Ona Wilhelmi (eow) <emma@mpi-hd.mpg.de> // Giovanna Pedaletti (gp) <gpedalet@aliga.ieec.uab.es>
 *   Institute: MPIk, Heidelberg                  
 *   NOTES: Valid for sources < 0.5 degree
 * =====================================================================================
 */
#include <makeCTAmaps_v6.h>
#define DISTMAX 4.5



/*************************************************************/
// BEGINNING USER SIMULATION OPTIONS
/*************************************************************/
/********************/
/* Creating a shape */
/********************/
bool create_excess=true;                       // This option allows you to create a morphology shape with a given spectrum (see below for options)
double user_Xpos=0;                            // For the ex. in the 2D histogram   -258.7;
double user_Ypos=0;                            // For the ex. in the 2D histogram -40.
int whichshape=0;                              // Shape to create 0: Gauss, 1: Shell-type, 2: Composite, 3: Cooling
double ExternalGaus=0.2;                       // For whichshape=1/2 define the inner and outer radius of the shell (substracting two gaussians functions)
double InternalGaus=0.18;
int whichspectrum=0;                           // 0: PL, 1: PL+Exp
double Gamma=2.7;                              // photon spectrum for a galactic source
double Ecut=1.;                                // Exponential Cutoff if selected
double FluxRatio=2;                            // Flux Ratio between the outer and inner emission region (for cases 2 and 3) Fshell/Fpoint_like or Fextended/Fpoint_like
/***************************************/
/* Reading a shape from a 2D histogram */
/***************************************/
TString fshapename = "./j1713_sim.root";       // Name of your root file
TString fhistoname = "Excess";                 // Name of your 2D histo (TH2D)

/**************************************/
/* To fit your photon map fitsrc=true */
/**************************************/
bool fitsrc=false;                                   

/***********************************************/
/* To obtain radial profile radialprofile=true */
/***********************************************/
bool radialprofile=false; 
bool user_xslice=true;                           // Select X or Y axis for your profile
double user_halfwidth=10.;                       // Size of your profile
double user_range1=0;                            // if range1 and 2 are set to 0, the whole map is used
double user_range2=0;
bool user_average=false;                         // Average bin content

/********************/
/* Output           */
/********************/
bool saveplots=true;                             // Save Plots
bool savefit=false;                              // Save Fit Results in a text file 
bool plotresults=true;                           // Plot results
bool Calc_Significance=true;

/********************/
/* General          */
/********************/
bool printonlyex=0;                                   // Include background fluctuations [0] or not [1]
bool Is_IFAEoffsetEA=true;                            // Include IFAE-Extended effective areas
TString configpath = "./config/";
bool differential=true;                               // Differential or Integral Energy Bins
const int Ebins=1;                                    // number of energy bins
double Eth=0.05;                                      // Energy threshold
double Emax=100; 
double observationTime=50*60*60;                       // observation time in seconds
double rangex=3;                                      // 3 degree in X axis
double rangey=3;                                      // 3 degree in Y axis
bool gaus_smooth=0;
float smooth=0.06;
/*************************************************************/
// END  -- USER SIMULATION OPTIONS
/*************************************************************/

double crabflux=2.86e-11; // TeV cm-2 s-1
const int NBINS_OFFSET=17;
float Src_Offset=0;                                   // Offset of the source with respect of the center of the image
float OffSet_Bins[NBINS_OFFSET]={0.,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0};


double EbinMin[Ebins];
double logEbinMin[Ebins];
double EbinMax[Ebins];
double logEbinMax[Ebins];
double Emin=Eth;
double logEmin=TMath::Log10(Emin);
double logEmax=TMath::Log10(Emax);
double Ebinsize=(logEmax-logEmin)/(Ebins);
float OffsetMin;
float OffsetMax;
int nbinsx=rangex/0.03;                               // Bin Size fixed to 0.01
int nbinsy=rangey/0.03;

double psfsigmaFit[Ebins];
TF2 * fGaussianFit[Ebins];
TH2D* fExcess[Ebins];  
TH1D* hradialprofile[Ebins];
TH2D* hSigma[Ebins];  // Significance maps
TGraphErrors *grradialprofile[Ebins]; // not used




// **************************************************************************************************************
//main function
// **************************************************************************************************************
void makeCTAmaps(TString subarray="I",
		 double FluxCrabUnits=1,                // per cent of flux compared to Crab at 1 TeV
		 double Gaussian_sigma=0.1,             // intrinsic extension of the source
		 TString basename="MySimulation"        // basename for output files
		 )
{

  /****************/
  /* Output Files */
  /****************/


  TString psf_outfile(basename);
  psf_outfile+="_FitResults.dat";
  TString plots_outfile(basename);
  plots_outfile+="_PlotsResults.root";

  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetPalette(60);

  /****************/
  /* Binning      */
  /****************/
  
  for(int ibin=0;ibin<Ebins;ibin++){
    logEbinMin[ibin]=logEmin+ibin*Ebinsize;
    logEbinMax[ibin]=logEmin+(ibin+1)*Ebinsize;

    EbinMin[ibin]=pow(10,logEbinMin[ibin]);
    EbinMax[ibin]=pow(10,logEbinMax[ibin]);
    
#ifdef DEBUG
    std::cout << " Ebin: ["<< EbinMin[ibin] << " , " << EbinMax[ibin] << " ] "  << std::endl; 
    std::cout << " logEbin: ["<< logEbinMin[ibin] << " , " << logEbinMax[ibin] << " ] " << " EbinSize " << Ebinsize <<  std::endl; 
#endif   
  }

  if(Is_IFAEoffsetEA ==1 && subarray.CompareTo("E") > 0){
    if(subarray.CompareTo("B") > 0)  {
      std::cout << "Sorry! Only configurations E and B are available for offset analysis and you select subarray " << subarray <<  std::endl;
      exit(0);
    }
  }



  /********************************************************************/
  /* Read the original skymap (TH2F or create it from function shape) */
  /********************************************************************/
  std::vector<TH2D*> fluxMapvector;
  if(create_excess)
    fluxMapvector = CreateShape(FluxCrabUnits,Gaussian_sigma,whichshape,whichspectrum,Ecut);
  else{
    std::cout<< "create a new skymap from a function with extension sigma" << std::endl;
    fluxMapvector = ReadFromFile(FluxCrabUnits,fshapename,fhistoname,whichspectrum);
  }


  if(Is_IFAEoffsetEA){
    std::cout << " Your source is offset w.r.t. centre of the camera a distance of " << Src_Offset << std::endl;
    OffsetMin=SetOffsetBins(Src_Offset); // Setting the bin number for the offset 
    OffsetMax=OffsetMin + 1;
#ifdef DEBUG_SHAPE
     std::cout << " The minimum offset :" <<  OffsetMin << " and the maximum :" << OffsetMax << std::endl;
#endif
  }


  /********************************************************************/
  /* Construct Photon Maps                                            */ 
  /********************************************************************/
  std::vector<TH2D*> photonMapVector;

  if(printonlyex)  {
    photonMapVector =   ConstructPhotonMap(fluxMapvector,subarray,printonlyex);
  }
  else{
    photonMapVector = ConstructPhotonMap(fluxMapvector,subarray,printonlyex);
  }

  /********************************************************************/
  /* (Simple) Analysis of the skymaps                                 */ 
  /********************************************************************/
  
  if(fitsrc)
    call_fit(photonMapVector, subarray, psf_outfile, Gaussian_sigma);
  
      
  if(radialprofile)
    call_radialprofile2(photonMapVector,user_xslice, user_halfwidth,user_range1,user_range2,user_average);
    //    call_radialprofile(photonMapVector);

  // Gaussian Smoothed
  TH2D * tmpgaus[Ebins];  
  for(int ibin=0;ibin<Ebins;ibin++){
    if(gaus_smooth){ 
      tmpgaus[ibin]=GausSmooth(photonMapVector[ibin],smooth);     
      TString hSmoothName("Smooth_Bin");
      hSmoothName+=ibin;
      tmpgaus[ibin]->SetName(hSmoothName.Data());
      tmpgaus[ibin]->SetTitle(hSmoothName.Data());
    }
  }
  


  /********************************************************************/
  /* Draw Results                                                     */ 
  /********************************************************************/

  TFile *fsaveplots=0;
  if(saveplots){
    fsaveplots = new TFile(plots_outfile.Data(), "RECREATE" );    
  }
  if(plotresults){
    




  for( size_t i=0; i!=photonMapVector.size(); ++i ) { // excess+bg
    if(Calc_Significance){
      TString cnames("Canvas_significance_");
      cnames+=i;
      TCanvas* cs = new TCanvas(cnames.Data(),"c",600,400);  
      cs->SetFillColor(10);
      cs->SetFillColor(10);
      cs->SetBorderMode(0);
      cs->SetBorderSize(0);
      
      cs->SetFrameFillColor(0);
      cs->SetFrameLineWidth(2);
      cs->SetFrameBorderMode(0);
      cs->SetFrameLineWidth(3);
      cs->SetFrameBorderMode(0);
      gPad->SetLeftMargin(0.12857);
      gPad->SetRightMargin(0.2);
      gPad->SetTopMargin(0.118217);
      gPad->SetBottomMargin(0.128915);
      //    hSigma[i]->SetTitle(titlestr1.Data());
      hSigma[i]->GetXaxis()->SetAxisColor(10);
      hSigma[i]->GetYaxis()->SetAxisColor(10);
      hSigma[i]->GetXaxis()->SetTitle("X");
      hSigma[i]->GetYaxis()->SetTitle("Y");
      hSigma[i]->GetZaxis()->SetTitle("Significance");
      hSigma[i]->GetZaxis()->SetTitleOffset(1.4);
      hSigma[i]->Draw("colz");
    }
    TString cname("Canvas_");  
    cname+=i;
    TString titlestr1("Flux Map : Energy ");
    TString titlestr2("Events Map : Energy ");
    TCanvas* c = new TCanvas(cname.Data(),"c",600,800);  
    c->SetFillColor(10);
    c->SetFillColor(10);
    c->SetBorderMode(0);
    c->SetBorderSize(0);
 
    c->SetFrameFillColor(0);
    c->SetFrameLineWidth(2);
    c->SetFrameBorderMode(0);
    c->SetFrameLineWidth(3);
    c->SetFrameBorderMode(0);

    c->Divide(1, 2);
    c->cd(1);
    gPad->SetLeftMargin(0.12857);
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.118217);
    gPad->SetBottomMargin(0.128915);

    
    titlestr1+=(float)EbinMin[i];
    titlestr1+=" TeV ";
    titlestr1+=(float)EbinMax[i];
    titlestr1+=" TeV ";
    fluxMapvector[i]->SetTitle(titlestr1.Data());
    fluxMapvector[i]->GetXaxis()->SetAxisColor(10);
    fluxMapvector[i]->GetYaxis()->SetAxisColor(10);
    fluxMapvector[i]->GetXaxis()->SetTitle("X");
    fluxMapvector[i]->GetYaxis()->SetTitle("Y");
    fluxMapvector[i]->GetZaxis()->SetTitle("Flux [TeV^{-1} cm^{-2} s^{-1}]");
    fluxMapvector[i]->GetZaxis()->SetTitleOffset(1.4);
    if(radialprofile){
      //      grradialprofile[i]->SetMarkerStyle(20);
      //      grradialprofile[i]->Draw("AP");	  
      TString slicetitle;
      if(user_xslice)
	slicetitle="x";
      else
	slicetitle="y";
      hradialprofile[i]->GetXaxis()->SetTitle(slicetitle.Data());
      hradialprofile[i]->GetYaxis()->SetTitle("# photons");
      hradialprofile[i]->GetYaxis()->SetTitleOffset(1.4);
      hradialprofile[i]->SetTitle("Profile");
      hradialprofile[i]->SetMaximum(1.5*hradialprofile[i]->GetMaximum());
      hradialprofile[i]->Draw();
      if(saveplots)
	hradialprofile[i]->Write();
      
    }else{
      //      fluxMapvector[i]->Draw("surf2");
      fluxMapvector[i]->Draw("colz");
    }
    gPad->Modified();
    gPad->Update();


    c->cd(2);
    gPad->SetLeftMargin(0.128571);
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.118217);
    gPad->SetBottomMargin(0.128915);
    titlestr2+=EbinMin[i];
    titlestr2+=" TeV ";
    titlestr2+=EbinMax[i];
    titlestr2+=" TeV ";

    photonMapVector[i]->SetTitle(titlestr2.Data());
    photonMapVector[i]->GetXaxis()->SetTitle("X");
    photonMapVector[i]->GetYaxis()->SetTitle("Y");
    photonMapVector[i]->GetZaxis()->SetTitle("Photons");
    photonMapVector[i]->GetXaxis()->SetAxisColor(10);
    photonMapVector[i]->GetYaxis()->SetAxisColor(10);
    photonMapVector[i]->GetZaxis()->SetTitleOffset(1.4);
    if(gaus_smooth){
      tmpgaus[i]->Draw("colz");     
    }else{
      photonMapVector[i]->Draw("colz");
    }
    //    photonMapVector[i]->Draw("surf2");
    if(fitsrc){
    double radius=fGaussianFit[i]->GetParameter(1);
    double eradius_max=radius+fGaussianFit[i]->GetParError(1);
    double eradius_min=radius-fGaussianFit[i]->GetParError(1);

    double psfradius=fGaussianFit[i]->GetParameter(2);
    double xfitted=fGaussianFit[i]->GetParameter(3);
    double yfitted=fGaussianFit[i]->GetParameter(4);
    
    TEllipse *fitcircle = new TEllipse(xfitted,yfitted,radius,radius,0,360,0);
    fitcircle->SetLineColor(kOrange-3);
    fitcircle->SetFillStyle(0);
    fitcircle->SetLineWidth(2);
    fitcircle->Draw();

    TEllipse *fitcircle_min = new TEllipse(xfitted,yfitted,eradius_min,eradius_min,0,360,0);
    fitcircle_min->SetLineColor(kOrange-3);
    fitcircle_min->SetFillStyle(0);
    fitcircle_min->SetLineWidth(2);
    fitcircle_min->Draw();

    TEllipse *fitcircle_max = new TEllipse(xfitted,yfitted,eradius_max,eradius_max,0,360,0);
    fitcircle_max->SetLineColor(kOrange-3);
    fitcircle_max->SetFillStyle(0);
    fitcircle_max->SetLineWidth(2);
    fitcircle_max->Draw();

    TEllipse *psfcircle = new TEllipse(xfitted,yfitted,psfradius,psfradius,0,360,0);
    psfcircle->SetLineColor(kYellow);
    psfcircle->SetFillStyle(0);
    psfcircle->SetLineWidth(2);
    psfcircle->SetLineStyle(2);
    psfcircle->Draw();
    }

    gPad->Modified();
    gPad->Update();
  
    //    if(fitsrc)
      //    fGaussianFit[i]->Draw("same");
    if(saveplots){
      fluxMapvector[i]->Write();
      photonMapVector[i]->Write();
//       if(radialprofile)
// 	grradialprofile[i]->Write();
      if(fitsrc)
 	fGaussianFit[i]->Write();
    }
  }
  }  
  if(saveplots)
    fsaveplots->Close();
  
}

// **************************************************************************************************************
// SkyMaps Radial Profile
// **************************************************************************************************************

void call_radialprofile(std::vector<TH2D*> photonMap )
{
  
  float rsize=0.1;
  float nr=1./rsize;
  if(nr<0.03){
    std::cout<< "Wedge in the radial profile is smaller than bin size, please increase the value " << std::endl;
    exit(0);
  }
  int inr=(int)nr;
#ifdef DEBUG_SHAPE
  std::cout<< " Step on the Radial Profile Histogram : " << inr << std::endl;
#endif
  const int iinr=inr;  
  double cont[iinr];
  double econt[iinr];
  double radiopoint[iinr];
  double radius_min[iinr];
  double radius_max[iinr];
  //  double rarea[iinr];
  for( size_t i=0; i!=photonMap.size(); ++i ) {
    
    TH2D* map2profile=(TH2D*)photonMap[i]->Clone("map2profile");
    
    for(int ir=0;ir<nr;ir++){
      radius_min[ir]=ir*rsize;
      radius_max[ir]=ir*rsize+rsize;
      //      rarea[ir]=TMath::Pi()*((radius_max[ir]*radius_max[ir])-(radius_min[ir]*radius_min[ir]));
      
    }
    for(int ir=0;ir<nr;ir++){
      cont[ir]=0;      
      econt[ir]=0;      
      radiopoint[ir]=radius_min[ir]+(radius_max[ir]-radius_min[ir])/2.;
    }
    
    for(int xbin=1;xbin<=map2profile->GetNbinsX();xbin++){
      for(int ybin=1;ybin<=map2profile->GetNbinsY();ybin++){
	      
	double xcenter=map2profile->GetXaxis()->GetBinCenter(xbin);
	double ycenter=map2profile->GetYaxis()->GetBinCenter(ybin);
	double r=0;
	double contbin=0;
	  r=sqrt((xcenter*xcenter)+(ycenter*ycenter));
	  contbin=map2profile->GetBinContent(xbin,ybin);

	for(int ir=0;ir<nr;ir++){
	  if(r>radius_min[ir] && r<=radius_max[ir]){
	    cont[ir]+=contbin;
	  }    
	}
      }
    }
    const int nmax=10000;
    int nrg=0;
    
    double radiopointg[nmax];
    double contg[nmax];
    
    for(int j=0;j<nr;j++){      
      if(cont[j]>0){
	radiopointg[nrg]=radiopoint[j];
	contg[nrg]=cont[j];
	econt[nrg]=sqrt(contg[j]);
	nrg++;
    }
      }


    grradialprofile[i]=new TGraphErrors(nrg,radiopointg,contg,NULL,econt);
  }

  std::cout << " Finishing radial profile " << std::endl;

}


// **************************************************************************************************************
// Fit SkyMaps
// **************************************************************************************************************
void call_fit(std::vector<TH2D*> photonMap, TString cnfname, TString psf_outfile, double Gaussian_sigma){
  
    FILE *fsavefit;
    fsavefit = fopen(psf_outfile.Data(), "a+" );
    
    TString fGaussianname("fGaussianname");
    for( size_t i=0; i!=photonMap.size(); ++i ) {
      fGaussianname+=i;
      fGaussianFit[i]=new TF2(fGaussianname.Data(),gausfunctionPSF,-rangex,rangex,-rangey,rangey,5);          
    
      TH2D* map2fit=(TH2D*)photonMap[i]->Clone("map2fit");
      map2fit->GetXaxis()->SetRangeUser(-Gaussian_sigma-0.1,Gaussian_sigma+0.1);
      map2fit->GetYaxis()->SetRangeUser(-Gaussian_sigma-0.1,Gaussian_sigma+0.1);
      psfsigmaFit[i] = GetPSF( logEbinMin[i], logEbinMax[i], cnfname);
      
#ifdef DEBUG
      std::cout << " the psf for [" << logEbinMin[i] << " , " << logEbinMax[i] << "] is " << psfsigmaFit[i] << std::endl; 
#endif
      TString fitname("gaussianfit_");
      fitname+=i;
      

      //      fGaussian[i]->SetNpx(50);
      fGaussianFit[i]->SetParNames("Constant","Sigma_Fitted","Sigma_PSF");
      fGaussianFit[i]->SetLineWidth(1);
      fGaussianFit[i]->SetLineColor(kBlack);
      fGaussianFit[i]->SetParLimits(0,1e-22,1000);
      fGaussianFit[i]->SetParLimits(1,0,Gaussian_sigma+2.5);
      fGaussianFit[i]->SetParameter(1,Gaussian_sigma);
      fGaussianFit[i]->FixParameter(2,psfsigmaFit[i]);
      
      if(create_excess){
	fGaussianFit[i]->SetParameter(4,0);
	fGaussianFit[i]->SetParameter(5,0);
      }
      else{
	fGaussianFit[i]->SetParLimits(4,map2fit->GetXaxis()->GetXmin(),map2fit->GetXaxis()->GetXmax());
	fGaussianFit[i]->SetParLimits(5,map2fit->GetYaxis()->GetXmin(),map2fit->GetYaxis()->GetXmax());
      }
      map2fit->Fit(fGaussianFit[i],"0Q");
      
      std::cout << "**************************************************************************************" << std::endl;
      std::cout << "Results of the Fit" << std::endl;
      std::cout << "**************************************************************************************" << std::endl;
      std::cout << "Norm :"  << fGaussianFit[i]->GetParameter(0) << std::endl;
      std::cout << "Sigma :"  << fGaussianFit[i]->GetParameter(1) << std::endl;
      std::cout << "SigmaPSF :"  << fGaussianFit[i]->GetParameter(2) << std::endl;
      std::cout << "Chi2/NDF :"  << fGaussianFit[i]->GetChisquare() << " / " << fGaussianFit[i]->GetNDF() << std::endl;
      
      if(savefit){
	//	fprintf(fsavefit,"Normalization, Sigma Gaussian, Sigma PSF, Chi2, NDF \n");
	fprintf(fsavefit,"%e %e %e %e %e %e %d \n",Gaussian_sigma,fGaussianFit[i]->GetParameter(0),fGaussianFit[i]->GetParameter(1),fGaussianFit[i]->GetParError(1),fGaussianFit[i]->GetParameter(2),fGaussianFit[i]->GetChisquare(),fGaussianFit[i]->GetNDF());
      }
    }
    
    fclose(fsavefit);
}

int SetOffsetBins(float selected_offset){

  for(int i=0;i<NBINS_OFFSET-1;i++){
    if(selected_offset >= OffSet_Bins[i] && selected_offset < OffSet_Bins[i+1]){
      return(OffSet_Bins[i]);
    }
  }
  if(selected_offset <= OffSet_Bins[0])
    return(OffSet_Bins[0]);
  if(selected_offset >= OffSet_Bins[NBINS_OFFSET-1])
    return(OffSet_Bins[NBINS_OFFSET-2]);
  
  return(-1);
  

}
// **************************************************************************************************************
// Create Gaussian Shape 
// **************************************************************************************************************

std::vector<TH2D*> CreateShape(double FluxCrabUnits,double Gaussian_sigma, int selected_shape, int selected_spectrum,double EExpcut)
{
  //  gRandom->SetSeed(0); 
  std::vector<TH2D*> photonMaps;
  double norm[Ebins];
  double normHi[Ebins];
  TH2D *Energy2DHisto[Ebins];
  TF1* iflux=0;
  TF1* ifluxHi=0;
  
  Src_Offset=sqrt(user_Xpos*user_Xpos + user_Ypos*user_Ypos); // distance to the center of the camera
  if((Src_Offset*Src_Offset)>DISTMAX*DISTMAX)
    {
      std::cout << "ERROR: source is too far away from the camera center, please use less than " << DISTMAX << " deg " << std::endl;
      exit(0);
    }

  if(selected_spectrum==0){
    iflux = new TF1("Differential Flux",pwlfunction,Eth,Emax,2); // pwl function from 0.01 to 100
    iflux->FixParameter(0,FluxCrabUnits*crabflux);
    iflux->FixParameter(1,Gamma); 

#ifdef DEBUG_SHAPE
    std::cout << "******************************************************************************************" << std::endl;
    std::cout << " You selected a power law spectrum with this parameters                                   " << std::endl; 
    std::cout << "% Crab " << FluxCrabUnits << " * " << crabflux << " : " << FluxCrabUnits*crabflux<<  std::endl
	      << " Spectral Index : " << Gamma << std::endl;
    std::cout << "******************************************************************************************" << std::endl;
#endif
  }else{
    iflux = new TF1("Differential Flux",expfunction,Eth,Emax,3); // pwl function from 0.01 to 100
    iflux->FixParameter(0,FluxCrabUnits*crabflux);
    iflux->FixParameter(1,Gamma);
    iflux->FixParameter(2,EExpcut);
  }
#ifdef DEBUG_SHAPE
  std::cout << "******************************************************************************************" << std::endl;
  std::cout << " You selected an exponential power law spectrum with this parameters                      " << std::endl; 
  std::cout << "% Crab " << FluxCrabUnits << " * " << crabflux << " : " << FluxCrabUnits*crabflux<<  std::endl
	    << " Spectral Index : " << Gamma << " EExpcut : " << EExpcut <<                           std::endl;
  std::cout << "******************************************************************************************" << std::endl;
#endif

  if(selected_shape==3){
    ifluxHi = new TF1("Differential FluxHi",pwlfunction,Eth,Emax,2); // pwl function from 0.01 to 100
    ifluxHi->FixParameter(0,FluxCrabUnits*crabflux);
    // A cooling effect of Gamma-2.7 for the high energy bin
    ifluxHi->FixParameter(1,Gamma-2.7); 
  }
  /*******************************************************/
  /* Set Differential or Integral Flux in each Energy Bin*/
  /*******************************************************/
    
  for(int ibin=0;ibin<Ebins;ibin++){
    normHi[ibin]=0;
    if(differential){
#ifdef DEBUG_SHAPE
  std::cout << "******************************************************************************************" << std::endl;
  std::cout << " You selected differential energy bin flux "                                                << std::endl;
  std::cout << "******************************************************************************************" << std::endl;
#endif
      norm[ibin]=iflux->Integral(EbinMin[ibin],EbinMax[ibin]);        
      if(selected_shape==3)
	normHi[ibin]=ifluxHi->Integral(EbinMin[ibin],EbinMax[ibin]);        
    }else{
#ifdef DEBUG_SHAPE
  std::cout << "******************************************************************************************" << std::endl;
  std::cout << " You selected integral energy bin flux "                                                << std::endl;
  std::cout << "******************************************************************************************" << std::endl;
#endif
      norm[ibin]=iflux->Integral(EbinMin[ibin],Emax);  
      if(selected_shape==3)
	normHi[ibin]=ifluxHi->Integral(EbinMin[ibin],Emax);        
    }
  }
 /*******************************************************/
  /* Define Functions for different shapes              */
  /******************************************************/
  TF2* fGaussian[Ebins];   // Gaussian-like Excess
  TF2* fRing[Ebins];       // Shell-like Excess
  TF2* fadd[Ebins];        // Shell+Gaussian Excess
  TF2* fadd2gaus[Ebins];   // 2 Gaussians  : hard-PL + soft-extended 
  //  TF2* fGaussianHi[Ebins];

  for(int ibin=0;ibin<Ebins;ibin++){

    TString gausname("fGaussian");
    gausname+=ibin;
    fGaussian[ibin] = new TF2(gausname.Data(),gausfunctionPSF,-rangex/2.,rangex/2.,-rangey/2.,rangey/2.,5);
    fGaussian[ibin]->FixParameter(2,0.0);   // PSF Gaussian =0
    fGaussian[ibin]->FixParameter(1,Gaussian_sigma);
    fGaussian[ibin]->FixParameter(4,user_Xpos);   // PSF Gaussian =0
    fGaussian[ibin]->FixParameter(5,user_Ypos);   // PSF Gaussian =0

    TString ringname("fRing");
    ringname+=ibin;
    fRing[ibin] =new TF2(ringname.Data(),ringfunction,-rangex/2.,rangex/2.,-rangey/2.,rangey/2.,6);
    fRing[ibin]->FixParameter(2,0.0);       // PSF Gaussian =0
    fRing[ibin]->FixParameter(1,ExternalGaus);
    fRing[ibin]->FixParameter(3,InternalGaus);
    fRing[ibin]->FixParameter(4,user_Xpos);
    fRing[ibin]->FixParameter(5,user_Ypos);
    

    TString faddname("faddname");
    faddname+=ibin;
    fadd[ibin]=new TF2(faddname.Data(),fadding,-rangex/2.,rangex/2.,-rangey/2.,rangey/2.,7);

    TString fadd2gausname("fadd2gausname");
    fadd2gausname+=ibin;
    fadd2gaus[ibin]=new TF2(fadd2gausname.Data(),fadding2gauss,-rangex/2.,rangex/2.,-rangey/2.,rangey/2.,6);
  }


  switch (selected_shape) 
    {
    case 0:                                 // ---------------------------------------------------- > Gaussian
      for(int ibin=0;ibin<Ebins;ibin++){
#ifdef DEBUG_SHAPE
	std::cout << "******************************************************************************************" << std::endl;
	std::cout << " YOU SELECTED A GAUSSIAN SHAPE "                                                            << std::endl;
	std::cout << "******************************************************************************************" << std::endl;
#endif
	fGaussian[ibin]->FixParameter(0,norm[ibin]);
	TString histoname("Energy2DHName_");	
      	histoname+=ibin;
	Energy2DHisto[ibin]=new TH2D(histoname.Data(),"Flux Histogram ",nbinsx,-rangex/2.,rangex/2.,nbinsy,-rangey/2.,rangey/2.);
	for( Int_t x=1; x<=Energy2DHisto[ibin]->GetNbinsX(); ++x ) {
	  for( Int_t y=1; y<=Energy2DHisto[ibin]->GetNbinsY(); ++y ) {	
	    Energy2DHisto[ibin]->Fill(Energy2DHisto[ibin]->GetXaxis()->GetBinCenter(x),Energy2DHisto[ibin]->GetYaxis()->GetBinCenter(y),fGaussian[ibin]->Eval(Energy2DHisto[ibin]->GetXaxis()->GetBinCenter(x),Energy2DHisto[ibin]->GetYaxis()->GetBinCenter(y)));
	  }
	}
      }
      break;
    case 1:                                  // ---------------------------------------------------- > Shell
#ifdef DEBUG_SHAPE
	std::cout << "******************************************************************************************" << std::endl;
	std::cout << " YOU SELECTED A SHELL SHAPE    "                                                            << std::endl;
	std::cout << "******************************************************************************************" << std::endl;
#endif
      for(int ibin=0;ibin<Ebins;ibin++){
	fRing[ibin]->FixParameter(0,norm[ibin]);
	TString histoname("Energy2DHName_");	
	histoname+=ibin;
	Energy2DHisto[ibin]=new TH2D(histoname.Data(),"Flux Histogram ",nbinsx,-rangex/2.,rangex/2.,nbinsy,-rangey/2.,rangey/2.);
	for( Int_t x=1; x<=Energy2DHisto[ibin]->GetNbinsX(); ++x ) {
	  for( Int_t y=1; y<=Energy2DHisto[ibin]->GetNbinsY(); ++y ) {	
	    Energy2DHisto[ibin]->Fill(Energy2DHisto[ibin]->GetXaxis()->GetBinCenter(x),Energy2DHisto[ibin]->GetYaxis()->GetBinCenter(y),fRing[ibin]->Eval(Energy2DHisto[ibin]->GetXaxis()->GetBinCenter(x),Energy2DHisto[ibin]->GetYaxis()->GetBinCenter(y)));
	  }
	}
      }
      break;
    case 2:                              // ---------------------------------------------------- > Shell+Gauss
#ifdef DEBUG_SHAPE
	std::cout << "******************************************************************************************" << std::endl;
	std::cout << " YOU SELECTED A SHELL+Gauss SHAPE  "                                                        << std::endl;
	std::cout << "******************************************************************************************" << std::endl;
#endif
      for(int ibin=0;ibin<Ebins;ibin++){

	fadd[ibin]->FixParameter(0,FluxRatio*norm[ibin]);
	fadd[ibin]->FixParameter(1,0.04); // Pl source
	//	fadd[ibin]->FixParameter(2,0.2); // inner radius
	//	fadd[ibin]->FixParameter(3,0.22); // outer radius
	fadd[ibin]->FixParameter(2,InternalGaus); // inner radius
	fadd[ibin]->FixParameter(3,ExternalGaus); // outer radius
	fadd[ibin]->FixParameter(4,norm[ibin]);
	fadd[ibin]->FixParameter(5,user_Xpos);
	fadd[ibin]->FixParameter(6,user_Ypos);

	TString histoname("Energy2DHName_");	
	histoname+=ibin;
	Energy2DHisto[ibin]=new TH2D(histoname.Data(),"Flux Histogram ",nbinsx,-rangex/2.,rangex/2.,nbinsy,-rangey/2.,rangey/2.);
	for( Int_t x=1; x<=Energy2DHisto[ibin]->GetNbinsX(); ++x ) {
	  for( Int_t y=1; y<=Energy2DHisto[ibin]->GetNbinsY(); ++y ) {	
	    Energy2DHisto[ibin]->Fill(Energy2DHisto[ibin]->GetXaxis()->GetBinCenter(x),Energy2DHisto[ibin]->GetYaxis()->GetBinCenter(y),fadd[ibin]->Eval(Energy2DHisto[ibin]->GetXaxis()->GetBinCenter(x),Energy2DHisto[ibin]->GetYaxis()->GetBinCenter(y)));
	  }
	}
      }
      break;
    case 3:                          // ---------------------------------------------------- > Cooling
	std::cout << "******************************************************************************************" << std::endl;
	std::cout << " YOU SELECTED A COOLING EFFECT     "                                                        << std::endl;
	std::cout << "******************************************************************************************" << std::endl;

	for(int ibin=0;ibin<Ebins;ibin++){
	  fadd2gaus[ibin]->FixParameter(0,FluxRatio*norm[ibin]); 
	  fadd2gaus[ibin]->FixParameter(1,Gaussian_sigma); // it has to be large enough (~1.7)
	  fadd2gaus[ibin]->FixParameter(2,0.1);
	  fadd2gaus[ibin]->FixParameter(3,normHi[ibin]);
	  fadd2gaus[ibin]->FixParameter(4,user_Xpos);
	  fadd2gaus[ibin]->FixParameter(5,user_Ypos);
#ifdef DEBUG_SHAPE	  
	  std::cout << " Norm 1 : "<<  norm[ibin] << " Norm 2 : " <<  normHi[ibin] << std::endl;
	  std::cout <<" Gaussian 1:"<< Gaussian_sigma << " Gaussian 2: " << 0.07 << std::endl;
#endif
	TString histoname("Energy2DHName_");	
	histoname+=ibin;
	Energy2DHisto[ibin]=new TH2D(histoname.Data(),"Flux Histogram ",nbinsx,-rangex,rangex,nbinsy,-rangey,rangey);
	for( Int_t x=1; x<=Energy2DHisto[ibin]->GetNbinsX(); ++x ) {
	  for( Int_t y=1; y<=Energy2DHisto[ibin]->GetNbinsY(); ++y ) {	
	    Energy2DHisto[ibin]->Fill(Energy2DHisto[ibin]->GetXaxis()->GetBinCenter(x),Energy2DHisto[ibin]->GetYaxis()->GetBinCenter(y),fadd2gaus[ibin]->Eval(Energy2DHisto[ibin]->GetXaxis()->GetBinCenter(x),Energy2DHisto[ibin]->GetYaxis()->GetBinCenter(y)));
	  }
	}
	}
	break;
	
    default:
      printf("Please select 0: Gaussian Shape, 1: Ring Shape, 2: Composite SNR, 3: Cooled Gaussian, \n");
      break;
    }
  for(int ibin=0;ibin<Ebins;ibin++){
    photonMaps.push_back(Energy2DHisto[ibin]);
    
  }
  return photonMaps;
}


// **************************************************************************************************************
// Read from file 
// **************************************************************************************************************
std::vector<TH2D*> ReadFromFile(double FluxCrabUnits,TString read_shapename,TString read_histoname, int selected_spectrum)
{

  std::vector<TH2D*> photonMaps;
  TFile* fShape = 0;
  fShape = new TFile(read_shapename.Data());
#ifdef DEBUG
  std::cout << " Using the file " << read_shapename.Data() << std::endl;
#endif
  if( !fShape->IsOpen() ) {
    std::cout << "Error: could not open shape file!" << std::endl;
    exit(0);
  }
  double norm[Ebins];
  TF1* iflux=0;
  if(selected_spectrum==0){
    iflux = new TF1("Differential Flux",pwlfunction,Eth,Emax,2); // pwl function from 0.01 to 100
    iflux->FixParameter(0,FluxCrabUnits*crabflux);
    iflux->FixParameter(1,Gamma); 
#ifdef DEBUG_SHAPE
    std::cout << "******************************************************************************************" << std::endl;
    std::cout << " You selected a power law spectrum with this parameters                                   " << std::endl; 
    std::cout << "% Crab " << FluxCrabUnits << " * " << crabflux << " : " << FluxCrabUnits*crabflux<<  std::endl
	      << " Spectral Index : " << Gamma << std::endl;
    std::cout << "******************************************************************************************" << std::endl;
#endif
  }else{
    iflux = new TF1("Differential Flux",expfunction,0.01,100,3); // pwl function from 0.01 to 100
    iflux->FixParameter(0,FluxCrabUnits*crabflux);
    iflux->FixParameter(1,Gamma);
    iflux->FixParameter(2,Ecut);
  }
#ifdef DEBUG_SHAPE
  std::cout << "******************************************************************************************" << std::endl;
  std::cout << " You selected an exponential power law spectrum with this parameters                      " << std::endl; 
  std::cout << "% Crab " << FluxCrabUnits << " * " << crabflux << " : " << FluxCrabUnits*crabflux<<  std::endl
	    << " Spectral Index : " << Gamma << " EExpcut : " << EExpcut <<                           std::endl;
  std::cout << "******************************************************************************************" << std::endl;
#endif
  

  for(int i=0; i<Ebins;i++){
    TString excess_name=read_histoname;    
    //excess_name+=i;
    fExcess[i]=(TH2D*)fShape->Get(excess_name.Data());
    if( fExcess[i] == 0 ) 
      std::cout << "Error: could not access Excess histogram" << std::endl;

    // Get distance to the camera //
    double xmin=fExcess[i]->GetXaxis()->GetXmin();
    double xmax=fExcess[i]->GetXaxis()->GetXmax();
    double ymin=fExcess[i]->GetYaxis()->GetXmin();
    double ymax=fExcess[i]->GetYaxis()->GetXmax();
    double xc = (xmin + xmax) /2.;
    double yc = (ymin + ymax) /2.;
    Src_Offset = sqrt((user_Xpos-xc)*(user_Xpos-xc) + (user_Ypos-yc)*(user_Ypos-yc));
    if((Src_Offset*Src_Offset)>DISTMAX*DISTMAX)
      {
	std::cout << "ERROR: source is too far away from the camera center, please use less than " << DISTMAX << " deg " << std::endl;
	exit(0);
      }
    // Normalize histo with user flux

    norm[i]=iflux->Integral(EbinMin[i],EbinMax[i]); 

    double scale=1./fExcess[i]->Integral();
    scale *=norm[i];

    std::cout << " integral " << fExcess[i]->Integral() << " " << norm[i] <<  endl;
    fExcess[i]->Scale(scale);


    photonMaps.push_back(fExcess[i]); // todo: energy dependence morphology
  }
  

  return photonMaps;    

 
  if( fShape->IsOpen() )
    fShape->Close();

  

}

// **************************************************************************************************************
// Construct Photon Maps with CTA effective Area and PSF 
// **************************************************************************************************************
std::vector<TH2D*> ConstructPhotonMap(std::vector<TH2D*> fluxMaps,TString subarray, bool IsEx)
{

  std::vector<TH2D*> photonMapVector;
  std::vector<TH2D*> protonMapVector;
  std::vector<TH2D*>::iterator it = fluxMaps.begin();

  TH2D *photonmap[Ebins];
  TH2D *bgmaps[Ebins];
  TH2D *bgmaps1[Ebins];
  TH2D *bgmaps2[Ebins];
  TH2D *bgmaps3[Ebins];
  TH2D *bgmaps4[Ebins];
  TH2D *hOnExcess[Ebins];

  double psfsigma=0;
  int j=0;

  for( ; it != fluxMaps.end(); ++it ) {
    gRandom->SetSeed(0);

    TString hExcessName("PhotonMap_Bin");
    hExcessName+=j;
    photonmap[j] = (TH2D*)(*it)->Clone(hExcessName.Data());
    photonmap[j]->SetTitle(hExcessName.Data());

    /*******************************************************/
    /* Get PSF for given configuration File               */
    /******************************************************/
    psfsigma = GetPSF( logEbinMin[j], logEbinMax[j], subarray);
    /*******************************************************/
    /* Modified to read the extended effective areas */
    /*******************************************************/
    //    if (Gaussian_sigma < psfsigma || Is_IFAEoffsetEA ) 
    //    if (Gaussian_sigma < psfsigma ) 
    //      efficgampsf=1.0;    
    //    else
    //      efficgampsf=1.2;
    //    if(Gaussian_sigma < psfsigma)
    //      Is_IFAEoffsetEA=0;
    //    else
    //      Is_IFAEoffsetEA=1;
	
    
    double effArea =  1e4 * GetEffectiveArea( logEbinMin[j], logEbinMax[j] , subarray); // cm^2
    //    double effArea = 1e4 * GetEffectiveArea( logEbinMin[j], logEbinMax[j] , subarray); // cm^2
    photonmap[j]->Scale( effArea );                        // Scale with effective Aera
    photonmap[j]->Scale( observationTime );                // Scale with observation time
    double binArea = GetBinArea( photonmap[j],photonmap[j]->GetNbinsX() / 2,
				 photonmap[j]->GetNbinsY() / 2 );
    
    binArea *= TMath::DegToRad() * TMath::DegToRad();
      
    photonmap[j]->Scale( binArea );                       // Scale with the bin area

#ifdef DEBUG
      std::cout << "logE: [" << logEbinMin[j] << ", " << logEbinMax[j] << "]"
		<< " effArea [cm^2]: " << effArea
		<< " observationTime [s]: " << observationTime
		<< " bin area [sr]: " << binArea
		<< std::endl;
#endif
      // todo: test random fluctuations
      bool createFluctuations=false;
      if(createFluctuations) {
	TRandom3* random = new TRandom3();
	TH2D* generator=(TH2D*)photonmap[j]->Clone("generator");
	for( Int_t x=1; x<=generator->GetNbinsX(); ++x ) {
	  for( Int_t y=1; y<=generator->GetNbinsY(); ++y ) {
	    int n = random->Poisson( generator->GetBinContent(x,y));
	    if( n >0)
	      std::cout << n << std::endl;
	  }
	}
      }
       
#ifdef DEBUG       
       std::cout << "Smooth skymap with " << psfsigma << " deg" << std::endl;
#endif

       TH2D* temp = 0;
       //       if(create_excess)
       temp = GausSmooth(photonmap[j],psfsigma);
       //       else
       //	 temp =(TH2D*)photonmap[j]->Clone("temp");
       TString hExcessSmoothName("Excess_Smooth_Bin");
       hExcessSmoothName+=j;
       temp->SetName(hExcessSmoothName.Data());
       temp->SetTitle(hExcessSmoothName.Data());


       photonMapVector.push_back( temp );


       /*******************************************************/
       /*  CREATE THE SKYMAP FOR THE BACKGROUND               */
       /*  ON OFF BACKGROUND MAKER - 3 OFF REGIONS            */
       /*******************************************************/
       TString hBgName("BgMap_Bin");
       hBgName+=j;
       bgmaps[j] = (TH2D*) (*it)->Clone(hBgName.Data());
       bgmaps[j]->SetTitle(hBgName.Data());
       bgmaps[j]->Reset("ICE");

       TString hBgName1("BgMap1_Bin");
       hBgName1+=j;
       bgmaps1[j] = (TH2D*) (*it)->Clone(hBgName1.Data());
       bgmaps1[j]->SetTitle(hBgName1.Data());
       bgmaps1[j]->Reset("ICE");

       TString hBgName2("BgMap2_Bin");
       hBgName2+=j;
       bgmaps2[j] = (TH2D*) (*it)->Clone(hBgName2.Data());
       bgmaps2[j]->SetTitle(hBgName2.Data());
       bgmaps2[j]->Reset("ICE");

       TString hBgName3("BgMap3_Bin");
       hBgName3+=j;
       bgmaps3[j] = (TH2D*) (*it)->Clone(hBgName3.Data());
       bgmaps3[j]->SetTitle(hBgName3.Data());
       bgmaps3[j]->Reset("ICE");

       TString hBgName4("BgMap4_Bin");
       hBgName4+=j;
       bgmaps4[j] = (TH2D*) (*it)->Clone(hBgName4.Data());
       bgmaps4[j]->SetTitle(hBgName4.Data());
       bgmaps4[j]->Reset("ICE");

       if(Calc_Significance){
       TString hSigmaName("Significance_Bin");
       hSigmaName+=j;
       //       hSigma[j] = (TH2D*) (*it)->Clone(hSigmaName.Data());
       //       hSigma[j]->Reset("ICE");
       int nbinsx_rb=(int)(*it)->GetNbinsX()/2;
       int nbinsy_rb=(int)(*it)->GetNbinsY()/2;
       int xmin_rb=(int)(*it)->GetXaxis()->GetXmin();
       int xmax_rb=(int)(*it)->GetXaxis()->GetXmax();
       int ymin_rb=(int)(*it)->GetYaxis()->GetXmin();
       int ymax_rb=(int)(*it)->GetYaxis()->GetXmax();

       hSigma[j] = new TH2D("Significance","Significance",nbinsx_rb,xmin_rb,xmax_rb,nbinsy_rb,ymin_rb,ymax_rb);
       }
       double rate = GetBgRate(logEbinMin[j],logEbinMax[j], subarray); // Hz/deg2
       rate *= observationTime*binArea;

       /*******************************************************/
       /*  Create Gaussian Distribution                       */
       /*******************************************************/
       TF1* fgenerator=new TF1("fgenerator","gaus(0)",rate-3*sqrt(rate),rate+3*sqrt(rate));
       fgenerator->FixParameter(0,1);
       fgenerator->FixParameter(1,rate);
       fgenerator->FixParameter(2,sqrt(rate));

       for (Int_t x=1; x<=bgmaps1[j]->GetNbinsX();x++) {
	   for( Int_t y=1; y<=bgmaps1[j]->GetNbinsY(); ++y ) {
	     double n=fgenerator->GetRandom();
	     double m=fgenerator->GetRandom();
	     double w=fgenerator->GetRandom();
	     double p=fgenerator->GetRandom();
#ifdef DEBUG
	     //	     //	     std::cout << "Generating bg number "  << n << " " << m << " " << w << std::endl;
#endif 

	     bgmaps1[j]->SetBinContent(x,y,n);
	     bgmaps2[j]->SetBinContent(x,y,m);
	     bgmaps3[j]->SetBinContent(x,y,w);
	     bgmaps4[j]->SetBinContent(x,y,p);
	   }
       }

       
       bgmaps[j]->Add(bgmaps1[j]);
       bgmaps[j]->Add(bgmaps2[j]);
       bgmaps[j]->Add(bgmaps3[j]);
       bgmaps[j]->Scale(1./3.);
       

       TH2D* tempBg = 0;
       TH2D* tempBgOn = 0;
       if(create_excess){
	 tempBg = GausSmooth(bgmaps[j],psfsigma);
	 tempBgOn = GausSmooth(bgmaps4[j],psfsigma);
       }       else{
	 tempBg =(TH2D*)bgmaps[j]->Clone("tempBg");
	 tempBgOn =(TH2D*)bgmaps4[j]->Clone("tempBgOn");
       }
       TString hBackgroundSmoothName("Background_Smooth_Bin");
       hBackgroundSmoothName+=j;
       tempBg->SetName(hBackgroundSmoothName.Data());
       tempBg->SetTitle(hBackgroundSmoothName.Data());

       TString hBackgroundOnSmoothName("BackgroundOn_Smooth_Bin");
       hBackgroundOnSmoothName+=j;
       tempBgOn->SetName(hBackgroundOnSmoothName.Data());
       tempBgOn->SetTitle(hBackgroundOnSmoothName.Data());
       

       TString hOnName("OnName_Bin");
       hOnName+=j;
       hOnExcess[j] = (TH2D*)tempBgOn->Clone(hOnName.Data());
       hOnExcess[j]->Add(temp);

 
       //       hOnExcess[j]->Add(bgmaps[j],-1);
       hOnExcess[j]->Add(tempBg,-1);

       // Calculate Significance Maps using Li&Ma approach

    if(Calc_Significance){
      TH2D* tempOn = 0;
      TH2D* tempOff = 0;
      tempOn = (TH2D*)hOnExcess[j]->Rebin2D();
      tempOff = (TH2D*)tempBg->Rebin2D();
      double alpha=3.;
      for (Int_t x=1; x<=tempOn->GetNbinsX();x++) {
	 for( Int_t y=1; y<=tempOn->GetNbinsY(); ++y ) {
	   double non=tempOn->GetBinContent(x,y);
	   double noff=tempOff->GetBinContent(x,y);
	   double sign=0;
	   if(non>0 && noff>0)
	     sign=LiMa(non,noff,alpha);
	   else{
	     sign=0;
	   }

	   hSigma[j]->SetBinContent(x,y,sign);
	 }
       }       

    }
       protonMapVector.push_back( hOnExcess[j] );
       
  

       
     j++;       
  }

  // todo: ring bg evaluation

  

  if(IsEx)
    return photonMapVector;
  else
    return protonMapVector;
  
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
  Double_t Xpos=par[5]; // position in camera coordinates
  Double_t Ypos=par[6]; 
  Double_t spl=par[1]; // point like source
  Double_t sinner=par[2];
  Double_t souter=par[3];
  Double_t r2 = pow(X-Xpos,2)+pow(Y-Ypos,2);
  Double_t z  = par[0]*(exp(-r2/(2*souter*souter))-exp(-r2/(2*sinner*sinner)))+par[4]*exp(-r2/(2*spl*spl));
  return z;
}
// **************************************************************************************************************
// Adding two Gaussian 
// **************************************************************************************************************
Double_t fadding2gauss(Double_t *x, Double_t *par)
{
  Double_t X=x[0];
  Double_t Y=x[1];
  Double_t Xpos=par[3]; // position in camera coordinates
  Double_t Ypos=par[4]; 
  Double_t s1=par[1];
  Double_t s2=par[2];
  Double_t r2 = pow(X-Xpos,2)+pow(Y-Ypos,2);
  Double_t z  = par[0]*(exp(-r2/(2*s1*s1))) + par[3]*exp(-r2/(2*s2*s2));
  return z;
}
// **************************************************************************************************************
// Shell-type function 
// **************************************************************************************************************
Double_t ringfunction(Double_t *x, Double_t *par)
{
  
  Double_t X=x[0];
  Double_t Y=x[1];
  Double_t Xpos=par[4]; // position in camera coordinates
  Double_t Ypos=par[5]; 
  Double_t s=par[1];
  Double_t sinner=par[3];
  Double_t spsf=par[2];
  Double_t r2 = pow(X-Xpos,2)+pow(Y-Ypos,2);

  Double_t z  = par[0]*(exp(-r2/(2*((s*s)+(spsf*spsf))))-exp(-r2/(2*((sinner*sinner)+(spsf*spsf)))));
    return z;
}
// **************************************************************************************************************
// Gaus Function and PSF 
// **************************************************************************************************************
Double_t gausfunctionPSF(Double_t *x, Double_t *par)
{
  
  Double_t X=x[0];
  Double_t Y=x[1];
  Double_t Xpos=par[3]; // position in camera coordinates
  Double_t Ypos=par[4]; 
  Double_t s=par[1];
  Double_t spsf=par[2];
  Double_t r2 = pow(X-Xpos,2)+pow(Y-Ypos,2);

  Double_t z  = par[0]*exp(-r2/(2*((s*s)+(spsf*spsf))));
    return z;
}

// **************************************************************************************************************
// Read Effective Area
// **************************************************************************************************************
double GetEffectiveArea( double logEMin, double logEMax, TString cnfname ) {
  
  static TH1F* effArea1D = 0; 
  static TH2F* effArea2D = 0; 
  TFile* fEffArea = 0;
  if( effArea1D == 0 || effArea2D ==0 ) {
    TString filename(configpath);
    if(Is_IFAEoffsetEA){
      filename+="Subarray";
      filename+=cnfname;
      filename+="_offaxis.root"; 
}    else{
      filename+="kb_";
      filename+=cnfname;
      filename+="_50h_20deg_v3.root";
    }
    fEffArea = new TFile(filename.Data());
    std::cout << " Using the file " << filename << " to read the effective area for configuration " << cnfname << std::endl;
    if( !fEffArea->IsOpen() ) {
      std::cerr << "Error: could not open effArea file!" << std::endl;
      exit(0);      
    }
    if(Is_IFAEoffsetEA){
       fEffArea->GetObject("EffectiveArea_offaxis", effArea2D);
       
   }  else{
      fEffArea->GetObject("EffectiveArea", effArea1D);
    }
    if( effArea1D == 0 && effArea2D == 0) {
      std::cerr << "Error: could not access effArea histogram" << std::endl;      
      exit(0);
    }
  }
  double area=0;
  if(Is_IFAEoffsetEA){
    int lowBinEnergy = effArea2D->GetXaxis()->FindBin(logEMin);//,OffsetMin); 
    int highBinEnergy = effArea2D->GetXaxis()->FindBin(logEMax);
    int lowBinOffset = effArea2D->GetYaxis()->FindBin(OffsetMin);
    int highBinOffset = effArea2D->GetYaxis()->FindBin(OffsetMax);
    area =effArea2D->Integral(lowBinEnergy, highBinEnergy,lowBinOffset,highBinOffset);
    area   /= ( (highBinOffset - lowBinOffset + 1) * (highBinEnergy - lowBinEnergy + 1) );
   }
   else{
     
     int lowBin  = effArea1D->FindBin(logEMin);
    int highBin = effArea1D->FindBin(logEMax);
    area = effArea1D->Integral(lowBin, highBin);
    area       /= highBin - lowBin + 1;
    
    //    fEffArea->Close();
   }

  if( fEffArea->IsOpen() ) 
    fEffArea->Close();
#ifdef DEBUG_SHAPE
  std::cout << "average effective area: " << area << std::endl;
#endif
  return area;
}

// **************************************************************************************************************
// Read Background Rate
// **************************************************************************************************************
double GetBgRate( double logEMin, double logEMax, TString cnfname ) {
  
  static TH1F* bgRate1D = 0; 
  static TH2F* bgRate2D = 0;
  TFile* fBgRate = 0;
  if( bgRate1D == 0 || bgRate2D ==0 ) {
    TString filename(configpath);
    if(Is_IFAEoffsetEA){
      filename+="Subarray";
      filename+=cnfname;
      filename+="_offaxis.root"; 
    }
    else{
      filename+="kb_";
      filename+=cnfname;
      filename+="_50h_20deg_v3.root";
    }
    
    fBgRate = new TFile(filename.Data());
    std::cout << " Using the file " << filename.Data() << " to read the Background rate for configuration " << cnfname <<std::endl;
    if( !fBgRate->IsOpen() ) {
      std::cerr << "Error: could not open BgRate file!" << std::endl;
      exit(0);
    }
    if(Is_IFAEoffsetEA)
      fBgRate->GetObject("BGRatePerSqDeg_offaxis", bgRate2D);
    else {
      fBgRate->GetObject("BGRatePerSqDeg", bgRate1D);
    }
    if( bgRate1D == 0 && bgRate2D ==0 ) {
      std::cerr << "Error: could not access background histogram" << std::endl;
      exit(0);
    }
    //    fBgRate->Close();
  }
  double rate = 0;
  if(Is_IFAEoffsetEA){
    int lowBinEnergy = bgRate2D->GetXaxis()->FindBin(logEMin);//,OffsetMin); 
    int highBinEnergy = bgRate2D->GetXaxis()->FindBin(logEMax);
    int lowBinOffset = bgRate2D->GetYaxis()->FindBin(OffsetMin);
    int highBinOffset = bgRate2D->GetYaxis()->FindBin(OffsetMax);
    rate =bgRate2D->Integral(lowBinEnergy, highBinEnergy,lowBinOffset,highBinOffset);
    rate   /= ( (highBinOffset - lowBinOffset + 1) * (highBinEnergy - lowBinEnergy + 1) );
  }else{
    int lowBin  = bgRate1D->FindBin(logEMin);
    int highBin = bgRate1D->FindBin(logEMax);
    rate = bgRate1D->Integral(lowBin, highBin);
    rate       /= highBin - lowBin + 1;
  }
#ifdef DEBUG
  std::cout << "average rate: " << rate << std::endl;
#endif
  if(fBgRate->IsOpen())
    fBgRate->Close();

  return rate;
}

// **************************************************************************************************************
// Read PSF
// **************************************************************************************************************
double GetPSF( double logEMin, double logEMax , TString cnfname) { 

  static TH1F* PSF1D = 0; 
  static TH2F* PSF2D = 0; 
  TFile* fPSF = 0;
  if( PSF1D == 0 && PSF2D == 0 ) {
    TString filename(configpath);
    if(Is_IFAEoffsetEA){
      filename+="Subarray";
      filename+=cnfname;
      filename+="_offaxis.root"; 
    }
    else{
    filename+="kb_";
    filename+=cnfname;
    filename+="_50h_20deg_v3.root";
    }
    fPSF = new TFile(filename.Data());

    if( !fPSF->IsOpen() ) {
      std::cerr << "Error: could not open PSF file!" << std::endl;
      exit(0);
    }
    if(Is_IFAEoffsetEA)
      fPSF->GetObject("AngRes80_offaxis", PSF2D);
    else{
      fPSF->GetObject("AngRes80",PSF1D);
    }
    if( PSF1D == 0 && PSF2D == 0 ) {
      std::cerr << "Error: could not access PSF histogram" << std::endl;
      exit(0);
    }
    //    fPSF->Close();
  }
  double value = 0;
  if(Is_IFAEoffsetEA){
    int lowBinEnergy = PSF2D->GetXaxis()->FindBin(logEMin);//,OffsetMin); 
    int highBinEnergy = PSF2D->GetXaxis()->FindBin(logEMax);
    int lowBinOffset = PSF2D->GetYaxis()->FindBin(OffsetMin);
    int highBinOffset = PSF2D->GetYaxis()->FindBin(OffsetMax);
    value = PSF2D->Integral(lowBinEnergy, highBinEnergy,lowBinOffset,highBinOffset);
    value  /= ( (highBinOffset - lowBinOffset + 1) * (highBinEnergy - lowBinEnergy + 1) );
  }else{
    int lowBin  = PSF1D->FindBin(logEMin);
    int highBin = PSF1D->FindBin(logEMax);
  value = PSF1D->Integral(lowBin, highBin);
  value       /= highBin - lowBin + 1;
  }  



  return value;
  if(fPSF->IsOpen())
    fPSF->Close();
}
// **************************************************************************************************************
float GetBinArea(TH2D *histo,int xbin, int ybin) 
{
  float xlow = histo->GetXaxis()->GetBinLowEdge(xbin);
  float xup  = histo->GetXaxis()->GetBinUpEdge(xbin);
  float ylow = histo->GetYaxis()->GetBinLowEdge(ybin);
  float yup  = histo->GetYaxis()->GetBinUpEdge(ybin);


  float dx = xup - xlow;
  float dy = yup - ylow;

  return dx*dy;

}
// **************************************************************************************************************
//  Produce a new, oversampled Histogram produced by
// filling each bin with the integrated contents within a defined
// minimum and maximum radius of the bin. The bin contents are
// weighted using a Gauss function.
//
// \param sigma - standard deviation (degrees)
//  HESS - HD team
// **************************************************************************************************************

TH2D* GausSmooth(TH2D *histo, float sigma)
{

  TString sname=histo->GetName();
  sname+="_smooth";
  TH2D *sksmooth = (TH2D*)histo->Clone("sname");

  
 sksmooth->Reset("ICE");
  //  std::cout << " The name of the smoothed histogram " << sname << std::endl;
  double s2 = pow(sigma,2);
  int range = static_cast<int>(30.0*sigma/histo->GetYaxis()->GetBinWidth(1));
  int startX = 1;
  int startY = 1;
  int binsX = histo->GetNbinsX();
  int binsY = histo->GetNbinsY();
  int pcount = 0;
  int nbins =0;

  std::vector<std::vector<double> > lookup(range+1); 
  for (int l=0; l<=range; ++l) {
    double dx = histo->GetXaxis()->GetBinWidth(1)*l;
    lookup[l] = std::vector<double>(range+1);
    for (int m=0; m<=range; ++m) {
      double dy = histo->GetYaxis()->GetBinWidth(1)*(double)(m);
      double dist = (dx*dx + dy*dy)/s2;
      lookup[l][m] = exp(-dist/2.);
    }
  }
  
  startX = histo->GetXaxis()->FindFixBin(user_Xpos-rangex/2.);
  binsX  = histo->GetXaxis()->FindFixBin(user_Xpos+range/2.);
  
  startY = histo->GetYaxis()->FindFixBin(user_Ypos-rangey/2.);
  binsY  = histo->GetYaxis()->FindFixBin(user_Ypos+rangey/2.);

  
  if (startX < 1) startX = 1;
  if (startY < 1) startY = 1;
  if (binsX > histo->GetXaxis()->GetNbins() ) binsX = histo->GetXaxis()->GetNbins();
  if (binsY > histo->GetYaxis()->GetNbins() ) binsY = histo->GetYaxis()->GetNbins();

#ifdef DEBUG
  std::cout << " binX : " << binsX << std::endl;
  std::cout << " binY : " << binsY << std::endl;
#endif
  for (int i = startX; i < binsX; ++i) {
    double bx = histo->GetXaxis()->GetBinCenter(i);
    for (int j = startY; j < binsY; ++j) {
      double by = histo->GetYaxis()->GetBinCenter(j);
      int bin = histo->FindBin(bx,by);
      double binContWobble = histo->GetBinContent(bin);
      double binErrorWobble = histo->GetBinError(bin);
      
      pcount++; 
      nbins++;
      if(pcount >= 80000) { 
	std::cout << "Smoothed " 
		  << (int)(nbins/((double)(binsX*binsY))*100)
		  << "%" << std::endl;
	pcount=0;
      }
      
      if ((fabs(binContWobble) > 0.0001)) {
	for (int l=i-range; l<=i+range; ++l) {
	  int l2=abs(l-i);
	  if (l>0 && l<=histo->GetXaxis()->GetNbins()) {
	    for (int m=j-range; m<=j+range; ++m) {
	      int m2=abs(m-j);
	      if (m > 0 && m<=histo->GetYaxis()->GetNbins()) {
		int bin2 = histo->GetBin(l,m);
		double binContWob2 = sksmooth->GetBinContent(bin2);
		double binError2 = sksmooth->GetBinError(bin2);
		binContWob2 += binContWobble * lookup[l2][m2];
		binError2 = sqrt(pow(binErrorWobble * lookup[l2][m2],2.) + pow(binError2,2.));
		sksmooth->SetBinContent(bin2,binContWob2);
		sksmooth->SetBinError(bin2,binError2);

	      }  
	    }
	  }
	}
      }
    }
  }
  

  return sksmooth;
}

TString sigmar(double gamma_zero, double radius, double index_change){

  TF1 * func = new TF1("func","[0]+x*[1]",0.0,2.);
  func->FixParameter(0,gamma_zero);
  func->FixParameter(1,index_change);

  TString specstrDeff;
  specstrDeff = Form("(%e)*pow((x/1.),-%.2f)",1.,func->Eval(radius));
  return specstrDeff;

}

/********** FROM HESS-HAP CODE : Heidelberg Team ********////////
/**
* Make slice along x or y axis.
* 
*/
void call_radialprofile2(std::vector<TH2D*> photonMap, bool xslice,double halfwidth,double range1,double range2,bool average)
{   
  TH1D* newhist = 0;
  for( size_t i=0; i!=photonMap.size(); ++i ) {
    Int_t  nbinsx_rp = photonMap[i]->GetXaxis()->GetNbins();
    Int_t  nbinsy_rp = photonMap[i]->GetYaxis()->GetNbins();  

    int nbinsSliceX = nbinsx_rp;
    int nbinsSliceY = nbinsy_rp;

    Float_t xmin  = 0.;
    Float_t xmax  = 0.;
    Float_t ymin  = 0.;
    Float_t ymax  = 0.;
    Float_t centre= 0.;
    if(xslice) {
      ymin  = photonMap[i]->GetYaxis()->GetXmin();
      ymax  = photonMap[i]->GetYaxis()->GetXmax();
      centre=(ymax+ymin)/2.;
      if(fabs(range1) > 1e-14 || fabs(range2) > 1e-14){
	xmin  = range1;
	xmax  = range2;
	int binmin = photonMap[i]->GetXaxis()->FindFixBin(xmin);
	int binmax = photonMap[i]->GetXaxis()->FindFixBin(xmax);
	nbinsSliceX = binmax - binmin + 1;


      } else {
	xmin  = photonMap[i]->GetXaxis()->GetXmin();
	xmax  = photonMap[i]->GetXaxis()->GetXmax();

      }
    } else {
      xmin  = photonMap[i]->GetXaxis()->GetXmin();
      xmax  = photonMap[i]->GetXaxis()->GetXmax();
      centre=(xmax+xmin)/2.;
      if(fabs(range1) > 1e-14 || fabs(range2) > 1e-14){
	ymin  = range1;
	ymax  = range2;
	int binmin = photonMap[i]->GetYaxis()->FindFixBin(ymin);
	int binmax = photonMap[i]->GetYaxis()->FindFixBin(ymax);
	nbinsSliceY = binmax - binmin + 1;

      } else {
	ymin  = photonMap[i]->GetYaxis()->GetXmin();
	ymax  = photonMap[i]->GetYaxis()->GetXmax();

      }
      
    }
#ifdef DEBUG
    std::cout << "Generating histogram with " << nbinsSliceX << " bins " << std::endl;
#endif

    if (xslice) {
      newhist = new TH1D("xslice","xslice",nbinsSliceX,xmin,xmax);
    } else {
      newhist = new TH1D("yslice","yslice",nbinsSliceY,ymin,ymax);
    }
    TH1D* binsFilled = (TH1D*)newhist->Clone("binsFilled");
    binsFilled->Sumw2();
    newhist->Sumw2();
  
    double maxZ = 0.0;
    int bincounter = 0;
  
    for (int biny=1;biny<=nbinsy_rp;biny++) {
      Float_t by  = photonMap[i]->GetYaxis()->GetBinCenter(biny);
      for (int binx=1;binx<=nbinsx_rp;binx++) {
      
	Float_t bx = photonMap[i]->GetXaxis()->GetBinCenter(binx);

	Int_t bin2d = photonMap[i]->GetBin(binx,biny);
	
	Int_t newbin = 0;
	if (xslice) {
	  newbin = newhist->FindBin(bx);
	} else {
	  newbin = newhist->FindBin(by);
	}

	if(newbin == 0 || newbin >= newhist->GetNbinsX()+1) continue;
	if ( ((xslice) && (fabs(by-centre)<halfwidth)) || 
	     ((!xslice) && (fabs(bx-centre)<halfwidth)) ) { 


	  //	if ( ((xslice) && (fabs(by)<halfwidth)) || 
	  //	     ((!xslice) && (fabs(bx)<halfwidth)) ) { 

	  ++bincounter;

	  
	  Double_t newcont = newhist->GetBinContent(newbin);
	  Double_t oldcont = photonMap[i]->GetBinContent(bin2d);
	  Double_t sum = newcont + oldcont;
	  if (sum > maxZ) maxZ = sum;
	  newhist->SetBinContent(newbin, sum);
	  double newbinError = 
	    sqrt((photonMap[i]->GetBinError(bin2d)*photonMap[i]->GetBinError(bin2d))+
		 (newhist->GetBinError(newbin)*newhist->GetBinError(newbin)));
	  newhist->SetBinError(newbin, newbinError);

	  double binsSoFar = binsFilled->GetBinContent(newbin);
	  binsSoFar += 1.;
	  binsFilled->SetBinContent(newbin,binsSoFar);
	  binsFilled->SetBinError(newbin,0.);
	  
	}
      }
    }

#ifdef DEBUG    
    std::cout << "MakeSlice number of bins: " << bincounter << std::endl;
#endif

    if(average) newhist->Divide(binsFilled);
    else newhist->SetMaximum(maxZ*1.1);
    delete binsFilled;


    TString hradialprofilename("hrpname_");
    hradialprofilename+=i;
    hradialprofile[i] = (TH1D*)newhist->Clone("hradialprofilename");
    delete newhist;
  }
  std::cout << " Radial Profile Created " << std::endl;
}

/******************************/

double LiMa(float non, float noff, double alpha) { 
   
  double excess = non - noff*alpha;
  double n1,n2,sign,a;

  if(alpha == 0) return sqrt(1.0*non);

  if (excess > 0) {
    n1 = non;
    n2 = noff;
    sign = 1.0;
    a = alpha;
  } else {
    n2 = non;
    n1 = noff;
    sign = -1.0;
    a = 1.0/alpha;
  }
  
  // A development rather gives this values for
  // non == 0 or noff == 0
  if(n2 == 0) 
    return sign * sqrt(2 * n1 * log((1+a)/a));
  if(n1 == 0)  
    return sign * sqrt(2 * n2 * log(1+a));
  
  
  double nt = n1+n2;
  double pa = 1 + a;

  double t1 = n1*log((pa/a)*(n1/nt));
  double t2 = n2*log(pa*(n2/nt));
  double sig = sqrt(2.)*sqrt(t1+t2);

  return sign*sig;  
}

/******************************/
void help(){
  
std::cout << "\n*********** USAGE ********" << std::endl;
 std::cout << "makeCTAmaps(\"Config - [E]\",  % Crab [1], Gaussian Sigma [0.1], Basename for output" << std::endl;
 // std::cout << "Option 2: Ring radius are fixed at 0.19 and 0.20 " << std::endl
 //	   << "Option 3: Ring radius are fixed at 0.2 and 0.22 and Gaussian source is PL sigma=0.04"<< std::endl
 //	   << "Option 4: Fist Gaussian with spectral index 2.8 and extension given by the use, Second Gaussian with spectral index 2.8-0.7 and PL (0.07)"
 //	   << std::endl;
 
 
 std::cout << std::endl << std::endl;
}
