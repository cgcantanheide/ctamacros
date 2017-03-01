// To be used with astroroot 
// To start astroroot
// source /cern/root5/setup.sh
// astroroot
// E. de Ona Wilhelmi
// Ex: test_fits("hess_rxj1713_2005.fits","cor_excess_smooth_zoom","outfilename.root")
void convertFITSto2D(TString file_name, TString image_name,TString outfile_name){

  TFFloatImg * img = TFReadImage(file_name.Data(),image_name.Data());
   if (img == NULL)
      {
      TFError::PrintErrors();
      return;
      }
   std::cout << " The header of your file is :" << endl;
   TFAttrIter i_attr = img->MakeAttrIterator();
   char hstr[500];
   while(i_attr.Next()){
     printf("%-16s : %20s %-5s | %s\n",i_attr->GetName(),i_attr->GetStringValue(hstr),i_attr->GetUnit(),i_attr->GetComment());
   }

   // create a TH2D histogram and display it
   TH2D * hist = img->MakeHisto();
   TH2D *hExcess = (TH2D*)hist->Clone("Excess");   

   Double_t CRVAL1 = dynamic_cast<TFDoubleAttr&>(img->GetAttribute("CRVAL1"));
   Double_t CRVAL2 = dynamic_cast<TFDoubleAttr&>(img->GetAttribute("CRVAL2"));
   Double_t CDELT1 = dynamic_cast<TFDoubleAttr&>(img->GetAttribute("CDELT1"));
   Double_t CDELT2 = dynamic_cast<TFDoubleAttr&>(img->GetAttribute("CDELT2"));

   Double_t CRPIX1 = dynamic_cast<TFDoubleAttr&>(img->GetAttribute("CRPIX1"));
   Double_t CRPIX2 = dynamic_cast<TFDoubleAttr&>(img->GetAttribute("CRPIX2"));

   Double_t nbinsx=2*CRPIX1 - 1.;
   Double_t nbinsy=2*CRPIX2 - 1.;

   double xmin=0;
   double xmax=0;
   double ymin=0;
   double ymax=0;
   if(CDELT1<0){
     xmin=-CRVAL1+(CRPIX1*CDELT1);
     xmax=-CRVAL1-((CRPIX1+1)*CDELT1);
   }else{
     xmin=CRVAL1-(CRPIX1*CDELT1);
     xmax=CRVAL1+((CRPIX1+1)*CDELT1);
   }
   ymin=CRVAL2-(CRPIX2*CDELT2);
   ymax=CRVAL2+((CRPIX2+1)*CDELT2);


   
   TH2D *hExcess=new TH2D("Excess","Excess",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
   for (int k=1; k<=hist->GetNbinsX(); k++) {
       for (int l=1; l<=hist->GetNbinsY(); l++) {
	 hExcess->SetBinContent(k,l,hist->GetBinContent(k,l));
       }
   }
   
   hExcess->Draw("colz");
   TFile *f=new TFile(outfile_name.Data(),"RECREATE");
   hExcess->Write();
   f->Close();
   
}
