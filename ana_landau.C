#include <iostream>
#include <fstream>


// Fit parameters for wide end landgau
// 200, 1174, 30000, 200
// Result of the fit:
// 54.86
// 1042
// 2.074e4
// 220.5
double langau(double *x, double *par) {

    //Fit parameters:
       //par[0]=Width (scale) parameter of Landau density
       //par[1]=Most Probable (MP, location) parameter of Landau density
       //par[2]=Total area (integral -inf to inf, normalization constant)
       //par[3]=Width (sigma) of convoluted Gaussian function
      
    // Numeric constants
          Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
          Double_t mpshift  = -0.22278298;       // Landau maximum location

          // Control constants
          Double_t np = 1000.0;      // number of convolution steps
          Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

          // Variables
          Double_t xx;
          Double_t mpc;
          Double_t fland;
          Double_t sum = 0.0;
          Double_t xlow,xupp;
          Double_t step;
          Double_t i;

          // MP shift correction
          mpc = par[1] - mpshift * par[0];

          // Range of convolution integral
          xlow = x[0] - sc * par[3];
          xupp = x[0] + sc * par[3];

          step = (xupp-xlow) / np;

          // Convolution integral of Landau and Gaussian by sum
          for(i=1.0; i<=np/2; i++) {
             xx = xlow + (i-.5) * step;
             fland = TMath::Landau(xx,mpc,par[0]) / par[0];
             sum += fland * TMath::Gaus(x[0],xx,par[3]);

             xx = xupp - (i-.5) * step;
             fland = TMath::Landau(xx,mpc,par[0]) / par[0];
             sum += fland * TMath::Gaus(x[0],xx,par[3]);
          }

          return (par[2]*step*sum*invsq2pi / par[3]);
}

struct landau_fit_result{
  double mu;
  double sigma_landau;
  double amplitude;
  double sigma_gauss;

  double mu_error;
  double sigma_landau_error;
  double amplitude_error;
  double sigma_gauss_error;
};

landau_fit_result ana_landau(std::string fname, double Q_to_phe, std::array<double,4> parameters, std::array<double,2> langau_function_range, std::array<double,2> shift_range, std::string plot_name, std::string h_title, bool do_fit = true, int rebin_value = -1){
  
  fstream ifile(fname);
  double xmin = 1e20, xmax = -1e20, delta; 


  double x, y;
  vector<double> yval; 
  vector<double> xval; 

  // rmeove useless first lines
  string dummy;
  for(int i=0; i<5; i++){
    getline(ifile, dummy);
  }
  
  // x = bin x, y = entries in i-th bin
  while(ifile>>x>>y){
    x = x * Q_to_phe;
    if(yval.size() == 0) delta = -x;
    if(yval.size() == 1) delta += x;
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    xval.push_back(x);
    yval.push_back(y);
  }


  TH1D* h_shift = new TH1D("h_shift", "h_shift", yval.size(), xmin-delta, xmax+delta);
  for(unsigned int i = 0; i < yval.size(); i++){
    double x = xval[i];
    double y = yval[i];
    for (unsigned int y_iter = 0; y_iter < y; y_iter++){
      h_shift->Fill(x);
    }
    // h->SetBinContent(i+1, yval[i]);
  }
  h_shift->GetXaxis()->SetRangeUser(shift_range[0], shift_range[1]);
  // double shift_mean = -300;
  double shift_mean = h_shift->GetMean();
  

  
  TH1D* h = new TH1D(Form("h_%s", plot_name.c_str()), h_title.c_str(), yval.size(), xmin-delta - shift_mean, xmax+delta - shift_mean);
  for(unsigned int i = 0; i < yval.size(); i++){
    double x = xval[i] - shift_mean;
    double y = yval[i];
    for (unsigned int y_iter = 0; y_iter < y; y_iter++){
      h->Fill(x);
    }
    // h->SetBinContent(i+1, yval[i]);
  }
  if(rebin_value > 0 ){
    h->Rebin(rebin_value);
  }

  std::cout << "xmin, xmax: " << xmin << ", " << xmax << std::endl;

  TF1* flangau = new TF1("flangau", langau, langau_function_range[0], langau_function_range[1], 4);
  flangau->SetParameters(parameters.data());
  flangau->SetParName(0, "#sigma_{L}");
  flangau->SetParName(1, "#mu");
  flangau->SetParName(2, "A");
  flangau->SetParName(3, "#sigma_{G}");
  // flangau->SetParLimits(1, 0, 1000); 

  TCanvas* c = new TCanvas(Form("c_%s", plot_name.c_str()), Form("c_%s", plot_name.c_str()), 800, 600);
  // c->SetLogy();
  h->Draw();
  if(do_fit){
    h->Fit(flangau, "RQ"); 
    flangau->Draw("same");
  }
  // c->SaveAs("test.pdf");



  // TCanvas* c = new TCanvas("c", "c", 800, 600);
  // // c->SetLogy();
  // h->Draw();
  // h->Fit(land_neg, "RQ"); 
  // land_neg->Draw("same");

  std::array<double, 4> fit_vals = {flangau->GetParameter(0),flangau->GetParameter(1),flangau->GetParameter(2),flangau->GetParameter(3)};
  // delete h;
  // delete c;
  // delete flangau;
  landau_fit_result lr;
  lr.mu = flangau->GetParameter(1);
  lr.sigma_landau = flangau->GetParameter(0);
  lr.amplitude = flangau->GetParameter(2);
  lr.sigma_gauss = flangau->GetParameter(3);

  lr.mu_error = flangau->GetParError(1);
  lr.sigma_landau_error = flangau->GetParError(0);
  lr.amplitude_error = flangau->GetParError(2);
  lr.sigma_gauss_error = flangau->GetParError(3);
  std::cout << "DEBUG" << std::endl << std::endl << "MU: " << lr.mu << std::endl;

  return lr;

}
