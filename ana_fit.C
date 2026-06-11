#include <iostream>
#include <fstream>

// Gain = 2
// gaussian fit:
/*
    NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 

   1  p0           1.53868e+02   1.74296e+00  -1.29523e-04   2.19794e-04
   2  p1           1.08227e-03   1.34069e-06   1.24279e-09  -2.07854e+02
   3  p2           8.59892e-05   1.05847e-06   1.38296e-09   1.16138e+03
   4  p3           6.69506e+01   1.11809e+00  -2.94614e-04   1.14015e-04
   5  p4           1.38493e-03   2.93738e-06   4.19477e-09  -2.44787e+02
   6  p5           1.06635e-04   2.46663e-06  -8.56876e-10  -3.47004e+01
   7  p6           8.79378e+00   4.21179e-01  -4.70300e-04   1.49476e-03
   8  p7           1.70000e-03   8.73416e-06  -4.95795e-04** at limit **
   9  p8           2.23329e-04   9.11291e-06   2.63978e-08   4.94797e+01

 */
// linear fit:
/*
p0                        =  0.000777317   +/-   3.45784e-06 
p1                        =  0.000304756   +/-   2.68352e-06 
 */

// Gain ~40
//Gaussian fit
/*
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           7.14115e+02   6.53338e+00   1.08288e-01  -1.29536e-04
   2  p1           9.84592e-03   6.01021e-06   1.15036e-07  -8.59999e+00
   3  p2           6.84463e-04   5.44414e-06   7.70628e-08  -3.52129e+01
   4  p3           7.95184e+01   2.32041e+00   3.98105e-02  -3.60837e-04
   5  p4           1.83608e-02   2.39102e-05   5.11208e-07   3.61036e+01
   6  p5           8.33767e-04   2.50960e-05   3.74746e-07  -4.20823e+00
   7  p6           7.63024e+01   9.96889e-01   1.44353e-02   7.99051e-04
   8  p7           1.23959e-02   6.34568e-05   9.60812e-07   1.72571e+01
   9  p8           4.94058e-03   4.58177e-05   6.05721e-07   1.26676e+01

 */

//linear fit:
/*
p0                        =   0.00133106   +/-   2.67617e-05 
p1                        =   0.00851486   +/-   2.4654e-05  

 */

// Finding Scale
// scale = p1_2 / p1_40 = 1/28.3

// Gain ~40 Integrated charge
// Gaussian fit:
/*
   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  p0           4.71210e+02   3.59339e+00  -1.16344e-03  -8.71978e-06
   2  p1          -4.09945e-12   1.73669e-12   9.93615e-16   1.45656e+07
   3  p2           2.49344e-10   1.48706e-12  -6.31406e-16  -8.27983e+07
   4  p3           1.95255e+02   2.21372e+00  -3.35568e-04  -1.03528e-04
   5  p4          -1.04134e-09   3.96207e-12   1.94864e-15  -2.03027e+07
   6  p5           2.88617e-10   4.03680e-12   1.30347e-15  -9.50525e+07
   7  p6           4.77218e+01   1.03146e+00  -7.04016e-05  -6.39189e-05
   8  p7          -2.12215e-09   1.27378e-11  -4.74156e-15  -7.08121e+06
   9  p8           3.93091e-10   1.26548e-11  -1.22771e-15  -1.16287e+07
*/
//Linear fit:
/*
p0                        =  1.04029e-09   +/-   4.64971e-12 
p1                        = -1.04382e-09   +/-   3.67718e-12 
 */
// Charge_2 = scale * Charge_40 = 

std::map<int,std::string> multigauss_map = {
  {1, "[0]*TMath::Gaus(x, [1], [2])"},
  {2, "[0]*TMath::Gaus(x, [1], [2])+[3]*TMath::Gaus(x, [4], [5])"},
  {3, "[0]*TMath::Gaus(x, [1], [2])+[3]*TMath::Gaus(x, [4], [5]) + [6] *TMath::Gaus(x, [7],[8])"},
  {4, "[0]*TMath::Gaus(x, [1], [2])+[3]*TMath::Gaus(x, [4], [5]) + [6] *TMath::Gaus(x, [7],[8])+ [9] *TMath::Gaus(x, [10],[11])"},
  {5, "[0]*TMath::Gaus(x, [1], [2])+[3]*TMath::Gaus(x, [4], [5]) + [6] *TMath::Gaus(x, [7],[8])+ [9] *TMath::Gaus(x, [10],[11])+ [12] *TMath::Gaus(x, [13],[14])"}
};

std::array<double,2> ana_fit(std::string fname, std::vector<double> fit_params, int n_points, int n_points_linear_fit, double* function_range, bool show_plots = false){
  
  //fstream ifile("./Test_crs/F1--phe-g2--00000.txt");
  fstream ifile(fname);

  double xmin = 1e20; // minimum x value
  double xmax = -1e20; // maximum x value
  double delta; // bin width

  double x, y; // temporary values to save x and y from file
  vector<double> yval; // array to save all values of y

  // remove useless first lines
  string dummy;
  for(int i=0; i<5; i++){
    getline(ifile, dummy);
  }
  
  while(ifile>>x>>y){ // while you read two values from file 
    // use first two entries to caluclate bin width
    if(yval.size() == 0) delta = -x;
    if(yval.size() == 1) delta += x;
    // update minimum and maximum x values
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    // add y value to list of y values
    yval.push_back(y);
  }
  std::cout << "yval size: " << yval.size() << std::endl;
  std::cout << "xmin: " << xmin << std::endl;
  std::cout << "xmax: " << xmax << std::endl;

 

  // initialize histogram. Parameters: name, title, number of bins, minimum x value, maximum x value
  TH1D* h = new TH1D("h", "h", yval.size(), xmin-delta, xmax+delta);
  for(unsigned int i = 0; i < yval.size(); i++){ // fill histogram with y values (bin content)
    h->SetBinContent(i+1, yval[i]);
  }

  /*
  // plot on canvas
  TCanvas* c = new TCanvas("c", "c", 800, 600);
  h->Draw("histo"); // Draw histogram with option "histo" to draw only the histogram without error bars

  // define a gaussian function. Parameters: name, formula, minimum x value, maximum x value
  TF1* fgaus = new TF1("f", "gaus", 0.002, 0.0026); 
  fgaus->SetParameters(1, 0.0023, 0.0001); // set initial parameters for fit (amplitude, mean, sigma)
  h->Fit(fgaus, "R"); // fit histogram with gaussian function in range of function
  

  TCanvas* cfit1 = new TCanvas("cfit1", "cfit1", 800, 600);
  h->Draw("histo");
  fgaus->Draw("same"); // draw fit function on same canvas
  */
 std::cout << "Function: " << multigauss_map[n_points] << std::endl;
 
  // define a function thatis a sum of two gaussians. Parameters: name, formula, minimum x value, maximum x value 
  TF1* fmultigaus = new TF1("fmulti", multigauss_map[n_points].c_str(), function_range[0], function_range[1]); // define a function with two gaussians
  
  //TF1* fmultigaus = new TF1("fmulti", "[0]*TMath::Gaus(x, [1], [2])+[3]*TMath::Gaus(x, [4], [5])", -1.5e-9, 1e-9); // define a function with two gaussians
  //TF1* fmultigaus = new TF1("fmulti", "[0]*TMath::Gaus(x, [1], [2])+[3]*TMath::Gaus(x, [4], [5])+[6]*TMath::Gaus(x, [7], [8])+[9]*TMath::Gaus(x, [10], [11])+[12]*TMath::Gaus(x, [13], [14])+[15]*TMath::Gaus(x, [16], [17])", 0.002, 0.006); // define a function with two gaussians
    
  //double params[9] = {480, 0, 11.9e-10, 195, -1.04e-9, 9.6e-10, 51, -2.05e-9, 11.9e-10};
  // double params[9] = {458, -2.9e-11, 2.5e-10, 211, -1.05e-9, 2.5e-10, 48.05, -2.15e-9, 4.02e-10};
  fmultigaus->SetParameters(fit_params.data()); // set initial parameters for fit (amplitude, mean, sigma for each gaussian)
  //fmultigaus->FixParameter(0,458);
  //fmultigaus->FixParameter(1, -2.9e-11);
  //fmultigaus->FixParameter(2, 2.5e-10);
  //fmultigaus->FixParameter(4,-1.17e-9);
  //fmultigaus->FixParameter(5, 2.5e-10);
  //fmultigaus->FixParameter(7, -2.15e-9);
  //fmultigaus->SetParLimits(7, .0017, .003);
  // fmultigaus->SetParLimits(3, 90, 135);
  // fmultigaus->SetParLimits(6, 70, 95);
  // fmultigaus->SetParLimits(9, 45, 65);
  
  // h->Fit(fmultigaus); // fit histogram with multi-gaussian function in range of function
  h->Fit(fmultigaus, "R0"); // fit histogram with multi-gaussian function in range of function

  // cfit2->SaveAs("gauss.pdf");

  TGraphErrors* gr = new TGraphErrors(); // initialize a graph to plot the phe amplitude vs the phe number  

for(int ii = 0; ii < n_points_linear_fit; ii++){
    double mu = fmultigaus->GetParameter(3*ii+1); // get mean of i-th gaussian
    double emu = fmultigaus->GetParError(3*ii+1); // get error on mean of i-th gaussian
    gr->SetPoint(ii, ii+1, mu); // set point on graph with x value equal to mean and y value equal to amplitude
    gr->SetPointError(ii, 0, emu); // set error on point with x error equal to 0 and y error equal to error on mean
  }


  TF1* flin = new TF1("flin", "pol1", 0.5, n_points_linear_fit + 0.5); // define a linear function to fit the graph

  
  gr->Fit(flin, "R"); // fit graph with linear function in range of function
  
  if(show_plots){
    TCanvas* cfit2 = new TCanvas("cfit2", "cfit2", 800, 600);
    h->Draw("histo");
    fmultigaus->Draw("same"); // draw fit function on same canvas
    TCanvas* cgr = new TCanvas("cgr", "cgr", 800, 600);
    gr->Draw("AP"); // draw graph with option "AP" to draw axes and points
    gr->SetMarkerStyle(20); // set marker style for points

  }

  // gStyle->SetOptFit(1111); // set option to display fit parameters on plot
  // cgr->SaveAs("./linear.pdf");
  std::array<double, 2> linear_fit_result = {flin->GetParameter(0), flin->GetParameter(1)};
  return linear_fit_result;
  

}
