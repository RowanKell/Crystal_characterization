#include <iostream>
#include <fstream>

void ana(){
  
  fstream ifile("./Test_crs/F2--large-end--00000.txt");

  double xmin = 1e20, xmax = -1e20, delta; 

  double x, y;
  vector<double> yval; 

  // remove useless first lines
  string dummy;
  for(int i=0; i<5; i++){
    getline(ifile, dummy);
  }
  
  // x = bin x, y = entries in i-th bin
  while(ifile>>x>>y){
    if(yval.size() == 0) delta = -x;
    if(yval.size() == 1) delta += x;
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    yval.push_back(y);
    cout << x << " " << y << endl;
  }

  cout << xmin << " " << xmax << " " << delta << endl;

  TH1D* h = new TH1D("h", "h", yval.size(), xmin-delta, xmax+delta);
  for(unsigned int i = 0; i < yval.size(); i++){
    h->SetBinContent(i+1, yval[i]);
  }

  TCanvas* c = new TCanvas("c", "c", 800, 600);
  h->Draw(); 

}
