void plot_MPVs(double energy_per_cm, std::vector<double> x_values, std::vector<double> x_err, std::vector<landau_fit_result> landau_fit_result_vector, std::vector<double> penetration_depths, bool use_amplitude){
    gStyle->SetOptFit(0);

    TGraphErrors* gr = new TGraphErrors();
    std::string using_amplitude_str;
    if(use_amplitude){
        using_amplitude_str = "amp";
    }
    else
    {
        using_amplitude_str = "charge";
    }
    gr->SetTitle(Form("MPV vs Distance from SiPM (using %s);Distance (cm);MPV (#phe/MeV)", using_amplitude_str.c_str()));
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(2);
    // gr->SetMarkerColor(kAzure + 2);
    // gr->SetLineColor(kAzure + 2);
    gr->SetLineWidth(4);

    for(int ii = 0; ii < x_values.size(); ii++){
        landau_fit_result lr = landau_fit_result_vector[ii];
        double x = x_values[ii];
        double y = lr.mu;
        double y_error = lr.mu_error;
        double x_error = x_err[ii];
        std::cout << "y_values"<<": "<<y<<", " << energy_per_cm*penetration_depths[ii]<<std::endl;
        y = y/energy_per_cm/penetration_depths[ii];
        gr->SetPoint(ii, x, y);
        y_error = y_error/energy_per_cm/penetration_depths[ii];
        gr->SetPointError(ii, x_error, y_error);
    }
    gr->Print();

    // TF1* flin = new TF1("flin", "pol1", x_values[0], x_values[x_values.size() - 1]);
    // gr->Fit(flin, "R0");
    TCanvas* cMPV = new TCanvas("cMPV", "cMPV", 800, 600);

    gr->Draw("AP");
    // flin->Draw("same"); 
    cMPV->SaveAs(Form("plots/MPV_%s_%dpoints.pdf",using_amplitude_str.c_str(), x_values.size()));
    

}