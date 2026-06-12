#include "ana_landau.C"
#include "ana_fit.C"
#include "ana_plot.C"

void ana_driver(bool use_amplitude = false){
    gStyle->SetOptFit(1111); // set option to display fit parameters on plot

    
    // Gain = 2, amplitude, gauss, SiPM measurement
    std::string fname_g2a = "./Test_crs/F1--phe-g2--00000.txt";
    std::string fname_plot_g2a =  "g2a";
    int n_points_gauss_g2a = 3;
    int n_points_linear_fit_g2a = 3;
    std::vector<double> gaussian_fit_params_g2a = {1.53868e+02,1.08227e-03,8.59892e-05,6.69506e+01,1.38493e-03,1.06635e-04,8.79378e+00,1.70000e-03,2.23329e-04};
    double gaussian_function_range_g2a[2] = {0.00091, 0.002}; // ????
    std::array<double, 2> fit_vals_g2a = ana_fit(fname_g2a, gaussian_fit_params_g2a, n_points_gauss_g2a, n_points_linear_fit_g2a, gaussian_function_range_g2a,fname_plot_g2a, true);

    double convert_to_phe;

    if(use_amplitude){
        convert_to_phe = 1/fit_vals_g2a[1];
    }
    else{
        // Gain = 40, amplitude, gauss, SiPM measurement
        std::string fname_g40a = "./Test_crs/F1--phe-g40--00000.txt";
    std::string fname_plot_g40a =  "g40a";
        int n_points_gauss_g40a = 3;
        int n_points_linear_fit_g40a = 2;
        std::vector<double> gaussian_fit_params_g40a = {7.14115e+02, 9.84592e-03, 6.84463e-04, 7.95184e+01, 1.83608e-02, 8.33767e-04, 7.63024e+01, 1.23959e-02, 4.94058e-03}; 
        double gaussian_function_range_g40a[2] = {0, 0.03};
        std::array<double, 2> fit_vals_g40a = ana_fit(fname_g40a, gaussian_fit_params_g40a, n_points_gauss_g40a, n_points_linear_fit_g40a, gaussian_function_range_g40a,fname_plot_g40a, true);
        std::cout << "gain = 40, linear fit vals: (" << fit_vals_g40a[0] << ", " << fit_vals_g40a[1] << std::endl;
        
        // Gain = 40, charge, gauss, SiPM measurement
        std::string fname_g40c = "./Test_crs/F2--phe-g40--00000.txt";
    std::string fname_plot_g40c =  "g40c";
        int n_points_gauss_g40c = 3;
        int n_points_linear_fit_g40c = 3;
        std::vector<double> gaussian_fit_params_g40c = {4.712e2, -4.099e-12, 2.493e-10, 1.953e2, -1.041e-9, 2.886e-10, 4.7721e1, -2.1222e-9, 3.9309e-10};
        double gaussian_function_range_g40c[2] = {-3e-9, 1e-9};
        std::array<double, 2> fit_vals_g40c = ana_fit(fname_g40c, gaussian_fit_params_g40c, n_points_gauss_g40c, n_points_linear_fit_g40c, gaussian_function_range_g40c,fname_plot_g40c, true);
        
        // Converts from Wb to Phe
        convert_to_phe = 1/(fit_vals_g40c[1]*fit_vals_g2a[1]/fit_vals_g40a[1]);

        double gain_40_real_gain = 2 * fit_vals_g40a[1] / fit_vals_g2a[1];
        std::cout << "Gain40 real gain: " << gain_40_real_gain << "\n";
        // double Q_to_phe = -1/3e-11;
    } 

    std::string fname_3p5;
    std::string fname_9;
    std::string fname_20;
    std::string fname_15;
    std::string fname_15_redo;
    std::string fname_6;
    std::array<double,4> landau_params_3p5;
    std::array<double,4> landau_params_9;
    std::array<double,4> landau_params_20;
    std::array<double,4> landau_params_15;
    std::array<double,4> landau_params_15_redo;
    std::array<double,4> landau_params_6;
    std::array<double,2> langau_function_range_3p5;
    std::array<double,2> langau_function_range_9;
    std::array<double,2> langau_function_range_20;
    std::array<double,2> langau_function_range_15;
    std::array<double,2> langau_function_range_15_redo;
    std::array<double,2> langau_function_range_6;
    std::array<double,2> landau_shift_range_3p5;
    std::array<double,2> landau_shift_range_9;
    std::array<double,2> landau_shift_range_20;
    std::array<double,2> landau_shift_range_15;
    std::array<double,2> landau_shift_range_15_redo;
    std::array<double,2> landau_shift_range_6;
    if(use_amplitude){
        fname_3p5 = "./Test_crs/F1--large-end--00000.txt";
        fname_9 = "./Test_crs/Middle_cosmics/F1--mid--00000.txt";
        fname_20 = "./Test_crs/small/F1--small--00000.txt";
        fname_15 = "./Test_crs/15cm/F1--15cm--00000.txt";
	fname_15_redo = "./Test_crs/15cm_again/F1--15cm_redo--00000.txt";
        fname_6 = "./Test_crs/6cm/F1--6cm--00000.txt";
        landau_params_3p5 = {30, 400, 5000, 30}; 
        landau_params_9 = {50, 350, 5000, 50}; 
        landau_params_20 = {20, 225, 3000, 20}; 
        landau_params_15 = {17, 212, 3000, 47};
        landau_params_15_redo = {17, 212, 3000, 47};
        landau_params_6 = {17, 212, 3000, 47};
        langau_function_range_3p5 = {250, 800};
        langau_function_range_9 = {200, 700};
        langau_function_range_20 =  {115, 500};
        langau_function_range_15 =  {128, 490};
        langau_function_range_15_redo =  {128, 490};
        langau_function_range_6 =  {128, 1200};
        landau_shift_range_3p5 = {10,20};
        landau_shift_range_9 = {1,10};
        landau_shift_range_20 = {2, 7};
        landau_shift_range_15 = {2.5, 10};
        landau_shift_range_15_redo = {-0.009, 0.009};
        landau_shift_range_6 = {2, 20};
    }
    else{
        fname_3p5 = "./Test_crs/F2--large-end--00000.txt";
        fname_9 = "./Test_crs/Middle_cosmics/F2--mid--00000.txt";
        fname_20 = "./Test_crs/small/F2--small--00000.txt";
        fname_15 = "./Test_crs/15cm/F2--15cm--00000.txt";
        fname_15_redo = "./Test_crs/15cm_again/F2--15cm_redo--00000.txt";
        fname_6 = "./Test_crs/6cm/F2--6cm--00000.txt";
        landau_params_3p5 = {200, 1174, 30000, 200}; 
        landau_params_9 = {450, 890, 25000, 450}; 
        landau_params_20 = {450, 890, 25000, 450};
        landau_params_15 = {100, 890, 15000, 100};
        landau_params_15_redo = {50, 803, 25000, 100};
        landau_params_6 = {64, 1100, 54000, 150};
        langau_function_range_3p5 = {0, 2000};
        langau_function_range_9 = {400, 2500};
        langau_function_range_20 =  {400, 1700};
        langau_function_range_15 =  {100, 1700};
        langau_function_range_15_redo =  {500, 2700};
        langau_function_range_6 =  {500, 3000};
        landau_shift_range_3p5 = {-350,0};
        landau_shift_range_9 = {-275,-150};
        landau_shift_range_20 = {-250,-220};
        landau_shift_range_15 = {-730,-325};
        landau_shift_range_15_redo = {-820,-650};
        landau_shift_range_6 = {-625, -525};
    }
    // First line x_values = distance from SiPM
    double x_hyp = 25.25;
    double x_adj = 25.15;
    double x_ruler_err = 0.1;
    double angle_diff = 2 * TMath::ACos(x_adj/x_hyp);
    double cos_angle = TMath::Cos(angle_diff);
    std::vector<double> x_values = {3.5,9,20,15,15,6};
    std::vector<double> x_err;
    for(int i = 0; i < x_values.size(); i++){
        double x_measured = x_values[i];
        double x_real = cos_angle * x_values[i];
        double x_angle_err = abs(x_real - x_measured);
        x_err.push_back(pow(x_angle_err * x_angle_err + x_ruler_err * x_ruler_err,0.5));
    }

    // Gain = 2, 3.5cm, charge, landau
    std::string plot_name_3p5 = "3.5";
    std::string histo_title_3p5 = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[0]);
    landau_fit_result fit_vals_3p5 = ana_landau(fname_3p5, convert_to_phe, landau_params_3p5, langau_function_range_3p5, landau_shift_range_3p5, plot_name_3p5, histo_title_3p5, true);

    // Gain = 2, 9cm, charge, landau
    std::string plot_name_9 = "9";
    std::string histo_title_9 = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[1]);
    landau_fit_result fit_vals_9 = ana_landau(fname_9, convert_to_phe, landau_params_9, langau_function_range_9, landau_shift_range_9, plot_name_9, histo_title_9, true);

    // Gain = 2, 20cm, charge, landau
    std::string plot_name_20 = "20";
    std::string histo_title_20 = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[2]);
    landau_fit_result fit_vals_20 = ana_landau(fname_20, convert_to_phe, landau_params_20, langau_function_range_20, landau_shift_range_20, plot_name_20, histo_title_20, true);
    
    // Gain = 2, 15cm, charge, landau
    std::string plot_name_15 = "15";
    std::string histo_title_15 = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[3]);
    landau_fit_result fit_vals_15 = ana_landau(fname_15, convert_to_phe, landau_params_15, langau_function_range_15, landau_shift_range_15, plot_name_15, histo_title_15, true, 4);

    // Gain = 2, 15cm REDO, charge, landau
    std::string plot_name_15_redo = "15_Redo";
    std::string histo_title_15_redo = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[3]);
    landau_fit_result fit_vals_15_redo = ana_landau(fname_15_redo, convert_to_phe, landau_params_15_redo, langau_function_range_15_redo, landau_shift_range_15_redo, plot_name_15_redo, histo_title_15_redo, true);

    // Gain = 2, 6cm, charge, landau
    std::string plot_name_6 = "6";
    std::string histo_title_6 = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[4]);
    landau_fit_result fit_vals_6 = ana_landau(fname_6, convert_to_phe, landau_params_6, langau_function_range_6, landau_shift_range_6, plot_name_6, histo_title_6, true);


    // Converts from Phe to Phe/MeV
    double density = 7.12; // g/cm3
    double energy_dep_per_gram_per_cm2 = 2; // MeV * cm2 / g
    double energy_per_cm = density * energy_dep_per_gram_per_cm2;

    // Second line penetration_depths is the thickness of scintillator
    std::vector<landau_fit_result> landau_fit_result_vector = {fit_vals_3p5, fit_vals_9, fit_vals_20, fit_vals_15, fit_vals_15_redo, fit_vals_6}; // landau fit result information
    std::vector<double> penetration_depths = {5.4, 4.2, 2.5, 3.5, 3.5, 4.5};


    plot_MPVs(energy_per_cm, x_values, x_err, landau_fit_result_vector, penetration_depths, use_amplitude); // E/cm, x_arr, y_arr, depth_arr
    
}
