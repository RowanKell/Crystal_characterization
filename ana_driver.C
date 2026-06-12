#include "ana_landau.C"
#include "ana_fit.C"
#include "ana_plot.C"

void ana_driver(bool use_amplitude = true){
    gStyle->SetOptFit(1111); // set option to display fit parameters on plot

    
    // Gain = 2, amplitude, gauss
    std::string fname_g2a = "./Test_crs/F1--phe-g2--00000.txt";
    int n_points_gauss_g2a = 3;
    int n_points_linear_fit_g2a = 3;
    std::vector<double> gaussian_fit_params_g2a = {1.53868e+02,1.08227e-03,8.59892e-05,6.69506e+01,1.38493e-03,1.06635e-04,8.79378e+00,1.70000e-03,2.23329e-04};
    double gaussian_function_range_g2a[2] = {0.00091, 0.002}; // ????
    std::array<double, 2> fit_vals_g2a = ana_fit(fname_g2a, gaussian_fit_params_g2a, n_points_gauss_g2a, n_points_linear_fit_g2a, gaussian_function_range_g2a, false);

    double convert_to_phe;

    if(use_amplitude){
        convert_to_phe = 1/fit_vals_g2a[1];
    }
    else{
        // Gain = 40, amplitude, gauss
        std::string fname_g40a = "./Test_crs/F1--phe-g40--00000.txt";
        int n_points_gauss_g40a = 3;
        int n_points_linear_fit_g40a = 2;
        std::vector<double> gaussian_fit_params_g40a = {7.14115e+02, 9.84592e-03, 6.84463e-04, 7.95184e+01, 1.83608e-02, 8.33767e-04, 7.63024e+01, 1.23959e-02, 4.94058e-03}; 
        double gaussian_function_range_g40a[2] = {0, 0.03};
        std::array<double, 2> fit_vals_g40a = ana_fit(fname_g40a, gaussian_fit_params_g40a, n_points_gauss_g40a, n_points_linear_fit_g40a, gaussian_function_range_g40a, false);
        std::cout << "gain = 40, linear fit vals: (" << fit_vals_g40a[0] << ", " << fit_vals_g40a[1] << std::endl;
        
        // Gain = 40, charge, gauss
        std::string fname_g40c = "./Test_crs/F2--phe-g40--00000.txt";
        int n_points_gauss_g40c = 3;
        int n_points_linear_fit_g40c = 3;
        std::vector<double> gaussian_fit_params_g40c = {4.712e2, -4.099e-12, 2.493e-10, 1.953e2, -1.041e-9, 2.886e-10, 4.7721e1, -2.1222e-9, 3.9309e-10};
        double gaussian_function_range_g40c[2] = {-3e-9, 1e-9};
        std::array<double, 2> fit_vals_g40c = ana_fit(fname_g40c, gaussian_fit_params_g40c, n_points_gauss_g40c, n_points_linear_fit_g40c, gaussian_function_range_g40c, false);
        
        // Converts from Wb to Phe
        convert_to_phe = 1/(fit_vals_g40c[1]*fit_vals_g2a[1]/fit_vals_g40a[1]);
        // double Q_to_phe = -1/3e-11;
    }
    std::string fname_large;
    std::string fname_middle;
    std::string fname_small;
    std::string fname_15;
    std::string fname_6;
    std::array<double,4> landau_params_large;
    std::array<double,4> landau_params_middle;
    std::array<double,4> landau_params_small;
    std::array<double,4> landau_params_15;
    std::array<double,4> landau_params_6;
    std::array<double,2> langau_function_range_large;
    std::array<double,2> langau_function_range_middle;
    std::array<double,2> langau_function_range_small;
    std::array<double,2> langau_function_range_15;
    std::array<double,2> langau_function_range_6;
    std::array<double,2> landau_shift_range_large;
    std::array<double,2> landau_shift_range_middle;
    std::array<double,2> landau_shift_range_small;
    std::array<double,2> landau_shift_range_15;
    std::array<double,2> landau_shift_range_6;
    if(use_amplitude){
        fname_large = "./Test_crs/F1--large-end--00000.txt";
        fname_middle = "./Test_crs/Middle_cosmics/F1--mid--00000.txt";
        fname_small = "./Test_crs/small/F1--small--00000.txt";
        fname_15 = "./Test_crs/15cm/F1--15cm--00000.txt";
        fname_6 = "./Test_crs/6cm/F1--6cm--00000.txt";
        landau_params_large = {30, 400, 5000, 30}; 
        landau_params_middle = {50, 350, 5000, 50}; 
        landau_params_small = {20, 225, 3000, 20}; 
        landau_params_15 = {17, 212, 3000, 47};
        landau_params_6 = {17, 212, 3000, 47};
        langau_function_range_large = {250, 800};
        langau_function_range_middle = {200, 700};
        langau_function_range_small =  {115, 500};
        langau_function_range_15 =  {128, 490};
        langau_function_range_6 =  {128, 1200};
        landau_shift_range_large = {10,20};
        landau_shift_range_middle = {1,10};
        landau_shift_range_small = {2, 7};
        landau_shift_range_15 = {2.5, 10};
        landau_shift_range_6 = {2, 20};
    }
    else{
        fname_large = "./Test_crs/F2--large-end--00000.txt";
        fname_middle = "./Test_crs/Middle_cosmics/F2--mid--00000.txt";
        fname_small = "./Test_crs/small/F2--small--00000.txt";
        fname_15 = "./Test_crs/15cm/F2--15cm--00000.txt";
        fname_6 = "./Test_crs/6cm/F2--6cm--00000.txt";
        landau_params_large = {200, 1174, 30000, 200}; 
        landau_params_middle = {450, 890, 25000, 450}; 
        landau_params_small = {450, 890, 25000, 450};
        landau_params_15 = {100, 890, 15000, 100};
        landau_params_6 = {64, 1100, 54000, 150};
        langau_function_range_large = {0, 2000};
        langau_function_range_middle = {400, 2500};
        langau_function_range_small =  {400, 1700};
        langau_function_range_15 =  {100, 1700};
        langau_function_range_6 =  {500, 3000};
        landau_shift_range_large = {-350,0};
        landau_shift_range_middle = {-275,-150};
        landau_shift_range_small = {-250,-220};
        landau_shift_range_15 = {-730,-325};
        landau_shift_range_6 = {-625, -525};
    }
    // First line x_values = distance from SiPM
    std::vector<double> x_values = {3.5,9,20,15, 6};

    // Gain = 2, large, charge, landau
    std::string plot_name_large = "Large";
    std::string histo_title_large = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[0]);
    landau_fit_result fit_vals_large = ana_landau(fname_large, convert_to_phe, landau_params_large, langau_function_range_large, landau_shift_range_large, plot_name_large, histo_title_large, true);

    // Gain = 2, middle, charge, landau
    std::string plot_name_middle = "Middle";
    std::string histo_title_middle = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[1]);
    landau_fit_result fit_vals_middle = ana_landau(fname_middle, convert_to_phe, landau_params_middle, langau_function_range_middle, landau_shift_range_middle, plot_name_middle, histo_title_middle, true);

    // Gain = 2, small, charge, landau
    std::string plot_name_small = "Small";
    std::string histo_title_small = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[2]);
    landau_fit_result fit_vals_small = ana_landau(fname_small, convert_to_phe, landau_params_small, langau_function_range_small, landau_shift_range_small, plot_name_small, histo_title_small, true);
    
    // Gain = 2, 15cm, charge, landau
    std::string plot_name_15 = "15";
    std::string histo_title_15 = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[3]);
    landau_fit_result fit_vals_15 = ana_landau(fname_15, convert_to_phe, landau_params_15, langau_function_range_15, landau_shift_range_15, plot_name_15, histo_title_15, true, 4);

    // Gain = 2, 6cm, charge, landau
    std::string plot_name_6 = "6";
    std::string histo_title_6 = Form("Number of Phe at X = %f cm;# Photo Electrons;Counts", x_values[4]);
    landau_fit_result fit_vals_6 = ana_landau(fname_6, convert_to_phe, landau_params_6, langau_function_range_6, landau_shift_range_6, plot_name_6, histo_title_6, true);


    // Converts from Phe to Phe/MeV
    double density = 7.12; // g/cm3
    double energy_dep_per_gram_per_cm2 = 2; // MeV * cm2 / g
    double energy_per_cm = density * energy_dep_per_gram_per_cm2;

    // Second line penetration_depths is the thickness of scintillator
    std::vector<landau_fit_result> landau_fit_result_vector = {fit_vals_large,fit_vals_middle, fit_vals_small, fit_vals_15, fit_vals_6}; // landau fit result information
    std::vector<double> penetration_depths = {5.4, 4.2, 2.5, 3.5, 4.5};


    plot_MPVs(energy_per_cm, x_values, landau_fit_result_vector, penetration_depths, use_amplitude); // E/cm, x_arr, y_arr, depth_arr
    
}
