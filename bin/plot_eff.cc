#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda
//-----------------------------------------------------------------

//input example: plot_eff --measure x_sec --channel 1 --bins pt_y --eff preeff
int main(int argc, char** argv)
{
  TString measure = "";
  int channel = 1;
  std::string bins = "pt";
  std::string eff = "";
  
  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--measure")
	{
	  convert << argv[++i];
	  convert >> measure;
	}
      if(argument == "--channel")
	{
	  convert << argv[++i];
	  convert >> channel;
	}
      if(argument == "--bins")
	{
	  convert << argv[++i];
	  convert >> bins;
	}
      if(argument == "--eff")
	{
	  convert << argv[++i];
	  convert >> eff;
	}
    }
  
  if(eff == "")
    std::cout << "ERROR: please choose the efficiency to plot with the --eff option" << std::endl;
  
  //to create the directories to save the .png files
  std::vector<std::string> dir_list;  
  
  dir_list.push_back(static_cast<const char*>(TString::Format(VERSION) +"/efficiencies"));
  create_dir(dir_list);
  
  //set up the vectors
  TString var1_name = "";
  TString var2_name = "";

  int n_var1_bins=1;
  int n_var2_bins=1;

  double* var1_bins = NULL;
  double* var2_bins = NULL;

  setup_bins(measure, channel, bins, &var1_name, &n_var1_bins, &var2_name, &n_var2_bins, &var1_bins, &var2_bins);
  
  //initialize arrays for the efficiencies
  double var1_bin_centre[n_var1_bins];
  double var1_bin_centre_lo[n_var1_bins];
  double var1_bin_centre_hi[n_var1_bins];
  
  double eff_array[n_var2_bins][n_var1_bins];
  double eff_err_lo[n_var2_bins][n_var1_bins];
  double eff_err_hi[n_var2_bins][n_var1_bins];

  TString x_axis_name = "";
  if(var1_name =="pt")
    x_axis_name = "p_{T}(B) [GeV]";
  else
    x_axis_name = "|y|(B)";
 
  TString b_title= "";
  switch(channel)
    {
    default:
    case 1:
      b_title = "B+";
      break;
    case 2:
      b_title = "B0";
      break;
    case 4:
      b_title = "Bs";
      break;
    }
  
  //calculate the mean values and bin size of var1 and var2
  for(int i=0; i<n_var1_bins; i++)
    {
      var1_bin_centre[i] = var1_bins[i] + (var1_bins[i+1] - var1_bins[i])/2;
      var1_bin_centre_lo[i] = var1_bin_centre[i] - var1_bins[i];
      var1_bin_centre_hi[i] = var1_bins[i+1] - var1_bin_centre[i];
    }
  
  //read efficiency
  read_vector(measure, channel, eff, var1_name , var2_name, n_var1_bins, n_var2_bins, var1_bins, var2_bins, eff_array[0], eff_err_lo[0], eff_err_hi[0]);
  
  //to calculate cross section
  for(int j=0; j<n_var2_bins; j++)
    {
      plot_eff(measure, eff, channel, n_var1_bins, var2_name, var2_bins[j], var2_bins[j+1], x_axis_name, b_title, var1_bin_centre, var1_bin_centre_lo, var1_bin_centre_hi, eff_array[j], eff_err_lo[j], eff_err_hi[j]);
    }
  
}//end of plot_eff
