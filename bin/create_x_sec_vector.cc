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

//input example: create_x_sec --channel 1 --bins pt_y
int main(int argc, char** argv)
{
  int channel = 1;
  std::string bins = "pt";

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

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
    }
  
  //create directories
  std::vector<std::string> dir_list;
  dir_list.push_back("x_sec/");
  create_dir(dir_list);
  
  //set up the vectors
  TString var1_name = "";
  TString var2_name = "";

  int n_var1_bins=1;
  int n_var2_bins=1;

  double* var1_bins = NULL;
  double* var2_bins = NULL;
  
  TString measure = "x_sec";
  
  setup_bins(measure, channel, bins, &var1_name, &n_var1_bins, &var2_name, &n_var2_bins, &var1_bins, &var2_bins);
  
  //initialize arrays for yield, efficiencies, etc 
  double var1_bin_size[n_var1_bins];
  double var2_bin_size[n_var2_bins];

  double yield_array[n_var2_bins][n_var1_bins];
  double yield_err_lo[n_var2_bins][n_var1_bins];
  double yield_err_hi[n_var2_bins][n_var1_bins];
  
  double total_eff[n_var2_bins][n_var1_bins];
  double total_eff_err_lo[n_var2_bins][n_var1_bins];
  double total_eff_err_hi[n_var2_bins][n_var1_bins];

  double x_sec[n_var2_bins][n_var1_bins];
  double x_sec_err_lo[n_var2_bins][n_var1_bins];
  double x_sec_err_hi[n_var2_bins][n_var1_bins];
    
  RooRealVar* branch = branching_fraction(channel);

  double bin_size;
  
  if(bins=="full")
    {
      var1_bin_size[0] = 1; //force the bin size equal to one when we want to calculate full dataset cross section.
      var2_bin_size[0] = 1;
    }
  else
    if(bins=="pt" || bins=="y")
      {
	var2_bin_size[0]=1;
      }
  
  for(int j=0; j<n_var2_bins; j++)
    {
      var2_bin_size[j] = var2_bins[j+1] - var2_bins[j];
    }

  for(int i=0; i<n_var1_bins; i++)
    {
      var1_bin_size[i] = var1_bins[i+1] - var1_bins[i];
    }

  //read yields
  read_vector(measure, channel, "yield", var1_name , var2_name, n_var1_bins, n_var2_bins, var1_bins, var2_bins, yield_array[0], yield_err_lo[0], yield_err_hi[0]);

  //read efficiency
  read_vector(measure, channel, "totaleff", var1_name , var2_name, n_var1_bins, n_var2_bins, var1_bins, var2_bins, total_eff[0], total_eff_err_lo[0], total_eff_err_hi[0]);
  
  //to calculate cross section
  for(int j=0; j<n_var2_bins; j++)
    {
      for(int i=0; i<n_var1_bins; i++)
	{
	  bin_size = var1_bin_size[i] * var2_bin_size[j];
	 
	  x_sec[j][i] = yield_array[j][i] / (2  * total_eff[j][i] * LUMINOSITY * bin_size * branch->getVal());
	  x_sec_err_lo[j][i] = x_sec[j][i] * (yield_err_lo[j][i] / yield_array[j][i]);
	  x_sec_err_hi[j][i] = x_sec[j][i] * (yield_err_hi[j][i] / yield_array[j][i]);
	}
    }
  
  //save x-sec array to file
  for(int j=0; j<n_var2_bins; j++)
    {
      TVectorD out_val(n_var1_bins);
      TVectorD out_err_lo(n_var1_bins);
      TVectorD out_err_hi(n_var1_bins);
      
      for(int i=0; i<n_var1_bins; i++)
	{
	  out_val[i] = x_sec[j][i];
	  out_err_lo[i] = x_sec_err_lo[j][i];
	  out_err_hi[i] = x_sec_err_hi[j][i];
	}
	  
      TString bins_str = "";
      TString out_file_name = "";

      if(var1_name == "pt")
	bins_str = TString::Format("%.2f_to_%.2f", var2_bins[j], var2_bins[j+1]);
      else
	bins_str = TString::Format("%d_to_%d", (int)var2_bins[j], (int)var2_bins[j+1]);

      out_file_name = "x_sec/x_sec_vector_" + channel_to_ntuple_name(channel) + "_" + var1_name + "_bins_" + var2_name + "_from_" + bins_str + "_" + TString::Format(VERSION) + ".root";

      if(bins == "full")
	out_file_name = "x_sec/x_sec_vector_" + channel_to_ntuple_name(channel) + "_full_bins_" + TString::Format(VERSION) + ".root";
	  
      TFile* fout = new TFile(out_file_name,"recreate");

      out_val.Write("val");
      out_err_lo.Write("err_lo");
      out_err_hi.Write("err_hi");

      fout->Close();
      delete fout;
    }
    
  //////////////////////////////////////////////////////////////////////////////////////////////
  //To show the values of cross section or signal yield and the errors at the end, like a table/
  //////////////////////////////////////////////////////////////////////////////////////////////

  //signal yield
  print_table("YIELDS", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bins, var2_bins, yield_array[0], yield_err_lo[0], yield_err_hi[0]);
  
  //total eff
  print_table("OVERALL EFFICIENCY", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bins, var2_bins, total_eff[0], total_eff_err_lo[0], total_eff_err_hi[0]);

  //cross section
  print_table("CROSS SECTION", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bins, var2_bins, x_sec[0], x_sec_err_lo[0], x_sec_err_hi[0]);

  //branching fraction
  std::cout << "branching fraction: " << branch->getVal() << std::endl;
  std::cout << std::endl;
  
}//end of create_x_sec_vector
