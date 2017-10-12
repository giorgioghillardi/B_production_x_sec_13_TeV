#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda
// channel = 7: Bc -> J/psi Pi
//-----------------------------------------------------------------

void plot_syst(TString measure, TString syst, int channel, int n_var1_bins, TString var2_name, double var2_min, double var2_max, TString x_axis_name, TString b_title, double* var1_bin_centre, double* var1_bin_centre_lo, double* var1_bin_centre_hi, double* syst_array);

//input example: plot_syst --measure x_sec --channel 1 --bins pt_y --syst signal_pdf
int main(int argc, char** argv)
{
  TString measure = "";
  int channel = 1;
  std::string bins = "pt";
  std::string syst = "";
  
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
      if(argument == "--syst")
	{
	  convert << argv[++i];
	  convert >> syst;
	}
    }
  
  if(syst == "")
    std::cout << "ERROR: please choose the syst error to plot with the --syst option" << std::endl;
  
  //to create the directories to save the .png files
  std::vector<std::string> dir_list;  
  
  dir_list.push_back(static_cast<const char*>(TString::Format(VERSION) + "/systematics"));
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
  
  double syst_array[n_var2_bins][n_var1_bins];
  double syst_err_lo[n_var2_bins][n_var1_bins];
  double syst_err_hi[n_var2_bins][n_var1_bins];

  TString x_axis_name = "";
  if(var1_name =="pt")
    x_axis_name = "p_{T}(B) [GeV]";
  else
    x_axis_name = "|y|(B)";
 
  TString b_title= "";
  switch(channel)
    {
    case 1:
      b_title = "B+";
      break;
    case 2:
      b_title = "B0";
      break;
    case 4:
      b_title = "Bs";
      break;
    case 7:
      b_title = "Bc";
      break;
    default:
      break;
    }
  
  //calculate the centre of var1 bins
  for(int i=0; i<n_var1_bins; i++)
    {
      var1_bin_centre[i] = var1_bins[i] + (var1_bins[i+1] - var1_bins[i])/2;
      var1_bin_centre_lo[i] = var1_bin_centre[i] - var1_bins[i];
      var1_bin_centre_hi[i] = var1_bins[i+1] - var1_bin_centre[i];
    }
  
  //read syst vector
  read_vector(channel, syst, var1_name , var2_name, n_var1_bins, n_var2_bins, var1_bins, var2_bins, syst_array[0], syst_err_lo[0], syst_err_hi[0]);
  
  //to plot the systematic errors
  for(int j=0; j<n_var2_bins; j++)
    {
      plot_syst(measure, syst, channel, n_var1_bins, var2_name, var2_bins[j], var2_bins[j+1], x_axis_name, b_title, var1_bin_centre, var1_bin_centre_lo, var1_bin_centre_hi, syst_err_hi[j]);
    }
  
}//end

void plot_syst(TString measure, TString syst, int channel, int n_var1_bins, TString var2_name, double var2_min, double var2_max, TString x_axis_name, TString b_title, double* var1_bin_centre, double* var1_bin_centre_lo, double* var1_bin_centre_hi, double* syst_array)
{
  TCanvas ce;
  TGraphAsymmErrors* graph_syst = new TGraphAsymmErrors(n_var1_bins, var1_bin_centre, syst_array, var1_bin_centre_lo, var1_bin_centre_hi);

  TString syst_title = measure + " " + syst + " " + "systematic error" ;

  graph_syst->SetTitle(syst_title);
  graph_syst->GetXaxis()->SetTitle(x_axis_name);
  graph_syst->Draw("AP");

  if(measure != "x_sec" && measure != "ratio")
    return;

  TString save_syst = TString::Format(VERSION) + "/systematics/" + measure + "_" + syst + "_" + channel_to_ntuple_name(channel) + "_" + var2_name + TString::Format("_from_%.2f_to_%.2f", var2_min, var2_max) + ".png";

  if(n_var1_bins == 1)
    save_syst = TString::Format(VERSION) + "/systematics/" + measure + "_" + syst + "_" + channel_to_ntuple_name(channel) + "_full_bins.png";

  ce.SaveAs(save_syst);
}
