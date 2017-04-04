#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/bins.h"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda
//-----------------------------------------------------------------

//input example: calculate_yields --channel 1 --yield 1 --preeff 0 --recoeff 0 --totaleff 0 --bins pt_y
int main(int argc, char** argv)
{
  int channel = 1;
  TString bins = "";
  int yield = 1;
  int preeff = 0;
  int recoeff = 0;
  int totaleff = 0;

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
      if(argument == "--yield")
	{
	  convert << argv[++i];
	  convert >> yield;
	}
      if(argument == "--preeff")
	{
	  convert << argv[++i];
	  convert >> preeff;
	}
      if(argument == "--recoeff")
	{
	  convert << argv[++i];
	  convert >> recoeff;
	}
      if(argument == "--totaleff")
	{
	  convert << argv[++i];
	  convert >> totaleff;
	}
    }

  if(bins == "")
    {
      std::cout << "Error: Enter the --bins option." << std::endl;
      return 0;
    }
  
  //to create the directories to save the files
  std::vector<std::string> dir_list;  
  
  dir_list.push_back(static_cast<const char*>("signal_yield_root/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION)));
  dir_list.push_back(static_cast<const char*>("mass_fits/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION)));
  create_dir(dir_list);

  //choose the bins to read

  int n_pt_bins=1;
  double* pt_bins=total_pt_bin_edges;

  int n_y_bins=1;
  double* y_bins=total_y_bin_edges;

  switch (channel)
    {
    default:
    case 1:
      pt_bins = ntkp_pt_bins;
      n_pt_bins = (sizeof(ntkp_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntkp_y_bins;
      n_y_bins = (sizeof(ntkp_y_bins) / sizeof(double)) - 1 ;
      break;

    case 2:
      pt_bins = ntkstar_pt_bins;
      n_pt_bins = (sizeof(ntkstar_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntkstar_y_bins;
      n_y_bins = (sizeof(ntkstar_y_bins) / sizeof(double)) - 1 ;
      break;

    case 3:
      pt_bins = ntks_pt_bins;
      n_pt_bins = (sizeof(ntks_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntks_y_bins;
      n_y_bins = (sizeof(ntks_y_bins) / sizeof(double)) - 1 ;
      break;

    case 4:
      pt_bins = ntphi_pt_bins;
      n_pt_bins = (sizeof(ntphi_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntphi_y_bins;
      n_y_bins = (sizeof(ntphi_y_bins) / sizeof(double)) - 1 ;
      break;

    case 5:
      pt_bins = ntmix_pt_bins;
      n_pt_bins = (sizeof(ntmix_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntmix_y_bins;
      n_y_bins = (sizeof(ntmix_y_bins) / sizeof(double)) - 1 ;
      break;

    case 6:
      pt_bins = ntlambda_pt_bins;
      n_pt_bins = (sizeof(ntlambda_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntlambda_y_bins;
      n_y_bins = (sizeof(ntlambda_y_bins) / sizeof(double)) - 1 ;
      break;
    }
  
  if(bins == "full")
    {
      n_pt_bins= 1;
      pt_bins = total_pt_bin_edges;
      
      n_y_bins= 1;
      y_bins = total_y_bin_edges;
    }
  else
    if(bins == "pt")
      {
	n_y_bins= 1;
	y_bins = total_y_bin_edges;
      }
    else
      if(bins == "y")
	{
	  n_pt_bins= 1;
	  pt_bins = total_pt_bin_edges;
	}
      else
	if(bins == "pt_y" || bins == "y_pt")
	  {
	    //the arrays are already in place in the switch above
	  }
	else
	  {
	    std::cout << "Error: Enter a valid --bins option. pt, y, pt_y, y_pt are implemented." << std::endl;
	    return 0;
	  }
	  
  TString line = "";
  TString command = "";
  TString opt = "";
  double pt_min;
  double pt_max;
  double y_min;
  double y_max;
  
  for(int j=0; j<n_y_bins; j++)
    {
      for(int i=0; i<n_pt_bins; i++)
        {
	  pt_min = pt_bins[i];
	  pt_max = pt_bins[i+1];
	  y_min = y_bins[j];
	  y_max = y_bins[j+1];
	  
	  opt = " --channel " + TString::Format("%d", channel) + " --ptmin " + TString::Format("%d", (int)pt_min) + " --ptmax " + TString::Format("%d", (int)pt_max) + " --ymin " + TString::Format("%.2f", y_min) + " --ymax " + TString::Format("%.2f", y_max); 
	  
	  if(yield)
	    {
	      command = "calculate_bin_yield";
	      line = command + opt;
	      gSystem->Exec(line);
	    }
	  if(preeff)
	    {
	      command = "calculate_bin_eff";
	      line = command + " --eff pre-filter" + opt;
	      gSystem->Exec(line);
	    }
	  if(recoeff)
	    {
	      command = "calculate_bin_eff";
	      line = command + " --eff reco" + opt;
	      gSystem->Exec(line);
	    }
	  if(totaleff)
	    {
	      command = "calculate_bin_eff";
	      line = command + " --eff total" + opt;
	      gSystem->Exec(line);
	    }
	}
    }
}//end of calculate_yields
