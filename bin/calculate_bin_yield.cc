#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda
// channel = 7: Bc -> J/psi Pi+
//-----------------------------------------------------------------

//input example: calculate_bin_yield --channel 1 --ptmin 10 --ptmax 20 --ymin 0.00 --ymax 0.50
int main(int argc, char** argv)
{
  int channel = 1;
  double pt_min = -1;
  double pt_max = -1;
  double y_min = -1;
  double y_max = -1;

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--channel")
	{
	  convert << argv[++i];
	  convert >> channel;
	}
      if(argument == "--ptmin")
	{
	  convert << argv[++i];
	  convert >> pt_min;
	}
      if(argument == "--ptmax")
	{
	  convert << argv[++i];
	  convert >> pt_max;
	}
      if(argument == "--ymin")
	{
	  convert << argv[++i];
	  convert >> y_min;
	}
      if(argument == "--ymax")
	{
	  convert << argv[++i];
	  convert >> y_max;
	}
    }
  
  if(pt_min == -1 || pt_max == -1 ||y_min == -1 ||y_max == -1)
    {
      std::cout << "Error: The bin was not well defined. Please enter pt and y bin." << std::endl;
      return 0;
    }
  
  //to create the directories to save the files
  std::vector<std::string> dir_list;  
  
  dir_list.push_back(static_cast<const char*>(TString::Format(VERSION) + "/signal_yield_root/" + channel_to_ntuple_name(channel)));
  dir_list.push_back(static_cast<const char*>(TString::Format(VERSION) + "/mass_fits/" + channel_to_ntuple_name(channel)));
  create_dir(dir_list);
  TString data_selection_input_file ="/lstore/cms/balves/Jobs/Full_Dataset_2016/2016_data_ntkp_with_cuts.root";
  //TString data_selection_input_file = TString::Format(BASE_DIR) + "selected_myloop_new_data_" + channel_to_ntuple_name(channel) + "_with_cuts.root";
  RooWorkspace* ws = new RooWorkspace("ws","Bmass");
  RooRealVar* signal_res; 
  
  //set up mass, pt and y variables inside ws  
  set_up_workspace_variables(*ws,channel);
  //read data from the selected data file, and import it as a dataset into the workspace.
  read_data(*ws, data_selection_input_file,channel);
  
  /////////////////////////////////////////////////////
  //calculate the signal yield for a bin of pt and y.//
  /////////////////////////////////////////////////////
  std::cout << "processing subsample: " << pt_min << " < " << "pt" << " < " << pt_max << " and " << y_min << " < " << "|y|" << " < " << y_max << std::endl;
  
  signal_res = bin_mass_fit(*ws, channel, pt_min, pt_max, y_min, y_max);

  //write signal_yield and statistical error to file
  TString yield_file_name = TString::Format(VERSION) + "/signal_yield_root/" + channel_to_ntuple_name(channel) + "/yield_" + channel_to_ntuple_name(channel) + "_pt_from_" + TString::Format("%d_to_%d", (int)pt_min, (int)pt_max) + "_y_from_" + TString::Format("%.2f_to_%.2f", y_min, y_max) + ".root";
  
  TFile* yield_file = new TFile(yield_file_name,"recreate");
  
  TVectorD signal_yield(1);
  TVectorD stat_err_lo(1);
  TVectorD stat_err_hi(1);

  signal_yield[0] = fabs(signal_res->getVal());
  stat_err_lo[0] = fabs(signal_res->getAsymErrorLo());
  stat_err_hi[0] = fabs(signal_res->getAsymErrorHi());
  
  signal_yield.Write("val");
  stat_err_lo.Write("err_lo");
  stat_err_hi.Write("err_hi");

  yield_file->Close();
  delete yield_file; 
}//end of calculate_bin_yield
