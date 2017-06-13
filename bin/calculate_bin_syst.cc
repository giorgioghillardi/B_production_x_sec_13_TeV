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

double pdf_syst(RooWorkspace& ws, int channel, double pt_min, double pt_max, double y_min, double y_max, double nominal_yield, TString syst);
double mass_window_syst(RooWorkspace& ws, int channel, double pt_min, double pt_max, double y_min, double y_max, double nominal_yield, TString input_file);

//input example: calculate_bin_syst --channel 1 --syst signal_pdf --ptmin 30 --ptmax 35 --ymin 0.00 --ymax 2.25
int main(int argc, char** argv)
{
  int channel = 1;
  TString syst = "";
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
      if(argument == "--syst")
	{
	  convert << argv[++i];
	  convert >> syst;
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
  
  dir_list.push_back(static_cast<const char*>("signal_yield_root/syst/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION)));
  dir_list.push_back(static_cast<const char*>("mass_fits/syst/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION)));

  create_dir(dir_list);

  TString data_selection_input_file = TString::Format(BASE_DIR) + "selected_myloop_new_data_" + channel_to_ntuple_name(channel) + "_with_cuts.root";
  RooWorkspace* ws = new RooWorkspace("ws","Bmass");
  
  //set up mass, pt and y variables inside ws  
  set_up_workspace_variables(*ws,channel);
  //read data from the selected data file, and import it as a dataset into the workspace.
  read_data(*ws, data_selection_input_file,channel);
  
  ///////////////////////////////////////////////////
  //calculate the syst error for a bin of pt and y.//
  ///////////////////////////////////////////////////
  std::cout << "processing subsample: " << pt_min << " < " << "pt" << " < " << pt_max << " and " << y_min << " < " << "|y|" << " < " << y_max << std::endl;
  
  //read nominal signal yield
  TString in_file_name = "signal_yield_root/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + "/yield_" + channel_to_ntuple_name(channel) + "_pt_from_" + TString::Format("%d_to_%d", (int)pt_min, (int)pt_max) + "_y_from_" + TString::Format("%.2f_to_%.2f", y_min, y_max) + "_" + TString::Format(VERSION) + ".root";

  std::cout << "read :" << in_file_name << std::endl;

  TFile* fin = new TFile(in_file_name);
  TVectorD *in_val = (TVectorD*)fin->Get("val");
  TVectorD *in_err_lo = (TVectorD*)fin->Get("err_lo");
  TVectorD *in_err_hi = (TVectorD*)fin->Get("err_hi");
  delete fin;

  RooRealVar nominal_yield("nominal_yield", "nominal_yield", in_val[0][0]);
  nominal_yield.setAsymError(-in_err_lo[0][0],in_err_hi[0][0]);

  //debug:
  std::cout << "nominal_yield: " << nominal_yield.getVal() << " err_lo: " << nominal_yield.getAsymErrorLo() << " err_hi: " << nominal_yield.getAsymErrorHi() << std::endl;

  //calculate syst yield
  double signal_res = 0;
  
  if(syst == "mass_window")
    signal_res = mass_window_syst(*ws, channel, pt_min, pt_max, y_min, y_max, nominal_yield.getVal(), data_selection_input_file);
  else
    if(syst == "signal_pdf" || syst == "cb_pdf")
      signal_res = pdf_syst(*ws, channel, pt_min, pt_max, y_min, y_max, nominal_yield.getVal(), syst);
    
  //calculate the diff in absolute value, i.e. from 0 to 1
  double absolute_syst_val;
  absolute_syst_val = fabs(nominal_yield.getVal() - signal_res)/nominal_yield.getVal();

  //debug:
  std::cout << "absolute_syst_val: " << absolute_syst_val << std::endl;

  //write to file with syst name
  TString yield_file_name = "signal_yield_root/syst/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + "/" + syst + "_" + channel_to_ntuple_name(channel) + "_pt_from_" + TString::Format("%d_to_%d", (int)pt_min, (int)pt_max) + "_y_from_" + TString::Format("%.2f_to_%.2f", y_min, y_max) + "_" + TString::Format(VERSION) + ".root";
  
  TFile* yield_file = new TFile(yield_file_name,"recreate");
  
  TVectorD val(1);
  TVectorD err_lo(1);
  TVectorD err_hi(1);

  val[0] = 1.00;
  err_lo[0] = fabs(absolute_syst_val);
  err_hi[0] = fabs(absolute_syst_val);
  
  val.Write("val");
  err_lo.Write("err_lo");
  err_hi.Write("err_hi");

  yield_file->Close();
  delete yield_file; 
}//end

double pdf_syst(RooWorkspace& ws, int channel, double pt_min, double pt_max, double y_min, double y_max, double nominal_yield, TString syst)
{
  //calculates the signal_yield using different pdfs for signal or backgrounds (seperatly), and returns the maximum difference from the nominal yield.

  RooRealVar* fit_res;
  
  std::vector<std::string> signal = {"1gauss"}; //,"crystal", "3gauss"};
  std::vector<std::string> background = {"bern"}; //, "2exp", "power"};
  
  std::vector<std::string> pdf;
  TString pdf_name = "";
  
  std::cout << syst << std::endl;
  
  //copy the names of the pdfs into the pdf vector.
  if(syst == "signal_pdf")
    {
      pdf_name = "signal";
      pdf.reserve((int)signal.size()); 
      
      for(int i=0; i<(int)signal.size(); i++)
	pdf.push_back(signal[i]);
    }
  else
    if(syst == "cb_pdf")
      {
	pdf_name = "background";
	pdf.reserve((int)background.size()); 
	
	for(int i=0; i<(int)background.size(); i++)
	  pdf.push_back(background[i]);
      }
  
  //to save the various yield results
  std::vector<double> yield_syst;
  yield_syst.reserve((int)pdf.size());
  
  //calculate systematics
  for(int i=0; i<(int)pdf.size(); i++)
    {
      fit_res = bin_mass_fit(ws, channel, pt_min, pt_max, y_min, y_max, pdf[i], pdf_name.Data());
      yield_syst.push_back((double)fit_res->getVal());
    }
  
  //print table at the end
  for(int i=0; i<(int)yield_syst.size(); i++)
    std::cout << pdf_name << "_syst[" << i << "]: " << pdf[i] << " : " << yield_syst[i] << std::endl;
    
  int i_max = 0;
  double max_diff = 0;
  
  for(int i=0; i<(int)yield_syst.size(); i++)
    {
      if(fabs(yield_syst[i] - nominal_yield) > max_diff)
	{
	  max_diff = fabs(yield_syst[i] - nominal_yield);
	  i_max = i;
	}
    }
  
  return yield_syst[i_max];
}

double mass_window_syst(RooWorkspace& ws, int channel, double pt_min, double pt_max, double y_min, double y_max, double nominal_yield, TString input_file)
{
  std::cout << "calculate mass_window_syst" << std::endl;

  RooRealVar mass = *(ws.var("mass"));
  RooRealVar* fit_res;

  std::vector<double> mass_min;
  std::vector<double> mass_max;

  mass_min.push_back(mass.getMin());
  mass_min.push_back(0.95 * mass.getMin());
  mass_max.push_back(mass.getMax());
  mass_max.push_back(1.05 * mass.getMax());

  std::vector<double> range_syst;

  range_syst.reserve((int)mass_min.size());

  //Mass Range Systematics
  for(int i=0; i<(int)mass_min.size(); i++)
    {                                                                                                                                                                                                               
      RooWorkspace* ws1 = new RooWorkspace("ws1","Bmass");
      set_up_workspace_variables(*ws1,channel,mass_min[i],mass_max[1-i]);
      read_data(*ws1, input_file, channel);
    
      fit_res = bin_mass_fit(*ws1, channel, pt_min, pt_max, y_min, y_max, "", "", mass_min[i], mass_max[1-i]);
      range_syst.push_back((double)fit_res->getVal());
    }
  
  for(int i=0; i<(int)range_syst.size(); i++)
    std::cout << "mass_window_syst[" << i << "]: " << range_syst[i] << std::endl;

  int i_max = 0;
  double max_diff = 0;
  
  for(int i=0; i<(int)range_syst.size(); i++)
    {
      if(fabs(range_syst[i] - nominal_yield) > max_diff)
	{
	  max_diff = fabs(range_syst[i] - nominal_yield);
	  i_max = i;
	}
    }
  
  return range_syst[i_max];
}
