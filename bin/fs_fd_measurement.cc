#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

using namespace RooFit;

#define LUMINOSITY          2.71

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 4: Bs -> J/psi phi

//input example: fs_fd_ratio --ratio fs_fu --bins pt_y --preeff 1 --recoeff 1 --mcstudy 0 --syst 0
int main(int argc, char** argv)
{
  std::string ratio = "fs_fu";
  std::string yield_sub_samples = "full";
  int calculate_pre_filter_eff = 0;
  int calculate_reco_eff = 0;
  int mcstudy = 0;
  int syst = 0;

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--ratio")
	{
	  convert << argv[++i];
	  convert >> ratio;
	}
      if(argument == "--bins")
	{
	  convert << argv[++i];
	  convert >> yield_sub_samples;
	}
      if(argument == "--preeff")
	{
	  convert << argv[++i];
	  convert >> calculate_pre_filter_eff;
	}
      if(argument == "--recoeff")
	{
	  convert << argv[++i];
	  convert >> calculate_reco_eff;
	}
      if(argument == "--mcstudy")
	{
	  convert << argv[++i];
	  convert >> mcstudy;
	}
      if(argument == "--syst")
	{
	  convert << argv[++i];
	  convert >> syst;
	}
    }

  //pt bins
  double pt_bins[]= {9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 90};
  double total_pt_bin_edges[]={0, 400};
  int number_of_pt_bins = (sizeof(pt_bins) / sizeof(double)) - 1; //if pt_bins is an empty array, then number_of_pt_bins is equal to 0
  
  //y bins
  double y_bins[]= {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2, 2.25};
  double total_y_bin_edges[]={0.0, 2.25};
  int number_of_y_bins = (sizeof(y_bins) / sizeof(double)) - 1;

  //var1 and var2 are associated to pt or y below, as needed
  int n_var1_bins =1;
  int n_var2_bins =1;
  std::string var1_name = "";
  std::string var2_name = "";
  double* var1_bin_edges = total_pt_bin_edges;
  double* var2_bin_edges = total_y_bin_edges;
  
  //set up of var1 and var2
  if(yield_sub_samples == "full")
    {
      n_var1_bins= 1;
      var1_name = "pt";
      var1_bin_edges = total_pt_bin_edges;
      
      n_var2_bins= 1;
      var2_name = "y";
      var2_bin_edges = total_y_bin_edges;
    }
  else
    if(yield_sub_samples == "pt")
      {
	n_var1_bins = number_of_pt_bins;
	var1_name = "pt";
	var1_bin_edges = pt_bins;
      
	n_var2_bins = 1;
	var2_name = "y";
	var2_bin_edges = total_y_bin_edges;
      }
    else
      if(yield_sub_samples == "pt_y")
	{
	  n_var1_bins = number_of_pt_bins;
	  var1_name = "pt";
	  var1_bin_edges = pt_bins;
      
	  n_var2_bins = number_of_y_bins;
	  var2_name = "y";
	  var2_bin_edges = y_bins;
	}
      else
	if(yield_sub_samples == "y")
	  {
	    n_var1_bins = number_of_y_bins;
	    var1_name = "y";
	    var1_bin_edges = y_bins;
      
	    n_var2_bins = 1;
	    var2_name = "pt";
	    var2_bin_edges = total_pt_bin_edges;
	  }
	else
	  if(yield_sub_samples == "y_pt")
	    {
	      n_var1_bins = number_of_y_bins;
	      var1_name = "y";
	      var1_bin_edges = y_bins;
      
	      n_var2_bins = number_of_pt_bins;
	      var2_name = "pt";
	      var2_bin_edges = pt_bins;
	    }
	  else
	    {
	      printf("ERROR: The variables you asked for are not deffined. Only full, pt, y, pt_y, y_pt are deffined. Please check in the code.\n");
	      return 0;
	    }
  
  double yield_array[2][n_var2_bins][n_var1_bins];
  double errLo_array[2][n_var2_bins][n_var1_bins];
  double errHi_array[2][n_var2_bins][n_var1_bins];
  double yield_syst_array[2][n_var2_bins][n_var1_bins];
 
  double var1_bin_centre[n_var1_bins];
  double var1_bin_centre_Lo[n_var1_bins];
  double var1_bin_centre_Hi[n_var1_bins];
      
  RooRealVar* pre_filter_eff;
      
  double pre_eff_array[2][n_var2_bins][n_var1_bins];
  double pre_eff_err_lo_array[2][n_var2_bins][n_var1_bins];
  double pre_eff_err_hi_array[2][n_var2_bins][n_var1_bins];

  RooRealVar* reco_eff;
      
  double reco_eff_array[2][n_var2_bins][n_var1_bins];
  double reco_eff_err_lo_array[2][n_var2_bins][n_var1_bins];
  double reco_eff_err_hi_array[2][n_var2_bins][n_var1_bins];
 
  double total_eff_array[2][n_var2_bins][n_var1_bins];
  double total_eff_err_lo_array[2][n_var2_bins][n_var1_bins];
  double total_eff_err_hi_array[2][n_var2_bins][n_var1_bins];
  
  double ratio_eff_array[n_var2_bins][n_var1_bins];
  double ratio_eff_err_lo_array[n_var2_bins][n_var1_bins];
  double ratio_eff_err_hi_array[n_var2_bins][n_var1_bins];

  double fs_fd_array[n_var2_bins][n_var1_bins];
  double fs_fd_errLo_array[n_var2_bins][n_var1_bins];
  double fs_fd_errHi_array[n_var2_bins][n_var1_bins];
  double fs_fd_syst_lo_array[n_var2_bins][n_var1_bins];
  double fs_fd_syst_hi_array[n_var2_bins][n_var1_bins];
  
  double b_fraction[2];
  double b_fraction_err[2];

  for(int i=0; i<n_var1_bins; i++)
    {
      var1_bin_centre[i] = var1_bin_edges[i] + (var1_bin_edges[i+1]-var1_bin_edges[i])/2;
      var1_bin_centre_Lo[i] = var1_bin_centre[i] - var1_bin_edges[i];
      var1_bin_centre_Hi[i] = var1_bin_edges[i+1] - var1_bin_centre[i];
    }
  
  for(int ch=0; ch<2; ch++)
    { 
      int channel = 1; //just the default value for the channel, it should be assigned below.
      
      if(ratio == "fs_fu")
	channel = 3*ch+1; //fs_fu: if ch=0 -> channel=1, if ch=1 -> channel=4
      else
	if(ratio == "fs_fd")
	  channel = 2*(ch+1); //fs_fd: if ch=0 -> channel=2, if ch=1 -> channel=4
	else
	  if(ratio == "fd_fu")
	    channel= ch+1; //fd_fu: if ch=0 -> channel=1, if ch=1 -> channel=2
	  else
	    {
	      printf("ERROR: The ratio you asked for is not deffined. Only fs_fu, fs_fd, fd_fu are deffined. Please check in the code.");
	      return 0;
	    }
      
      //to create the directories to save the .png files
      std::vector<std::string> dir_list;
      dir_list.push_back(static_cast<const char*>("mass_fits/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION)));
      dir_list.push_back("fs_fd");
      dir_list.push_back("efficiencies");
      
      create_dir(dir_list);
      
      //------------data input---------------------
      //TString data_selection_input_file = TString::Format(BASE_DIR) + "selected_data_" + channel_to_ntuple_name(channel) + ".root";
      TString data_selection_input_file = TString::Format(BASE_DIR) + "selected_data_with_cuts.root";
      RooWorkspace* ws = new RooWorkspace("ws","Bmass");
      RooRealVar* signal_res; 
  
      //set up mass, pt and y variables inside ws  
      set_up_workspace_variables(*ws,channel);
      //read data from the selected data file, and import it as a dataset into the workspace.
      read_data(*ws, data_selection_input_file,channel);

      ws->Print();
      //------------------------------------------
            
      RooRealVar* branch = branching_fraction(channel);
      b_fraction[ch] = branch->getVal();
      b_fraction_err[ch] = branch->getError();

      for(int j=0; j<n_var2_bins; j++)
	{ 
	  std::cout << "processing subsample: " << var2_bin_edges[j] << " < " << var2_name << " < " << var2_bin_edges[j+1] << std::endl;

	  for(int i=0; i<n_var1_bins; i++)
	    {
	      std::cout << "processing subsample: " << var1_bin_edges[i] << " < " << var1_name << " < " << var1_bin_edges[i+1] << std::endl;
	      
	      //calculate the signal yield for a bin of pt and y.
	      if(var1_name == "pt")
		{
		  signal_res = bin_mass_fit(*ws,channel,var1_bin_edges[i],var1_bin_edges[i+1], var2_bin_edges[j], var2_bin_edges[j+1]);
		  yield_syst_array[ch][j][i] = bin_systematics(*ws, channel, var1_bin_edges[i], var1_bin_edges[i+1], var2_bin_edges[j], var2_bin_edges[j+1],signal_res->getVal(), data_selection_input_file, syst);
		}
	      else
		{
		  signal_res = bin_mass_fit(*ws,channel,var2_bin_edges[j],var2_bin_edges[j+1], var1_bin_edges[i], var1_bin_edges[i+1]);
		  yield_syst_array[ch][j][i] = bin_systematics(*ws, channel, var2_bin_edges[j], var2_bin_edges[j+1], var1_bin_edges[i], var1_bin_edges[i+1],signal_res->getVal(), data_selection_input_file, syst);
		}
	      
	      yield_array[ch][j][i] = signal_res->getVal();
	      errLo_array[ch][j][i] = -(signal_res->getAsymErrorLo());
	      errHi_array[ch][j][i] = signal_res->getAsymErrorHi();
	      
	      //MC study. only for --ratio full or pt or pt_y
	      if(mcstudy)
		{
		  std::cout << "MC study of bin: " << (int)var1_bin_edges[i] << " < pt < " << (int)var1_bin_edges[i+1] << std::endl;

		  mc_study(*ws,channel,var1_bin_edges[i],var1_bin_edges[i+1], var2_bin_edges[j], var2_bin_edges[j+1]);
		}
	    }
	}
      
      //to calculate pre-filter efficiency
      if(calculate_pre_filter_eff)
	{
	  for(int j=0; j<n_var2_bins; j++)
	    { 
	      std::cout << "processing subsample: " << var2_bin_edges[j] << " < " << var2_name << " < " << var2_bin_edges[j+1] << std::endl;
	      for(int i=0; i<n_var1_bins; i++)
		{
		  std::cout <<"calculating pre-filter efficiency: " << var1_bin_edges[i] << " < " << var1_name << " < " << var1_bin_edges[i+1] << std::endl;
		  		  	  
		  if(var1_name == "pt")
		    pre_filter_eff = prefilter_efficiency(channel,var1_bin_edges[i],var1_bin_edges[i+1],var2_bin_edges[j],var2_bin_edges[j+1]);
		  else
		    pre_filter_eff = prefilter_efficiency(channel,var2_bin_edges[j],var2_bin_edges[j+1],var1_bin_edges[i],var1_bin_edges[i+1]);
		  
		  pre_eff_array[ch][j][i] = pre_filter_eff->getVal();
		  pre_eff_err_lo_array[ch][j][i] = -(pre_filter_eff->getAsymErrorLo());
		  pre_eff_err_hi_array[ch][j][i] = pre_filter_eff->getAsymErrorHi();
		}

	      //plot of the pre-filter efficiency as a function of var1
	      TCanvas ce;
	      TGraphAsymmErrors* graph_pre_eff = new TGraphAsymmErrors(n_var1_bins, var1_bin_centre, pre_eff_array[ch][j], var1_bin_centre_Lo, var1_bin_centre_Hi, pre_eff_err_lo_array[ch][j], pre_eff_err_hi_array[ch][j]);
	      graph_pre_eff->SetTitle("Pre-filter efficiency");
	      graph_pre_eff->SetMarkerColor(4);
	      graph_pre_eff->SetMarkerStyle(21);
	      
	      TString x_axis_name = "";
	      if(var1_name =="pt")
	        x_axis_name = "p_{T}(B) [GeV]";
	      else
		x_axis_name = "|y|(B)";
	      
	      graph_pre_eff->GetXaxis()->SetTitle(x_axis_name);
	      graph_pre_eff->Draw("AP");
	      
	      TString eff1_name = "";
	      eff1_name = "efficiencies/pre_filter_efficiency_" + channel_to_ntuple_name(channel) + "_" + var2_name + TString::Format("_from_%.2f_to_%.2f",var2_bin_edges[j],var2_bin_edges[j+1]) + ".png";
	      ce.SaveAs(eff1_name);
	    }
	}

      //to calculate reconstruction efficiency
      if(calculate_reco_eff)
	{
	  for(int j=0; j<n_var2_bins; j++)
	    { 
	      std::cout << "processing subsample: " << var2_bin_edges[j] << " < " << var2_name << " < " << var2_bin_edges[j+1] << std::endl;
	      for(int i=0; i<n_var1_bins; i++)
		{
		  std::cout << "calculating reconstruction efficiency: " << var1_bin_edges[i] << " < " << var1_name << " < " << var1_bin_edges[i+1] << std::endl;
		   
		  if(var1_name =="pt")
		    reco_eff = reco_efficiency(channel,var1_bin_edges[i],var1_bin_edges[i+1],var2_bin_edges[j],var2_bin_edges[j+1]);
		  else
		    reco_eff = reco_efficiency(channel,var2_bin_edges[j],var2_bin_edges[j+1],var1_bin_edges[i],var1_bin_edges[i+1]);
		  
		  reco_eff_array[ch][j][i] = reco_eff->getVal();
		  reco_eff_err_lo_array[ch][j][i] = -(reco_eff->getAsymErrorLo());
		  reco_eff_err_hi_array[ch][j][i] = reco_eff->getAsymErrorHi();
		}

		  //plot of the reco efficiency as a function of var1
		  TCanvas cp;
		  TGraphAsymmErrors* graph_reco_eff = new TGraphAsymmErrors(n_var1_bins, var1_bin_centre, reco_eff_array[ch][j], var1_bin_centre_Lo, var1_bin_centre_Hi, reco_eff_err_lo_array[ch][j], reco_eff_err_hi_array[ch][j]);
		  graph_reco_eff->SetTitle("Reconstruction efficiency");
		  graph_reco_eff->SetMarkerColor(4);
		  graph_reco_eff->SetMarkerStyle(21);

		  TString x_axis_name = "";
		  if(var1_name =="pt")
		    x_axis_name = "p_{T}(B) [GeV]";
		  else
		    x_axis_name = "|y|(B)";
		  
		  graph_reco_eff->GetXaxis()->SetTitle(x_axis_name);
		  graph_reco_eff->Draw("AP");
		  
		  TString eff2_name = "";
		  eff2_name = "efficiencies/reco_efficiency_" + channel_to_ntuple_name(channel) + "_" + var2_name + TString::Format("_from_%.2f_to_%.2f",var2_bin_edges[j],var2_bin_edges[j+1]) + ".png";
		  cp.SaveAs(eff2_name);
	    }
	}
      
      if(calculate_pre_filter_eff && calculate_reco_eff)
	{
	  for(int j=0; j<n_var2_bins; j++)
	    { 
	      std::cout << "processing subsample: " << var2_bin_edges[j] << " < " << var2_name << " < " << var2_bin_edges[j+1] << std::endl;
	      for(int i=0; i<n_var1_bins; i++)
		{
		  std::cout << "calculating total efficiency: " << var1_bin_edges[i] << " < " << var1_name << " < " << var1_bin_edges[i+1] << std::endl;
		   
		  total_eff_array[ch][j][i] = pre_eff_array[ch][j][i] * reco_eff_array[ch][j][i];
		  total_eff_err_lo_array[ch][j][i] = total_eff_array[ch][j][i] * sqrt(pow(pre_eff_err_lo_array[ch][j][i]/pre_eff_array[ch][j][i],2) + pow(reco_eff_err_lo_array[ch][j][i]/reco_eff_array[ch][j][i],2));
		  total_eff_err_hi_array[ch][j][i] = total_eff_array[ch][j][i] * sqrt(pow(pre_eff_err_hi_array[ch][j][i]/pre_eff_array[ch][j][i],2) + pow(reco_eff_err_hi_array[ch][j][i]/reco_eff_array[ch][j][i],2));
		}

		  //plot of the total efficiency as a function of var1
		  TCanvas cp;
		  TGraphAsymmErrors* graph_total_eff = new TGraphAsymmErrors(n_var1_bins, var1_bin_centre, total_eff_array[ch][j], var1_bin_centre_Lo, var1_bin_centre_Hi, total_eff_err_lo_array[ch][j], total_eff_err_hi_array[ch][j]);
		  graph_total_eff->SetTitle("Overall efficiency");
		  graph_total_eff->SetMarkerColor(4);
		  graph_total_eff->SetMarkerStyle(21);

		  TString x_axis_name = "";
		  if(var1_name =="pt")
		    x_axis_name = "p_{T}(B) [GeV]";
		  else
		    x_axis_name = "|y|(B)";
		  
		  graph_total_eff->GetXaxis()->SetTitle(x_axis_name);
		  graph_total_eff->Draw("AP");
		  
		  TString eff3_name = "";
		  eff3_name = "efficiencies/total_efficiency_" + channel_to_ntuple_name(channel) + "_" + var2_name + TString::Format("_from_%.2f_to_%.2f",var2_bin_edges[j],var2_bin_edges[j+1]) + ".png";
		  cp.SaveAs(eff3_name);
	    }
	}
      
      if(calculate_pre_filter_eff && calculate_reco_eff)
	{
	  for(int j=0; j<n_var2_bins; j++)
	    { 
	      std::cout << "processing subsample: " << var2_bin_edges[j] << " < " << var2_name << " < " << var2_bin_edges[j+1] << std::endl;
	      for(int i=0; i<n_var1_bins; i++)
		{
		  std::cout << "calculating efficiency ratio: " << var1_bin_edges[i] << " < " << var1_name << " < " << var1_bin_edges[i+1] << std::endl;
		   
		  ratio_eff_array[j][i] = total_eff_array[0][j][i] * total_eff_array[1][j][i];
		  ratio_eff_err_lo_array[j][i] = ratio_eff_array[j][i] * sqrt(pow(total_eff_err_lo_array[0][j][i]/total_eff_array[0][j][i],2) + pow(total_eff_err_lo_array[1][j][i]/total_eff_array[1][j][i],2));
		  ratio_eff_err_hi_array[j][i] = ratio_eff_array[j][i] * sqrt(pow(total_eff_err_hi_array[0][j][i]/total_eff_array[0][j][i],2) + pow(total_eff_err_hi_array[1][j][i]/total_eff_array[1][j][i],2));
		}

		  //plot of the total efficiency as a function of var1
		  TCanvas cp;
		  TGraphAsymmErrors* graph_ratio_eff = new TGraphAsymmErrors(n_var1_bins, var1_bin_centre, ratio_eff_array[j], var1_bin_centre_Lo, var1_bin_centre_Hi, ratio_eff_err_lo_array[j], ratio_eff_err_hi_array[j]);
		  graph_ratio_eff->SetTitle("Efficiency ratio");
		  graph_ratio_eff->SetMarkerColor(4);
		  graph_ratio_eff->SetMarkerStyle(21);

		  TString x_axis_name = "";
		  if(var1_name =="pt")
		    x_axis_name = "p_{T}(B) [GeV]";
		  else
		    x_axis_name = "|y|(B)";
		  
		  graph_ratio_eff->GetXaxis()->SetTitle(x_axis_name);
		  graph_ratio_eff->Draw("AP");
		  
		  TString eff4_name = "";
		  eff4_name = "efficiencies/ratio_efficiency_" + ratio + "_" + var2_name + TString::Format("_from_%.2f_to_%.2f",var2_bin_edges[j],var2_bin_edges[j+1]) + ".png";
		  cp.SaveAs(eff4_name);
	    }
	}
      
      delete ws;
    }//end of channel cicle
  
  //to calculate the fragmentation fraction ratio
  for(int j=0; j<n_var2_bins; j++)
    {
      for(int i=0; i<n_var1_bins; i++)
	{
	  if(calculate_reco_eff && calculate_pre_filter_eff)
	    {
	      fs_fd_array[j][i] = (yield_array[1][j][i]/yield_array[0][j][i])* ((pre_eff_array[0][j][i]*reco_eff_array[0][j][i])/(pre_eff_array[1][j][i]*reco_eff_array[1][j][i])) * (b_fraction[0]/b_fraction[1]) * pow(10,j);
	      fs_fd_syst_lo_array[j][i]  = fs_fd_array[j][i] * sqrt(pow(yield_syst_array[0][j][i]/yield_array[0][j][i],2) + pow(yield_syst_array[1][j][i]/yield_array[1][j][i],2) + pow(pre_eff_err_lo_array[0][j][i]/pre_eff_array[0][j][i],2) + pow(reco_eff_err_lo_array[0][j][i]/reco_eff_array[0][j][i],2) + pow(pre_eff_err_lo_array[1][j][i]/pre_eff_array[1][j][i],2) + pow(reco_eff_err_lo_array[1][j][i]/reco_eff_array[1][j][i],2) + pow(b_fraction_err[0]/b_fraction[0],2) + pow(b_fraction_err[1]/b_fraction[1],2));
	      
	      fs_fd_syst_hi_array[j][i]  = fs_fd_array[j][i] * sqrt(pow(yield_syst_array[0][j][i]/yield_array[0][j][i],2) + pow(yield_syst_array[1][j][i]/yield_array[1][j][i],2) + pow(pre_eff_err_hi_array[0][j][i]/pre_eff_array[0][j][i],2) + pow(reco_eff_err_hi_array[0][j][i]/reco_eff_array[0][j][i],2) + pow(pre_eff_err_hi_array[1][j][i]/pre_eff_array[1][j][i],2) + pow(reco_eff_err_hi_array[1][j][i]/reco_eff_array[1][j][i],2) + pow(b_fraction_err[0]/b_fraction[0],2) + pow(b_fraction_err[1]/b_fraction[1],2));
	    }
	  else
	    {
	      fs_fd_array[j][i] = (yield_array[1][j][i]/yield_array[0][j][i]) * pow(10,j);
	
	      fs_fd_syst_lo_array[j][i]  = fs_fd_array[j][i] * sqrt( pow(yield_syst_array[0][j][i]/yield_array[0][j][i],2) + pow(yield_syst_array[1][j][i]/yield_array[1][j][i],2) + pow(b_fraction_err[0]/b_fraction[0],2) + pow(b_fraction_err[1]/b_fraction[1],2));
	      fs_fd_syst_hi_array[j][i]  = fs_fd_syst_lo_array[j][i]; //the syst errors without the efficiencies are symmetrical
	    }
	  
	  fs_fd_errLo_array[j][i] = fs_fd_array[j][i] * sqrt( pow(errLo_array[0][j][i]/yield_array[0][j][i],2) + pow(errLo_array[1][j][i]/yield_array[1][j][i],2) );
	  fs_fd_errHi_array[j][i] = fs_fd_array[j][i] * sqrt( pow(errHi_array[0][j][i]/yield_array[0][j][i],2) + pow(errHi_array[1][j][i]/yield_array[1][j][i],2) );
	}
    }
  
  ////////////////////////////////
  //plot fragmentation fraction//
  //////////////////////////////
  TCanvas cz;
  TPad *pad = new TPad("pad", "pad", 0.05, 0.05, 0.99, 0.99);
  pad->Draw();
  
  double leg_x2 = 0.98;
  double leg_x1 = leg_x2 - 0.25;
  double leg_y2 = 0.98;
  double leg_y1 = leg_y2 - (0.05 * n_var2_bins);  

  TLegend *leg = new TLegend (leg_x1, leg_y1, leg_x2, leg_y2);
  
  TLatex * tex = new TLatex(0.69,0.91,TString::Format("%.2f fb^{-1} (13 TeV)",LUMINOSITY));
  tex->SetNDC(kTRUE);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.04);
  tex->Draw();
  tex = new TLatex(0.11,0.91,"CMS Preliminary");
  tex->SetNDC(kTRUE);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  
  for(int j=0; j<n_var2_bins; j++)
    {
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(n_var1_bins, var1_bin_centre, fs_fd_array[j], var1_bin_centre_Lo, var1_bin_centre_Hi, fs_fd_errLo_array[j], fs_fd_errHi_array[j]);
      graph->SetTitle("Fragmentation fraction");
      graph->SetMarkerColor(2+j);
      graph->SetMarkerSize(0.5);
      graph->SetMarkerStyle(20+j);
      
      //draw this for the first bin of var2, or in case there is only one bin.
      if(j==0)
	{
	  TString x_axis_name = "";
	  if(var1_name == "pt")
	    x_axis_name = "p_{T}(B) [GeV]";
	  else
	    x_axis_name = "|y|(B)";
	  
	  graph->GetXaxis()->SetTitle(x_axis_name);
	  
	  //to set the range of the plot.
	  if(n_var2_bins > 1)
	    graph->GetYaxis()->SetRangeUser(1e-2, 10*fs_fd_array[n_var2_bins-1][0]);
	  else
	    graph->GetYaxis()->SetRangeUser(0, 2*fs_fd_array[n_var2_bins-1][0]);
	    
	  graph->Draw("ap same");
	}
      else //for the rest of the var2 bins.
	{
	  graph->Draw("p same");
	}
      
      //print the var2 bin label
      TString label = "";
      
      if(j!=0)
	label += TString::Format("(#times 10^{%d}) ",j);
      
      if(var1_name == "pt")
	label += TString::Format("%.2f",var2_bin_edges[j]) + " < |y| < " + TString::Format("%.2f",var2_bin_edges[j+1]);
      else
	label += TString::Format("%d",(int)var2_bin_edges[j]) + " < pt < " + TString::Format("%d",(int)var2_bin_edges[j+1]);
      
      leg->AddEntry(graph, label, "lp");
      
      //systematic errors, show the statistical err of eff as syst always
      if(calculate_pre_filter_eff && calculate_reco_eff)
	{
	  TGraphAsymmErrors* graph_syst = new TGraphAsymmErrors(n_var1_bins, var1_bin_centre, fs_fd_array[j], var1_bin_centre_Lo, var1_bin_centre_Hi, fs_fd_syst_lo_array[j], fs_fd_syst_hi_array[j]);
	  graph_syst->SetFillColor(j+2);
	  graph_syst->SetFillStyle(3001);
	  graph_syst->Draw("2 same");
	}
    }//end of the var2 for
  
  leg->Draw("same");
  cz.Update();
  
  if(n_var2_bins > 1)
    cz.SetLogy();
      
  TString systematic = "";
      
  if(syst)
    systematic = "_syst";
      
  if(calculate_pre_filter_eff && calculate_reco_eff)
    cz.SaveAs("fs_fd/" + ratio + "_" + yield_sub_samples + "_bins" + systematic + "_" + TString::Format(VERSION) + ".png");
  else
    cz.SaveAs("fs_fd/" + ratio + "_yield_ratio_" + yield_sub_samples + "_bins" + systematic + "_" + TString::Format(VERSION) + ".png");
  
  /////////////////////////////////////////////
  //To show the values at the end like tables//
  /////////////////////////////////////////////
  
  //signal yield
  std::cout << "SIGNAL YIELD:" << std::endl;
  for(int ch=0; ch<2; ch++)
    {
      std::cout << std::endl;
      std::cout << "CHANNEL: " << 2*(ch+1) << std::endl;
      
      for(int j=0; j<n_var2_bins; j++)
	{
	  std::cout << "BIN: " << var2_name << " " << var2_bin_edges[j] << " to " << var2_bin_edges[j+1] << " : " << std::endl;
	  
	  for(int i=0; i<n_var1_bins; i++)
	    {
	      std::cout << "BIN: " << var1_name << " " << var1_bin_edges[i] << " to " << var1_bin_edges[i+1] << " : " <<  yield_array[ch][j][i] << " +" << errHi_array[ch][j][i] << " -"<< errLo_array[ch][j][i] << " +-" << yield_syst_array[ch][j][i] << std::endl;
	    }
	  
	  std::cout << std::endl;
	}
    }//end of channel cicle

  //RATIO EFFICIENCY
  std::cout << "EFFICIENCY RATIO:" << std::endl;
  for(int j=0; j<n_var2_bins; j++)
    {
      std::cout << "BIN: " << var2_name << " " << var2_bin_edges[j] << " to " << var2_bin_edges[j+1] << " : " << std::endl;
      
      for(int i=0; i<n_var1_bins; i++)
	{
	  std::cout << "BIN: " << var1_name << " " << var1_bin_edges[i] << " to " << var1_bin_edges[i+1] << " : " <<  ratio_eff_array[j][i] << " +" << ratio_eff_err_hi_array[j][i] << " -"<< ratio_eff_err_lo_array[j][i] << std::endl;
	}
      
      std::cout << std::endl;
    }

  //Fracmentation fraction
  std::cout << "FRAGMENTATION FRACTION:" << ratio << std::endl;
  for(int j=0; j<n_var2_bins; j++)
    {
      std::cout << "BIN: " << var2_name << " " << var2_bin_edges[j] << " to " << var2_bin_edges[j+1] << " : " << std::endl;
      
      for(int i=0; i<n_var1_bins; i++)
	{
	  std::cout << "BIN: " << var1_name << " " << var1_bin_edges[i] << " to " << var1_bin_edges[i+1] << " : " <<  fs_fd_array[j][i] << " +" << fs_fd_errHi_array[j][i] << " -"<< fs_fd_errLo_array[j][i] << " +" << fs_fd_syst_hi_array[j][i] << " -" << fs_fd_syst_lo_array[j][i] << std::endl;
	}
      
      std::cout << std::endl;
    }

}//end of main
