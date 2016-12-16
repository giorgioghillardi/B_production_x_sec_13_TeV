#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

using namespace RooFit;

#define LUMINOSITY          2.71

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 4: Bs -> J/psi phi

//RooRealVar* branching_fraction(int channel);

//input example: fs_fd_ratio --ratio fs_fd --bins pt/y --preeff 1 --recoeff 1 --mc 0 --syst 0
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
      if(argument == "--mc")
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
  double pt_bins[]= {10, 15, 20, 25, 30, 35, 40, 60, 100}; //{9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 90};
  double total_pt_bin_edges[]={0, 400};
  int nptbins=1;
  double* pt_bin_edges=total_pt_bin_edges;

  //y bins
  double y_bins[]={0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2, 2.25};
  double total_y_bin_edges[]={0.0, 2.25};
  int nybins=1;
  double* y_bin_edges=total_y_bin_edges;

  //set the pt and y bin edges
  if(yield_sub_samples=="pt")
    {
	pt_bin_edges = pt_bins;
	nptbins = (sizeof(pt_bins) / sizeof(double)) - 1; //if pt_bin_edges is an empty array, then nptbins is equal to 0
    }
  if(yield_sub_samples=="y")
    {    
	y_bin_edges = y_bins;
	nybins = (sizeof(y_bins) / sizeof(double)) - 1;
    }
  if(yield_sub_samples=="pt_y" || yield_sub_samples=="y_pt")
    {
	pt_bin_edges = pt_bins;
	nptbins = (sizeof(pt_bins) / sizeof(double)) - 1;
	
	y_bin_edges = y_bins;
	nybins = (sizeof(y_bins) / sizeof(double)) - 1;
    }

  double yield_array[2][nybins][nptbins];
  double errLo_array[2][nybins][nptbins];
  double errHi_array[2][nybins][nptbins];
  double yield_syst_array[2][nybins][nptbins];
 
  double pt_bin_centre[nptbins];
  double pt_bin_centre_Lo[nptbins];
  double pt_bin_centre_Hi[nptbins];

  double y_bin_centre[nybins-1];
  double y_bin_centre_Lo[nybins-1];
  double y_bin_centre_Hi[nybins-1];
      
  RooRealVar* pre_filter_eff;
      
  double pre_eff_array[2][nybins][nptbins];
  double pre_eff_err_lo_array[2][nybins][nptbins];
  double pre_eff_err_hi_array[2][nybins][nptbins];

  double y_pre_eff_array[2][nybins];
  double y_pre_eff_err_lo_array[2][nybins];
  double y_pre_eff_err_hi_array[2][nybins];
      
  RooRealVar* reco_eff;
      
  double reco_eff_array[2][nybins][nptbins];
  double reco_eff_err_lo_array[2][nybins][nptbins];
  double reco_eff_err_hi_array[2][nybins][nptbins];
 
  double y_reco_eff_array[2][nybins];
  double y_reco_eff_err_lo_array[2][nybins];
  double y_reco_eff_err_hi_array[2][nybins];

  double eff_ratio_array[nybins][nptbins];
  double eff_ratio_err_lo_array[nybins][nptbins];
  double eff_ratio_err_hi_array[nybins][nptbins];

  double fs_fd_array[nybins][nptbins];
  double fs_fd_errLo_array[nybins][nptbins];
  double fs_fd_errHi_array[nybins][nptbins];
  double fs_fd_syst_lo_array[nybins][nptbins];
  double fs_fd_syst_hi_array[nybins][nptbins];
  
  double b_fraction[2];
  double b_fraction_err[2];

  for(int i=0; i<nptbins; i++)
    {
      pt_bin_centre[i] = pt_bin_edges[i] + (pt_bin_edges[i+1]-pt_bin_edges[i])/2;
      pt_bin_centre_Lo[i] = pt_bin_centre[i] - pt_bin_edges[i];
      pt_bin_centre_Hi[i] = pt_bin_edges[i+1] - pt_bin_centre[i];
    }
  
  for(int i=0; i<nybins; i++)
    {
      y_bin_centre[i] = y_bin_edges[i] + (y_bin_edges[i+1]-y_bin_edges[i])/2;
      y_bin_centre_Lo[i] = y_bin_centre[i] - y_bin_edges[i];
      y_bin_centre_Hi[i] = y_bin_edges[i+1] - y_bin_centre[i];	  
    }
  
  for(int ch=0; ch<2; ch++)
    { 
      //CHOOSE THE CHANNEL HERE!!
      int channel= 1; //just the default value for the channel, it should be assigned below.
      
      if(ratio == "fs_fu")
	channel = 3*ch+1; //fs_fu: if ch=0 -> channel=1, if ch=1 -> channel=4
      else
	if(ratio == "fs_fd")
	  channel = 2*(ch+1); //fs_fd: if ch=0 -> channel=2, if ch=1 -> channel=4
	else
	  if(ratio == "fd_fu")
	    {
	      channel= ch+1;
	      printf("THIS IS STILL NOT TESTED: TAKE IT WITH A GRAIN OF SALT...");
	    }
	  else
	    {
	      printf("ERROR: The ratio you asked for is not deffined. Please check in the code.");
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
      TString data_selection_input_file = "/exper-sw/cmst3/cmssw/users/brunogal/CMSSW_7_4_15/src/UserCode/B_production_x_sec_13_TeV/selected_data_with_cuts.root";
      RooWorkspace* ws = new RooWorkspace("ws","Bmass");
      RooRealVar* signal_res; 
  
      //set up mass, pt and y variables inside ws  
      set_up_workspace_variables(*ws,channel);
      //read data from the selected data file, and import it as a dataset into the workspace.
      read_data(*ws, data_selection_input_file,channel);

      ws->Print();
      //------------------------------------------
            
      RooRealVar* branch = branching_fraction(channel); //need to correct the function  
      b_fraction[ch] = branch->getVal();
      b_fraction_err[ch] = branch->getError();

      for(int j=0; j<nybins; j++)
	{ 
	  std::cout << "processing subsample: " << y_bin_edges[j] << " < |y| < " << y_bin_edges[j+1] << std::endl;

	  for(int i=0; i<nptbins; i++)
	    {
	      std::cout << "processing subsample: " << (int)pt_bin_edges[i] << " < pt < " << (int)pt_bin_edges[i+1] << std::endl;
	      
	      //calculate the signal yield for a bin of pt and y.
	      signal_res = bin_mass_fit(*ws,channel,pt_bin_edges[i],pt_bin_edges[i+1], y_bin_edges[j], y_bin_edges[j+1], mcstudy);
	      
	      yield_array[ch][j][i] = signal_res->getVal();
	      errLo_array[ch][j][i] = -(signal_res->getAsymErrorLo());
	      errHi_array[ch][j][i] = signal_res->getAsymErrorHi();
	      yield_syst_array[ch][j][i] = bin_systematics(*ws, channel, pt_bin_edges[i], pt_bin_edges[i+1], y_bin_edges[j], y_bin_edges[j+1],signal_res->getVal(), data_selection_input_file, syst);
	    }
	}
      
      //to calculate pre-filter efficiency
      if(calculate_pre_filter_eff)
	{
	  for(int c=0; c<nybins; c++)
	    { 
	      std::cout << "processing subsample: " << y_bin_edges[c] << " < |y| < " << y_bin_edges[c+1] << std::endl;
	      for(int i=0; i<nptbins; i++)
		{
		  std::cout <<"calculating pre-filter efficiency: "<< (int)pt_bin_edges[i] <<" < pt < "<< (int)pt_bin_edges[i+1] << std::endl;
		  		  	  
		  //pre-filter efficiency
		  pre_filter_eff = prefilter_efficiency(channel,pt_bin_edges[i],pt_bin_edges[i+1],y_bin_edges[c],y_bin_edges[c+1]);
		  
		  pre_eff_array[ch][c][i] = pre_filter_eff->getVal();
		  pre_eff_err_lo_array[ch][c][i] = -(pre_filter_eff->getAsymErrorLo());
		  pre_eff_err_hi_array[ch][c][i] = pre_filter_eff->getAsymErrorHi();
		}

	      if(yield_sub_samples=="pt_y" || yield_sub_samples=="pt")
		{
		  //plot of the pre-filter efficiency as a function of pT for each y bin
		  TCanvas ce;
		  TGraphAsymmErrors* graph_pre_eff = new TGraphAsymmErrors(nptbins, pt_bin_centre, pre_eff_array[ch][c], pt_bin_centre_Lo, pt_bin_centre_Hi, pre_eff_err_lo_array[ch][c], pre_eff_err_hi_array[ch][c]);
		  graph_pre_eff->SetTitle("pre-filter efficiency");
		  graph_pre_eff->SetMarkerColor(4);
		  graph_pre_eff->SetMarkerStyle(21);
		  graph_pre_eff->Draw("AP");
		  TString eff1_name = "";
		  eff1_name = "efficiencies/pre_filter_efficiency_" + channel_to_ntuple_name(channel) + "_" + TString::Format("y_from_%.1f_to_%.1f",y_bin_edges[c],y_bin_edges[c+1]) + ".png";
		  ce.SaveAs(eff1_name);
		}
	      if(yield_sub_samples=="y")
		{
		  //transpose the arrays to plot as function of y
		  y_pre_eff_array[ch][c] = pre_eff_array[ch][c][0];
		  y_pre_eff_err_lo_array[ch][c] = pre_eff_err_lo_array[ch][c][0];
		  y_pre_eff_err_hi_array[ch][c] = pre_eff_err_hi_array[ch][c][0];
		}
	    }
	  if(yield_sub_samples=="y")
	    {
	      //plot of the pre-filter efficiency as a function of y
	      TCanvas cj;
	      TGraphAsymmErrors* y_graph_pre_eff = new TGraphAsymmErrors(nybins, y_bin_centre, y_pre_eff_array[ch], y_bin_centre_Lo, y_bin_centre_Hi, y_pre_eff_err_lo_array[ch], y_pre_eff_err_hi_array[ch]);
	      y_graph_pre_eff->SetTitle("pre-filter efficiency");
	      y_graph_pre_eff->SetMarkerColor(4);
	      y_graph_pre_eff->SetMarkerStyle(21);
	      y_graph_pre_eff->Draw("AP");
	      TString eff1_name = "";
	      eff1_name = "efficiencies/y_pre_filter_efficiency_" + channel_to_ntuple_name(channel) + ".png";
	      cj.SaveAs(eff1_name);
	    }
	}

      //to calculate reconstruction efficiency
      if(calculate_reco_eff)
	{
	  for(int c=0; c<nybins; c++)
	    { 
	      std::cout << "processing subsample: " << y_bin_edges[c] << " < |y| < " << y_bin_edges[c+1] << std::endl;
	      for(int i=0; i<nptbins; i++)
		{
		  std::cout << "calculating reconstruction efficiency: " << (int)pt_bin_edges[i] << " < pt < " << (int)pt_bin_edges[i+1] << std::endl;
		   
		  //reco efficiency
		  reco_eff = reco_efficiency(channel,pt_bin_edges[i],pt_bin_edges[i+1],y_bin_edges[c],y_bin_edges[c+1]);
		 
		  reco_eff_array[ch][c][i] = reco_eff->getVal();
		  reco_eff_err_lo_array[ch][c][i] = -(reco_eff->getAsymErrorLo());
		  reco_eff_err_hi_array[ch][c][i] = reco_eff->getAsymErrorHi();
		}

	      if(yield_sub_samples=="pt_y" || yield_sub_samples=="pt")
		{
		  //plot of the reco efficiency as a function of pT for each y bin
		  TCanvas cp;
		  TGraphAsymmErrors* graph_reco_eff = new TGraphAsymmErrors(nptbins, pt_bin_centre, reco_eff_array[ch][c], pt_bin_centre_Lo, pt_bin_centre_Hi, reco_eff_err_lo_array[ch][c], reco_eff_err_hi_array[ch][c]);
		  graph_reco_eff->SetTitle("reconstruction efficiency");
		  graph_reco_eff->SetMarkerColor(4);
		  graph_reco_eff->SetMarkerStyle(21);
		  graph_reco_eff->Draw("AP");
		  TString eff2_name = "";
		  eff2_name = "efficiencies/reco_efficiency_" + channel_to_ntuple_name(channel) + "_" + TString::Format("y_from_%.1f_to_%.1f",y_bin_edges[c],y_bin_edges[c+1]) + ".png";
		  cp.SaveAs(eff2_name);
		}
	      if(yield_sub_samples=="y")
		{
		  //transpose the arrays to plot as function of y
		  y_reco_eff_array[ch][c] = reco_eff_array[ch][c][0];
		  y_reco_eff_err_lo_array[ch][c] = reco_eff_err_lo_array[ch][c][0];
		  y_reco_eff_err_hi_array[ch][c] = reco_eff_err_hi_array[ch][c][0];
		}
	    }
	  if(yield_sub_samples=="y")
	    {
	      //plot of the reco efficiency as a function of y
	      TCanvas cq;
	      TGraphAsymmErrors* y_graph_reco_eff = new TGraphAsymmErrors(nybins, y_bin_centre, y_reco_eff_array[ch], y_bin_centre_Lo, y_bin_centre_Hi, y_reco_eff_err_lo_array[ch], y_reco_eff_err_hi_array[ch]);
	      y_graph_reco_eff->SetTitle("reconstruction efficiency");
	      y_graph_reco_eff->SetMarkerColor(4);
	      y_graph_reco_eff->SetMarkerStyle(21);
	      y_graph_reco_eff->Draw("AP");
	      TString eff2_name = "";
	      eff2_name = "efficiencies/y_reco_efficiency_" + channel_to_ntuple_name(channel) + ".png";
	      cq.SaveAs(eff2_name);
	    }
	}
      delete ws;
    }//end of channel cicle
  
  //To show the values of signal yield and the errors at the end, like a table
  for(int ch=0; ch<2; ch++)
    {
      std::cout << std::endl;
      std::cout << "CHANNEL: " << 2*(ch+1) << std::endl;
      
      for(int j=0; j<nybins; j++)
	{
	  std::cout << "BIN y: " << y_bin_edges[j] << " to " << y_bin_edges[j+1] << " : " << std::endl;
	  
	  for(int i=0; i<nptbins; i++)
	    {
	      std::cout << "BIN pt: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  yield_array[ch][j][i] << " +" << errHi_array[ch][j][i] << " -"<< errLo_array[ch][j][i] << " +-" << yield_syst_array[ch][j][i] << std::endl;
	    }
	  
	  std::cout << std::endl;
	}
    }//end of channel cicle
  
  //to calculate the fd/fs
  for(int j=0; j<nybins; j++)
    {
      for(int i=0; i<nptbins; i++)
	{
	  if(calculate_reco_eff && calculate_pre_filter_eff)
	    {
	      if(yield_sub_samples=="y")
		fs_fd_array[j][i] =(yield_array[1][j][i]/yield_array[0][j][i])* ((pre_eff_array[0][j][i]*reco_eff_array[0][j][i])/(pre_eff_array[1][j][i]*reco_eff_array[1][j][i])) * (b_fraction[0]/b_fraction[1]);
	      else
		fs_fd_array[j][i] = (yield_array[1][j][i]/yield_array[0][j][i])* ((pre_eff_array[0][j][i]*reco_eff_array[0][j][i])/(pre_eff_array[1][j][i]*reco_eff_array[1][j][i])) * (b_fraction[0]/b_fraction[1]) * pow(10,j);

	      //calculate eff ratio, only implemented here, to extract the value and the stat error.	      
	      eff_ratio_array[j][i] = (pre_eff_array[0][j][i]*reco_eff_array[0][j][i])/(pre_eff_array[1][j][i]*reco_eff_array[1][j][i]);

	      eff_ratio_err_lo_array[j][i] = eff_ratio_array[j][i] * sqrt( pow(pre_eff_err_lo_array[0][j][i]/pre_eff_array[0][j][i],2) + pow(reco_eff_err_lo_array[0][j][i]/reco_eff_array[0][j][i],2) + pow(pre_eff_err_lo_array[1][j][i]/pre_eff_array[1][j][i],2) + pow(reco_eff_err_lo_array[1][j][i]/reco_eff_array[1][j][i],2));

	      eff_ratio_err_hi_array[j][i] = eff_ratio_array[j][i] * sqrt(pow(pre_eff_err_hi_array[0][j][i]/pre_eff_array[0][j][i],2) + pow(reco_eff_err_hi_array[0][j][i]/reco_eff_array[0][j][i],2) + pow(pre_eff_err_hi_array[1][j][i]/pre_eff_array[1][j][i],2) + pow(reco_eff_err_hi_array[1][j][i]/reco_eff_array[1][j][i],2));
	      //----------------------------------------------------------------------------------------------------

	      fs_fd_syst_lo_array[j][i]  = fs_fd_array[j][i] * sqrt(pow(yield_syst_array[0][j][i]/yield_array[0][j][i],2) + pow(yield_syst_array[1][j][i]/yield_array[1][j][i],2) + pow(pre_eff_err_lo_array[0][j][i]/pre_eff_array[0][j][i],2) + pow(reco_eff_err_lo_array[0][j][i]/reco_eff_array[0][j][i],2) + pow(pre_eff_err_lo_array[1][j][i]/pre_eff_array[1][j][i],2) + pow(reco_eff_err_lo_array[1][j][i]/reco_eff_array[1][j][i],2) + pow(b_fraction_err[0]/b_fraction[0],2) + pow(b_fraction_err[1]/b_fraction[0],2));

	      fs_fd_syst_hi_array[j][i]  = fs_fd_array[j][i] * sqrt(pow(yield_syst_array[0][j][i]/yield_array[0][j][i],2) + pow(yield_syst_array[1][j][i]/yield_array[1][j][i],2) + pow(pre_eff_err_hi_array[0][j][i]/pre_eff_array[0][j][i],2) + pow(reco_eff_err_hi_array[0][j][i]/reco_eff_array[0][j][i],2) + pow(pre_eff_err_hi_array[1][j][i]/pre_eff_array[1][j][i],2) + pow(reco_eff_err_hi_array[1][j][i]/reco_eff_array[1][j][i],2) + pow(b_fraction_err[0]/b_fraction[0],2) + pow(b_fraction_err[1]/b_fraction[0],2));
	    }
	  else
	    {
	      if(yield_sub_samples=="y")
		fs_fd_array[j][i] = (yield_array[1][j][i]/yield_array[0][j][i]);
	      else
		{
		  printf("debug: signal yield ratio %d \n", j);
		  fs_fd_array[j][i] = (yield_array[1][j][i]/yield_array[0][j][i]) * pow(10,j);
		}

	      fs_fd_syst_lo_array[j][i]  = fs_fd_array[j][i] * sqrt( pow(yield_syst_array[0][j][i]/yield_array[0][j][i],2) + pow(yield_syst_array[1][j][i]/yield_array[1][j][i],2) );
	      fs_fd_syst_hi_array[j][i]  = fs_fd_syst_lo_array[j][i];
	    }
	  fs_fd_errLo_array[j][i] = fs_fd_array[j][i] * sqrt( pow(errLo_array[0][j][i]/yield_array[0][j][i],2) + pow(errLo_array[1][j][i]/yield_array[1][j][i],2) );
	  fs_fd_errHi_array[j][i] = fs_fd_array[j][i] * sqrt( pow(errHi_array[0][j][i]/yield_array[0][j][i],2) + pow(errHi_array[1][j][i]/yield_array[1][j][i],2) );
	}
    }

  //To show the values of fs/fd and the errors at the end, like a table
  for(int j=0; j<nybins; j++)
    {
      std::cout << "BIN y: " << y_bin_edges[j] << " to " << y_bin_edges[j+1] << " : " << std::endl;
      
      for(int i=0; i<nptbins; i++)
	{
	  std::cout << "BIN pt: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  fs_fd_array[j][i] << " +" << fs_fd_errHi_array[j][i] << " -"<< fs_fd_errLo_array[j][i] << " +" << fs_fd_syst_hi_array[j][i] << " -" << fs_fd_syst_lo_array[j][i] << std::endl;
	}
      
      std::cout << std::endl;
    }

    //To show the values of pre eff at the end, like a table
  Printf("=====================DEBUG: Bs pre eff values==========================");
  for(int j=0; j<nybins; j++)
    {
      std::cout << "BIN y: " << y_bin_edges[j] << " to " << y_bin_edges[j+1] << " : " << std::endl;
      
      for(int i=0; i<nptbins; i++)
	{
	  std::cout << "BIN pt: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  pre_eff_array[1][j][i] << " +" << pre_eff_err_hi_array[1][j][i] << " -"<< pre_eff_err_lo_array[1][j][i] << std::endl;
	}
      
      std::cout << std::endl;
    }

  //To show the values of reco eff at the end, like a table
  Printf("=====================DEBUG: Bs reco eff values==========================");
  for(int j=0; j<nybins; j++)
    {
      std::cout << "BIN y: " << y_bin_edges[j] << " to " << y_bin_edges[j+1] << " : " << std::endl;
      
      for(int i=0; i<nptbins; i++)
	{
	  std::cout << "BIN pt: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  reco_eff_array[1][j][i] << " +" << reco_eff_err_hi_array[1][j][i] << " -"<< reco_eff_err_lo_array[1][j][i] << std::endl;
	}
      
      std::cout << std::endl;
    }

  //To show the values of eff ratio at the end, like a table
  Printf("=====================DEBUG: eff ratio values==========================");
  for(int j=0; j<nybins; j++)
    {
      std::cout << "BIN y: " << y_bin_edges[j] << " to " << y_bin_edges[j+1] << " : " << std::endl;
      
      for(int i=0; i<nptbins; i++)
	{
	  std::cout << "BIN pt: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  eff_ratio_array[j][i] << " +" << eff_ratio_err_hi_array[j][i] << " -"<< eff_ratio_err_lo_array[j][i] << std::endl;
	}
      
      std::cout << std::endl;
    }
  
  //plot fs/fd
  TCanvas cz;
  TPad *pad = new TPad("pad", "pad", 0.05, 0.05, 0.99, 0.99);
  pad->Draw();      
  TLegend *leg = new TLegend (0.70, 0.70, 0.90, 0.90);
      
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
  
  if(yield_sub_samples=="pt_y" || yield_sub_samples=="pt") //plot as a funtion of pt, with one or more y bins
    {
  for(int i=0; i<nybins; i++)
    {
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(nptbins, pt_bin_centre, fs_fd_array[i], pt_bin_centre_Lo, pt_bin_centre_Hi, fs_fd_errLo_array[i], fs_fd_errHi_array[i]);
      graph->SetTitle("Fragmentation fraction");
      graph->SetMarkerColor(i+2);
      graph->SetMarkerSize(0.5);
      graph->SetMarkerStyle(20+i);
	  
      //draw this for the first rapidity bin, or in case there is only one rapidity bin.
      if(i==0) 
	{
	  graph->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
	  //to set the range of the plot.
	  graph->GetYaxis()->SetRangeUser(0,2*fs_fd_array[nybins-1][0]);
	  graph->Draw("ap same");
	  //print the rapidity bin
	  leg->AddEntry(graph, TString::Format("%.1f < |y| < %.1f",y_bin_edges[i],y_bin_edges[i+1]), "lp");
	}
      else //for the rest of the rapidity bins.
	{     
	  graph->Draw("p same");
	  leg->AddEntry(graph, TString::Format("(#times 10^{%d}) %.1f < |y| < %.1f",i,y_bin_edges[i],y_bin_edges[i+1]), "lp");
	}

      //systematic errors, show the statistical err of eff as syst always
      if(calculate_pre_filter_eff && calculate_reco_eff)
	{
	  TGraphAsymmErrors* graph_syst = new TGraphAsymmErrors(nptbins, pt_bin_centre, fs_fd_array[i], pt_bin_centre_Lo, pt_bin_centre_Hi, fs_fd_syst_lo_array[i], fs_fd_syst_hi_array[i]);
	  graph_syst->SetFillColor(i+2);
	  graph_syst->SetFillStyle(3001);
	  graph_syst->Draw("2 same");
	}
    }//end of the ybins for

  leg->Draw("same");
  cz.Update();
  //cz.SetLogy();
    }
  
  if(yield_sub_samples=="y") //plot as a function of y, only one pt bin for now
    {
      //to transpose the arrays fd_fd_
      double y_fs_fd_array[nybins];
      double y_fs_fd_errLo_array[nybins];
      double y_fs_fd_errHi_array[nybins];
      double y_fs_fd_syst_lo_array[nybins];
      double y_fs_fd_syst_hi_array[nybins];

      for(int i=0; i<nybins; i++)
	{
	  y_fs_fd_array[i] = fs_fd_array[i][0];
	  y_fs_fd_errLo_array[i] = fs_fd_errLo_array[i][0];
	  y_fs_fd_errHi_array[i] = fs_fd_errHi_array[i][0];
	  y_fs_fd_syst_lo_array[i] = fs_fd_syst_lo_array[i][0];
	  y_fs_fd_syst_hi_array[i] = fs_fd_syst_hi_array[i][0];
	}

      TGraphAsymmErrors* graph = new TGraphAsymmErrors(nybins, y_bin_centre, y_fs_fd_array, y_bin_centre_Lo, y_bin_centre_Hi, y_fs_fd_errLo_array, y_fs_fd_errHi_array);
      graph->SetTitle("Fragmentation fraction");
      graph->SetMarkerColor(2);
      graph->SetMarkerSize(0.5);
      graph->SetMarkerStyle(20);
      
      //draw this for the first rapidity bin, or in case there is only one rapidity bin.
      graph->GetXaxis()->SetTitle("y(B)");
      //to set the range of the plot.
      graph->GetYaxis()->SetRangeUser(0.1*y_fs_fd_array[0],2*y_fs_fd_array[0]);
      graph->Draw("ap same");
      
      //systematic errors, show the statistical err of eff as syst always
      if(calculate_pre_filter_eff && calculate_reco_eff)
	{
	  TGraphAsymmErrors* graph_syst = new TGraphAsymmErrors(nybins, y_bin_centre, y_fs_fd_array, y_bin_centre_Lo, y_bin_centre_Hi, y_fs_fd_syst_lo_array, y_fs_fd_syst_hi_array);
	  graph_syst->SetFillColor(2);
	  graph_syst->SetFillStyle(3001);
	  graph_syst->Draw("2 same");
	}
      
      cz.Update();
      //cz.SetLogy();
    }
      
  TString systematic = "";
      
  if(syst)
    systematic = "_syst";
      
  if(calculate_pre_filter_eff && calculate_reco_eff)
    {
      cz.SaveAs("fs_fd/" + ratio + "_" + yield_sub_samples + "_bins" + systematic + "_" + TString::Format(VERSION) + ".png");
    }
  else
    {
      cz.SaveAs("fs_fd/" + ratio + "_yield_ratio_" + yield_sub_samples + "_bins" + systematic + "_" + TString::Format(VERSION) + ".png");
    }

  /* 
  Printf("+++++++++++++++Debug: table++++++++++++++++++");

  std::string table = "table_test";

  std::vector<std::string> col_name;
  col_name.push_back("col1");
  col_name.push_back("col2");

  std::vector<std::string> labels;
  labels.push_back("line1");
  labels.push_back("line2");
  labels.push_back("line3");

  std::vector<std::vector<double> > numbers = {{1,2} , {3,4} , {5,6}};
  std::string caption = "caption_test";

  latex_table(table, 2, 3, col_name, labels, numbers, caption);
  */

}//end of main
