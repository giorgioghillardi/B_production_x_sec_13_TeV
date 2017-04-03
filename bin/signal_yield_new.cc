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

//input example: signal_yield_new --channel 1 --bins pt_y --preeff 1 --recoeff 1 --mcstudy 0 --syst 0
int main(int argc, char** argv)
{
  int channel = 0;
  std::string yield_sub_samples = "full";
  int calculate_pre_filter_eff = 0;
  int calculate_reco_eff = 0;
  int mcstudy = 0;
  int syst = 0;

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

  if(channel==0)
    {
      std::cout << "Warning: No channel was provided as input. Please use --channel. Example: signal_yield_new --channel 1" << std::endl;
      return 0;
    }
  
  //to create the directories to save the .png files
  std::vector<std::string> dir_list;  
  
  dir_list.push_back(static_cast<const char*>("mass_fits/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION)));
  dir_list.push_back(static_cast<const char*>("mc_study/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION)));
  dir_list.push_back("signal_yield");
  dir_list.push_back("efficiencies");
  
  create_dir(dir_list);

  TString data_selection_input_file = TString::Format(BASE_DIR) + "selected_myloop_new_data_" + channel_to_ntuple_name(channel) + "_with_cuts.root";
  RooWorkspace* ws = new RooWorkspace("ws","Bmass");
  RooRealVar* signal_res; 
  
  //set up mass, pt and y variables inside ws  
  set_up_workspace_variables(*ws,channel);
  //read data from the selected data file, and import it as a dataset into the workspace.
  read_data(*ws, data_selection_input_file,channel);
  
  ws->Print();
  
  //pt bins
  double ntkp_pt_bins[]={9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 80, 90, 120};
  double ntkstar_pt_bins[]={9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 90};
  double ntks_pt_bins[]={0};
  double ntphi_pt_bins[]={9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 90};
  double ntmix_pt_bins[]={0};
  double ntlambda_pt_bins[]={0};
  
  int number_of_pt_bins=1;
  double total_pt_bin_edges[]={10, 150};
  double* pt_bins=total_pt_bin_edges;

  //y bins
  double ntkp_y_bins[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntkstar_y_bins[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntks_y_bins[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntphi_y_bins[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntmix_y_bins[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntlambda_y_bins[]={0.0, 0.5, 1.0, 1.5, 2.25};

  int number_of_y_bins=1;  
  double total_y_bin_edges[]={0.0, 2.25};
  double* y_bins=total_y_bin_edges;
  
  //var1 and var2 are associated to pt or y below, as needed
  int n_var1_bins =1;
  int n_var2_bins =1;
  std::string var1_name = "";
  std::string var2_name = "";
  double* var1_bin_edges = total_pt_bin_edges;
  double* var2_bin_edges = total_y_bin_edges;

  switch (channel)
    {
    default:
    case 1:
      pt_bins = ntkp_pt_bins;
      number_of_pt_bins = (sizeof(ntkp_pt_bins) / sizeof(double)) - 1 ; //if pt_bin_edges is an empty array, then nptbins is equal to 0 
      y_bins = ntkp_y_bins;
      number_of_y_bins = (sizeof(ntkp_y_bins) / sizeof(double)) - 1 ; //if y_bin_edges is an empty array, then nptbins is equal to 0
      break;
      
    case 2:
      pt_bins = ntkstar_pt_bins;
      number_of_pt_bins = (sizeof(ntkstar_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntkstar_y_bins;	
      number_of_y_bins = (sizeof(ntkstar_y_bins) / sizeof(double)) - 1 ;
      break;
      
    case 3:
      pt_bins = ntks_pt_bins;
      number_of_pt_bins = (sizeof(ntks_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntks_y_bins;
      number_of_y_bins = (sizeof(ntks_y_bins) / sizeof(double)) - 1 ;
      break;
      
    case 4:
      pt_bins = ntphi_pt_bins;
      number_of_pt_bins = (sizeof(ntphi_pt_bins) / sizeof(double)) - 1 ;	  
      y_bins = ntphi_y_bins;
      number_of_y_bins = (sizeof(ntphi_y_bins) / sizeof(double)) - 1 ;
      break;
      
    case 5:
      pt_bins = ntmix_pt_bins;
      number_of_pt_bins = (sizeof(ntmix_pt_bins) / sizeof(double)) - 1 ;	  
      y_bins = ntmix_y_bins;
      number_of_y_bins = (sizeof(ntmix_y_bins) / sizeof(double)) - 1 ;
      break;
      
    case 6:
      pt_bins = ntlambda_pt_bins;
      number_of_pt_bins = (sizeof(ntlambda_pt_bins) / sizeof(double)) - 1 ;	  
      y_bins = ntlambda_y_bins;
      number_of_y_bins = (sizeof(ntlambda_y_bins) / sizeof(double)) - 1 ;
      break;
    }
  
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

  double var1_bin_size[n_var1_bins];
  double var2_bin_size[n_var2_bins];

  double var1_bin_means[n_var1_bins];
  double var1_bin_edges_lo[n_var1_bins];
  double var1_bin_edges_hi[n_var1_bins];

  double var1_bin_centre[n_var1_bins];
  double var1_bin_centre_lo[n_var1_bins];
  double var1_bin_centre_hi[n_var1_bins];

  double yield_array[n_var2_bins][n_var1_bins];
  double errLo_array[n_var2_bins][n_var1_bins];
  double errHi_array[n_var2_bins][n_var1_bins];
  double yield_syst_array[n_var2_bins][n_var1_bins];

  double pre_eff_array[n_var2_bins][n_var1_bins];
  double pre_eff_err_lo_array[n_var2_bins][n_var1_bins];
  double pre_eff_err_hi_array[n_var2_bins][n_var1_bins];
    
  double reco_eff_array[n_var2_bins][n_var1_bins];
  double reco_eff_err_lo_array[n_var2_bins][n_var1_bins];
  double reco_eff_err_hi_array[n_var2_bins][n_var1_bins];
    
  double total_eff_array[n_var2_bins][n_var1_bins];
  double total_eff_err_lo_array[n_var2_bins][n_var1_bins];
  double total_eff_err_hi_array[n_var2_bins][n_var1_bins];

  double x_sec_array[n_var2_bins][n_var1_bins];
  double x_sec_errLo_array[n_var2_bins][n_var1_bins];
  double x_sec_errHi_array[n_var2_bins][n_var1_bins];
  double x_sec_syst_lo_array[n_var2_bins][n_var1_bins];
  double x_sec_syst_hi_array[n_var2_bins][n_var1_bins];
    
  RooRealVar* branch = branching_fraction(channel);

  double bin_size;

  TString eff_title = "";
  TString eff_name = "";

  TString sample = "";
  if(yield_sub_samples == "full")
    sample = "_" + yield_sub_samples;

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
  
  //calculate the centres, mean values, and bin size of var1
  for(int i=0; i<n_var1_bins; i++)
    {
      var1_bin_size[i] = var1_bin_edges[i+1] - var1_bin_edges[i];

      var1_bin_means[i] = var_mean_value(*ws, var1_name, var1_bin_edges[i], var1_bin_edges[i+1]);
      var1_bin_edges_lo[i] = var1_bin_means[i] - var1_bin_edges[i];
      var1_bin_edges_hi[i] = var1_bin_edges[i+1] - var1_bin_means[i];
		
      var1_bin_centre[i] = var1_bin_edges[i] + (var1_bin_edges[i+1]-var1_bin_edges[i])/2;
      var1_bin_centre_lo[i] = var1_bin_centre[i] - var1_bin_edges[i];
      var1_bin_centre_hi[i] = var1_bin_edges[i+1] - var1_bin_centre[i];
    }
  
  //calculate the size of the var2 bins
  for(int j=0; j<n_var2_bins; j++)
    {
      var2_bin_size[j] = var2_bin_edges[j+1] - var2_bin_edges[j];
    }
  
  if(yield_sub_samples=="full")
    {
      var1_bin_size[0] = 1; //force the bin size equal to one when we want to calculate full dataset cross section.
      var2_bin_size[0] = 1;
    }
  else
    if(yield_sub_samples=="pt" || yield_sub_samples=="y")
      {
	var2_bin_size[0]=1;
      }
  
  TString yield_file_name = channel_to_ntuple_name(channel) + "_yield.root";
  TFile* yield_file = new TFile(yield_file_name,"recreate");
  TVectorD signal_yield(n_var1_bins);
  
  for(int j=0; j<n_var2_bins; j++)
    { 
      std::cout << "processing subsample: " << var2_bin_edges[j] << " < " << var2_name << " < " << var2_bin_edges[j+1] << std::endl;
      
      for(int i=0; i<n_var1_bins; i++)
	{
	  std::cout << "processing subsample: " << var1_bin_edges[i] << " < " << var1_name << " < " << var1_bin_edges[i+1] << std::endl;
	    
	  //calculate the signal yield for a bin of pt and y.
	  if(var1_name == "pt")
	    signal_res = bin_mass_fit(*ws,channel,var1_bin_edges[i],var1_bin_edges[i+1], var2_bin_edges[j], var2_bin_edges[j+1]);
	  else
	    signal_res = bin_mass_fit(*ws,channel,var2_bin_edges[j],var2_bin_edges[j+1], var1_bin_edges[i], var1_bin_edges[i+1]);
	  
	  bin_size = var1_bin_size[i]*var2_bin_size[j];

	  yield_array[j][i] = signal_res->getVal()/bin_size;
	  errLo_array[j][i] = -signal_res->getAsymErrorLo()/bin_size;
	  errHi_array[j][i] = signal_res->getAsymErrorHi()/bin_size;
	  	    
	  //write to vector
	  signal_yield[i] = signal_res->getVal();

	}//end of var1 cicle
    }//end of var2 cicle
  
  signal_yield.Write("signal_yield");
  yield_file->Close();
  delete yield_file;

  TFile* f = new TFile(yield_file_name);
  TVectorD *yield = (TVectorD*)f->Get("signal_yield");
  
  yield->Print();
  delete f;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /*
    for(int j=0; j<n_var2_bins; j++)
    { 
      std::cout << "processing subsample: " << var2_bin_edges[j] << " < " << var2_name << " < " << var2_bin_edges[j+1] << std::endl;
      
      for(int i=0; i<n_var1_bins; i++)
	{
	  std::cout << "processing subsample: " << var1_bin_edges[i] << " < " << var1_name << " < " << var1_bin_edges[i+1] << std::endl;
	    
	  //calculate the signal yield for a bin of pt and y.
	  if(var1_name == "pt")
	    yield_syst_array[j][i] = bin_systematics(*ws, channel, var1_bin_edges[i], var1_bin_edges[i+1], var2_bin_edges[j], var2_bin_edges[j+1],signal_res->getVal(), data_selection_input_file, syst);
	  else
	    yield_syst_array[j][i] = bin_systematics(*ws, channel, var2_bin_edges[j], var2_bin_edges[j+1], var1_bin_edges[i], var1_bin_edges[i+1],signal_res->getVal(), data_selection_input_file, syst);
	  
	  bin_size = var1_bin_size[i]*var2_bin_size[j];

	  yield_syst_array[j][i] /= bin_size;
	    
	}//end of var1 cicle
    }//end of var2 cicle
  */
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //to calculate and plot pre-filter efficiency
  if(calculate_pre_filter_eff)
    {
      eff_name = "pre_eff";
      calculate_eff(eff_name, channel, n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, pre_eff_array[0], pre_eff_err_lo_array[0], pre_eff_err_hi_array[0]);
      
      for(int j=0; j<n_var2_bins; j++)
	{
	  plot_eff(eff_name, channel, n_var1_bins, var2_name, var2_bin_edges[j], var2_bin_edges[j+1], x_axis_name, b_title, var1_bin_centre, var1_bin_centre_lo, var1_bin_centre_hi, pre_eff_array[j], pre_eff_err_lo_array[j], pre_eff_err_hi_array[j]);
	}
    }
  
  //to calculate and plot reco efficiency
  if(calculate_reco_eff)
    {
      eff_name = "reco_eff";
      calculate_eff(eff_name, channel, n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, reco_eff_array[0], reco_eff_err_lo_array[0], reco_eff_err_hi_array[0]);

      for(int j=0; j<n_var2_bins; j++)
	{
	  plot_eff(eff_name, channel, n_var1_bins, var2_name, var2_bin_edges[j], var2_bin_edges[j+1], x_axis_name, b_title, var1_bin_centre, var1_bin_centre_lo, var1_bin_centre_hi, reco_eff_array[j], reco_eff_err_lo_array[j], reco_eff_err_hi_array[j]);
	}
    }

  //to calculate total efficiency
  if(calculate_pre_filter_eff && calculate_reco_eff)
    {
      eff_name = "total";
      
      for(int j=0; j<n_var2_bins; j++)
	{
	  std::cout << "processing bin: " << var2_bin_edges[j] << " < " << var2_name << " < " << var2_bin_edges[j+1] << std::endl;
	  for(int i=0; i<n_var1_bins; i++)
	    {
	      std::cout << "calculating total efficiency: " << var1_bin_edges[i] << " < " << var1_name << " < " << var1_bin_edges[i+1] << std::endl;

	      total_eff_array[j][i] = pre_eff_array[j][i] * reco_eff_array[j][i];
	      total_eff_err_lo_array[j][i] = total_eff_array[j][i] * sqrt(pow(pre_eff_err_lo_array[j][i]/pre_eff_array[j][i],2) + pow(reco_eff_err_lo_array[j][i]/reco_eff_array[j][i],2));
	      total_eff_err_hi_array[j][i] = total_eff_array[j][i] * sqrt(pow(pre_eff_err_hi_array[j][i]/pre_eff_array[j][i],2) + pow(reco_eff_err_hi_array[j][i]/reco_eff_array[j][i],2));
	    }

	  plot_eff(eff_name, channel, n_var1_bins, var2_name, var2_bin_edges[j], var2_bin_edges[j+1], x_axis_name, b_title, var1_bin_centre, var1_bin_centre_lo, var1_bin_centre_hi, total_eff_array[j], total_eff_err_lo_array[j], total_eff_err_hi_array[j]);
	}
    }

  //to calculate cross section
  for(int j=0; j<n_var2_bins; j++)
    {
      for(int i=0; i<n_var1_bins; i++)
	{
	  if(calculate_reco_eff && calculate_pre_filter_eff)
	    {
	      x_sec_array[j][i] = (yield_array[j][i] / (2  * reco_eff_array[j][i] * pre_eff_array[j][i] * LUMINOSITY * branch->getVal())) * (1e-9) * pow(10,j);
	      x_sec_errLo_array[j][i] = ((x_sec_array[j][i] * errLo_array[j][i]) / yield_array[j][i]);
	      x_sec_errHi_array[j][i] = ((x_sec_array[j][i] * errHi_array[j][i]) / yield_array[j][i]);

	      x_sec_syst_lo_array[j][i] = x_sec_array[j][i] * sqrt( pow(yield_syst_array[j][i]/yield_array[j][i],2) + pow(branch->getError()/branch->getVal(),2) + pow(pre_eff_err_lo_array[j][i]/pre_eff_array[j][i],2) + pow(reco_eff_err_lo_array[j][i]/reco_eff_array[j][i],2) + pow(0.04,2) );
	      x_sec_syst_hi_array[j][i] = x_sec_array[j][i] * sqrt( pow(yield_syst_array[j][i]/yield_array[j][i],2) + pow(branch->getError()/branch->getVal(),2) + pow(pre_eff_err_hi_array[j][i]/pre_eff_array[j][i],2) + pow(reco_eff_err_hi_array[j][i]/reco_eff_array[j][i],2) + pow(0.04,2) );
	    }
	  else
	    {
	      x_sec_array[j][i] = yield_array[j][i] * pow(10,j);
	      x_sec_errLo_array[j][i] = errLo_array[j][i] * pow(10,j);
	      x_sec_errHi_array[j][i] = errHi_array[j][i] * pow(10,j);
	      x_sec_syst_lo_array[j][i] = yield_syst_array[j][i] * pow(10,j);
	      x_sec_syst_hi_array[j][i] = yield_syst_array[j][i] * pow(10,j);
	    }
	}
    }

  //plot of the cross section
  TCanvas cz;
  
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
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(n_var1_bins, var1_bin_means, x_sec_array[j], var1_bin_edges_lo, var1_bin_edges_hi, x_sec_errLo_array[j], x_sec_errHi_array[j]);

      TString x_sec_title = b_title;
      if(calculate_pre_filter_eff && calculate_reco_eff)
        x_sec_title += " differential cross section";
      else
        x_sec_title += " signal yield";

      graph->SetTitle(x_sec_title);
      
      graph->SetMarkerColor(j+2);
      graph->SetMarkerSize(0.5);
      graph->SetMarkerStyle(20+j);
	  
      //draw this for the first var2 bin, or in case there is only one bin.
      if(j==0) 
	{
	  graph->GetXaxis()->SetTitle(x_axis_name);
	  graph->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
	  
	  //to set the range of the plot, it takes the min and max value of cross section.
	  if(n_var2_bins > 1)
            graph->GetYaxis()->SetRangeUser(0.1*x_sec_array[0][n_var1_bins-1], 10*x_sec_array[n_var2_bins-1][0]);
          else
	    graph->GetYaxis()->SetRangeUser(0.5*x_sec_array[0][n_var1_bins-1], 2*x_sec_array[0][0]);
	  
	  graph->Draw("ap same");
	}
      else //for the rest of the rapidity bins.
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
	  TGraphAsymmErrors* graph_syst = new TGraphAsymmErrors(n_var1_bins, var1_bin_means, x_sec_array[j], var1_bin_edges_lo, var1_bin_edges_hi, x_sec_syst_lo_array[j], x_sec_syst_hi_array[j]);
	  graph_syst->SetFillColor(j+2);
	  graph_syst->SetFillStyle(3001);
	  graph_syst->Draw("2 same");
	}
    }//end of var2 cicle
  
  leg->Draw("same");
  cz.Update();
  cz.SetLogy();

  TString systematic = "";

  if(syst)
    systematic = "_syst";

  if(calculate_pre_filter_eff && calculate_reco_eff)
    cz.SaveAs("signal_yield/x_sec_" + yield_sub_samples + "_bins_" + channel_to_ntuple_name(channel) + systematic + "_" + TString::Format(VERSION) + ".png");
  else
    cz.SaveAs("signal_yield/signal_yield_" + yield_sub_samples + "_bins_" + channel_to_ntuple_name(channel) + systematic + "_" + TString::Format(VERSION) + ".png");
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  //To show the values of cross section or signal yield and the errors at the end, like a table/
  //////////////////////////////////////////////////////////////////////////////////////////////

  //signal yield
  print_table("YIELDS", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, yield_array[0], errLo_array[0], errHi_array[0], yield_syst_array[0], yield_syst_array[0]);
  
  //pre-filter eff
  print_table("PRE-FILTER EFFICIENCY", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, pre_eff_array[0], pre_eff_err_lo_array[0], pre_eff_err_hi_array[0]);
  
  //reco eff
  print_table("RECO EFFICIENCY", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, reco_eff_array[0], reco_eff_err_lo_array[0], reco_eff_err_hi_array[0]);

  //total eff
  print_table("OVERALL EFFICIENCY", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, total_eff_array[0], total_eff_err_lo_array[0], total_eff_err_hi_array[0]);

  //cross section
  print_table("CROSS SECTION", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, x_sec_array[0], x_sec_errLo_array[0], x_sec_errHi_array[0], x_sec_syst_lo_array[0], x_sec_syst_hi_array[0]);

  //branching fraction
  std::cout << "branching fraction: " << branch->getVal() << std::endl;
  std::cout << std::endl;
  
}//end of signal_yield_new
