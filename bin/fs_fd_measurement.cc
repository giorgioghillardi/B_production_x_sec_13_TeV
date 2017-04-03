#include "TF1.h"
#include "TPaveStats.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K* -> J/psi K+pi-
// channel = 4: Bs -> J/psi phi -> J/psi K+K-

//input example: fs_fd_ratio --ratio fs_fu --bins pt_y --preeff 1 --recoeff 1 --poly 1 --syst 0 --mcstudy 0
int main(int argc, char** argv)
{
  std::string ratio = "fs_fu";
  std::string yield_sub_samples = "full";
  int calculate_pre_filter_eff = 0;
  int calculate_reco_eff = 0;
  int poly = 1;
  int syst = 0;
  int mcstudy = 0;

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
  double pt_bins[] = {9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 90};
  double total_pt_bin_edges[]={10, 150};
  int number_of_pt_bins = (sizeof(pt_bins) / sizeof(double)) - 1; //if pt_bins is an empty array, then number_of_pt_bins is equal to 0
  
  //y bins
  double y_bins[] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2, 2.25};
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
  
  double var1_bin_centre[n_var1_bins];
  double var1_bin_centre_lo[n_var1_bins];
  double var1_bin_centre_hi[n_var1_bins];
  
  double yield_array[2][n_var2_bins][n_var1_bins];
  double errLo_array[2][n_var2_bins][n_var1_bins];
  double errHi_array[2][n_var2_bins][n_var1_bins];
  double yield_syst_array[2][n_var2_bins][n_var1_bins];
        
  double pre_eff_array[2][n_var2_bins][n_var1_bins];
  double pre_eff_err_lo_array[2][n_var2_bins][n_var1_bins];
  double pre_eff_err_hi_array[2][n_var2_bins][n_var1_bins];
      
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
  double fs_fd_syst_lo_sqrt[n_var2_bins][n_var1_bins];
  double fs_fd_syst_hi_sqrt[n_var2_bins][n_var1_bins];
  double fs_fd_syst_lo_array[n_var2_bins][n_var1_bins];
  double fs_fd_syst_hi_array[n_var2_bins][n_var1_bins];
  
  double b_fraction[2];
  double b_fraction_err[2];

  //used later on in the eff plots
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
  
  for(int i=0; i<n_var1_bins; i++)
    {
      var1_bin_centre[i] = var1_bin_edges[i] + (var1_bin_edges[i+1]-var1_bin_edges[i])/2;
      var1_bin_centre_lo[i] = var1_bin_centre[i] - var1_bin_edges[i];
      var1_bin_centre_hi[i] = var1_bin_edges[i+1] - var1_bin_centre[i];
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

      //to create the directories to save the .png files
      std::vector<std::string> dir_list;
      dir_list.push_back(static_cast<const char*>("mass_fits/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION)));
      dir_list.push_back("fs_fd");
      dir_list.push_back("efficiencies");
      
      create_dir(dir_list); 
      
      //------------data input---------------------
      TString data_selection_input_file = TString::Format(BASE_DIR) + "selected_myloop_new_data_" + channel_to_ntuple_name(channel) + "_with_cuts.root";
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
	      if(var1_name == "pt") //to make sure that the function bin_mass_fit gets the right input.
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
		  std::cout << "MC study of bin: " << var1_bin_edges[i] << " < " << var1_name << " < " << var1_bin_edges[i+1] << std::endl;

		  mc_study(*ws,channel,var1_bin_edges[i],var1_bin_edges[i+1], var2_bin_edges[j], var2_bin_edges[j+1]);
		}
	    }
	}
      
      //to calculate pre-filter efficiency
      if(calculate_pre_filter_eff)
	{
	  eff_name = "pre_eff";
	  calculate_eff(eff_name, channel, n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, pre_eff_array[ch][0], pre_eff_err_lo_array[ch][0], pre_eff_err_hi_array[ch][0]);

	  for(int j=0; j<n_var2_bins; j++)
	    {
	      plot_eff(eff_name, channel, n_var1_bins, var2_name, var2_bin_edges[j], var2_bin_edges[j+1], x_axis_name, b_title, var1_bin_centre, var1_bin_centre_lo, var1_bin_centre_hi, pre_eff_array[ch][j], pre_eff_err_lo_array[ch][j], pre_eff_err_hi_array[ch][j]);
	    }
	}

      //to calculate reconstruction efficiency
      if(calculate_reco_eff)
	{
	  eff_name = "reco_eff";
	  calculate_eff(eff_name, channel, n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, reco_eff_array[ch][0], reco_eff_err_lo_array[ch][0], reco_eff_err_hi_array[ch][0]);

	  for(int j=0; j<n_var2_bins; j++)
	    {
	      plot_eff(eff_name, channel, n_var1_bins, var2_name, var2_bin_edges[j], var2_bin_edges[j+1], x_axis_name, b_title, var1_bin_centre, var1_bin_centre_lo, var1_bin_centre_hi, reco_eff_array[ch][j], reco_eff_err_lo_array[ch][j], reco_eff_err_hi_array[ch][j]);
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

		  total_eff_array[ch][j][i] = pre_eff_array[ch][j][i] * reco_eff_array[ch][j][i];
		  total_eff_err_lo_array[ch][j][i] = total_eff_array[ch][j][i] * sqrt(pow(pre_eff_err_lo_array[ch][j][i]/pre_eff_array[ch][j][i],2) + pow(reco_eff_err_lo_array[ch][j][i]/reco_eff_array[ch][j][i],2));
		  total_eff_err_hi_array[ch][j][i] = total_eff_array[ch][j][i] * sqrt(pow(pre_eff_err_hi_array[ch][j][i]/pre_eff_array[ch][j][i],2) + pow(reco_eff_err_hi_array[ch][j][i]/reco_eff_array[ch][j][i],2));
		}

	      plot_eff(eff_name, channel, n_var1_bins, var2_name, var2_bin_edges[j], var2_bin_edges[j+1], x_axis_name, b_title, var1_bin_centre, var1_bin_centre_lo, var1_bin_centre_hi, total_eff_array[ch][j], total_eff_err_lo_array[ch][j], total_eff_err_hi_array[ch][j]);
	    }
	}

      //to calculate efficiency ratio
      if(calculate_pre_filter_eff && calculate_reco_eff)
	{
	  eff_name = "ratio";

	  for(int j=0; j<n_var2_bins; j++)
	    { 
	      std::cout << "processing bin: " << var2_bin_edges[j] << " < " << var2_name << " < " << var2_bin_edges[j+1] << std::endl;
	      for(int i=0; i<n_var1_bins; i++)
		{
		  std::cout << "calculating efficiency ratio: " << var1_bin_edges[i] << " < " << var1_name << " < " << var1_bin_edges[i+1] << std::endl;
		   
		  ratio_eff_array[j][i] = total_eff_array[0][j][i] / total_eff_array[1][j][i];
		  ratio_eff_err_lo_array[j][i] = ratio_eff_array[j][i] * sqrt(pow(total_eff_err_lo_array[0][j][i]/total_eff_array[0][j][i],2) + pow(total_eff_err_lo_array[1][j][i]/total_eff_array[1][j][i],2));
		  ratio_eff_err_hi_array[j][i] = ratio_eff_array[j][i] * sqrt(pow(total_eff_err_hi_array[0][j][i]/total_eff_array[0][j][i],2) + pow(total_eff_err_hi_array[1][j][i]/total_eff_array[1][j][i],2));
		}

	      plot_eff(eff_name, channel, n_var1_bins, var2_name, var2_bin_edges[j], var2_bin_edges[j+1], x_axis_name, ratio, var1_bin_centre, var1_bin_centre_lo, var1_bin_centre_hi, ratio_eff_array[j], ratio_eff_err_lo_array[j], ratio_eff_err_hi_array[j]);
	    }
	}
      
      delete ws;
    }//end of channel cicle
  
  //to calculate the fragmentation fraction ratio
  for(int j=0; j<n_var2_bins; j++)
    {
      for(int i=0; i<n_var1_bins; i++)
	{
	  fs_fd_array[j][i] = (yield_array[1][j][i]/yield_array[0][j][i]) * pow(10,j);
	  
	  fs_fd_syst_lo_sqrt[j][i] = pow(yield_syst_array[0][j][i]/yield_array[0][j][i],2) + pow(yield_syst_array[1][j][i]/yield_array[1][j][i],2);
	  fs_fd_syst_hi_sqrt[j][i] = fs_fd_syst_lo_sqrt[j][i]; //the syst errors without the efficiencies are symmetrical.
	  
	  if(calculate_reco_eff && calculate_pre_filter_eff) //yield ratio and syst corrected by efficiency ratio and BF ratio.
	    {
	      fs_fd_array[j][i] *= ratio_eff_array[j][i] * (b_fraction[0]/b_fraction[1]);
	      fs_fd_syst_lo_sqrt[j][i]  += pow(ratio_eff_err_lo_array[j][i]/ratio_eff_array[j][i],2) + pow(b_fraction_err[0]/b_fraction[0],2) + pow(b_fraction_err[1]/b_fraction[1],2);
	      fs_fd_syst_hi_sqrt[j][i] += pow(ratio_eff_err_hi_array[j][i]/ratio_eff_array[j][i],2) + pow(b_fraction_err[0]/b_fraction[0],2) + pow(b_fraction_err[1]/b_fraction[1],2);
	    }
	  
	  fs_fd_syst_lo_array[j][i]  = fs_fd_array[j][i] * sqrt(fs_fd_syst_lo_sqrt[j][i]);
	  fs_fd_syst_hi_array[j][i]  = fs_fd_array[j][i] * sqrt(fs_fd_syst_hi_sqrt[j][i]);

	  fs_fd_errLo_array[j][i] = fs_fd_array[j][i] * sqrt( pow(errLo_array[0][j][i]/yield_array[0][j][i],2) + pow(errLo_array[1][j][i]/yield_array[1][j][i],2) );
	  fs_fd_errHi_array[j][i] = fs_fd_array[j][i] * sqrt( pow(errHi_array[0][j][i]/yield_array[0][j][i],2) + pow(errHi_array[1][j][i]/yield_array[1][j][i],2) );
	}
    }
  
  /////////////////////////////////////
  //plot fragmentation fraction ratio//
  ////////////////////////////////////
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
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(n_var1_bins, var1_bin_centre, fs_fd_array[j], var1_bin_centre_lo, var1_bin_centre_hi, fs_fd_errLo_array[j], fs_fd_errHi_array[j]);
      
      TString fragmentation_title = ratio;
      if(calculate_pre_filter_eff && calculate_reco_eff)
	fragmentation_title += " fragmentation fraction ratio";
      else
	fragmentation_title += " yield ratio";

      graph->SetTitle(fragmentation_title);
      graph->SetMarkerColor(2+j);
      graph->SetMarkerSize(0.5);
      graph->SetMarkerStyle(20+j);      
 
      //draw this for the first bin of var2, or in case there is only one bin.
      if(j==0)
	{
	  graph->GetXaxis()->SetTitle(x_axis_name);
	  
	  //to set the range of the plot.
	  if(n_var2_bins > 1)
	    graph->GetYaxis()->SetRangeUser(1e-2, 10*fs_fd_array[n_var2_bins-1][0]);
	  else
	    {
	      graph->GetYaxis()->SetRangeUser(0, 2*fs_fd_array[n_var2_bins-1][0]);
	      
	      if(poly)
		{
		  //fit the ratio with a polynomial function
		  graph->Fit("pol1","W","");
		  graph->GetFunction("pol1")->SetLineColor(4);
		  gStyle->SetOptStat(1111);
		  gStyle->SetOptFit(11); //show p0 and p1 parameters from the pol1 fit
		  
		  gStyle->SetStatX(0.9);
		  gStyle->SetStatY(0.9);
		}
	      
	      cz.Update();
	    }
	  
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
	  TGraphAsymmErrors* graph_syst = new TGraphAsymmErrors(n_var1_bins, var1_bin_centre, fs_fd_array[j], var1_bin_centre_lo, var1_bin_centre_hi, fs_fd_syst_lo_array[j], fs_fd_syst_hi_array[j]);
	  graph_syst->SetFillColor(j+2);
	  graph_syst->SetFillStyle(3001);
	  graph_syst->Draw("2 same");
	}
    }//end of var2 cicle
  
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

  for(int ch=0; ch<2; ch++)
    {
      //signal yield
      print_table("YIELDS", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, yield_array[ch][0], errLo_array[ch][0], errHi_array[ch][0], yield_syst_array[ch][0], yield_syst_array[ch][0]);

      //pre-filter eff
      print_table("PRE-FILTER EFFICIENCY", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, pre_eff_array[ch][0], pre_eff_err_lo_array[ch][0], pre_eff_err_hi_array[ch][0]);

      //reco eff
      print_table("RECO EFFICIENCY", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, reco_eff_array[ch][0], reco_eff_err_lo_array[ch][0], reco_eff_err_hi_array[ch][0]);
  
      //total eff
      print_table("OVERALL EFFICIENCY", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, total_eff_array[ch][0], total_eff_err_lo_array[ch][0], total_eff_err_hi_array[ch][0]);

      //branching fraction
      std::cout << "BRANCHING FRACTION : " << b_fraction[ch] << " +- " << b_fraction_err[ch] << std::endl;
      std::cout << std::endl;
    }

  //ratio eff
  print_table("EFFICIENCY RATIO", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, ratio_eff_array[0], ratio_eff_err_lo_array[0], ratio_eff_err_hi_array[0]);

  //Fragmentation fraction
  print_table("FRAGMENTATION FRACTION", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bin_edges, var2_bin_edges, fs_fd_array[0], fs_fd_errLo_array[0], fs_fd_errHi_array[0], fs_fd_syst_lo_array[0], fs_fd_syst_hi_array[0]);

}//end of main
