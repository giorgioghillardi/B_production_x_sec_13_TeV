#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

using namespace RooFit;

#define LUMINOSITY          2.71

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda
//-----------------------------------------------------------------

//input example: signal_yield_new --channel 1 --bins pt/y --preeff 1 --recoeff 1 --mcstudy 0 --syst 0
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

  //pt bins
  int nptbins=1;

  double ntkp_pt_bin_edges[]={9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 80, 90, 120};
  double ntkstar_pt_bin_edges[]={9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 90};
  double ntks_pt_bin_edges[]={0};
  double ntphi_pt_bin_edges[]={9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 90};
  double ntmix_pt_bin_edges[]={0};
  double ntlambda_pt_bin_edges[]={0};
  
  double total_pt_bin_edges[]={0, 400};
  double* pt_bin_edges=total_pt_bin_edges;

  //y bins
  int nybins=1;

  double ntkp_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntkstar_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntks_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntphi_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntmix_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntlambda_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  
  double total_y_bin_edges[]={0.0, 2.25};
  double* y_bin_edges=total_y_bin_edges;
    
  TString data_selection_input_file = TString::Format(BASE_DIR) + "selected_myloop_new_data_" + channel_to_ntuple_name(channel) + "_with_cuts.root";
  RooWorkspace* ws = new RooWorkspace("ws","Bmass");
  RooRealVar* signal_res; 

  TString pt_dist_directory="";
  TString y_dist_directory="";
  TString mass_fit_directory="";

  //set up mass, pt and y variables inside ws  
  set_up_workspace_variables(*ws,channel);
  //read data from the selected data file, and import it as a dataset into the workspace.
  read_data(*ws, data_selection_input_file,channel);

  ws->Print();

  if(yield_sub_samples=="full")
    {    
      pt_bin_edges = total_pt_bin_edges;
      nptbins = 1 ;
      y_bin_edges = total_y_bin_edges;
      nybins = 1;
    }  
  if(yield_sub_samples=="pt")
    {    
      switch (channel) {
      default:
      case 1:
	pt_bin_edges = ntkp_pt_bin_edges;
	nptbins = (sizeof(ntkp_pt_bin_edges) / sizeof(double)) - 1 ; //if pt_bin_edges is an empty array, then nptbins is equal to 0
	break;
      case 2:
	pt_bin_edges = ntkstar_pt_bin_edges;
	nptbins = (sizeof(ntkstar_pt_bin_edges) / sizeof(double)) - 1 ;
	break;
      case 3:
	pt_bin_edges = ntks_pt_bin_edges;
	nptbins = (sizeof(ntks_pt_bin_edges) / sizeof(double)) - 1 ;
	break;
      case 4:
	pt_bin_edges = ntphi_pt_bin_edges;
	nptbins = (sizeof(ntphi_pt_bin_edges) / sizeof(double)) - 1 ;
	break;
      case 5:
	pt_bin_edges = ntmix_pt_bin_edges;
	nptbins = (sizeof(ntmix_pt_bin_edges) / sizeof(double)) - 1 ;
	break;
      case 6:
	pt_bin_edges = ntlambda_pt_bin_edges;
	nptbins = (sizeof(ntlambda_pt_bin_edges) / sizeof(double)) - 1 ;
	break;
      }
    }
  if(yield_sub_samples=="pt/y")
    {
      switch (channel) {
      default:
      case 1:
	pt_bin_edges = ntkp_pt_bin_edges;
	nptbins = (sizeof(ntkp_pt_bin_edges) / sizeof(double)) - 1 ; //if pt_bin_edges is an empty array, then nptbins is equal to 0
	
	y_bin_edges = ntkp_y_bin_edges;
	nybins = (sizeof(ntkp_y_bin_edges) / sizeof(double)) - 1 ; //if y_bin_edges is an empty array, then nptbins is equal to 0
	break;
      case 2:
	pt_bin_edges = ntkstar_pt_bin_edges;
	nptbins = (sizeof(ntkstar_pt_bin_edges) / sizeof(double)) - 1 ;

	y_bin_edges = ntkstar_y_bin_edges;	
	nybins = (sizeof(ntkstar_y_bin_edges) / sizeof(double)) - 1 ;
	break;
      case 3:
	pt_bin_edges = ntks_pt_bin_edges;
	nptbins = (sizeof(ntks_pt_bin_edges) / sizeof(double)) - 1 ;
	    
	y_bin_edges = ntks_y_bin_edges;
	nybins = (sizeof(ntks_y_bin_edges) / sizeof(double)) - 1 ;
	break;
      case 4:
	pt_bin_edges = ntphi_pt_bin_edges;
	nptbins = (sizeof(ntphi_pt_bin_edges) / sizeof(double)) - 1 ;
	    
	y_bin_edges = ntphi_y_bin_edges;
	nybins = (sizeof(ntphi_y_bin_edges) / sizeof(double)) - 1 ;
	break;
      case 5:
	pt_bin_edges = ntmix_pt_bin_edges;
	nptbins = (sizeof(ntmix_pt_bin_edges) / sizeof(double)) - 1 ;

	y_bin_edges = ntmix_y_bin_edges;
	nybins = (sizeof(ntmix_y_bin_edges) / sizeof(double)) - 1 ;
	break;
      case 6:
	pt_bin_edges = ntlambda_pt_bin_edges;
	nptbins = (sizeof(ntlambda_pt_bin_edges) / sizeof(double)) - 1 ;
	    
	y_bin_edges = ntlambda_y_bin_edges;
	nybins = (sizeof(ntlambda_y_bin_edges) / sizeof(double)) - 1 ;
	break;
      }
    }

  double pt_bin_size[nptbins];
  double y_bin_size[nybins];
  double pt_bin_means[nptbins];
  double pt_bin_edges_Lo[nptbins];
  double pt_bin_edges_Hi[nptbins];

  double yield_array[nybins][nptbins];
  double errLo_array[nybins][nptbins];
  double errHi_array[nybins][nptbins];
  double yield_syst_array[nybins][nptbins];
 
  double pt_bin_centres_eff[nptbins];
  double pt_bin_edges_eff_Lo[nptbins];
  double pt_bin_edges_eff_Hi[nptbins];
      
  RooRealVar* pre_filter_eff;
      
  double pre_eff_array[nybins][nptbins];
  double pre_eff_err_lo_array[nybins][nptbins];
  double pre_eff_err_hi_array[nybins][nptbins];

  RooRealVar* reco_eff;
      
  double reco_eff_array[nybins][nptbins];
  double reco_eff_err_lo_array[nybins][nptbins];
  double reco_eff_err_hi_array[nybins][nptbins];

  double x_sec_array[nybins][nptbins];
  double x_sec_errLo_array[nybins][nptbins];
  double x_sec_errHi_array[nybins][nptbins];
  double x_sec_syst_lo_array[nybins][nptbins];
  double x_sec_syst_hi_array[nybins][nptbins];
      
  RooRealVar* branch = branching_fraction(channel);
  
  for(int c=0; c<nybins; c++)
    { 
      std::cout << "processing subsample: " << y_bin_edges[c] << " < |y| < " << y_bin_edges[c+1] << std::endl;

      //calculate the size of the y bins
      y_bin_size[c] = y_bin_edges[c+1]-y_bin_edges[c];

      for(int i=0; i<nptbins; i++)
	{
	  std::cout << "processing subsample: " << (int)pt_bin_edges[i] << " < pt < " << (int)pt_bin_edges[i+1] << std::endl;
	      
	  //calculate the size of the pt bins
	  pt_bin_size[i] = pt_bin_edges[i+1]-pt_bin_edges[i];
	      
	  //calculate the mean of the pt bin, to plot the cross section at this point
	  pt_bin_means[i] = var_mean_value(*ws,"pt",pt_bin_edges[i],pt_bin_edges[i+1]);
 
	  //calculate the new edges of the bin, since the centre is the mean. This is just for plotting.
	  pt_bin_edges_Lo[i] = pt_bin_means[i] - pt_bin_edges[i];
	  pt_bin_edges_Hi[i] = pt_bin_edges[i+1] - pt_bin_means[i];
	      
	  //calculate the signal yield for a bin of pt and y.
	  signal_res = bin_mass_fit(*ws, channel, pt_bin_edges[i], pt_bin_edges[i+1], y_bin_edges[c], y_bin_edges[c+1]);
	  
	  //MC study
	  if(mcstudy)
	    {
	      std::cout << "MC study of bin: " << (int)pt_bin_edges[i] << " < pt < " << (int)pt_bin_edges[i+1] << std::endl;

	      mc_study(*ws, channel, pt_bin_edges[i], pt_bin_edges[i+1], y_bin_edges[c], y_bin_edges[c+1]);
	    }
	  
	  if(yield_sub_samples=="full")
	    {
	      pt_bin_size[0] = 1; //force the bin size equal to one when we want to calculate full dataset cross section.
	      y_bin_size[0] = 1;
	    }
	  if(yield_sub_samples=="pt")
	    {
	      y_bin_size[0]=1;
	    }

	  yield_array[c][i] = (signal_res->getVal())/(pt_bin_size[i]*y_bin_size[c]);
	  errLo_array[c][i] = -(signal_res->getAsymErrorLo())/(pt_bin_size[i]*y_bin_size[c]);
	  errHi_array[c][i] = (signal_res->getAsymErrorHi())/(pt_bin_size[i]*y_bin_size[c]);
	  yield_syst_array[c][i] = bin_systematics(*ws, channel, pt_bin_edges[i], pt_bin_edges[i+1], y_bin_edges[c], y_bin_edges[c+1],signal_res->getVal(), data_selection_input_file, syst)/(pt_bin_size[i]*y_bin_size[c]);
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
		  
	      pt_bin_centres_eff[i] = pt_bin_edges[i] + (pt_bin_edges[i+1]-pt_bin_edges[i])/2;
	      pt_bin_edges_eff_Lo[i] = pt_bin_centres_eff[i] - pt_bin_edges[i];
	      pt_bin_edges_eff_Hi[i] = pt_bin_edges[i+1] - pt_bin_centres_eff[i];
		  
	      //pre-filter efficiency
	      pre_filter_eff = prefilter_efficiency(channel,pt_bin_edges[i],pt_bin_edges[i+1],y_bin_edges[c],y_bin_edges[c+1]);
	      //pre_filter_eff = prefilter_efficiency(channel,pt_bin_edges[i],pt_bin_edges[i+1],y_bin_edges[0],y_bin_edges[nybins]);
		 
	      pre_eff_array[c][i] = pre_filter_eff->getVal();
	      pre_eff_err_lo_array[c][i] = -(pre_filter_eff->getAsymErrorLo());
	      pre_eff_err_hi_array[c][i] = pre_filter_eff->getAsymErrorHi();
	    }
	  //plot of the pre-filter efficiency as a function of pT for each y bin
	  TCanvas ce;
	  TGraphAsymmErrors* graph_pre_eff = new TGraphAsymmErrors(nptbins, pt_bin_centres_eff, pre_eff_array[c], pt_bin_edges_eff_Lo, pt_bin_edges_eff_Hi, pre_eff_err_lo_array[c], pre_eff_err_hi_array[c]);
	  graph_pre_eff->SetTitle("pre-filter efficiency");
	  graph_pre_eff->SetMarkerColor(4);
	  graph_pre_eff->SetMarkerStyle(21);
	  graph_pre_eff->Draw("AP");
	  TString eff1_name = "";
	  eff1_name = "efficiencies/pre_filter_efficiency_" + channel_to_ntuple_name(channel) + "_" + TString::Format("y_from_%.1f_to_%.1f",y_bin_edges[c],y_bin_edges[c+1]) + ".png";
	  ce.SaveAs(eff1_name);
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
		  
	      pt_bin_centres_eff[i] = pt_bin_edges[i] + (pt_bin_edges[i+1]-pt_bin_edges[i])/2;
	      pt_bin_edges_eff_Lo[i] = pt_bin_centres_eff[i] - pt_bin_edges[i];
	      pt_bin_edges_eff_Hi[i] = pt_bin_edges[i+1] - pt_bin_centres_eff[i];
		  
	      //reco efficiency
	      reco_eff = reco_efficiency(channel,pt_bin_edges[i],pt_bin_edges[i+1],y_bin_edges[c],y_bin_edges[c+1]);
		 
	      reco_eff_array[c][i] = reco_eff->getVal();
	      reco_eff_err_lo_array[c][i] = -(reco_eff->getAsymErrorLo());
	      reco_eff_err_hi_array[c][i] = reco_eff->getAsymErrorHi();
	    }
	  //plot of the reco efficiency as a function of pT for each y bin
	  TCanvas cp;
	  TGraphAsymmErrors* graph_reco_eff = new TGraphAsymmErrors(nptbins, pt_bin_centres_eff, reco_eff_array[c], pt_bin_edges_eff_Lo, pt_bin_edges_eff_Hi, reco_eff_err_lo_array[c], reco_eff_err_hi_array[c]);
	  graph_reco_eff->SetTitle("reconstruction efficiency");
	  graph_reco_eff->SetMarkerColor(4);
	  graph_reco_eff->SetMarkerStyle(21);
	  graph_reco_eff->Draw("AP");
	  TString eff2_name = "";
	  eff2_name = "efficiencies/reco_efficiency_" + channel_to_ntuple_name(channel) + "_" + TString::Format("y_from_%.1f_to_%.1f",y_bin_edges[c],y_bin_edges[c+1]) + ".png";
	  cp.SaveAs(eff2_name);
	}
    }
      
  //to calculate cross section
  for(int j=0; j<nybins; j++)
    {
      for(int i=0; i<nptbins; i++)
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

  //To show the values of cross section or signal yield and the errors at the end, like a table
  std::cout << "CROSS SECTION" << std::endl;

  for(int c=0; c<nybins; c++)
    {
      std::cout << "BIN y: " << y_bin_edges[c] << " to " << y_bin_edges[c+1] << " : " << std::endl;
	  
      for(int i=0; i<nptbins; i++)
	{
	  std::cout << "BIN pt: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  x_sec_array[c][i] << " +" << x_sec_errHi_array[c][i] << " -"<< x_sec_errLo_array[c][i] << " -" << x_sec_syst_lo_array[c][i] << " +" << x_sec_syst_hi_array[c][i] << std::endl;
	}

      std::cout << std::endl;
    }

  if(calculate_pre_filter_eff)
    {
      //to show the values of pre-filter efficiency like a table.
      std::cout << "PRE FILTER EFFICIENCY" << std::endl;

      for(int c=0; c<nybins; c++)
	{
	  std::cout << "BIN y: " << y_bin_edges[c] << " to " << y_bin_edges[c+1] << " : " << std::endl;
	  
	  for(int i=0; i<nptbins; i++)
	    {
	      std::cout << "BIN pt: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  pre_eff_array[c][i] << " -" << pre_eff_err_lo_array[c][i] << " +"<< pre_eff_err_hi_array[c][i] << std::endl;
	    }

	  std::cout << std::endl;
	}
    }
  
  //to show the values of pre-filter efficiency like a table.
  if(calculate_reco_eff)
    {
      std::cout << "RECONSTRUCTION EFFICIENCY" << std::endl;
      
      for(int c=0; c<nybins; c++)
	{
	  std::cout << "BIN y: " << y_bin_edges[c] << " to " << y_bin_edges[c+1] << " : " << std::endl;
	  
	  for(int i=0; i<nptbins; i++)
	    {
	      std::cout << "BIN pt: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  reco_eff_array[c][i] << " -" << reco_eff_err_lo_array[c][i] << " +"<< reco_eff_err_hi_array[c][i] << std::endl;
	    }

	  std::cout << std::endl;
	}
    }
      
  //print the branching fraction
  std::cout << "branching fraction: " << branch->getVal() << std::endl;

  //plot of the cross section
  TCanvas cz;
  TPad *pad = new TPad("pad", "pad", 0.05, 0.05, 0.99, 0.99);
  pad->Draw();      
  TLegend *leg = new TLegend (0.65, 0.65, 0.85, 0.85);
      
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
      
  for(int i=0; i<nybins; i++)
    {
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(nptbins, pt_bin_means, x_sec_array[i], pt_bin_edges_Lo, pt_bin_edges_Hi, x_sec_errLo_array[i], x_sec_errHi_array[i]);
      graph->SetTitle("Differential cross section");
      graph->SetMarkerColor(i+2);
      graph->SetMarkerSize(0.5);
      graph->SetMarkerStyle(20+i);
	  
      //draw this for the first rapidity bin, or in case there is only one rapidity bin.
      if(i==0) 
	{
	  graph->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
	  graph->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
	  //to set the range of the plot, it takes the min and max value of cross section.
	  graph->GetYaxis()->SetRangeUser(1e-6/*0.1*x_sec_array[0][0]*/,10*x_sec_array[nybins-1][0]);
	  graph->Draw("ap same");
	  //print the rapidity bin
	  leg->AddEntry(graph, TString::Format("%.1f < |y| < %.1f",y_bin_edges[i],y_bin_edges[i+1]), "lp");
	}
      else //for the rest of the rapidity bins.
	{     
	  graph->Draw("p same");
	  leg->AddEntry(graph, TString::Format("(#times 10^{%d}) %.1f < |y| < %.1f",i,y_bin_edges[i],y_bin_edges[i+1]), "lp");
	}
      
      //systematic errors
      TGraphAsymmErrors* graph_syst = new TGraphAsymmErrors(nptbins, pt_bin_means, x_sec_array[i], pt_bin_edges_Lo, pt_bin_edges_Hi, x_sec_syst_lo_array[i], x_sec_syst_hi_array[i]);
      graph_syst->SetFillColor(i+2);
      graph_syst->SetFillStyle(3001);
      graph_syst->Draw("2 same");
      
    }//end of the ybins for
  
  if(yield_sub_samples=="pt/y") yield_sub_samples="pt_y";
      
  leg->Draw("same");
  cz.Update();
  cz.SetLogy();
      
  TString systematic = "";
      
  if(syst)
    systematic = "_syst";
      
  if(calculate_pre_filter_eff && calculate_reco_eff)
    {
      cz.SaveAs("signal_yield/x_sec_" + yield_sub_samples + "_bins_" + channel_to_ntuple_name(channel) + systematic + "_" + TString::Format(VERSION) + ".png");
      cz.SaveAs("signal_yield/x_sec_" + yield_sub_samples + "_bins_" + channel_to_ntuple_name(channel) + systematic + "_" + TString::Format(VERSION) + ".C");
    }
  else
    {
      cz.SaveAs("signal_yield/signal_yield_" + yield_sub_samples + "_bins_" + channel_to_ntuple_name(channel) + systematic + "_" + TString::Format(VERSION) + ".png");
      cz.SaveAs("signal_yield/signal_yield_" + yield_sub_samples + "_bins_" + channel_to_ntuple_name(channel) + systematic + "_" + TString::Format(VERSION) + ".C");
    }
  //    }//end of else
}//end of signal_yield_new
