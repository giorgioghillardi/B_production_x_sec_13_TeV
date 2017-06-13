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

//input example: plot_x_sec --channel 1 --bins pt_y --eff 1 --syst 0
int main(int argc, char** argv)
{
  int channel = 1;
  std::string bins = "pt";
  int eff = 1;
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
	  convert >> bins;
	}
      if(argument == "--eff")
	{
	  convert << argv[++i];
	  convert >> eff;
	}
      if(argument == "--syst")
	{
	  convert << argv[++i];
	  convert >> syst;
	}
    }
  
  //to create the directories to save the .png files
  std::vector<std::string> dir_list;  
  dir_list.push_back("x_sec");
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
  double bin_size;

  double var1_bin_means[n_var2_bins][n_var1_bins];
  double var1_bin_edges_lo[n_var2_bins][n_var1_bins];
  double var1_bin_edges_hi[n_var2_bins][n_var1_bins];
  
  double x_sec[n_var2_bins][n_var1_bins];
  double x_sec_err_lo[n_var2_bins][n_var1_bins];
  double x_sec_err_hi[n_var2_bins][n_var1_bins];

  double totaleff[n_var2_bins][n_var1_bins];
  double totaleff_err_lo[n_var2_bins][n_var1_bins];
  double totaleff_err_hi[n_var2_bins][n_var1_bins];

  double yield_syst[n_var2_bins][n_var1_bins];
  double yield_syst_lo[n_var2_bins][n_var1_bins];
  double yield_syst_hi[n_var2_bins][n_var1_bins];
 
  double x_sec_syst_sqrt_lo[n_var2_bins][n_var1_bins];
  double x_sec_syst_sqrt_hi[n_var2_bins][n_var1_bins];
  
  double x_sec_syst_lo[n_var2_bins][n_var1_bins];
  double x_sec_syst_hi[n_var2_bins][n_var1_bins];
    
  RooRealVar* branch = branching_fraction(channel);

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
  
  //read data to calculate bin mean
  TString data_selection_input_file = TString::Format(BASE_DIR) + "selected_myloop_new_data_" + channel_to_ntuple_name(channel) + "_with_cuts.root";
  RooWorkspace* ws = new RooWorkspace("ws","Bmass");

  set_up_workspace_variables(*ws,channel);
  read_data(*ws, data_selection_input_file,channel);

  //calculate the mean values and bin size of var1 and var2
  for(int j=0; j<n_var2_bins; j++)
    {
      var2_bin_size[j] = var2_bins[j+1] - var2_bins[j];
      
      for(int i=0; i<n_var1_bins; i++)
	{
	  var1_bin_size[i] = var1_bins[i+1] - var1_bins[i];
	  
	  var1_bin_means[j][i] = var_mean_value(*ws, var1_name, var1_bins[i], var1_bins[i+1], var2_name, var2_bins[j], var2_bins[j+1]);
	  var1_bin_edges_lo[j][i] = var1_bin_means[j][i] - var1_bins[i];
	  var1_bin_edges_hi[j][i] = var1_bins[i+1] - var1_bin_means[j][i];
	}
    }
  
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
  
  TString vector = "";
  
  if(eff)
    vector = "x_sec";
  else
    vector = "yield";
  
  read_vector(measure, channel, vector, var1_name , var2_name, n_var1_bins, n_var2_bins, var1_bins, var2_bins, x_sec[0], x_sec_err_lo[0], x_sec_err_hi[0]);

  if(eff)
    read_vector(measure, channel, "totaleff", var1_name , var2_name, n_var1_bins, n_var2_bins, var1_bins, var2_bins, totaleff[0], totaleff_err_lo[0], totaleff_err_hi[0]);

  //read syst
  if(syst)
    {
      read_vector(measure, channel, "combined_syst", var1_name , var2_name, n_var1_bins, n_var2_bins, var1_bins, var2_bins, yield_syst[0], yield_syst_lo[0], yield_syst_hi[0]);
    }
  else
    for(int j=0; j<n_var2_bins; j++)
	{
	  for(int i=0; i<n_var1_bins; i++)
	    {
	      yield_syst[j][i] = 0;
	      yield_syst_lo[j][i] = 0;
	      yield_syst_hi[j][i] = 0;
	    }
	}
  
  
  //to correct the yield for bin size
  if(!eff)
    {
      for(int j=0; j<n_var2_bins; j++)
	{
	  for(int i=0; i<n_var1_bins; i++)
	    {
	      bin_size = var1_bin_size[i] * var2_bin_size[j];
	  
	      x_sec[j][i] = (x_sec[j][i] / bin_size) * pow(10,j);
	      x_sec_err_lo[j][i] = (x_sec_err_lo[j][i] / bin_size) * pow(10,j);
	      x_sec_err_hi[j][i] = (x_sec_err_hi[j][i] / bin_size) * pow(10,j);
	      
	      x_sec_syst_lo[j][i] = (yield_syst_lo[j][i] / bin_size) * pow(10,j);
	      x_sec_syst_hi[j][i] = (yield_syst_hi[j][i] / bin_size) * pow(10,j);
	    }
	}
    }
  else
    {
      for(int j=0; j<n_var2_bins; j++)
	{
	  for(int i=0; i<n_var1_bins; i++)
	    {
	      x_sec[j][i] *=  (1e-9) * pow(10,j);
	      x_sec_err_lo[j][i] *= (1e-9) * pow(10,j);
	      x_sec_err_hi[j][i] *= (1e-9) * pow(10,j);

	      x_sec_syst_sqrt_lo[j][i] = pow(branch->getError()/branch->getVal(),2) + pow(totaleff_err_lo[j][i]/totaleff[j][i],2) + pow(0.04,2) + pow(yield_syst[j][i],2);
	      x_sec_syst_sqrt_hi[j][i] = pow(branch->getError()/branch->getVal(),2) + pow(totaleff_err_hi[j][i]/totaleff[j][i],2) + pow(0.04,2) + pow(yield_syst[j][i],2);
	      
	      x_sec_syst_lo[j][i] = x_sec[j][i] * sqrt(x_sec_syst_sqrt_lo[j][i]);
	      x_sec_syst_hi[j][i] = x_sec[j][i] * sqrt(x_sec_syst_sqrt_hi[j][i]);
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
      
  //find x_sec max and min
  double x_sec_max = x_sec[n_var2_bins-1][0] + x_sec_syst_hi[n_var2_bins-1][0];
  double x_sec_min = x_sec[0][n_var1_bins-1] - x_sec_syst_lo[0][n_var1_bins-1];

  for(int i=0; i<n_var1_bins; i++)
    {
      if((x_sec[n_var2_bins-1][i]+x_sec_syst_hi[n_var2_bins-1][i] ) > x_sec_max)
	x_sec_max = x_sec[n_var2_bins-1][i]+x_sec_syst_hi[n_var2_bins-1][i];
    }

  for(int i=0; i<n_var1_bins; i++)
    {
      if((x_sec[0][i] - x_sec_syst_lo[0][i]) < x_sec_min)
	x_sec_min = x_sec[0][i] - x_sec_syst_lo[0][i];
    }

  for(int j=0; j<n_var2_bins; j++)
    {
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(n_var1_bins, var1_bin_means[j], x_sec[j], var1_bin_edges_lo[j], var1_bin_edges_hi[j], x_sec_err_lo[j], x_sec_err_hi[j]);

      TString x_sec_title = b_title;
      if(eff)
        x_sec_title += " differential cross section";
      else
        x_sec_title += " signal yield";

      graph->SetTitle(x_sec_title);
      
      graph->SetMarkerColor(j+2);
      graph->SetMarkerSize(0.2);
      graph->SetMarkerStyle(20+j);
	  
      //draw this for the first var2 bin, or in case there is only one bin.
      if(j==0) 
	{
	  graph->GetXaxis()->SetTitle(x_axis_name);
	  graph->GetYaxis()->SetTitle("d#sigma/dp_{T} [#mub/GeV]");
	  
	  //to set the range of the plot, it takes the min and max value of cross section.
	  if(n_var2_bins > 1)
            graph->GetYaxis()->SetRangeUser(0.1*x_sec_min, 10*x_sec_max);
          else
	    graph->GetYaxis()->SetRangeUser(0.5*x_sec_min, 2*x_sec_max);
	  
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
        label += TString::Format("%.2f",var2_bins[j]) + " < |y| < " + TString::Format("%.2f",var2_bins[j+1]);
      else
        label += TString::Format("%d",(int)var2_bins[j]) + " < pt < " + TString::Format("%d",(int)var2_bins[j+1]);
      
      leg->AddEntry(graph, label, "lp");
      
      //systematic errors
      if(eff)
        {
	  TGraphAsymmErrors* graph_syst = new TGraphAsymmErrors(n_var1_bins, var1_bin_means[j], x_sec[j], var1_bin_edges_lo[j], var1_bin_edges_hi[j], x_sec_syst_lo[j], x_sec_syst_hi[j]);
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

  if(eff)
    cz.SaveAs("x_sec/x_sec_" + bins + "_bins_" + channel_to_ntuple_name(channel) + systematic + "_" + TString::Format(VERSION) + ".png");
  else
    cz.SaveAs("x_sec/signal_yield_" + bins + "_bins_" + channel_to_ntuple_name(channel) + systematic + "_" + TString::Format(VERSION) + ".png");
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  //To show the values of cross section or signal yield and the errors at the end, like a table/
  //////////////////////////////////////////////////////////////////////////////////////////////

  //cross section
  print_table("CROSS SECTION", n_var1_bins, n_var2_bins, var1_name, var2_name, var1_bins, var2_bins, x_sec[0], x_sec_err_lo[0], x_sec_err_hi[0], x_sec_syst_lo[0], x_sec_syst_hi[0]);

  //branching fraction
  std::cout << "branching fraction: " << branch->getVal() << std::endl;
  std::cout << std::endl;
  
}//end of calculate_x_sec
