#include <RooHist.h>
#include <TSystem.h>
#include <sstream>
#include <vector>
#include <string>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TString.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooWorkspace.h>
#include <RooAddPdf.h>
#include <RooGenericPdf.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <RooMCStudy.h>
#include <RooPlot.h>
#include <RooPlotable.h>
#include <RooThresholdCategory.h>
#include <Roo1DTable.h>
#include "TMath.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/myloop.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/plotDressing.h"
#include <TLegend.h>
using namespace RooFit;

// General fitting options
#define NUMBER_OF_CPU       1

#define VERSION             "v7"
#define BASE_DIR            "/lstore/cms/brunogal/input_for_B_production_x_sec_13_TeV/"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda

void create_dir(std::vector<std::string> list);
void plot_pt_dist(RooWorkspace& w, int channel, TString directory);
void plot_mass_fit(RooWorkspace& w, int channel, TString directory,int pt_high, int pt_low);
void plot_mass_fit(RooWorkspace& w, int channel, TString directory);
RooRealVar* bin_mass_fit(RooWorkspace& w, int channel, double pt_min, double pt_max, double y_min, double y_max);
double pt_bin_mean(RooWorkspace& w, double pt_min, double pt_max);
RooRealVar* pre_filter_efficiency(int channel, double pt_min, double pt_max);

void build_pdf(RooWorkspace& w, int channel);
void read_data(RooWorkspace& w, TString filename,int channel);
void read_data_cut(RooWorkspace& w, RooAbsData* data);
void set_up_workspace_variables(RooWorkspace& w, int channel);

TString channel_to_ntuple_name(int channel);
TString channel_to_xaxis_title(int channel);
int channel_to_nbins(int channel);

//input example: signal_yield_new --channel 1 --bins pt/y --eff 1 --mc 1
int main(int argc, char** argv)
{
  int channel = 0;
  std::string yield_sub_samples = "0";
  int calculate_efficiency = 0;
  int mcstudy = 0;
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
      if(argument == "--eff")
	{
	  convert << argv[++i];
	  convert >> calculate_efficiency;
	} 
      if(argument == "--mc")
	{
	  convert << argv[++i];
	  convert >> mcstudy;
	}
    }


  if(channel==0)
    {
      std::cout << "No channel was provided as input. Please use --channel. Example: signal_yield_new --channel 1" << std::endl;
      return 0;
    }
  
  //to create the directories to save the .png files
  std::vector<std::string> dir_list;
  dir_list.push_back("full_dataset_mass_fit");
  dir_list.push_back("full_dataset_mass_pt_histo");
  dir_list.push_back(static_cast<const char*>("pt_bin_mass_fit/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION)));
  dir_list.push_back("signal_yield");
  
  create_dir(dir_list);

  //pt bins
  double ntkp_pt_bin_edges[]={10,20,30,40,50,60,70,80,90,100,120,150};
  double ntkstar_pt_bin_edges[]={9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 90};
  double ntks_pt_bin_edges[]={10,20,30,40,50,60,70};
  double ntphi_pt_bin_edges[]={9, 13, 16, 20, 25, 30, 35, 42, 50, 60, 70, 90};
  double ntmix_pt_bin_edges[]={0};
  double ntlambda_pt_bin_edges[]={0};
  double total_pt_bin_edges[]={0, 400};
  double* pt_bin_edges=total_pt_bin_edges;
  
  int nptbins=0;

  //y bins
  double ntkp_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntkstar_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntks_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntphi_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntmix_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double ntlambda_y_bin_edges[]={0.0, 0.5, 1.0, 1.5, 2.25};
  double total_y_bin_edges[]={0.0, 4.0};
  double* y_bin_edges=total_y_bin_edges;
  
  int nybins=0;

  
  TString data_selection_input_file = "selected_data_" + channel_to_ntuple_name(channel) + ".root";
  
  RooWorkspace* ws = new RooWorkspace("ws","Bmass");
  RooAbsData* data;
  RooAbsPdf* model;
  RooFitResult* fit_res;
  RooRealVar* signal_res;
  

  TString pt_dist_directory="";
  TString y_dist_directory="";
  TString mass_fit_directory="";

  //set up mass, pt and y variables inside ws  
  set_up_workspace_variables(*ws,channel);

  //read data from the selected data file, and import it as a dataset into the workspace.
  read_data(*ws, data_selection_input_file,channel);
  ws->Print();
  
  if(yield_sub_samples=="0") //mass fit and plot the full dataset
    { 
      
      //build the pdf for the channel selected above, it uses the dataset which is saved in ws. need to change the dataset to change the pdf.
      RooRealVar* mass = ws->var("mass"); 
      build_pdf(*ws,channel);     
      
      data = ws->data("data");
      model = ws->pdf("model");     
      
      model->fitTo(*data,Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));

      signal_res = ws->var("n_signal");
      /* RooRealVar* mean = ws->var("m_mean");
	 RooRealVar* sigma1 = ws->var("m_sigma1");
	 RooRealVar* sigma2 = ws->var("m_sigma2");
	 RooRealVar* lambda = ws->var("m_exp");*/
  
      std::cout <<"SIGNAL: "<< signal_res->getVal() << " " << signal_res->getAsymErrorLo() << " +" << signal_res->getAsymErrorHi() << std::endl;
      
      mass_fit_directory = "full_dataset_mass_fit/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION);
      pt_dist_directory = "full_dataset_mass_pt_histo/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION);

      plot_mass_fit(*ws,channel,mass_fit_directory);
      plot_pt_dist(*ws,channel,pt_dist_directory);
    

      if(mcstudy)
	{

	  RooMCStudy* mctoy = new RooMCStudy (*model, *model, *mass, "", "mhv"); 
	  mctoy->generateAndFit(1000, data->sumEntries());

	  RooPlot* f_pull_signal = mctoy->plotPull(*signal_res, FitGauss(kTRUE));
	  RooPlot* f_param_signal = mctoy->plotParam(*signal_res);
	  RooPlot* f_error_signal = mctoy->plotError(*signal_res);
	  RooPlot* f_nll = mctoy->plotNLL();
	  //RooPlot* f_pull_sigma1 = mctoy->plotPull(sigma1, FitGauss(kTRUE));
	  //RooPlot* f_pull_mean = mctoy->plotPull(mean, FitGauss(kTRUE));
	  //RooPlot* f_pull_mean = mctoy->plotPull(mean, FitGauss(kTRUE));
     
	  TCanvas c;
	  f_param_signal->Draw();
	  c.SaveAs("param_signal.png");

	  TCanvas c2;
	  f_error_signal->Draw();
	  c2.SaveAs("error_signal.png");
	  
	  TCanvas c3;
	  f_pull_signal->Draw();
	  c3.SaveAs("pull_signal.png");

	  TCanvas c4;
	  f_nll->Draw();
	  c4.SaveAs("nll_signal.png");
	}

    }
  else
    {
      if(yield_sub_samples=="pt")
	{
	  y_bin_edges = total_y_bin_edges;
	  nybins=1;

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
      else if(yield_sub_samples=="y")
	{
	  pt_bin_edges = total_pt_bin_edges;
	  nptbins=1;

	  switch (channel) {
	  default:
	  case 1:
	    y_bin_edges = ntkp_y_bin_edges;
	    nybins = (sizeof(ntkp_y_bin_edges) / sizeof(double)) - 1 ; //if y_bin_edges is an empty array, then nptbins is equal to 0
	    break;
	  case 2:
	    y_bin_edges = ntkstar_y_bin_edges;	
	    nybins = (sizeof(ntkstar_y_bin_edges) / sizeof(double)) - 1 ;
	    break;
	  case 3:
	    y_bin_edges = ntks_y_bin_edges;
	    nybins = (sizeof(ntks_y_bin_edges) / sizeof(double)) - 1 ;
	    break;
	  case 4:
	    y_bin_edges = ntphi_y_bin_edges;
	    nybins = (sizeof(ntphi_y_bin_edges) / sizeof(double)) - 1 ;
	    break;
	  case 5:
	    y_bin_edges = ntmix_y_bin_edges;
	    nybins = (sizeof(ntmix_y_bin_edges) / sizeof(double)) - 1 ;
	    break;
	  case 6:
	    y_bin_edges = ntlambda_y_bin_edges;
	    nybins = (sizeof(ntlambda_y_bin_edges) / sizeof(double)) - 1 ;
	    break;
	  }

	}
      else if(yield_sub_samples=="pt/y")
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

      std::cout << "nybins: " << nybins << std::endl;
      std::cout << "y_bin_edges[0]: " << y_bin_edges[0] << std::endl;
      std::cout << "y_bin_edges[1]: " << y_bin_edges[1] << std::endl;

      double pt_bin_size[nptbins];
      double pt_bin_means[nptbins];
      double pt_bin_edges_Lo[nptbins];
      double pt_bin_edges_Hi[nptbins];

      double yield_array[nybins][nptbins];
      double errLo_array[nybins][nptbins];
      double errHi_array[nybins][nptbins];
 
      double pt_bin_centres_eff[nptbins];
      double pt_bin_edges_eff_Lo[nptbins];
      double pt_bin_edges_eff_Hi[nptbins];
      double eff_array[nptbins];
      double effLo_array[nptbins];
      double effHi_array[nptbins];

      RooRealVar* pre_filter_eff;
      
      for(int c=0; c<nybins; c++)
	{ 
	  std::cout << "processing subsample: " << y_bin_edges[c] << " < |y| < " << y_bin_edges[c+1] << std::endl;

	  for(int i=0; i<nptbins; i++)
	    {
	      std::cout << "processing subsample: " << (int)pt_bin_edges[i] << " < pt < " << (int)pt_bin_edges[i+1] << std::endl;
	      
	      pt_bin_size[i] = pt_bin_edges[i+1]-pt_bin_edges[i];
	      std::cout << pt_bin_size[i] << std::endl;
	      
	      pt_bin_means[i] = pt_bin_mean(*ws,pt_bin_edges[i],pt_bin_edges[i+1]);
	      pt_bin_edges_Lo[i] = pt_bin_means[i] - pt_bin_edges[i];
	      pt_bin_edges_Hi[i] = pt_bin_edges[i+1] - pt_bin_means[i];

	      std::cout << pt_bin_edges_Lo[i] << pt_bin_edges_Hi[i] << std::endl;
	      
	      signal_res = bin_mass_fit(*ws,channel,pt_bin_edges[i],pt_bin_edges[i+1], y_bin_edges[c], y_bin_edges[c+1]);
	      
	      yield_array[c][i] = (signal_res->getVal())/(pt_bin_size[i]*(y_bin_edges[c+1]-y_bin_edges[c]))*pow(10,nybins-c);
	      errLo_array[c][i] = -(signal_res->getAsymErrorLo())/(pt_bin_size[i]*(y_bin_edges[c+1]-y_bin_edges[c]))*pow(10,nybins-c);
	      errHi_array[c][i] = (signal_res->getAsymErrorHi())/(pt_bin_size[i]*(y_bin_edges[c+1]-y_bin_edges[c]))*pow(10,nybins-c);
	    }
	}

      //to show the values of signal_yield and the errors at the end, like a table
      for(int c=0; c<nybins; c++)
	{
	  std::cout << "BIN y: " << y_bin_edges[c] << " to " << y_bin_edges[c+1] << " : " << std::endl;
	  for(int i=0; i<nptbins; i++)
	    {
	      std::cout << "BIN pt: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  yield_array[c][i] << " +" << errHi_array[c][i] << " -"<< errLo_array[c][i] << std::endl;
	    }
	  std::cout << std::endl;
	}

      if(calculate_efficiency)
	{
	  for(int i=0; i<nptbins; i++)
	    {
	      std::cout << "calculating pre-filter efficiency: " << (int)pt_bin_edges[i] << " < pt < " << (int)pt_bin_edges[i+1] << std::endl;
	 
	      pt_bin_centres_eff[i] = pt_bin_edges[i] + (pt_bin_edges[i+1]-pt_bin_edges[i])/2;
	      pt_bin_edges_eff_Lo[i] = pt_bin_centres_eff[i] - pt_bin_edges[i];
	      pt_bin_edges_eff_Hi[i] = pt_bin_edges[i+1] - pt_bin_centres_eff[i];
	 
	      pre_filter_eff = pre_filter_efficiency(channel,pt_bin_edges[i],pt_bin_edges[i+1]);
	 
	      eff_array[i] = pre_filter_eff->getVal();
	      effLo_array[i] = -pre_filter_eff->getAsymErrorLo();
	      effHi_array[i] = pre_filter_eff->getAsymErrorHi();
	    }
	  //plot of the pre-filter efficiency as a function of pT
	  TCanvas ce;
	  TGraphAsymmErrors* graph_eff = new TGraphAsymmErrors(nptbins, pt_bin_centres_eff, eff_array, pt_bin_edges_eff_Lo, pt_bin_edges_eff_Hi,effLo_array,effHi_array);
	  graph_eff->SetTitle("pre filter efficiency");
	  graph_eff->SetMarkerColor(4);
	  graph_eff->SetMarkerStyle(21);
	  graph_eff->Draw("AP");
    
	  ce.SaveAs("pre_filter_efficiency_err.png");
	}
      
      //plot of the signal_yield as a function of pt, in the future should be the x-sec corrected by efficiency and other factors
      TCanvas cz;
      TPad *pad = new TPad("pad", "pad", 0.05, 0.05, 0.99, 0.99);
      //      pad->cd();
      pad->Draw();

      TH1D* empty = new TH1D("Raw Signal Yield in p_{T} Bins", "Raw Signal Yield in p_{T} Bins; p_{T} [GeV]; Signal Yield", nptbins, 0, 100);
      /*      empty->SetMinimum(1);

      if(yield_sub_samples=="pt"||yield_sub_samples=="y")
	empty->SetMaximum(4e10);
      */ 
      TLegend *leg = new TLegend (0.65, 0.65, 0.85, 0.85);

      TGraphAsymmErrors* graph = new TGraphAsymmErrors(nptbins, pt_bin_means, yield_array[0], pt_bin_edges_Lo, pt_bin_edges_Hi, errLo_array[0], errHi_array[0]);
      graph->SetTitle("Raw signal yield in Pt bins");
      graph->SetMarkerColor(1);
      graph->SetMarkerSize(0.5);
      graph->SetMarkerStyle(20);
      empty->SetMinimum(graph->GetMinimum());
      if(nybins<=1)
	{
	  std::cout << "ENTRAS??" << std::endl << "graph->GetMaximum(): " << graph->GetMaximum() << std::endl << "graph->GetMinimum(): " << graph->GetMaximum() << std::endl;
	  empty->SetMaximum(graph->GetMaximum());
	  empty->Draw("hist");
	}
      //      graph->SetFillColor(2);
      //graph->SetFillStyle(3001);
      //      graph->Draw("a");
      graph->Draw("p same");/*
      graph->GetYaxis()->SetMinimum(0.);
      graph->GetYaxis()->SetMaximum(3000000.);
			*/
      leg->AddEntry(graph, "(#times 10^{3}) 0<|y|<0.5", "lp");
           
      for(int i=1; i<nybins; i++)
	{
	  TGraphAsymmErrors* graph2 = new TGraphAsymmErrors(nptbins, pt_bin_means, yield_array[i], pt_bin_edges_Lo, pt_bin_edges_Hi, errLo_array[i], errHi_array[i]);
	  graph2->SetTitle("Raw signal yield in Pt bins");
	  graph2->SetMarkerColor(i+1);
	  graph2->SetMarkerSize(0.5);
	  graph2->SetMarkerStyle(20);
	  //graph2->Draw("a2 same");
	  graph2->Draw("p same");

	  if(i==1)
	    leg->AddEntry(graph2,"(#times 10^{2}) 0.5<|y|<1", "lp");

	  if(i==2)
	    leg->AddEntry(graph2,"(#times 10) 1<|y|<1.5", "lp");

	  if(i==3)
	    {
	      leg->AddEntry(graph2," 1.5<|y|<2.25", "lp");
	      empty->SetMaximum(graph2->GetMaximum());
	    }
	  /*    graph2->GetYaxis()->SetMinimum(0.);
		graph2->GetYaxis()->SetMaximum(3000000.);*/
	}

      leg->Draw("same");
      cz.Update();
      cz.SetLogy();
      cz.SaveAs("signal_yield/signal_yield_" + yield_sub_samples + "_bins_" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + ".png");
      
    }//end of else
}//end of signal_yield_new

RooRealVar* pre_filter_efficiency(int channel, double pt_min, double pt_max)
{
  ReducedGenBranches gen;
  TString mc_gen_input_file = TString(BASE_DIR) + "myloop_gen_bfilter.root";
  TFile *fin = new TFile(mc_gen_input_file);
  
  TString ntuple = channel_to_ntuple_name(channel) + "_gen";
  TTree *tin = (TTree*)fin->Get(ntuple);
  gen.setbranchadd(tin);

  //use histograms to count the events, and TEfficiency for efficiency, because it takes care of the errors and propagation
  TH1D* hist_tot = new TH1D("hist_tot","hist_tot",1,pt_min,pt_max);
  TH1D* hist_passed = new TH1D("hist_passed","hist_passed",1,pt_min,pt_max);

  for (int evt=0;evt<tin->GetEntries();evt++)
    {
      tin->GetEntry(evt);
      
      if (fabs(gen.eta) > 2.4) continue; //B mesons inside the detector region eta < 2.4
      //if (gen.pt<pt_min || gen.pt>=pt_max) continue; // within the -gen- pt binning
           
      hist_tot->Fill(gen.pt);
      
      bool muon1Filter = fabs(gen.mu1eta)<2.4 && gen.mu1pt>2.8;
      bool muon2Filter = fabs(gen.mu2eta)<2.4 && gen.mu2pt>2.8;

      if (muon1Filter && muon2Filter) hist_passed->Fill(gen.pt);//count only the events with the muon selection above
    }
  /*
    TCanvas ch1;
    hist_tot->Draw();
    ch1.SaveAs(TString::Format("hist_tot_from_%d_to_%d.png",(int)pt_min,(int)pt_max));  
    TCanvas ch2;
    hist_passed->Draw();
    ch2.SaveAs(TString::Format("hist_passed_from_%d_to_%d.png",(int)pt_min,(int)pt_max));
  */

  //calculates the efficiency by dividing the histograms
  TEfficiency* efficiency = new TEfficiency(*hist_passed, *hist_tot);
  
  double eff;
  double eff_lo;
  double eff_hi;

  eff = efficiency->GetEfficiency(1);
  eff_lo = efficiency->GetEfficiencyErrorLow(1);
  eff_hi = efficiency->GetEfficiencyErrorUp(1);
  
  RooRealVar* eff1 = new RooRealVar("eff1","eff1",eff);
  eff1->setAsymError(-eff_lo,eff_hi);

  fin->Close();
  delete fin;

  return eff1;
}

RooRealVar* bin_mass_fit(RooWorkspace& w, int channel, double pt_min, double pt_max, double y_min, double y_max)
{
  RooRealVar pt = *(w.var("pt"));
  RooRealVar pt_low("pt_low","pt_low",pt_min);
  RooRealVar pt_high("pt_high","pt_high",pt_max);
  RooRealVar y = *(w.var("y"));
  RooRealVar y_low("y_low","y_low",y_min);
  RooRealVar y_high("y_high","y_high",y_max);
  RooAbsData* data_original;
  RooAbsData* data_cut;
  RooWorkspace ws_cut;
  RooAbsPdf* model_cut;
  RooRealVar* signal_res;

  data_original = w.data("data");
  
  set_up_workspace_variables(ws_cut,channel);
   
  std::cout << "y_low.getVal(): " << y_low.getVal() << std::endl;  
  std::cout << "y_high.getVal(): " << y_high.getVal() << std::endl;  
  std::cout << "y_min: " << y_min << std::endl;  
  std::cout << "y_max: " << y_max << std::endl;  

  RooFormulaVar cut("cut","pt>pt_low && pt<pt_high && ((y>y_low && y<y_high) || (y>-y_high && y<-y_low))",
		    RooArgList(pt,pt_low,pt_high,y,y_low,y_high));
  
  data_cut = data_original->reduce(cut);
  read_data_cut(ws_cut,data_cut);

  build_pdf(ws_cut,channel);
  
  model_cut = ws_cut.pdf("model");
  //ws_cut.Print();
  
  model_cut->fitTo(*data_cut,Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));

  TString dir = "";
  dir = "pt_bin_mass_fit/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + "/" + channel_to_ntuple_name(channel) + "mass_fit_" + TString::Format("pt_from_%d_to_%d_y_from_%lf_to_%lf",(int)pt_min,(int)pt_max,y_min,y_max);
  
  plot_mass_fit(ws_cut,channel,dir, (int) pt_max, (int) pt_min);

  //how to put the legend indicating each pt bin ??
  //change the plot_mass_fit to output a TCanvas, and write the legend on top after, and then have a function just to save the plots.
  //Legend(channel,(int)pt_bin_lo,(int)pt_bin_hi,1);
  
  signal_res = ws_cut.var("n_signal");
  
  return signal_res;
}

double pt_bin_mean(RooWorkspace& w, double pt_min, double pt_max)
{
  RooRealVar pt = *(w.var("pt"));
  RooRealVar pt_low("pt_low","pt_low",pt_min);
  RooRealVar pt_high("pt_high","pt_high",pt_max);
  RooAbsData* data_original;
  RooAbsData* data_cut;
  double centre;

  data_original = w.data("data");

  RooFormulaVar ptcut("pt_cut","pt>pt_low && pt<pt_high",RooArgList(pt,pt_low,pt_high));
  data_cut = data_original->reduce(ptcut);

  centre = (double) data_cut->meanVar(pt)->getVal();
  
  return centre;
}

void plot_mass_fit(RooWorkspace& w, int channel, TString directory, int pt_high, int pt_low)
{
  RooRealVar mass = *(w.var("mass"));
  RooAbsData* data = w.data("data");
  RooAbsPdf* model = w.pdf("model");
  RooRealVar lambda = *(w.var("m_exp"));
  RooRealVar mean = *(w.var("m_mean"));
  RooRealVar sigma1 = *(w.var("m_sigma1"));
  RooRealVar sigma2 = *(w.var("m_sigma2"));
  RooRealVar n_signal = *(w.var("n_signal"));
  RooRealVar n_back = *(w.var("n_combinatorial"));    
  RooPlot* frame_m = mass.frame();
  
  TH1D* histo_data = (TH1D*)data->createHistogram("histo_data", mass, Binning(channel_to_nbins(channel), mass.getMin(), mass.getMax() ));
  histo_data->Sumw2(false);
  histo_data->SetBinErrorOption(TH1::kPoisson);
  histo_data->SetMarkerStyle(20);
  histo_data->SetMarkerSize(0.8);
  histo_data->SetLineColor(kBlack);
  
  for (int i=1; i<=channel_to_nbins(channel); i++)
    if (histo_data->GetBinContent(i)==0) histo_data->SetBinError(i,0.);
  
  data->plotOn(frame_m,Name("theData"),Binning(channel_to_nbins(channel)),Invisible());
  
  model->plotOn(frame_m,Name("thePdf"),Precision(2E-4));
  
  //model->paramOn(frame_m); //show all the parameters of the fit in the plot.
  
  model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_signal"),LineColor(8),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(8), VLines(), DrawOption("F"));
  
  if (channel==1 || channel==2 || channel==3 || channel==4 || channel==6)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_combinatorial_exp"),LineColor(9),LineWidth(2),LineStyle(2));
  else
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_combinatorial_bern"),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
  
  if (channel==1 || channel==3)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_jpsix"),LineColor(kViolet),LineWidth(2),LineStyle(7));
  if (channel==5)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_x3872"),LineColor(kOrange),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kOrange), VLines(), DrawOption("F"));
  
  frame_m->SetTitle("");
  //frame_m->GetXaxis()->SetTitle(channel_to_xaxis_title(channel));
  //frame_m->GetXaxis()->SetLabelFont(42);
  //frame_m->GetXaxis()->SetLabelOffset(0.01);
  //frame_m->GetXaxis()->SetTitleSize(0.06);
  //frame_m->GetXaxis()->SetTitleOffset(1.09);
  //frame_m->GetXaxis()->SetLabelFont(42);
  //frame_m->GetXaxis()->SetLabelSize(0.055);
  //frame_m->GetXaxis()->SetTitleFont(42);
  frame_m->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass.getMax()-mass.getMin())*1000./channel_to_nbins(channel)));
  frame_m->GetYaxis()->SetLabelFont(42);
  frame_m->GetYaxis()->SetLabelOffset(0.01);
  frame_m->GetYaxis()->SetTitleOffset(0.8);
  frame_m->GetYaxis()->SetTitleSize(0.05);
  frame_m->GetYaxis()->SetTitleFont(42);
  frame_m->GetYaxis()->SetLabelFont(42);
  frame_m->GetYaxis()->SetLabelSize(0.055);
  
  RooHist* pull_hist = frame_m->pullHist("theData","thePdf");
  
  RooPlot* pull_plot = mass.frame();
  
  RooGenericPdf* line_ref = new RooGenericPdf("ref_0", "ref_0", RooConst(0.));
  line_ref->plotOn(pull_plot, LineStyle(7), LineColor(13), LineWidth(2));  


  pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
  pull_plot->SetTitle("");
  pull_plot->GetXaxis()->SetTitle(channel_to_xaxis_title(channel));
  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelOffset(0.01);
  pull_plot->GetXaxis()->SetTitleSize(0.17);
  pull_plot->GetXaxis()->SetTitleOffset(1.09);
  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelSize(0.17);
  pull_plot->GetXaxis()->SetTitleFont(42);
  pull_plot->GetXaxis()->SetTickLength(0.15);
  //  pull_plot->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass.getMax()-mass.getMin())*1000./channel_to_nbins(channel)));
  pull_plot->GetYaxis()->SetLabelFont(42);
  pull_plot->GetYaxis()->SetLabelOffset(0.01);
  pull_plot->GetYaxis()->SetLabelSize(0.17);
  pull_plot->GetYaxis()->SetTitleOffset(1.6);
  pull_plot->GetYaxis()->SetTitleSize(0.17);
  pull_plot->GetYaxis()->SetTitleFont(42);
  pull_plot->GetYaxis()->SetNdivisions(305);
  TCanvas *c1 = canvasDressing("c1"); c1->cd();
  
  // TPad *p1 = new TPad("p1","p1",0,0,1,1);
  //p1->Draw();
   
  TPad *p1 = new TPad("p1","p1",0.05,0.27,0.99,0.99);
  // TPad *p1 = new TPad("p1","p1",0.05,0.05,0.99,0.99);
  p1->SetBorderMode(0); 
  p1->SetFrameBorderMode(0); 
  p1->SetBorderSize(2);
  p1->SetBottomMargin(0.0);
  p1->Draw(); 
  
    
  TPad *p2 = new TPad("p2","p2",0.05,0.075,0.99,0.27); 
  p2->SetTopMargin(0.);    
  p2->SetBorderMode(0);
  p2->SetBorderSize(2); 
  p2->SetFrameBorderMode(0); 
  //p2->SetTicks(1,2); 
  p2->Draw();
//  p2->SetGridy(true);  
  RooAbsReal* nll = model->createNLL(*data);
  double log_likelihood= nll->getVal();
  std::stringstream ll_str;
  ll_str >> log_likelihood;
  double chis = frame_m->chiSquare();
  double lambda_exp = lambda.getVal();
  double lambda_exp_err = lambda.getError();
  double mean_gauss = mean.getVal();
  double mean_gauss_err = mean.getError();
  double sigma1_gauss = sigma1.getVal();
  double sigma1_gauss_err = sigma1.getError();
  double sigma2_gauss = sigma2.getVal();
  double sigma2_gauss_err = sigma2.getError();
  double signal_yield = n_signal.getVal();
  double signal_yield_err = n_signal.getError();
  double back_yield = n_back.getVal();
  double back_yield_err = n_back.getError();

  TLatex* tex1 = new TLatex(0.165, 0.88, Form("#lambda_{exp} = %.3lf+/-%.3lf",lambda_exp,lambda_exp_err));
  tex1->SetNDC(kTRUE);
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.03);
  tex1->Draw();  

  TLatex* tex2 = new TLatex(0.165, 0.84, Form("#mu_{gauss} = %.5lf+/-%.5lf",mean_gauss,mean_gauss_err));
  tex2->SetNDC(kTRUE);
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.03);
  tex2->Draw();  

  TLatex* tex3 = new TLatex(0.165, 0.80, Form("#sigma_{gauss1} = %.5lf+/-%.5lf",sigma1_gauss,sigma1_gauss_err));
  tex3->SetNDC(kTRUE);
  tex3->SetTextFont(42);
  tex3->SetTextSize(0.03);
  tex3->Draw();  

  TLatex* tex4 = new TLatex(0.165, 0.76, Form("#sigma_{gauss2} = %.5lf+/-%.5lf",sigma2_gauss,sigma2_gauss_err));
  tex4->SetNDC(kTRUE);
  tex4->SetTextFont(42);
  tex4->SetTextSize(0.03);
  if(data->sumEntries()>250){
    tex4->Draw();  
  }

  TLatex* tex5 = new TLatex(0.165, 0.70, Form("Signal = %.0lf+/-%.0lf",signal_yield,signal_yield_err));
  tex5->SetNDC(kTRUE);
  tex5->SetTextFont(42);
  tex5->SetTextSize(0.03);
  tex5->Draw();  

  TLatex* tex6 = new TLatex(0.165, 0.66, Form("Background = %.0lf+/-%.0lf",back_yield,back_yield_err));
  tex6->SetNDC(kTRUE);
  tex6->SetTextFont(42);
  tex6->SetTextSize(0.03);
  tex6->Draw();  

  TLatex* tex7 = new TLatex(0.165, 0.60, Form("lnL = %.3lf", log_likelihood));
  tex7->SetNDC(kTRUE);
  tex7->SetTextFont(42);
  tex7->SetTextSize(0.03);
  tex7->Draw();  

  TLatex* tex8 = new TLatex(0.165, 0.56, Form("#chi^{2} = %.3lf", chis));
  tex8->SetNDC(kTRUE);
  tex8->SetTextFont(42);
  tex8->SetTextSize(0.03);
  tex8->Draw();  
 
  p1->cd();
  frame_m->Draw();
  histo_data->Draw("Esame");
  Legend(channel,pt_low,pt_high,1);
  p2->cd();
  pull_plot->Draw();
  
  c1->SaveAs(directory + ".root");
  c1->SaveAs(directory + ".png"); }

void plot_mass_fit(RooWorkspace& w, int channel, TString directory)
{
  RooRealVar mass = *(w.var("mass"));
  RooAbsData* data = w.data("data");
  RooAbsPdf* model = w.pdf("model");
  RooRealVar lambda = *(w.var("m_exp"));
  RooRealVar mean = *(w.var("m_mean"));
  RooRealVar sigma1 = *(w.var("m_sigma1"));
  RooRealVar sigma2 = *(w.var("m_sigma2"));
  RooRealVar n_signal = *(w.var("n_signal"));
  RooRealVar n_back = *(w.var("n_combinatorial"));    
  RooPlot* frame_m = mass.frame();
  
  TH1D* histo_data = (TH1D*)data->createHistogram("histo_data", mass, Binning(channel_to_nbins(channel), mass.getMin(), mass.getMax() ));
  histo_data->Sumw2(false);
  histo_data->SetBinErrorOption(TH1::kPoisson);
  histo_data->SetMarkerStyle(20);
  histo_data->SetMarkerSize(0.8);
  histo_data->SetLineColor(kBlack);
  
  for (int i=1; i<=channel_to_nbins(channel); i++)
    if (histo_data->GetBinContent(i)==0) histo_data->SetBinError(i,0.);
  
  data->plotOn(frame_m,Name("theData"),Binning(channel_to_nbins(channel)),Invisible());
  
  model->plotOn(frame_m,Name("thePdf"),Precision(2E-4));
  
  //model->paramOn(frame_m); //show all the parameters of the fit in the plot.
  
  model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_signal"),LineColor(8),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(8), VLines(), DrawOption("F"));
  
  if (channel==1 || channel==2 || channel==3 || channel==4 || channel==6)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_combinatorial_exp"),LineColor(9),LineWidth(2),LineStyle(2));
  else
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_combinatorial_bern"),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
  
  if (channel==1 || channel==3)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_jpsix"),LineColor(kViolet),LineWidth(2),LineStyle(7));
  if (channel==5)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_x3872"),LineColor(kOrange),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kOrange), VLines(), DrawOption("F"));
  
  frame_m->SetTitle("");
  //frame_m->GetXaxis()->SetTitle(channel_to_xaxis_title(channel));
  //frame_m->GetXaxis()->SetLabelFont(42);
  //frame_m->GetXaxis()->SetLabelOffset(0.01);
  //frame_m->GetXaxis()->SetTitleSize(0.06);
  //frame_m->GetXaxis()->SetTitleOffset(1.09);
  //frame_m->GetXaxis()->SetLabelFont(42);
  //frame_m->GetXaxis()->SetLabelSize(0.055);
  //frame_m->GetXaxis()->SetTitleFont(42);
  frame_m->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass.getMax()-mass.getMin())*1000./channel_to_nbins(channel)));
  frame_m->GetYaxis()->SetLabelFont(42);
  frame_m->GetYaxis()->SetLabelOffset(0.01);
  frame_m->GetYaxis()->SetTitleOffset(1.4);
  frame_m->GetYaxis()->SetTitleSize(0.05);
  frame_m->GetYaxis()->SetTitleFont(42);
  frame_m->GetYaxis()->SetLabelFont(42);
  frame_m->GetYaxis()->SetLabelSize(0.055);
  
  RooHist* pull_hist = frame_m->pullHist("theData","thePdf");
  
  RooPlot* pull_plot = mass.frame();
  
  RooGenericPdf* line_ref = new RooGenericPdf("ref_0", "ref_0", RooConst(0.));
  line_ref->plotOn(pull_plot, LineStyle(7), LineColor(13), LineWidth(2));  


  pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
  pull_plot->SetTitle("");
  pull_plot->GetXaxis()->SetTitle(channel_to_xaxis_title(channel));
  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelOffset(0.01);
  pull_plot->GetXaxis()->SetTitleSize(0.17);
  pull_plot->GetXaxis()->SetTitleOffset(1.09);
  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelSize(0.17);
  pull_plot->GetXaxis()->SetTitleFont(42);
  pull_plot->GetXaxis()->SetTickLength(0.15);
  //  pull_plot->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass.getMax()-mass.getMin())*1000./channel_to_nbins(channel)));
  pull_plot->GetYaxis()->SetLabelFont(42);
  pull_plot->GetYaxis()->SetLabelOffset(0.01);
  pull_plot->GetYaxis()->SetLabelSize(0.17);
  pull_plot->GetYaxis()->SetNdivisions(305);
  // pull_plot->GetYaxis()->SetTitleOffset(1.6);
  //pull_plot->GetYaxis()->SetTitleSize(0.17);
  //pull_plot->GetYaxis()->SetTitleFont(42);
  TCanvas *c1 = canvasDressing("c1"); c1->cd();
  
  // TPad *p1 = new TPad("p1","p1",0,0,1,1);
  //p1->Draw();
   
  TPad *p1 = new TPad("p1","p1",0.05,0.27,0.99,0.99);
  // TPad *p1 = new TPad("p1","p1",0.05,0.05,0.99,0.99);
  p1->SetBorderMode(0); 
  p1->SetFrameBorderMode(0); 
  p1->SetBorderSize(2);
  p1->SetBottomMargin(0.0);
  p1->Draw(); 
  
    
  TPad *p2 = new TPad("p2","p2",0.05,0.075,0.99,0.27); 
  p2->SetTopMargin(0.);    
  p2->SetBorderMode(0);
  p2->SetBorderSize(2); 
  p2->SetFrameBorderMode(0); 
  //p2->SetTicks(1,2); 
  p2->Draw();
  
  RooAbsReal* nll = model->createNLL(*data);
  double log_likelihood= nll->getVal();
  std::stringstream ll_str;
  ll_str >> log_likelihood;
  double chis = frame_m->chiSquare();
  double lambda_exp = lambda.getVal();
  double lambda_exp_err = lambda.getError();
  double mean_gauss = mean.getVal();
  double mean_gauss_err = mean.getError();
  double sigma1_gauss = sigma1.getVal();
  double sigma1_gauss_err = sigma1.getError();
  double sigma2_gauss = sigma2.getVal();
  double sigma2_gauss_err = sigma2.getError();
  double signal_yield = n_signal.getVal();
  double signal_yield_err = n_signal.getError();
  double back_yield = n_back.getVal();
  double back_yield_err = n_back.getError();

  TLatex* tex1 = new TLatex(0.165, 0.88, Form("#lambda_{exp} = %.3lf+/-%.3lf",lambda_exp,lambda_exp_err));
  tex1->SetNDC(kTRUE);
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.03);
  tex1->Draw();  

  TLatex* tex2 = new TLatex(0.165, 0.84, Form("#mu_{gauss} = %.5lf+/-%.5lf",mean_gauss,mean_gauss_err));
  tex2->SetNDC(kTRUE);
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.03);
  tex2->Draw();  

  TLatex* tex3 = new TLatex(0.165, 0.80, Form("#sigma_{gauss1} = %.5lf+/-%.5lf",sigma1_gauss,sigma1_gauss_err));
  tex3->SetNDC(kTRUE);
  tex3->SetTextFont(42);
  tex3->SetTextSize(0.03);
  tex3->Draw();  

  TLatex* tex4 = new TLatex(0.165, 0.76, Form("#sigma_{gauss2} = %.5lf+/-%.5lf",sigma2_gauss, sigma2_gauss_err));
  tex4->SetNDC(kTRUE);
  tex4->SetTextFont(42);
  tex4->SetTextSize(0.03);
  if(data->sumEntries()>250){
    tex4->Draw();  
  }

  TLatex* tex5 = new TLatex(0.165, 0.70, Form("Signal = %.0lf+/-%.0lf",signal_yield, signal_yield_err));
  tex5->SetNDC(kTRUE);
  tex5->SetTextFont(42);
  tex5->SetTextSize(0.03);
  tex5->Draw();  

  TLatex* tex6 = new TLatex(0.165, 0.66, Form("Background = %.0lf+/-%.0lf",back_yield, back_yield_err));
  tex6->SetNDC(kTRUE);
  tex6->SetTextFont(42);
  tex6->SetTextSize(0.03);
  tex6->Draw();  

  TLatex* tex7 = new TLatex(0.165, 0.60, Form("lnL = %.3lf", log_likelihood));
  tex7->SetNDC(kTRUE);
  tex7->SetTextFont(42);
  tex7->SetTextSize(0.03);
  tex7->Draw();  

  TLatex* tex8 = new TLatex(0.165, 0.56, Form("#chi^{2} = %.3lf", chis));
  tex8->SetNDC(kTRUE);
  tex8->SetTextFont(42);
  tex8->SetTextSize(0.03);
  tex8->Draw();  
 
  p1->cd();
  frame_m->Draw();
  histo_data->Draw("Esame");
  Legend(channel,0,0,0);
  p2->cd();
  pull_plot->Draw();
  
  c1->SaveAs(directory + ".root");
  c1->SaveAs(directory + ".png"); }

void plot_pt_dist(RooWorkspace& w, int channel, TString directory)
{
  //full dataset pt distribution
  RooRealVar pt = *(w.var("pt"));
  RooAbsData* data = w.data("data");

  TCanvas c2;
  TH1D* pt_dist = (TH1D*)data->createHistogram("pt_dist",pt);
  pt_dist->Draw();
  c2.SetLogy();
  //c2.SaveAs(directory + ".root");
  c2.SaveAs(directory + ".png");
}

void build_pdf(RooWorkspace& w, int channel)
{
  double mass_peak;

  RooRealVar mass = *(w.var("mass"));
  RooRealVar pt = *(w.var("pt"));  
  RooAbsData* data = w.data("data");

  switch (channel) {
  default:
  case 1:
    mass_peak = BP_MASS;
    break;
  case 2:
    mass_peak = B0_MASS;
    break;
  case 3:
    mass_peak = B0_MASS;
    break;
  case 4:
    mass_peak = BS_MASS;
    break;
  case 5:
    mass_peak = PSI2S_MASS;
    break;
  case 6:
    mass_peak = LAMBDAB_MASS;
    break;
  }
  
  double n_signal_initial = data->sumEntries(TString::Format("abs(mass-%g)<0.015",mass_peak))
    - data->sumEntries(TString::Format("abs(mass-%g)<0.030&&abs(mass-%g)>0.015",mass_peak,mass_peak));
  
  if(n_signal_initial<0)
  n_signal_initial=1;

  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  
  //-----------------------------------------------------------------
  // signal PDF 
  RooRealVar m_mean("m_mean","m_mean",mass_peak,mass_peak-0.09,mass_peak+0.09);
  RooRealVar m_sigma1("m_sigma1","m_sigma1",0.015,0.005,0.07);
  RooRealVar m_sigma2("m_sigma2","m_sigma2",0.030,0.001,0.100);
  RooRealVar m_fraction("m_fraction","m_fraction",0.5);
  RooGaussian m_gaussian1("m_gaussian1","m_gaussian1",mass,m_mean,m_sigma1);
  RooGaussian m_gaussian2("m_gaussian2","m_gaussian2",mass,m_mean,m_sigma2);
  RooAddPdf pdf_m_signal("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));
  
  // use single Gaussian for J/psi Ks and J/psi Lambda due to low statistics
  if (channel==3 || channel==6 || data->sumEntries()<250) {
    m_sigma2.setConstant(kTRUE);
    m_fraction.setVal(1.);
  }
  
  //-----------------------------------------------------------------
  // combinatorial background PDF (exponential or bernstean poly.)
  
  RooRealVar m_exp("m_exp","m_exp",-0.3,-4.,0.);
  RooExponential pdf_m_combinatorial_exp("pdf_m_combinatorial_exp","pdf_m_combinatorial_exp",mass,m_exp);
  
  RooRealVar m_par1("m_par1","m_par2",1.,0,+10.);
  RooRealVar m_par2("m_par2","m_par3",1.,0,+10.);
  RooRealVar m_par3("m_par3","m_par3",1.,0,+10.);
  
  RooBernstein pdf_m_combinatorial_bern("pdf_m_combinatorial_bern","pdf_m_combinatorial_bern",mass,RooArgList(RooConst(1.),m_par1,m_par2,m_par3));
  //erfc component on channel 1 and 3
  RooFormulaVar pdf_m_jpsix("pdf_m_jpsix","2.7*erfc((mass-5.14)/(0.5*0.08))",{mass});
  
  //-----------------------------------------------------------------
  // X(3872) PDF, only for J/psi pipi fit
  
  RooRealVar m_x3872_mean("m_x3872_mean","m_x3872_mean",3.872,3.7,3.9);
  RooRealVar m_x3872_sigma("m_x3872_sigma","m_x3872_sigma",0.01,0.001,0.010);
  RooGaussian pdf_m_x3872("pdf_m_x3872","pdf_m_x3872",mass,m_x3872_mean,m_x3872_sigma);
  
  //-----------------------------------------------------------------
  // full model
  
  RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
  RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
  RooRealVar n_x3872("n_x3872","n_x3872",200.,0.,data->sumEntries());
  
  RooRealVar n_jpsix("n_jpsix","n_jpsix",data->sumEntries(TString::Format("mass>4.9&&mass<5.14")),data->sumEntries(TString::Format("mass>4.9&&mass<5.14")),data->sumEntries());
  
  RooAddPdf* model;

  switch(channel)
    {
    default:
    case 1:// B+ -> J/psi K+
    case 3://B0 -> J/psi Ks
      model = new RooAddPdf("model","model",
			    RooArgList(pdf_m_signal, pdf_m_combinatorial_exp, pdf_m_jpsix),
			    RooArgList(n_signal, n_combinatorial, n_jpsix));
      break;
    case 2:// B0 -> J/psi K* 
    case 4://Bs -> J/psi phi
    case 6://Lambda_b -> J/psi Lambda
      model = new RooAddPdf("model","model",
			    RooArgList(pdf_m_signal, pdf_m_combinatorial_exp),
			    RooArgList(n_signal, n_combinatorial));
      break;
    case 5:// J/psi pipi
      model = new RooAddPdf("model","model",
			    RooArgList(pdf_m_signal, pdf_m_combinatorial_bern, pdf_m_x3872),
			    RooArgList(n_signal, n_combinatorial, n_x3872));
      break;
    }

  w.import(*model);
}

void read_data(RooWorkspace& w, TString filename,int channel)
{
  TFile* f = new TFile(filename);
  TNtupleD* _nt = (TNtupleD*)f->Get(channel_to_ntuple_name(channel));
 
  RooDataSet* data = new RooDataSet("data","data",_nt,RooArgSet( *(w.var("mass")) , *(w.var("pt")) , *(w.var("y")) ));
  
  w.import(*data);
}

void read_data_cut(RooWorkspace& w, RooAbsData* data)
{
  w.import(*data);
}

void set_up_workspace_variables(RooWorkspace& w, int channel)
{
  double mass_min, mass_max;
  double pt_min, pt_max;
  double y_min, y_max;

  pt_min=0;
  pt_max=400;

  y_min=-2.4;
  y_max=2.4;
  
  switch (channel) {
  default: 
  case 1:
    mass_min = 5.0; mass_max = 6.0;
    break;
  case 2:
    mass_min = 5.0; mass_max = 5.6;
    break;
  case 3:
    mass_min = 5.0; mass_max = 6.0;
    break;
  case 4:
    mass_min = 5.0; mass_max = 6.0;
    break;
  case 5:
    mass_min = 3.6; mass_max = 4.0;
    break;
  case 6:
    mass_min = 5.3; mass_max = 6.3;
    break;
  }

  RooRealVar mass("mass","mass",mass_min,mass_max);
  RooRealVar pt("pt","pt",pt_min,pt_max);
  RooRealVar y("y", "y", y_min, y_max);

  w.import(mass);
  w.import(pt);
  w.import(y);
}

TString channel_to_ntuple_name(int channel)
{
  //returns a TString with the ntuple name corresponding to the channel. It can be used to find the data on each channel saved in a file. or to write the name of a directory

  TString ntuple_name = "";

  switch(channel){
  default:
  case 1:
    ntuple_name="ntkp";
    break;
  case 2:
    ntuple_name="ntkstar";
    break;
  case 3:
    ntuple_name="ntks";
    break;
  case 4:
    ntuple_name="ntphi";
    break;
  case 5:
    ntuple_name="ntmix";
    break;
  case 6:
    ntuple_name="ntlambda";
    break;
  }
  return ntuple_name;
}

TString channel_to_xaxis_title(int channel)
{
  TString xaxis_title = "";

  switch (channel) {
  default:
  case 1:
    xaxis_title = "M_{J/#psi K^{#pm}} [GeV]";
    break;
  case 2:
    xaxis_title = "M_{J/#psi K^{#pm}#pi^{#mp}} [GeV]";
    break;
  case 3:
    xaxis_title = "M_{J/#psi K^{0}_{S}} [GeV]";
    break;
  case 4:
    xaxis_title = "M_{J/#psi K^{#pm}K^{#mp}} [GeV]";
    break;
  case 5:
    xaxis_title = "M_{J/#psi #pi^{#pm}#pi^{#mp}} [GeV]";
    break;
  case 6:
    xaxis_title = "M_{J/#psi #Lambda} [GeV]";
    break;
  }
  return xaxis_title;
}

int channel_to_nbins(int channel)
{
  int nbins;

  switch (channel) {
  default:
  case 1:
    nbins = 50;
    break;
  case 2:
    nbins = 50;
    break;
  case 3:
    nbins = 50;
    break;
  case 4:
    nbins = 50;
    break;
  case 5:
    nbins = 80;
    break;
  case 6:
    nbins = 50;
    break;
  }
  return nbins;
}

void create_dir(std::vector<std::string> list)
{
  //to create the directories needed to save the output files, like .png and .root
  for(size_t i=0 ; i< list.size() ; ++i)
    {
      gSystem->Exec(("mkdir -p " + list[i]).c_str());
    }
}

/*
  switch (channel) {
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
      
  RooDataSet* data_original  = new RooDataSet("data_original", "data_original", *(ws->data("data")->get()),Import( *(dynamic_cast<RooDataSet *>(ws->data("data"))) ));
      
  RooRealVar pt = *(ws->var("pt"));
      
  RooThresholdCategory ptRegion("ptRegion", "region of pt", pt);
  ptRegion.addThreshold(*(pt_bin_edges),"below 1st bin");

  for(int i=0; i<nptbins; i++)
  {
  TString reg = TString::Format("PtBin%d",i+1);
  ptRegion.addThreshold(*(pt_bin_edges+i+1),reg);
  }
  data_original->addColumn(ptRegion);
      
  Roo1DTable * tab = data_original->table(ptRegion);
  tab->Print("v");
  delete tab;
      
  //to produce and process each pt subsample                                                                                                    
  RooDataSet *data_cut;
  RooWorkspace* ws_cut;
  RooAbsPdf* model_cut;
  RooRealVar* pt_mean;
  TString directory="";
  double pt_bin_centre[nptbins];
  double pt_bin_edges_Lo[nptbins];
  double pt_bin_edges_Hi[nptbins];
  double yield_array[nptbins];
  double errLo_array[nptbins];
  double errHi_array[nptbins];
      
  for(int i=0; i<nptbins; i++)
  {
  cout << "processing subsample pt: " << i+1 << std::endl;

  TString ptcut(TString::Format("(ptRegion==ptRegion::PtBin%d)", i+1));

  data_cut = new RooDataSet("data", "data", *(data_original->get()),Import(*data_original), Cut(ptcut));

  // TString ptcut(TString::Format("(pt>(pt_bin_edges+%d))&&(pt<(pt_bin_edges+%d+1))",i));
  //RooFormulaVar ptcut("pt_cut","pt_cut","pt>*(pt_bin_edges+i) && pt<*(pt_bin_edges+i+1)",RooArgList(pt_bin_edges,i));
  //data_cut=data_original->reduce(Cut(ptcut));

  ws_cut = new RooWorkspace("ws_cut","Bmass_cut");
  set_up_workspace_variables(*ws_cut,channel);
  read_data_cut(*ws_cut,data_cut);
  build_pdf(*ws_cut,channel);

  model_cut = ws_cut->pdf("model");
  ws_cut->Print();
	
  pt_mean = data_cut->meanVar(pt); //older way: pt_bin_centre[i] = *(pt_bin_edges+i) + (*(pt_bin_edges+i+1)-*(pt_bin_edges+i))/2;
  pt_bin_centre[i] = (double) pt_mean->getVal();
  pt_bin_edges_Lo[i] = pt_bin_centre[i] - *(pt_bin_edges+i);
  pt_bin_edges_Hi[i] = *(pt_bin_edges+i+1) - pt_bin_centre[i];

  fit_res = model_cut->fitTo(*data_cut,Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));
  signal_res = ws_cut->var("n_signal");

  yield_array[i] = signal_res->getVal();
  errLo_array[i] = -signal_res->getAsymErrorLo();
  errHi_array[i] = signal_res->getAsymErrorHi();

  directory = "pt_bin_mass_fit/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + "/" + "mass_fit_" + channel_to_ntuple_name(channel) + TString::Format("_bin_%d_%d", (int)*(pt_bin_edges+i), (int)*(pt_bin_edges+i+1));
	  
  plot_mass_fit(*ws_cut,channel,directory);
	  
  //how to put the legend indicating each pt bin ??
  //change the plot_mass_fit to output a TCanvas, and write the legend on top after, and then have a function just to save the plots.
  //  Legend(channel,(int)pt_bin_lo,(int)pt_bin_hi,1);
  }
  for(int i=0; i<nptbins; i++)
  {
  std::cout << "BIN: "<< (int) *(pt_bin_edges+i) << " to " << (int) *(pt_bin_edges+i+1) << " : " <<  yield_array[i] << " +" << errHi_array[i] << " -"<< errLo_array[i] << std::endl;
  }

  TCanvas cz;
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(nptbins, pt_bin_centre, yield_array, pt_bin_edges_Lo, pt_bin_edges_Hi, errLo_array, errHi_array);
  graph->SetTitle("Raw signal yield in Pt bins");
  graph->SetFillColor(2);
  graph->SetFillStyle(3001);
  graph->Draw("a2");
  graph->Draw("p");
  cz.SetLogy();
  // cz.SaveAs("signal_yield/signal_yield_" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + ".root");
  cz.SaveAs("signal_yield/signal_yield_" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + ".png");     
*/
