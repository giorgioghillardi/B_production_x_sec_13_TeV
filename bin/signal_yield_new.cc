#include <RooPlotable.h>
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
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <RooPlot.h>
#include "TMath.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/myloop.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/plotDressing.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/channel.h"

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
void plot_mass_fit(RooWorkspace& w, int channel, TString directory);

RooRealVar* bin_mass_fit(RooWorkspace& w, int channel, double pt_min, double pt_max);
double pt_bin_mean(RooWorkspace& w, double pt_min, double pt_max);
RooRealVar* prefilter_efficiency(int channel, double pt_min, double pt_max);
RooRealVar* reco_efficiency(RooWorkspace& ws_mc, int channel, double pt_min, double pt_max);

void build_pdf(RooWorkspace& w, int channel, int mc);
void read_data(RooWorkspace& w, TString filename,int channel,double pt_min, double pt_max);
void read_data_cut(RooWorkspace& w, RooAbsData* data);
void set_up_workspace_variables(RooWorkspace& w, int channel);

//input example: signal_yield_new --channel 1 --bins 1 --preeff 1 --recoeff 1

int main(int argc, char** argv)
{
  int channel = 0;
  int yield_sub_samples = 0;
  int calculate_pre_filter_eff = 0;
  int calculate_reco_eff = 0;

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
  dir_list.push_back(static_cast<const char*>("pt_bin_mass_fit/" + channel_to_ntuple_name(channel) + "_mc_" + TString::Format(VERSION)));
  dir_list.push_back("signal_yield");
  
  create_dir(dir_list);

  double ntkp_pt_bin_edges[]={9,13,16,20,25,30,35,42,50,60,70,80,90,100,120,150};
  double ntkstar_pt_bin_edges[]={9,13,16,20,25,30,35,42,50,60,70,90};
  double ntks_pt_bin_edges[]={0};
  double ntphi_pt_bin_edges[]={9,13,16,20,25,30,35,42,50,60,70,90};
  double ntmix_pt_bin_edges[]={0};
  double ntlambda_pt_bin_edges[]={0};
  double* pt_bin_edges;
  
  int nptbins;
  
  //to read the data into the workspace
  TString data_selection_input_file = TString::Format(BASE_DIR) + "selected_data_" + channel_to_ntuple_name(channel) + ".root";
  RooRealVar* signal_res; 
  RooWorkspace* ws;

  //to read MC into a different workspace
  //TString mc_input_file = TString::Format(BASE_DIR) + "reduced_mc_tree_" + channel_to_ntuple_name(channel) + ".root";
  RooWorkspace* ws_mc;
  
  TString pt_dist_directory="";
  TString mass_fit_directory="";

  if(!yield_sub_samples) //mass fit and plot the full dataset
    { 
      //NOT WORKING. NEED TO BE CORRECTED.
      /*
      //read data from the selected data file, and import it as a dataset into the workspace.
      read_data(*ws, data_selection_input_file,channel);
      
      //build the pdf for the channel selected above, it uses the dataset which is saved in ws. need to change the dataset to change the pdf.
      build_pdf(*ws,channel,0);
      data = ws->data("data");
      model = ws->pdf("model");     
      
      model->fitTo(*data,Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));

      signal_res = ws->var("n_signal");
      
      std::cout <<"SIGNAL: "<< signal_res->getVal() << " " << signal_res->getAsymErrorLo() << " +" << signal_res->getAsymErrorHi() << std::endl;
      
      mass_fit_directory = "full_dataset_mass_fit/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION);
      pt_dist_directory = "full_dataset_mass_pt_histo/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION);

      plot_mass_fit(*ws,channel,mass_fit_directory);
      plot_pt_dist(*ws,channel,pt_dist_directory);
      */
      return 0;
    }
  else
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
      double pt_bin_size[nptbins];
      double pt_bin_means[nptbins];
      double pt_bin_edges_Lo[nptbins];
      double pt_bin_edges_Hi[nptbins];
      
      double yield_array[nptbins];
      double errLo_array[nptbins];
      double errHi_array[nptbins];
      
      double pt_bin_centres_eff[nptbins];
      double pt_bin_edges_eff_Lo[nptbins];
      double pt_bin_edges_eff_Hi[nptbins];
      
      RooRealVar* pre_filter_eff;
      
      double pre_eff_array[nptbins];
      double pre_effLo_array[nptbins];
      double pre_effHi_array[nptbins];
      
      RooRealVar* reco_eff;
      
      double reco_eff_array[nptbins];
      double reco_effLo_array[nptbins];
      double reco_effHi_array[nptbins];
      
      double x_sec_array[nptbins];
      double x_sec_errLo_array[nptbins];
      double x_sec_errHi_array[nptbins];

      //to extract the signal_yield using the fit
      for(int i=0; i<nptbins; i++)
	{
	  std::cout << "processing subsample: " << (int)pt_bin_edges[i] << " < pt < " << (int)pt_bin_edges[i+1] << std::endl;
	  
	   ws = new RooWorkspace("ws","Bmass");

	  //read the data
	  set_up_workspace_variables(*ws,channel);
	  read_data(*ws, data_selection_input_file, channel, pt_bin_edges[i], pt_bin_edges[i+1]);

	  pt_bin_size[i] = pt_bin_edges[i+1]-pt_bin_edges[i];
	  
	  pt_bin_means[i] = pt_bin_mean(*ws,pt_bin_edges[i],pt_bin_edges[i+1]);
	  pt_bin_edges_Lo[i] = pt_bin_means[i] - pt_bin_edges[i];
	  pt_bin_edges_Hi[i] = pt_bin_edges[i+1] - pt_bin_means[i];
	  
	  signal_res = bin_mass_fit(*ws,channel,pt_bin_edges[i],pt_bin_edges[i+1]);
              
	  yield_array[i] = signal_res->getVal();
	  errLo_array[i] = -signal_res->getAsymErrorLo();
	  errHi_array[i] = signal_res->getAsymErrorHi();
	  
	  delete ws;
	}
      
      //to show the values of signal_yield, like a table
      for(int i=0; i<nptbins; i++)
	{
	  std::cout << "BIN: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  yield_array[i] << " +" << errHi_array[i] << " -"<< errLo_array[i] << std::endl;
	}
      
      //to calculate the prefilter efficiency
      if(calculate_pre_filter_eff)
	{
	  for(int i=0; i<nptbins; i++)
	    {
	      std::cout << "calculating pre-filter efficiency: " << (int)pt_bin_edges[i] << " < pt < " << (int)pt_bin_edges[i+1] << std::endl;
	      
	      pt_bin_centres_eff[i] = pt_bin_edges[i] + (pt_bin_edges[i+1]-pt_bin_edges[i])/2;
	      pt_bin_edges_eff_Lo[i] = pt_bin_centres_eff[i] - pt_bin_edges[i];
	      pt_bin_edges_eff_Hi[i] = pt_bin_edges[i+1] - pt_bin_centres_eff[i];
	      
	      //pre-filter efficiency
	      pre_filter_eff = prefilter_efficiency(channel,pt_bin_edges[i],pt_bin_edges[i+1]);
	      
	      pre_eff_array[i] = pre_filter_eff->getVal();
	      pre_effLo_array[i] = -pre_filter_eff->getAsymErrorLo();
	      pre_effHi_array[i] = pre_filter_eff->getAsymErrorHi();
	    }
	  //plot of the pre-filter efficiency as a function of pT
	  TCanvas ce;
	  TGraphAsymmErrors* graph_pre_eff = new TGraphAsymmErrors(nptbins, pt_bin_centres_eff, pre_eff_array, pt_bin_edges_eff_Lo, pt_bin_edges_eff_Hi, pre_effLo_array, pre_effHi_array);
	  graph_pre_eff->SetTitle("efficiency");
	  graph_pre_eff->SetMarkerColor(4);
	  graph_pre_eff->SetMarkerStyle(21);
	  graph_pre_eff->Draw("AP");
	  TString eff1_name = "";
	  eff1_name = "pre_filter_efficiency_" + channel_to_ntuple_name(channel) + ".png";
	  ce.SaveAs(eff1_name);
	}
      //to calculate the reco efficiency
      if(calculate_reco_eff)
	{
	  for(int i=0; i<nptbins; i++)
	    {
	      std::cout << "calculating reconstruction efficiency: " << (int)pt_bin_edges[i] << " < pt < " << (int)pt_bin_edges[i+1] << std::endl;
	      ws_mc = new RooWorkspace("ws_mc","Bmass_mc");

	      //read the data
	      set_up_workspace_variables(*ws_mc,channel);
	      //read_data(*ws_mc, mc_input_file,channel, pt_bin_edges[i], pt_bin_edges[i+1]);

	      pt_bin_centres_eff[i] = pt_bin_edges[i] + (pt_bin_edges[i+1]-pt_bin_edges[i])/2;
	      pt_bin_edges_eff_Lo[i] = pt_bin_centres_eff[i] - pt_bin_edges[i];
	      pt_bin_edges_eff_Hi[i] = pt_bin_edges[i+1] - pt_bin_centres_eff[i];
	      
	      //reconstruction efficiency
	      reco_eff = reco_efficiency(*ws_mc,channel, pt_bin_edges[i], pt_bin_edges[i+1]);
	      
	      reco_eff_array[i] = reco_eff->getVal();
	      reco_effLo_array[i] = -reco_eff->getAsymErrorLo();
	      reco_effHi_array[i] = reco_eff->getAsymErrorHi();
	    
	      delete ws_mc;
	    }
	  
	  //plot of the reconstruction efficiency as a function of pT
	  TCanvas cp;
	  TGraphAsymmErrors* graph_reco_eff = new TGraphAsymmErrors(nptbins, pt_bin_centres_eff, reco_eff_array, pt_bin_edges_eff_Lo, pt_bin_edges_eff_Hi, reco_effLo_array, reco_effHi_array);
	  graph_reco_eff->SetTitle("efficiency");
	  graph_reco_eff->SetMarkerColor(4);
	  graph_reco_eff->SetMarkerStyle(21);
	  graph_reco_eff->Draw("AP");
	  TString eff2_name = "";
	  eff2_name = "reco_efficiency_" + channel_to_ntuple_name(channel) + ".png";
	  cp.SaveAs(eff2_name);
	}
      
      //plot x-sec corrected by efficiency and other factors
      for(int i=0; i<nptbins; i++)
	{
	  if(calculate_reco_eff && calculate_pre_filter_eff)
	    {
	      x_sec_array[i] = (yield_array[i] / (2 * pt_bin_size[i] * reco_eff_array[i] * pre_eff_array[i] * (2.71) * (1.014e-3) * (5.93e-2)) )* (1e-9);
	      x_sec_errLo_array[i] = (errLo_array[i] / (2 * pt_bin_size[i] * reco_eff_array[i] * pre_eff_array[i] * (2.71) * (1.014e-3) * (5.93e-2)) ) * (1e-9);
	      x_sec_errHi_array[i] = (errHi_array[i] / (2 * pt_bin_size[i] * reco_eff_array[i] * pre_eff_array[i] * (2.71) * (1.014e-3) * (5.93e-2)) ) * (1e-9);
	    }
	  else
	    {
	      x_sec_array[i] = yield_array[i] / pt_bin_size[i];
	      x_sec_errLo_array[i] = errLo_array[i] / pt_bin_size[i];
	      x_sec_errHi_array[i] = errHi_array[i] / pt_bin_size[i];
	    }
	}

      //to show the values of signal_yield, like a table
      for(int i=0; i<nptbins; i++)
	{
	  std::cout << "BIN: "<< (int) pt_bin_edges[i] << " to " << (int) pt_bin_edges[i+1] << " : " <<  x_sec_array[i] << " +" << x_sec_errHi_array[i] << " -"<< x_sec_errLo_array[i] << std::endl;
	}

      TCanvas cz;
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(nptbins, pt_bin_means, x_sec_array, pt_bin_edges_Lo, pt_bin_edges_Hi, x_sec_errLo_array, x_sec_errHi_array);
      graph->SetTitle("Cross section in Pt bins");
      graph->SetFillColor(2);
      graph->SetFillStyle(3001);
      graph->Draw("a2");
      graph->Draw("p");
      cz.SetLogy();
      cz.SaveAs("signal_yield/signal_yield_" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + ".png");
      
    }//end of else
}//end of signal_yield_new

//to caltulate the prefilter efficiency. it counts the evevnts that pass the pre-filter cut.
//the input file must be produced with myloop_gen.cc to have the gen info. otherwise the signal needs to be extracted using a fit.
RooRealVar* prefilter_efficiency(int channel, double pt_min, double pt_max)
{
  TString mc_gen_input_file = TString::Format(BASE_DIR) + "myloop_gen_" + channel_to_ntuple_name(channel) + "_bfilter.root";
  TFile *fin = new TFile(mc_gen_input_file);
  
  TString ntuple_name = channel_to_ntuple_name(channel) + "_gen";
  TTree *tin = (TTree*)fin->Get(ntuple_name);

  ReducedGenBranches gen;
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

RooRealVar* reco_efficiency(RooWorkspace& w,int channel, double pt_min, double pt_max)
{
  //read monte carlo without cuts
  TString mc_input_no_cuts = TString::Format(BASE_DIR) + "reduced_" + channel_to_ntuple_name(channel) + "_bmuonfilter_no_cuts.root";
  TFile *fin_no_cuts = new TFile(mc_input_no_cuts);
  TTree *tin_no_cuts = (TTree*)fin_no_cuts->Get(channel_to_ntuple_name(channel));
  ReducedBranches mc_1;
  mc_1.setbranchadd(tin_no_cuts);

  //read monte carlo with cuts
  // TString mc_input_with_cuts = TString::Format(BASE_DIR) + "reduced_" + channel_to_ntuple_name(channel) + "_bmuonfilter_with_cuts.root";
  TString mc_input_with_cuts = TString::Format(BASE_DIR) + "selected_mc_" + channel_to_ntuple_name(channel) + "_bmuonfilter_with_cuts.root";
  TFile *fin_with_cuts = new TFile(mc_input_with_cuts);
  TTree *tin_with_cuts = (TTree*)fin_with_cuts->Get(channel_to_ntuple_name(channel));
  ReducedBranches mc_2;
  mc_2.setbranchadd(tin_with_cuts);

  //use histograms to count the events, and TEfficiency for efficiency, because it takes care of the errors and propagation
  TH1D* hist_tot = new TH1D("hist_tot","hist_tot",1,pt_min,pt_max);
  TH1D* hist_passed = new TH1D("hist_passed","hist_passed",1,pt_min,pt_max);

  Printf("counting no_cuts \n");

  for (int evt=0;evt<tin_no_cuts->GetEntries();evt++)
    {
      tin_no_cuts->GetEntry(evt);
      
      if (fabs(mc_1.eta) > 2.4) continue; //B mesons inside the detector region eta < 2.4
      if (fabs(mc_1.mass) < 5 || mc_1.mass > 5.5) continue;
      //if (gen.pt<pt_min || gen.pt>=pt_max) continue; // within the -gen- pt binning
           
      bool muon1Filter = fabs(mc_1.mu1eta)<2.4 && mc_1.mu1pt>2.8;
      bool muon2Filter = fabs(mc_1.mu2eta)<2.4 && mc_1.mu2pt>2.8;
 
      if (muon1Filter && muon2Filter) hist_tot->Fill(mc_1.pt);//count only the events with the muon selection above
      
      //hist_tot->Fill(mc_1.pt);
    }
  
    Printf("counting with_cuts \n");

  for (int evt=0;evt<tin_with_cuts->GetEntries();evt++)
    {
      tin_with_cuts->GetEntry(evt);
      
      if (fabs(mc_2.eta) > 2.4) continue; //B mesons inside the detector region eta < 2.4
      if (fabs(mc_2.mass) < 5 || mc_2.mass > 5.5) continue;
      
      bool muon1Filter = fabs(mc_2.mu1eta)<2.4 && mc_2.mu1pt>2.8;
      bool muon2Filter = fabs(mc_2.mu2eta)<2.4 && mc_2.mu2pt>2.8;
 
      if (muon1Filter && muon2Filter) hist_passed->Fill(mc_2.pt);//count only the events with the muon selection above
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
  
  RooRealVar* eff2 = new RooRealVar("eff2","eff2",eff);
  eff2->setAsymError(-eff_lo,eff_hi);

  fin_no_cuts->Close();
  delete fin_no_cuts;

  fin_with_cuts->Close();
  delete fin_with_cuts;

  return eff2;
  
  /*
  RooAbsData* data;
  RooAbsPdf* model;
  RooRealVar* signal_res;

  data = w.data("data");
  build_pdf(w,channel,1);
  model = w.pdf("model");
  model->fitTo(*data,Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));

  TString dir = "";
  dir = "pt_bin_mass_fit/" + channel_to_ntuple_name(channel) + "_mc_" + TString::Format(VERSION) + "/" + channel_to_ntuple_name(channel) + "_mc_mass_fit_" + TString::Format("pt_from_%d_to_%d",(int)pt_min,(int)pt_max);
  
  plot_mass_fit(w,channel,dir);
  
  signal_res = w.var("n_signal");

  std::cout <<"SIGNAL_MC: "<< signal_res->getVal() << " " << signal_res->getAsymErrorLo() << " +" << signal_res->getAsymErrorHi() << std::endl;
  */
}

RooRealVar* bin_mass_fit(RooWorkspace& w, int channel, double pt_min, double pt_max)
{
  RooAbsData* data;
  RooAbsPdf* model;
  RooRealVar* signal_res;

  data = w.data("data");
  build_pdf(w,channel,0);
  model = w.pdf("model");
  model->fitTo(*data,Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));

  TString dir = "";
  dir = "pt_bin_mass_fit/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + "/" + channel_to_ntuple_name(channel) + "_mass_fit_" + TString::Format("pt_from_%d_to_%d",(int)pt_min,(int)pt_max);
  
  plot_mass_fit(w,channel,dir);

  //how to put the legend indicating each pt bin ??
  //change the plot_mass_fit to output a TCanvas, and write the legend on top after, and then have a function just to save the plots.
  //Legend(channel,(int)pt_bin_lo,(int)pt_bin_hi,1);
  
  signal_res = w.var("n_signal");
  
  return signal_res;
}

double pt_bin_mean(RooWorkspace& w, double pt_min, double pt_max)
{
  Printf("calculating bin mean \n");

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

void plot_mass_fit(RooWorkspace& w, int channel, TString directory)
{
  //plots the fit of the mass spectrum in the workspace, and saves it in directory, in .png

  RooRealVar mass = *(w.var("mass"));
  RooAbsData* data = w.data("data");
  RooAbsPdf* model = w.pdf("model");
  
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
  
  model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_signal"),LineColor(kRed),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kRed), VLines(), DrawOption("F"));
  
  if (channel==1 || channel==2 || channel==3 || channel==4 || channel==6)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_combinatorial_exp"),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
  else
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_combinatorial_bern"),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
  
  if (channel==1 || channel==3)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_jpsix"),LineColor(kViolet),LineWidth(2),LineStyle(7));
  if (channel==5)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_x3872"),LineColor(kOrange),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kOrange), VLines(), DrawOption("F"));
  
  frame_m->SetTitle("");
  frame_m->GetXaxis()->SetTitle(channel_to_xaxis_title(channel));
  frame_m->GetXaxis()->SetLabelFont(42);
  frame_m->GetXaxis()->SetLabelOffset(0.01);
  frame_m->GetXaxis()->SetTitleSize(0.06);
  frame_m->GetXaxis()->SetTitleOffset(1.09);
  frame_m->GetXaxis()->SetLabelFont(42);
  frame_m->GetXaxis()->SetLabelSize(0.055);
  frame_m->GetXaxis()->SetTitleFont(42);
  frame_m->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass.getMax()-mass.getMin())*1000./channel_to_nbins(channel)));
  frame_m->GetYaxis()->SetLabelFont(42);
  frame_m->GetYaxis()->SetLabelOffset(0.01);
  frame_m->GetYaxis()->SetTitleOffset(1.6);
  frame_m->GetYaxis()->SetTitleSize(0.05);
  frame_m->GetYaxis()->SetTitleFont(42);
  frame_m->GetYaxis()->SetLabelFont(42);
  frame_m->GetYaxis()->SetLabelSize(0.055);
  
  RooHist* pull_hist = frame_m->pullHist("theData","thePdf");
  
  RooPlot* pull_plot = mass.frame();
  pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
  pull_plot->SetTitle("");
  pull_plot->GetXaxis()->SetTitle(channel_to_xaxis_title(channel));
  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelOffset(0.01);
  pull_plot->GetXaxis()->SetTitleSize(0.06);
  pull_plot->GetXaxis()->SetTitleOffset(1.09);
  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelSize(0.055);
  pull_plot->GetXaxis()->SetTitleFont(42);
  pull_plot->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass.getMax()-mass.getMin())*1000./channel_to_nbins(channel)));
  pull_plot->GetYaxis()->SetLabelFont(42);
  pull_plot->GetYaxis()->SetLabelOffset(0.01);
  pull_plot->GetYaxis()->SetTitleOffset(1.14);
  pull_plot->GetYaxis()->SetTitleSize(0.06);
  pull_plot->GetYaxis()->SetTitleFont(42);
  pull_plot->GetYaxis()->SetLabelFont(42);
  pull_plot->GetYaxis()->SetLabelSize(0.055);

  TCanvas *c1 = canvasDressing("c1"); c1->cd();
  
  // TPad *p1 = new TPad("p1","p1",0,0,1,1);
  //p1->Draw();
   
  //  TPad *p1 = new TPad("p1","p1",0.03,0.27,0.99,0.99);
  TPad *p1 = new TPad("p1","p1",0.05,0.05,0.99,0.99);
  p1->SetBorderMode(0); 
  p1->Draw(); 
  
  /*
  TPad *p2 = new TPad("p2","p2",0.03,0.01,0.99,0.26); 
  p2->SetTopMargin(0.);    
  p2->SetBorderMode(0); 
  p2->SetTicks(1,2); 
  p2->Draw();
  */
  p1->cd(); frame_m->Draw(); histo_data->Draw("Esame"); Legend(channel,0,0,0);
  //p2->cd(); pull_plot->Draw();
  
  // c1->SaveAs(directory + ".root");
  c1->SaveAs(directory + ".png"); 
}

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

void build_pdf(RooWorkspace& w, int channel,int mc)
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
  
  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  
  //-----------------------------------------------------------------
  // signal PDF 
  RooRealVar m_mean("m_mean","m_mean",mass_peak,mass.getMin(),mass.getMax());
  RooRealVar m_sigma1("m_sigma1","m_sigma1",0.015,0.001,0.050);
  RooRealVar m_sigma2("m_sigma2","m_sigma2",0.030,0.001,0.100);
  RooRealVar m_fraction("m_fraction","m_fraction",0.5);
  RooGaussian m_gaussian1("m_gaussian1","m_gaussian1",mass,m_mean,m_sigma1);
  RooGaussian m_gaussian2("m_gaussian2","m_gaussian2",mass,m_mean,m_sigma2);
  RooAddPdf pdf_m_signal("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));
  
  // use single Gaussian for J/psi Ks and J/psi Lambda due to low statistics
  if (channel==3 || channel==6) {
    m_sigma2.setConstant(kTRUE);
    m_fraction.setVal(1.);
  }
  
  //-----------------------------------------------------------------
  // combinatorial background PDF (exponential or bernstean poly.)
  
  RooRealVar m_exp("m_exp","m_exp",-0.3,-4.,+4.);
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
      if(!mc)
	{
	  model = new RooAddPdf("model","model", RooArgList(pdf_m_signal, pdf_m_combinatorial_exp, pdf_m_jpsix), RooArgList(n_signal, n_combinatorial, n_jpsix));
	}
      else
	{
	  model = new RooAddPdf("model","model", RooArgList(pdf_m_signal, pdf_m_combinatorial_exp), RooArgList(n_signal, n_combinatorial));
	}
      break;
    case 2:// B0 -> J/psi K* 
    case 4://Bs -> J/psi phi
    case 6://Lambda_b -> J/psi Lambda
      model = new RooAddPdf("model","model", RooArgList(pdf_m_signal, pdf_m_combinatorial_exp), RooArgList(n_signal, n_combinatorial));
      break;
    case 5:// J/psi pipi
      model = new RooAddPdf("model","model", RooArgList(pdf_m_signal, pdf_m_combinatorial_bern, pdf_m_x3872), RooArgList(n_signal, n_combinatorial, n_x3872));
      break;
    }

  w.import(*model);
}

void read_data(RooWorkspace& w, TString filename,int channel,double pt_min, double pt_max)
{
  Printf("-----------importing data-----------\n");

  TFile* f = new TFile(filename);
  TNtupleD* _nt = (TNtupleD*)f->Get(channel_to_ntuple_name(channel));

  RooRealVar pt = *(w.var("pt"));
  RooRealVar pt_low("pt_low","pt_low",pt_min);
  RooRealVar pt_high("pt_high","pt_high",pt_max);
  RooFormulaVar ptcut("pt_cut","pt>pt_low && pt<pt_high",RooArgList(pt,pt_low,pt_high));
  
  RooDataSet* data = new RooDataSet("data","data",_nt,RooArgSet( *(w.var("mass")) , *(w.var("pt")), *(w.var("eta")) ), ptcut);
  
  w.import(*data);
  Printf("----------data imported----------\n");
}

void read_data_cut(RooWorkspace& w, RooAbsData* data)
{
  w.import(*data);
}

void set_up_workspace_variables(RooWorkspace& w, int channel)
{
  double mass_min, mass_max;
  double pt_min, pt_max;
  double eta_min, eta_max;

  pt_min=0;
  pt_max=400;

  eta_min=-5;
  eta_max=5;
  
  switch (channel) {
  default: 
  case 1:
    mass_min = 5.0; mass_max = 6.0;
    break;
  case 2:
    mass_min = 5.0; mass_max = 6.0;
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
  RooRealVar eta("eta","eta",eta_min,eta_max);

  w.import(mass);
  w.import(pt);
  w.import(eta);
}

void create_dir(std::vector<std::string> list)
{
  //to create the directories needed to save the output files, like .png and .root
  for(size_t i=0 ; i< list.size() ; ++i)
    {
      gSystem->Exec(("mkdir -p " + list[i]).c_str());
    }
}
