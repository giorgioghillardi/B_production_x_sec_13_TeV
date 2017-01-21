#include <RooHist.h>
#include <TSystem.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
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
#include <TChain.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <RooRealVar.h>
#include <RooProduct.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooBifurGauss.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooWorkspace.h>
#include <RooAddPdf.h>
#include <RooGenericPdf.h>
#include <RooCBShape.h>
#include <RooArgusBG.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <RooMCStudy.h>
#include <RooPlot.h>
#include <RooPlotable.h>
#include <RooThresholdCategory.h>
#include <Roo1DTable.h>
#include "TMath.h"
#include <TLegend.h>
#include "UserCode/B_production_x_sec_13_TeV/interface/myloop.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/format.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/plotDressing.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/channel.h"

using namespace RooFit;

// General fitting options
#define NUMBER_OF_CPU       1
#define VERSION             "v5"
#define BASE_DIR            "/lstore/cms/brunogal/input_for_B_production_x_sec_13_TeV/"

//////////////////////////////////////////////
// Definition of channel #                  //
// channel = 1: B+ -> J/psi K+              //
// channel = 2: B0 -> J/psi K*              //
// channel = 2: B0 -> J/psi Ks              //
// channel = 4: Bs -> J/psi phi             //
// channel = 5: Jpsi + pipi                 //
// channel = 6: Lambda_b -> Jpsi + Lambda   //
//////////////////////////////////////////////

void create_dir(std::vector<std::string> list);

void set_up_workspace_variables(RooWorkspace& w, int channel, double mass_min = 0.0 , double mass_max = 0.0);
void read_data(RooWorkspace& w, TString filename,int channel);
void read_data_cut(RooWorkspace& w, RooAbsData* data);
void build_pdf(RooWorkspace& w, int channel, std::string choice = "", std::string choice2 = "");

double var_mean_value(RooWorkspace& w, std::string var_name, double var_min, double var_max);
void plot_var_dist(RooWorkspace& w, std::string var_name, int channel, TString directory);
void plot_mass_fit(RooWorkspace& w, int channel, TString directory,int pt_high, int pt_low, double y_min, double y_max);

RooRealVar* bin_mass_fit(RooWorkspace& w, int channel, double pt_min, double pt_max, double y_min, double y_max, std::string choice = "", std::string choice2 = "", double mass_min = 0.0, double mass_max = 0.0);
double bin_systematics(RooWorkspace& ws, int channel, double pt_min, double pt_max, double y_min, double y_max, double signal_res, TString data_selection_input_file, int syst);

RooRealVar* prefilter_efficiency(int channel, double pt_min, double pt_max, double y_min, double y_max);
RooRealVar* reco_efficiency(int channel, double pt_min, double pt_max, double y_min, double y_max);
RooRealVar* branching_fraction(int channel);

void mc_study(RooWorkspace& w, int channel, double pt_min, double pt_max, double y_min, double y_max);
void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, std::vector<std::vector<double> > numbers, std::string caption);

//////////////////////////////////////////FUNCIONS////////////////////////////////////////////////////
void mc_study(RooWorkspace& w, int channel, double pt_min, double pt_max, double y_min, double y_max)
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
  
  RooFormulaVar cut("cut","pt>pt_low && pt<pt_high && ((y>y_low && y<y_high) || (y>-y_high && y<-y_low))", RooArgList(pt,pt_low,pt_high,y,y_low,y_high));
  
  data_cut = data_original->reduce(cut);
  read_data_cut(ws_cut,data_cut);

  build_pdf(ws_cut,channel);
  model_cut = ws_cut.pdf("model");
  model_cut->fitTo(*data_cut,Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));

  signal_res = ws_cut.var("n_signal");
  
  RooRealVar* mass = w.var("mass");
  
  RooMCStudy* mctoy = new RooMCStudy (*model_cut, *model_cut, *mass, "", "mhv"); 
  mctoy->generateAndFit(5000, data_cut->sumEntries());
  
  RooPlot* f_pull_signal = mctoy->plotPull(*signal_res, FitGauss(kTRUE));
  RooPlot* f_param_signal = mctoy->plotParam(*signal_res);
  RooPlot* f_error_signal = mctoy->plotError(*signal_res);
  RooPlot* f_nll = mctoy->plotNLL();
  
  TString dir_mc = "";
  dir_mc = "mc_study/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + "/" + channel_to_ntuple_name(channel) + "_mc_study_" + TString::Format("pt_from_%d_to_%d_y_from_%.2f_to_%.2f",(int)pt_min,(int)pt_max,y_min,y_max);
  
  TCanvas c;
  f_param_signal->Draw();
  c.SaveAs(dir_mc + "_param_signal.png");
  
  TCanvas c2;
  f_error_signal->Draw();
  c2.SaveAs(dir_mc + "_error_signal.png");
  
  TCanvas c3;
  f_pull_signal->Draw();
  c3.SaveAs(dir_mc + "_pull_signal.png");
  
  TCanvas c4;
  f_nll->Draw();
  c4.SaveAs(dir_mc + "_nll_signal.png");
}

RooRealVar* bin_mass_fit(RooWorkspace& w, int channel, double pt_min, double pt_max, double y_min, double y_max, std::string choice, std::string choice2, double mass_min, double mass_max)
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
  
  //if mass_min and mass_max are provided as input, it sets that mass windos in the workspace ws_cut.
  //if mass_min or mass_max are not provided it uses mass window from workspace w.
  if(mass_min!=0.0 && mass_max!=0.0)
    set_up_workspace_variables(ws_cut,channel,mass_min,mass_max);

  set_up_workspace_variables(ws_cut,channel);
  
  RooFormulaVar cut("cut","pt>pt_low && pt<pt_high && ((y>y_low && y<y_high) || (y>-y_high && y<-y_low))", RooArgList(pt,pt_low,pt_high,y,y_low,y_high));
  
  data_cut = data_original->reduce(cut);
  read_data_cut(ws_cut,data_cut);

  build_pdf(ws_cut,channel, choice, choice2);
  
  model_cut = ws_cut.pdf("model");
   
  model_cut->fitTo(*data_cut,Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));
  
  TString mass_info = "";
  if(mass_min!=0.0 && mass_max!=0.0)
    mass_info = TString::Format("_mass_from_%.2f_to_%.2f",mass_min,mass_max);
  
  TString syst_info = "";
  if(choice != "" && choice2 != "")
    syst_info = "_syst_" + choice + "_" + choice2;

  TString dir = "";
  
  dir = "mass_fits/" + channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + "/" + channel_to_ntuple_name(channel) + syst_info + "_mass_fit_" + TString::Format("pt_from_%d_to_%d_y_from_%.2f_to_%.2f",(int)pt_min,(int)pt_max,y_min,y_max) + mass_info;
  
  plot_mass_fit(ws_cut, channel, dir, (int) pt_max, (int) pt_min, y_min, y_max);

  signal_res = ws_cut.var("n_signal");
  return signal_res;
}

double bin_systematics(RooWorkspace& ws, int channel, double pt_min, double pt_max, double y_min, double y_max,double signal_res, TString data_selection_input_file, int syst)
{
  if(syst==0) return 0.;

  RooRealVar* fit_res;

  std::vector<std::string> background = {"2exp", "bern", "power"};
  std::vector<std::string> signal = {"crystal", "1gauss", "3gauss"};

  std::vector<double> mass_min(2);
  std::vector<double> mass_max(2);
  std::vector<std::string> mass_min_str(2);
  std::vector<std::string> mass_max_str(2);

  switch(channel)
    {
    case 1:
      mass_min[0] = 5.1;
      mass_min[1] = 5.0;
      mass_max[0] = 5.6;
      mass_max[1] = 6.0;

      mass_min_str[0] = "5.1";
      mass_min_str[1] = "5.0";
      mass_max_str[0] = "5.6";
      mass_max_str[1] = "6.0";
      break;
    case 2:
      mass_min[0] = 4.65;
      mass_min[1] = 5.05;
      mass_max[0] = 5.5;
      mass_max[1] = 6.0;

      mass_min_str[0] = "4.65";
      mass_min_str[1] = "5.05";
      mass_max_str[0] = "5.5";
      mass_max_str[1] = "6.0";
      break;

    case 4:
      mass_min[0] = 4.65;
      mass_min[1] = 5.11;
      mass_max[0] = 5.75;
      mass_max[1] = 6.2;

      mass_min_str[0] = "4.65";
      mass_min_str[1] = "5.11";
      mass_max_str[0] = "5.75";
      mass_max_str[1] = "6.2";
      break;
    }

  std::vector<double> signal_syst;
  std::vector<double> back_syst;
  std::vector<double> range_syst;
  signal_syst.reserve(4);
  back_syst.reserve(4);
  range_syst.reserve(3);
  signal_syst.push_back(signal_res);
  back_syst.push_back(signal_res);
  range_syst.push_back(signal_res);

  std::cout << std::endl << std::endl << std::endl << " signal_syst[0]: " << signal_syst[0] << std::endl << std::endl << std::endl;

  //Background Systematics
  for(int i=0; i<3; i++)
    {
      fit_res = bin_mass_fit(ws, channel, pt_min, pt_max, y_min, y_max, background[i], "background");
      back_syst.push_back((double)fit_res->getVal());
      std::cout << std::endl << std::endl << std::endl << " back_syst[" << i+1 << "]: " << back_syst[i+1] << std::endl << std::endl << std::endl;
    }

  //Signal Systematics
  for(int i=0; i<3; i++)
    {
      fit_res = bin_mass_fit(ws, channel, pt_min, pt_max, y_min, y_max, signal[i], "signal");
      signal_syst.push_back((double)fit_res->getVal());
      std::cout << std::endl << std::endl << std::endl << " signal_syst[" << i+1 << "]: " << signal_syst[i+1] << std::endl << std::endl << std::endl;
    }

  //Mass Range Systematics
  for(int i=0; i<2; i++)
    {
      RooWorkspace* ws1 = new RooWorkspace("ws1","Bmass");

      set_up_workspace_variables(*ws1,channel,mass_min[i],mass_max[1-i]);
      read_data(*ws1, data_selection_input_file,channel);
      
      fit_res = bin_mass_fit(*ws1, channel, pt_min, pt_max, y_min, y_max, "", "", mass_min[i], mass_max[1-i]);
      range_syst.push_back((double)fit_res->getVal());
      std::cout << std::endl << std::endl << std::endl << " range_syst[" << i+1 << "]: " << range_syst[i+1] << std::endl << std::endl << std::endl;
    }

  for(unsigned int i=0; i<signal_syst.size(); i++)
    std::cout << "signal_syst[" << i << "]: " << signal_syst[i] << std::endl;

  for(unsigned int i=0; i<back_syst.size(); i++)
    std::cout << "back_syst[" << i << "]: " << back_syst[i] << std::endl;
  for(unsigned int i=0; i<range_syst.size(); i++)
    std::cout << "range_syst[" << i << "]: " << range_syst[i] << std::endl;

  std::vector<double> deviation_signal;
  std::vector<double> deviation_back;
  std::vector<double> deviation_range;

  for(unsigned int i=0; i<signal_syst.size(); i++)
    deviation_signal.push_back(abs(signal_syst[i] - signal_syst[0]));   

  for(unsigned int i=0; i<back_syst.size(); i++)
    deviation_back.push_back(abs(back_syst[i] - back_syst[0]));   
 
  for(unsigned int i=0; i<range_syst.size(); i++)
    deviation_range.push_back(abs(range_syst[i] - range_syst[0]));    
  
  //find the max difference from the initial result, for each systematic error
  double signal_diff = *max_element(deviation_signal.begin(), deviation_signal.end());
  double bkg_diff = *max_element(deviation_back.begin(), deviation_back.end());
  double range_diff = *max_element(deviation_range.begin(), deviation_range.end());
  
  //sum the different systematic errors in quadrature to get the overall systematic error
  double overall_syst = sqrt(pow(signal_diff,2) + pow(bkg_diff,2) + pow(range_diff,2));
  
  return overall_syst;
}

double var_mean_value(RooWorkspace& w, std::string var_name, double var_min, double var_max)
{
  const char* var_name_str = var_name.c_str();

  RooRealVar var = *(w.var(var_name_str));
  RooRealVar var_low("var_low","var_low",var_min);
  RooRealVar var_high("var_high","var_high",var_max);
  RooAbsData* data_original;
  RooAbsData* data_cut;
  double mean_value;

  data_original = w.data("data");

  TString cut_formula = var_name + ">var_low && " + var_name + "<var_high";
  
  RooFormulaVar var_cut("var_cut", cut_formula, RooArgList(var,var_low,var_high));
  
  data_cut = data_original->reduce(var_cut);

  mean_value = (double) data_cut->meanVar(var)->getVal();
  
  return mean_value;
}

void plot_mass_fit(RooWorkspace& w, int channel, TString directory, int pt_high, int pt_low, double y_low, double y_high)
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
  
  model->paramOn(frame_m,Layout(0.505,0.89,0.57)); // Layout(0.45,0.85,0.57) //show all the parameters of the fit in the plot.
  
  model->plotOn(frame_m,Name("signal"),Precision(2E-4),Components("pdf_m_signal"),LineColor(8),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(8), VLines(), DrawOption("F"));
  
  if(channel==1 || channel==2 || channel==3 || channel==4 || channel==6)
    model->plotOn(frame_m,Name("combinatorial"),Precision(2E-4),Components("pdf_m_combinatorial_exp"),LineColor(9),LineWidth(2),LineStyle(2));
  else
    model->plotOn(frame_m,Name("combinatorial"),Precision(2E-4),Components("pdf_m_combinatorial_bern"),LineColor(kCyan+1),LineWidth(2),LineStyle(2));

  if(channel==2) // k pi swap component
    model->plotOn(frame_m,Name("kpiswap"),Precision(2E-4),Components("k_pi_swap"),LineColor(kViolet),LineWidth(2),LineStyle(7),FillStyle(3008),FillColor(kViolet), VLines(), DrawOption("F"));

  if(channel==1 || channel==3)
    model->plotOn(frame_m,Name("nonprompt"),Precision(2E-4),Components("pdf_m_nonprompt_erf"),LineColor(kViolet),LineWidth(2),LineStyle(2));
  
  if(channel==1)
    model->plotOn(frame_m,Name("jpsipi"),Precision(2E-4),Components("pdf_m_jpsipi"),LineColor(kRed),LineWidth(2),LineStyle(7),FillStyle(3008),FillColor(kRed), VLines(), DrawOption("F"));
  
  if(channel==5)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_x3872"),LineColor(kOrange),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kOrange), VLines(), DrawOption("F"));
  
  frame_m->SetTitle("");
  frame_m->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass.getMax()-mass.getMin())*1000./channel_to_nbins(channel)));
  frame_m->GetYaxis()->SetLabelFont(42);
  frame_m->GetYaxis()->SetLabelOffset(0.01);
  frame_m->GetYaxis()->SetTitleOffset(0.8);
  frame_m->GetYaxis()->SetTitleSize(0.06);
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
  pull_plot->GetYaxis()->SetLabelFont(42);
  pull_plot->GetYaxis()->SetLabelOffset(0.01);
  pull_plot->GetYaxis()->SetLabelSize(0.17);
  pull_plot->GetYaxis()->SetTitleOffset(1.6);
  pull_plot->GetYaxis()->SetTitleSize(0.17);
  pull_plot->GetYaxis()->SetTitleFont(42);
  pull_plot->GetYaxis()->SetNdivisions(305);
  
  TCanvas *c1 = canvasDressing("c1"); c1->cd();
  
  TPad *p1 = new TPad("p1","p1",0.05,0.27,0.99,0.99);
  //p1->SetLogy();
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
  p2->Draw();
  
  double chi_square = frame_m->chiSquare("thePdf","theData");

  RooAbsReal* nll = model->createNLL(*data);
  double log_likelihood= nll->getVal();
  std::stringstream ll_str;
  ll_str >> log_likelihood;
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

  bool show_legend = 0;

  if(show_legend)
    {
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
    
      if(data->sumEntries()>250)
	{
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
    }
 
  TLatex* tex8 = new TLatex(0.17, 0.5, Form("#chi^{2} = %.3lf", chi_square));
  tex8->SetNDC(kTRUE);
  tex8->SetTextFont(42);
  tex8->SetTextSize(0.035);
  tex8->Draw();

  p1->cd();
  frame_m->Draw();
  histo_data->Draw("Esame");
  Legend(channel, pt_low, pt_high, y_low, y_high, 1);

  //////////////////////////////////////////
  double x_1 = 0.505; //0.45
  double x_2 = 0.89; //0.85
  double y_1;
  double y_2= 0.88;
  double y_space = 0.05;

  int nitems = 6;
  y_1 = y_2 - y_space*nitems;

  TLegend *leg = new TLegend(x_1, y_1, x_2, y_2);
  leg->SetTextSize(0.04);
  leg->AddEntry(histo_data,"Data", "EPL");
  leg->AddEntry("thePdf","Total Fit", "L");
  
  switch(channel)
    {
    case 1:
      leg->AddEntry("signal","B^{#pm} #rightarrow J/#psi K^{#pm} Signal", "F");
      leg->AddEntry("combinatorial","Combinatorial background", "L");
      leg->AddEntry("jpsipi","B^{#pm} #rightarrow J/#psi #pi^{#pm} background", "F");
      leg->AddEntry("nonprompt","B^{#pm} #rightarrow J/#psi + hadrons background", "L");
      break;
    case 2:
      leg->AddEntry("signal","B^{0} #rightarrow J/#psi K^{*0} Signal", "F");
      leg->AddEntry("combinatorial","Combinatorial background", "L");
      leg->AddEntry("kpiswap","Swapped K^{#pm} #pi^{#mp} background", "F");
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
    case 6:
      break;      
    }

  leg->Draw();
  //////////////////////////////////////////
  
  p2->cd();
  pull_plot->Draw();
  
  c1->SaveAs(directory + ".png");
}

void plot_var_dist(RooWorkspace& w, std::string var_name, int channel, TString directory)
{
  const char* var_name_str = var_name.c_str();

  RooRealVar var = *(w.var(var_name_str));
  RooAbsData* data = w.data("data");

  TCanvas c2;
  TString hist_name = var_name + "_dist";
  TH1D* var_dist = (TH1D*)data->createHistogram(hist_name,var);
  
  var_dist->Draw();
  c2.SetLogy();
  c2.SaveAs(directory + ".png");
}

void build_pdf(RooWorkspace& w, int channel, std::string choice, std::string choice2)
{
  //choice is either signal or background
  //choice2 is the funntion to describe signal or background
  //if choice and choice 2 are not provided, the nominal fit is used.

  double mass_peak;

  RooRealVar mass = *(w.var("mass"));
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

  //Two Gaussians
  RooRealVar m_mean("m_mean","m_mean",mass_peak,mass_peak-0.09,mass_peak+0.09);
  RooRealVar m_sigma1("m_sigma1","m_sigma1",0.015,0.005,0.07);
  RooRealVar m_sigma2("m_sigma2","m_sigma2",0.030,0.001,0.100);
  RooRealVar m_fraction("m_fraction","m_fraction", 0.5, 0, 1);
  RooGaussian m_gaussian1("m_gaussian1","m_gaussian1",mass,m_mean,m_sigma1);
  RooGaussian m_gaussian2("m_gaussian2","m_gaussian2",mass,m_mean,m_sigma2);

  //Crystal Ball
  RooRealVar m_alpha("m_alpha", "m_alpha", mass_peak-0.015/2, mass_peak-0.08, mass_peak-0.003);
  RooRealVar m_n("m_n", "m_n", 2.7, 1, 7);
  RooCBShape m_crystal("m_crystal", "m_crystal", mass, m_mean, m_sigma1, m_alpha, m_n);

  //Three Gaussians
  RooRealVar m_sigma3("m_sigma3","m_sigma3",0.030,0.001,0.100);
  RooGaussian m_gaussian3("m_gaussian3","m_gaussian3",mass,m_mean,m_sigma3);
  RooRealVar m_fraction2("m_fraction2","m_fraction2",0.5);

  RooAddPdf* pdf_m_signal;
      
  // use single Gaussian for low statistics
  if(data->sumEntries()<250)
  {
    m_sigma2.setConstant(kTRUE);
    m_sigma3.setConstant(kTRUE);
    m_fraction.setVal(1.);
    m_fraction2.setVal(1.);
  }
  
  if(choice2=="signal" && choice=="crystal")
    {
      pdf_m_signal = new RooAddPdf("pdf_m_signal", "pdf_m_signal", RooArgList(m_crystal,m_gaussian2), RooArgList(m_fraction));
      m_sigma2.setConstant(kTRUE);
      m_fraction.setVal(1.);
    }
  else if(choice2=="signal" && choice=="1gauss")
    {
      pdf_m_signal = new RooAddPdf("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));
      m_sigma2.setConstant(kTRUE);
      m_fraction.setVal(1.);
    }
  else if(choice2=="signal" && choice=="3gauss")
    pdf_m_signal = new RooAddPdf("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2,m_gaussian3),RooArgList(m_fraction,m_fraction2));
  else //this is the nominal signal
    pdf_m_signal = new RooAddPdf("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));

  //-----------------------------------------------------------------
  // combinatorial background PDF
  
  //One Exponential
  RooRealVar m_exp("m_exp","m_exp",-0.3,-4.,0.);
  RooExponential pdf_m_combinatorial_exp("pdf_m_combinatorial_exp","pdf_m_combinatorial_exp",mass,m_exp);
  
  //Two Exponentials
  RooRealVar m_exp2("m_exp2","m_exp2",-0.3,-4.,0.);
  RooExponential pdf_m_combinatorial_exp2("pdf_m_combinatorial_exp2","pdf_m_combinatorial_exp2",mass,m_exp2);
  RooRealVar m_fraction_exp("m_fraction_exp", "m_fraction_exp", 0.5);

  //Bernstein
  RooRealVar m_par1("m_par1","m_par2",1.,0,+10.);
  RooRealVar m_par2("m_par2","m_par3",1.,0,+10.);
  RooRealVar m_par3("m_par3","m_par3",1.,0,+10.);
  
  RooBernstein pdf_m_combinatorial_bern("pdf_m_combinatorial_bern","pdf_m_combinatorial_bern",mass,RooArgList(RooConst(1.),m_par1,m_par2,m_par3));

  //Power Law (doesn't work)
  RooRealVar m_k("m_k", "m_k", -3., -1000., 0.);
  RooGenericPdf pdf_m_power("pdf_m_power", "pdf_m_power", "pow(mass, m_k)", RooArgSet(mass,m_k));

  //Argus
  RooRealVar m_arg1("m_arg1","m_arg1", -20., -100., -1.);
  RooRealVar m_arg2("m_arg2","m_arg2", -20., -100., -1.);
  RooArgusBG pdf_m_argus("pdf_m_argus", "pdf_m_argus", mass, m_arg1, m_arg2);

  RooAddPdf* pdf_m_combinatorial;

  if(choice2=="background" && choice=="2exp")
    pdf_m_combinatorial=new RooAddPdf("pdf_m_combinatorial","pdf_m_combinatorial",RooArgList(pdf_m_combinatorial_exp,pdf_m_combinatorial_exp2),RooArgList(m_fraction_exp));
  else if(choice2=="background" && choice=="bern")
    {
      pdf_m_combinatorial=new RooAddPdf("pdf_m_combinatorial","pdf_m_combinatorial",RooArgList(pdf_m_combinatorial_bern,pdf_m_combinatorial_exp),RooArgList(m_fraction_exp));
      m_exp.setConstant(kTRUE);
      m_fraction_exp.setVal(1.);    
    }
  else if(choice2=="background" && choice=="power")
    {
      pdf_m_combinatorial=new RooAddPdf("pdf_m_combinatorial","pdf_m_combinatorial",RooArgList(pdf_m_power,pdf_m_combinatorial_exp),RooArgList(m_fraction_exp));
      m_exp.setConstant(kTRUE);
      m_fraction_exp.setVal(1.);    
    }
  else //this is the nominal bkg
    {
      pdf_m_combinatorial=new RooAddPdf("pdf_m_combinatorial","pdf_m_combinatorial",RooArgList(pdf_m_combinatorial_exp,pdf_m_combinatorial_exp2),RooArgList(m_fraction_exp));
      m_exp2.setConstant(kTRUE);
      m_fraction_exp.setVal(1.);    
    }
  ////////////////////////////////////////////////////////////////////////////////////////////
  //The components below have no systematic variation yet, they are part of the nominal fit.//
  ////////////////////////////////////////////////////////////////////////////////////////////

  //K pi swap component, for channel 2. B0->jpsi K*0  
  RooRealVar sigma_swapped1("sigma_swapped1","sigma_swapped1", 0.0419);//0.178);
  RooRealVar sigma_swapped2("sigma_swapped2","sigma_swapped2", 0.1138);//0.047);
  RooRealVar mean("mean","mean", B0_MASS);
  RooGaussian swapped1("swapped1","swapped1",mass, mean, sigma_swapped1);
  RooGaussian swapped2("swapped2","swapped2",mass, mean, sigma_swapped2);
  RooRealVar r12("r12","r12", 0.655);//0.223);
  RooAddPdf k_pi_swap("k_pi_swap","k_pi_swap",RooArgSet(swapped1,swapped2),r12);
  //--------------------------------------------------------------------

  //jpsi_pi component, for channel 1.
  RooRealVar m_jpsipi_mean1("m_jpsipi_mean1","m_jpsipi_mean1",5.34693e+00,mass.getAsymErrorLo(),mass.getAsymErrorHi());
  RooRealVar m_jpsipi_mean2("m_jpsipi_mean2","m_jpsipi_mean2",5.46876e+00,mass.getAsymErrorLo(),mass.getAsymErrorHi());
  RooRealVar m_jpsipi_mean3("m_jpsipi_mean3","m_jpsipi_mean3",5.48073e+00,mass.getAsymErrorLo(),mass.getAsymErrorHi());
  RooRealVar m_jpsipi_sigma1l("m_jpsipi_sigma1l","m_jpsipi_sigma1l",2.90762e-02,0.010,0.150);
  RooRealVar m_jpsipi_sigma1r("m_jpsipi_sigma1r","m_jpsipi_sigma1r",6.52519e-02,0.010,0.150);
  RooRealVar m_jpsipi_sigma2("m_jpsipi_sigma2","m_jpsipi_sigma2",9.94712e-02,0.020,0.500);
  RooRealVar m_jpsipi_sigma3("m_jpsipi_sigma3","m_jpsipi_sigma3",3.30152e-01,0.020,0.500);
  RooRealVar m_jpsipi_fraction2("m_jpsipi_fraction2","m_jpsipi_fraction2",2.34646e-01,0.0,1.0);
  RooRealVar m_jpsipi_fraction3("m_jpsipi_fraction3","m_jpsipi_fraction3",1.14338e-01,0.0,1.0);

  m_jpsipi_mean1.setConstant(kTRUE);
  m_jpsipi_mean2.setConstant(kTRUE);
  m_jpsipi_mean3.setConstant(kTRUE);
  m_jpsipi_sigma1l.setConstant(kTRUE);
  m_jpsipi_sigma1r.setConstant(kTRUE);
  m_jpsipi_sigma2.setConstant(kTRUE);
  m_jpsipi_sigma3.setConstant(kTRUE);
  m_jpsipi_fraction2.setConstant(kTRUE);
  m_jpsipi_fraction3.setConstant(kTRUE);

  RooBifurGauss m_jpsipi_gaussian1("m_jpsipi_gaussian1","m_jpsipi_gaussian1",mass,m_jpsipi_mean1,m_jpsipi_sigma1l,m_jpsipi_sigma1r);
  RooGaussian m_jpsipi_gaussian2("m_jpsipi_gaussian2","m_jpsipi_gaussian2",mass,m_jpsipi_mean2,m_jpsipi_sigma2);
  RooGaussian m_jpsipi_gaussian3("m_jpsipi_gaussian3","m_jpsipi_gaussian3",mass,m_jpsipi_mean3,m_jpsipi_sigma3);

  RooAddPdf pdf_m_jpsipi("pdf_m_jpsipi","pdf_m_jpsipi",RooArgList(m_jpsipi_gaussian3,m_jpsipi_gaussian2,m_jpsipi_gaussian1),RooArgList(m_jpsipi_fraction3,m_jpsipi_fraction2));
  //--------------------------------------------------------------------

  //erfc component on channel 1 and 3
  //RooFormulaVar pdf_m_jpsix("pdf_m_jpsix","2.7*erfc((mass-5.14)/(0.5*0.08))",{mass});
  
  RooRealVar m_nonprompt_scale("m_nonprompt_scale","m_nonprompt_scale",4.74168e-02); //1.93204e-02, 0.001, 0.3);
  RooRealVar m_nonprompt_shift("m_nonprompt_shift","m_nonprompt_shift",5.14425); //5.14357e+00,5.12,5.16);
  RooGenericPdf pdf_m_nonprompt_erf("pdf_m_nonprompt_erf","pdf_m_nonprompt_erf","TMath::Erfc((mass-m_nonprompt_shift)/m_nonprompt_scale)", RooArgList(mass,m_nonprompt_scale,m_nonprompt_shift));
  
  m_nonprompt_scale.setConstant(kTRUE);
  m_nonprompt_shift.setConstant(kTRUE);
  
  //-------------------------------------------------------------------

  // X(3872) PDF, only for J/psi pipi fit
  RooRealVar m_x3872_mean("m_x3872_mean","m_x3872_mean",3.872,3.7,3.9);
  RooRealVar m_x3872_sigma("m_x3872_sigma","m_x3872_sigma",0.01,0.001,0.010);
  RooGaussian pdf_m_x3872("pdf_m_x3872","pdf_m_x3872",mass,m_x3872_mean,m_x3872_sigma);  
  //-----------------------------------------------------------------

  // full model  
  RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
  RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
  RooRealVar n_x3872("n_x3872","n_x3872",200.,0.,data->sumEntries());

  //RooRealVar n_jpsix("n_jpsix","n_jpsix",data->sumEntries(TString::Format("mass>4.9&&mass<5.14")),data->sumEntries(TString::Format("mass>4.9&&mass<5.14")),data->sumEntries());

  RooRealVar f_swap("f_swap","f_swap", 0.136765); //0.182273); //,0.140618); //for the k pi swap component of channel 2
  RooProduct n_swap("n_swap","n_swap",RooArgList(n_signal,f_swap));

  RooRealVar f_jpsipi("f_jpsipi","f_jpsipi",4.1E-5/1.026E-3,0.,0.1); //BF(jpsi_pi) = (4.1+-0.4)*10^-5 / BF(jpsi K) = (1.026+-0.031)*10^-3
  RooProduct n_jpsipi("n_jpsipi","n_jpsipi",RooArgList(n_signal,f_jpsipi));
  
  RooRealVar f_nonprompt("f_nonprompt","f_nonprompt",2.50259e-01,0.0,0.3);
  RooProduct n_nonprompt("n_nonprompt","n_nonprompt",RooArgList(n_signal,f_nonprompt));

  f_swap.setConstant(kTRUE);
  f_jpsipi.setConstant(kTRUE);
  f_nonprompt.setConstant(kTRUE);
  
  RooAddPdf* model;

  switch(channel)
    {
    default:
    case 1:// B+ -> J/psi K+
      model = new RooAddPdf("model","model", RooArgList(*pdf_m_signal, pdf_m_combinatorial_exp, pdf_m_nonprompt_erf, pdf_m_jpsipi),RooArgList(n_signal, n_combinatorial, n_nonprompt, n_jpsipi));
      break;
    case 2:// B0 -> J/psi K*
      model = new RooAddPdf("model","model", RooArgList(*pdf_m_signal, k_pi_swap, *pdf_m_combinatorial), RooArgList(n_signal, n_swap, n_combinatorial));
      break;
    case 3://B0 -> J/psi Ks
      model = new RooAddPdf("model","model", RooArgList(*pdf_m_signal, pdf_m_combinatorial_exp, pdf_m_nonprompt_erf),RooArgList(n_signal, n_combinatorial, n_nonprompt));
      break;
    case 4://Bs -> J/psi phi
    case 6://Lambda_b -> J/psi Lambda
      model = new RooAddPdf("model","model", RooArgList(*pdf_m_signal, *pdf_m_combinatorial), RooArgList(n_signal, n_combinatorial));
      break;
    case 5:// J/psi pipi
      model = new RooAddPdf("model","model", RooArgList(*pdf_m_signal, pdf_m_combinatorial_bern, pdf_m_x3872), RooArgList(n_signal, n_combinatorial, n_x3872));
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

void set_up_workspace_variables(RooWorkspace& w, int channel, double mass_min, double mass_max)
{
  double pt_min, pt_max;
  double y_min, y_max;

  pt_min=0;
  pt_max=400;

  y_min=-2.4;
  y_max=2.4;

  if(mass_min == 0.0 || mass_max == 0.0) //Default value, when no mass_min or mass_max are provided. Check the function declaration above.
    {
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
    }

  RooRealVar mass("mass","mass",mass_min,mass_max);
  RooRealVar pt("pt","pt",pt_min,pt_max);
  RooRealVar y("y", "y", y_min, y_max);

  w.import(mass);
  w.import(pt);
  w.import(y);
}

void create_dir(std::vector<std::string> list)
{
  //to create the directories needed to save the output files, like .png and .root
  for(size_t i=0 ; i< list.size() ; ++i)
    {
      gSystem->Exec(("mkdir -p " + list[i]).c_str());
    }
}

void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, std::vector<std::vector<double> > numbers, std::string caption)
{
  std::ofstream file;

  //Begin Document                                                                                                                               
  file.open(filename + ".tex");

  file << "\\documentclass{article}" << std::endl;
  file << "\\usepackage{cancel}" << std::endl;
  file << "\\usepackage{geometry}" << std::endl;
  file << "\\usepackage{booktabs}" << std::endl;
  file << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}" << std::endl;
  file << "\\title{B production at 13 TeV}" << std::endl;
  file << "\\begin{document}" << std::endl;
  file << "\\maketitle" << std::endl;

  // Create table                                                                                                                                
  file << "\\begin{table}[!h]" << std::endl;

  //setup table size                                                                                                                             
  std::string col="c";

  for(int i=1; i<n_col; i++)
    col+="|c";

  file << "\\begin{tabular}{"+col+"}" << std::endl;
  file << "\\toprule" << std::endl;

  for(int c=0; c<n_col-1; c++)
    file << col_name[c] << " & ";

  file << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;

  for(int i=1; i<n_lin; i++)
    {
      file << labels[i-1] << " & ";

      for(int c=1; c<n_col-1; c++)
	file << numbers[c-1][i-1] << " & ";

      file << numbers[n_col-2][i-1] << " \\\\" << std::endl;
    }

  file << "\\bottomrule" << std::endl;

  //End Table                                                                                                                                
  file << "\\end{tabular}" << std::endl;
  file << "\\caption{"+caption+"}" << std::endl;
  file << "\\end{table}" << std::endl;

  //End document
  file << "\\end{document}" << std::endl;

  system(("pdflatex " + filename + ".tex").c_str());
  system(("gnome-open " + filename + ".pdf").c_str());
}

//the input file must be produced with myloop_gen.cc to have the gen info. otherwise the signal needs to be extracted using a fit.
RooRealVar* prefilter_efficiency(int channel, double pt_min, double pt_max, double y_min, double y_max)
{
  TString mc_gen_input_file = TString::Format(BASE_DIR) + "reduced_myloop_gen_" + channel_to_ntuple_name(channel) + "_bfilter.root";
  TFile *fin = new TFile(mc_gen_input_file);
  
  TString ntuple_name = channel_to_ntuple_name(channel);
  TTree *tin = (TTree*)fin->Get(ntuple_name);
  
  //set up the variables needed
  double pt_b, eta_b, y_b, pt_mu1, pt_mu2, eta_mu1, eta_mu2;
  
  //read the ntuple from selected_data
  tin->SetBranchAddress("eta", &eta_b);
  tin->SetBranchAddress("y", &y_b);
  tin->SetBranchAddress("pt", &pt_b);
  tin->SetBranchAddress("mu1pt", &pt_mu1);
  tin->SetBranchAddress("mu2pt", &pt_mu2);
  tin->SetBranchAddress("mu1eta", &eta_mu1);
  tin->SetBranchAddress("mu2eta", &eta_mu2);
  
  //use histograms to count the events, and TEfficiency for efficiency, because it takes care of the errors and propagation
  TH1D* hist_tot = new TH1D("hist_tot","hist_tot",1,pt_min,pt_max);
  TH1D* hist_passed = new TH1D("hist_passed","hist_passed",1,pt_min,pt_max);

  //DEBUG: count using numbers and calculate the error with low statistics formula
  double total = 0;
  double passed = 0;

  for (int evt=0;evt<tin->GetEntries();evt++)
    {
      tin->GetEntry(evt);
      
      if (fabs(eta_b) > 2.4) continue; //B mesons inside the detector region eta < 2.4
      if (fabs(y_b)<y_min || fabs(y_b)>y_max) continue; // within the y binning
      if (pt_b<pt_min || pt_b>pt_max) continue; //within the pt bin
      
      hist_tot->Fill(pt_b);
      total ++;
      
      bool muon1Filter = fabs(eta_mu1) < 2.4 && pt_mu1>2.8;
      bool muon2Filter = fabs(eta_mu2) < 2.4 && pt_mu2>2.8;
      
      if (muon1Filter && muon2Filter) 
	{
	  hist_passed->Fill(pt_b);//count only the events with the muon selection above
	  passed ++;
	}
    }
  
  //debug
  //std::cout << "debug: passed: " << hist_passed->GetBinContent(1) << std::endl;
  //std::cout << "debug: total: " << hist_tot->GetBinContent(1) << std::endl;
  //std::cout << "debug: total number: " << total << std::endl;
  //std::cout << "debug: passed number: " << passed << std::endl;
  //----------------------------------
  
  //calculates the efficiency by dividing the histograms
  TEfficiency* efficiency = new TEfficiency(*hist_passed, *hist_tot);
  
  double eff;
  double eff_lo;
  double eff_hi;

  //DEBUG: either use Tefficiency, or use the formula.
  
  eff = efficiency->GetEfficiency(1);
  eff_lo = -(efficiency->GetEfficiencyErrorLow(1));
  eff_hi = efficiency->GetEfficiencyErrorUp(1);
  
  
  //eff = passed/total;
  //eff_lo = -eff * sqrt(((passed+1)*(total-passed+1))/((total+3)*(total+2)*(total+2)));
  //eff_hi =  eff * sqrt(((passed+1)*(total-passed+1))/((total+3)*(total+2)*(total+2)));
  

  RooRealVar* eff1 = new RooRealVar("eff1","eff1",eff);
  eff1->setAsymError(eff_lo,eff_hi);

  fin->Close();
  delete fin;

  return eff1; 
}

RooRealVar* reco_efficiency(int channel, double pt_min, double pt_max, double y_min, double y_max)
{
  //------------read monte carlo gen without cuts-----------------------------
  TString mc_input_no_cuts = TString::Format(BASE_DIR) + "reduced_myloop_gen_" + channel_to_ntuple_name(channel) + "_bmuonfilter.root";
  TFile *fin_no_cuts = new TFile(mc_input_no_cuts);
  TString ntuple_name = channel_to_ntuple_name(channel);
  TTree *tin_no_cuts = (TTree*)fin_no_cuts->Get(ntuple_name);
  
  //set up the variables needed
  double pt_b, eta_b, y_b, pt_mu1, pt_mu2, eta_mu1, eta_mu2;
  
  //read the ntuple from selected_data
  tin_no_cuts->SetBranchAddress("eta", &eta_b);
  tin_no_cuts->SetBranchAddress("y", &y_b);
  tin_no_cuts->SetBranchAddress("pt", &pt_b);
  tin_no_cuts->SetBranchAddress("mu1pt", &pt_mu1);
  tin_no_cuts->SetBranchAddress("mu2pt", &pt_mu2);
  tin_no_cuts->SetBranchAddress("mu1eta", &eta_mu1);
  tin_no_cuts->SetBranchAddress("mu2eta", &eta_mu2);
 
  //use histograms to count the events, and TEfficiency for efficiency, because it takes care of the errors and propagation
  TH1D* hist_tot = new TH1D("hist_tot","hist_tot",1,pt_min,pt_max);
  
  for (int evt=0; evt < tin_no_cuts->GetEntries(); evt++)
    {
      tin_no_cuts->GetEntry(evt);
      
      if (fabs(y_b)<y_min || fabs(y_b)>y_max) continue; // within the y binning
      
      if (fabs(eta_b) > 2.4) continue; //B mesons inside the detector region eta < 2.4
                 
      bool muon1Filter = fabs(eta_mu1) < 2.4 && pt_mu1 > 2.8;
      bool muon2Filter = fabs(eta_mu2) < 2.4 && pt_mu2 > 2.8;
 
      if (muon1Filter && muon2Filter) hist_tot->Fill(pt_b);//count only the events with the muon selection above
    }
      
    //--------------------------------read monte carlo with cuts------------------------
    TString mc_input_with_cuts = TString::Format(BASE_DIR) + "reduced_myloop_new_mc_truth_" + channel_to_ntuple_name(channel) + "_with_cuts.root";
    TFile *fin_with_cuts = new TFile(mc_input_with_cuts);
    TTree *tin_with_cuts = (TTree*)fin_with_cuts->Get(channel_to_ntuple_name(channel));
   
    //read the ntuple
    tin_with_cuts->SetBranchAddress("eta", &eta_b);
    tin_with_cuts->SetBranchAddress("y", &y_b);
    tin_with_cuts->SetBranchAddress("pt", &pt_b);
    tin_with_cuts->SetBranchAddress("mu1pt", &pt_mu1);
    tin_with_cuts->SetBranchAddress("mu2pt", &pt_mu2);
    tin_with_cuts->SetBranchAddress("mu1eta", &eta_mu1);
    tin_with_cuts->SetBranchAddress("mu2eta", &eta_mu2);
 
    TH1D* hist_passed = new TH1D("hist_passed","hist_passed",1,pt_min,pt_max);

    for (int evt=0; evt < tin_with_cuts->GetEntries(); evt++)
      {
	tin_with_cuts->GetEntry(evt);
	
	if (fabs(y_b)<y_min || fabs(y_b)>y_max) continue; // within the y binning
	
	if (fabs(eta_b) > 2.4) continue; //B mesons inside the detector region eta < 2.4
	
	bool muon1Filter = fabs(eta_mu1) < 2.4 && pt_mu1 > 2.8;
	bool muon2Filter = fabs(eta_mu2) < 2.4 && pt_mu2 > 2.8;
	
	if (muon1Filter && muon2Filter) hist_passed->Fill(pt_b);//count only the events with the muon selection above
      }
    
    //calculates the efficiency by dividing the histograms
    TEfficiency* efficiency = new TEfficiency(*hist_passed, *hist_tot);
    
    double eff;
    double eff_lo;
    double eff_hi;
    
    eff = efficiency->GetEfficiency(1);
    eff_lo = -(efficiency->GetEfficiencyErrorLow(1));
    eff_hi = efficiency->GetEfficiencyErrorUp(1);
        
    RooRealVar* eff2 = new RooRealVar("eff2","eff2",eff);
    eff2->setAsymError(eff_lo,eff_hi);
    
    fin_no_cuts->Close();
    delete fin_no_cuts;
    
    fin_with_cuts->Close();
    delete fin_with_cuts;
    
    return eff2;
}

RooRealVar* branching_fraction(int channel)
{
  RooRealVar* b_fraction = new RooRealVar("b_fraction","b_fraction",1);
  b_fraction->setError(1);
  
  RooRealVar* jpsi_to_mu_mu = new RooRealVar("jpsi","jpsi",5.93e-2);
  jpsi_to_mu_mu->setError(0.06e-2);

  RooRealVar* bu_to_jpsi_ku = new RooRealVar("bu","bu",1.026e-3);
  bu_to_jpsi_ku->setError(0.031e-3);

  RooRealVar* bd_to_jpsi_kstar = new RooRealVar("bd","bd",1.32e-3);
  bd_to_jpsi_kstar->setError(6e-5);
  
  RooRealVar* kstar_to_k_pi = new RooRealVar("kstar","kstar",0.99901);
  kstar_to_k_pi->setError(9e-5);

  RooRealVar* bs_to_jpsi_phi = new RooRealVar("bs","bs",1.08e-3);
  bs_to_jpsi_phi->setError(9e-5);

  RooRealVar* phi_to_k_k = new RooRealVar("phi","phi",48.9e-2);
  phi_to_k_k->setError(5e-3);

  double err =1;
  
  switch (channel) 
    {
    default:
    case 1:
      b_fraction->setVal(bu_to_jpsi_ku->getVal() * jpsi_to_mu_mu->getVal());
      err = b_fraction->getVal()*sqrt( pow(bu_to_jpsi_ku->getError()/bu_to_jpsi_ku->getVal(),2) + pow(jpsi_to_mu_mu->getError()/jpsi_to_mu_mu->getVal(),2) );
      break;      
    case 2:
      b_fraction->setVal( bd_to_jpsi_kstar->getVal() * kstar_to_k_pi->getVal() * jpsi_to_mu_mu->getVal());
      err = b_fraction->getVal() * sqrt( pow(bd_to_jpsi_kstar->getError()/bd_to_jpsi_kstar->getVal(),2) + pow(kstar_to_k_pi->getError()/kstar_to_k_pi->getVal(),2) + pow(jpsi_to_mu_mu->getError()/jpsi_to_mu_mu->getVal(),2) );
      break;
    case 4:
      b_fraction->setVal( bs_to_jpsi_phi->getVal() * phi_to_k_k->getVal() * jpsi_to_mu_mu->getVal());
      err = b_fraction->getVal() * sqrt(pow(bs_to_jpsi_phi->getError()/bs_to_jpsi_phi->getVal(),2) + pow(phi_to_k_k->getError()/phi_to_k_k->getVal(),2) + pow(jpsi_to_mu_mu->getError()/jpsi_to_mu_mu->getVal(),2));
      break;
    }

  b_fraction->setError(err);
  
  return b_fraction;
}
