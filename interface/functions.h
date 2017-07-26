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
#include "TVectorD.h"
#include <TLegend.h>
#include "TF1.h"
#include "TPaveStats.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/myloop.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/format.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/plotDressing.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/channel.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/bins.h"

using namespace RooFit;

#define LUMINOSITY          2.71
#define NUMBER_OF_CPU       1
#define VERSION             "v15"
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
void setup_bins(TString measure, int channel, TString bins, TString* var1_name, int* n_var1_bins, TString* var2_name, int* n_var2_bins, double** var1_bins, double** var2_bins);

double var_mean_value(RooWorkspace& w, TString var1_name, double var1_min, double var1_max, TString var2_name, double var2_min, double var2_max);
void plot_var_dist(RooWorkspace& w, std::string var_name, int channel, TString directory);
void plot_mass_fit(RooWorkspace& w, int channel, TString directory,int pt_high, int pt_low, double y_min, double y_max);
void plot_eff(TString measure, TString eff_name, int channel, int n_var1_bins, TString var2_name, double var2_min, double var2_max, TString x_axis_name, TString b_title, double* var1_bin_centre, double* var1_bin_centre_lo, double* var1_bin_centre_hi, double* eff_array, double* eff_err_lo_array, double* eff_err_hi_array);

RooRealVar* bin_mass_fit(RooWorkspace& w, int channel, double pt_min, double pt_max, double y_min, double y_max, std::string choice = "", std::string choice2 = "", double mass_min = 0.0, double mass_max = 0.0, Bool_t verb = kFALSE);

RooRealVar* prefilter_efficiency(int channel, double pt_min, double pt_max, double y_min, double y_max);
RooRealVar* reco_efficiency(int channel, double pt_min, double pt_max, double y_min, double y_max);
RooRealVar* branching_fraction(int channel);

//void read_vector(TString measure, int channel, TString vector, TString var1_name , TString var2_name, int n_var1_bins, int n_var2_bins,  double* var1_bins, double* var2_bins, double* array, double* err_lo = NULL, double* err_hi = NULL);
void read_vector(int channel, TString vector, TString var1_name , TString var2_name, int n_var1_bins, int n_var2_bins,  double* var1_bins, double* var2_bins, double* array, double* err_lo = NULL, double* err_hi = NULL);
void print_table(TString title, int n_var1_bins, int n_var2_bins, TString var1_name, TString var2_name, double* var1_bin_edges, double* var2_bin_edges, double* array, double* stat_err_lo, double* stat_err_hi, double* syst_err_lo = NULL, double* syst_err_hi = NULL);
void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, std::vector<std::vector<double> > numbers, std::string caption);

void mc_study(RooWorkspace& w, int channel, double pt_min, double pt_max, double y_min, double y_max);

//////////////////////////////////////////FUNCIONS////////////////////////////////////////////////////
void create_dir(std::vector<std::string> list)
{
  //to create the directories needed to save the output files, like .png and .root
  for(size_t i=0 ; i< list.size() ; ++i)
    {
      gSystem->Exec(("mkdir -p " + list[i]).c_str());
    }
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
      case 2:
	mass_min = 5.05; mass_max = 5.5;
	break;
      case 3:
	mass_min = 5.0; mass_max = 6.0;
	break;
      case 4:
	mass_min = 5.2; mass_max = 5.52;
	break;
      case 5:
	mass_min = 3.6; mass_max = 4.0;
	break;
      case 6:
	mass_min = 5.3; mass_max = 6.3;
	break;
      case 7:
	mass_min = 5.9; mass_max = 6.9;
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

void build_pdf(RooWorkspace& w, int channel, std::string choice, std::string choice2)
{
  //choice is either signal or background
  //choice2 is the function to describe signal or background
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
  case 7:
    mass_peak = BC_MASS;
    break;

  }
  
  double n_signal_initial = data->sumEntries(TString::Format("abs(mass-%g)<0.015",mass_peak)) - data->sumEntries(TString::Format("abs(mass-%g)<0.030&&abs(mass-%g)>0.015",mass_peak,mass_peak));
  
  if(n_signal_initial<0)
    n_signal_initial=1;

  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  
  //-----------------------------------------------------------------
  // signal PDF 

  //Two Gaussians
  RooRealVar m_mean("m_mean","m_mean",mass_peak,mass_peak-0.09,mass_peak+0.09);
  RooRealVar m_sigma1("m_sigma1","m_sigma1",0.010,0.009,0.200);
  RooRealVar m_sigma2("m_sigma2","m_sigma2",0.005,0.004,0.100);
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
  if(n_signal_initial < 500)
  {
    //m_sigma2.setConstant(kTRUE);
    //m_sigma3.setConstant(kTRUE);
    m_fraction.setVal(1.);
    m_fraction2.setVal(1.);
  }
  
  if(choice2=="signal" && choice=="crystal")
    {
      pdf_m_signal = new RooAddPdf("pdf_m_signal", "pdf_m_signal", RooArgList(m_crystal,m_gaussian2), RooArgList(m_fraction));
      m_sigma2.setConstant(kTRUE);
      m_fraction.setVal(1.);
    }
  else 
    if(choice2=="signal" && choice=="1gauss")
      {
	pdf_m_signal = new RooAddPdf("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));
	m_sigma2.setConstant(kTRUE);
	m_fraction.setVal(1.);
      }
    else
      if(choice2=="signal" && choice=="3gauss")
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

  RooAddPdf* pdf_m_combinatorial;

  if(choice2=="background" && choice=="2exp")
    pdf_m_combinatorial=new RooAddPdf("pdf_m_combinatorial","pdf_m_combinatorial",RooArgList(pdf_m_combinatorial_exp,pdf_m_combinatorial_exp2),RooArgList(m_fraction_exp));
  else
    if(choice2=="background" && choice=="bern")
      {
	pdf_m_combinatorial=new RooAddPdf("pdf_m_combinatorial","pdf_m_combinatorial",RooArgList(pdf_m_combinatorial_bern,pdf_m_combinatorial_exp),RooArgList(m_fraction_exp));
	m_exp.setConstant(kTRUE);
	m_fraction_exp.setVal(1.);    
      }
    else 
      if(choice2=="background" && choice=="power")
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
  RooRealVar sigma_swapped1("sigma_swapped1","sigma_swapped1", 0.0419);
  RooRealVar sigma_swapped2("sigma_swapped2","sigma_swapped2", 0.1138);
  
  RooGaussian swapped1("swapped1","swapped1",mass, m_mean, sigma_swapped1);
  RooGaussian swapped2("swapped2","swapped2",mass, m_mean, sigma_swapped2);
  RooRealVar r12("r12","r12", 0.655);
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
  RooRealVar m_nonprompt_scale("m_nonprompt_scale","m_nonprompt_scale",4.74168e-02, 0, 1); //1.93204e-02, 0.001, 0.3);
  RooRealVar m_nonprompt_shift("m_nonprompt_shift","m_nonprompt_shift",5.14425); //5.14357e+00,5.12,5.16);
  RooGenericPdf pdf_m_nonprompt_erf("pdf_m_nonprompt_erf","pdf_m_nonprompt_erf","TMath::Erfc((mass-m_nonprompt_shift)/m_nonprompt_scale)", RooArgList(mass,m_nonprompt_scale,m_nonprompt_shift));
  
  m_nonprompt_shift.setConstant(kTRUE);
  m_nonprompt_scale.setConstant(kTRUE);

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

  RooRealVar f_swap("f_swap","f_swap", 0.136765); //for the k pi swap component of channel 2
  //set n_swap like n_jpsipi in case we want to count the k_pi_swap component as background.
  
  RooRealVar f_jpsipi("f_jpsipi","f_jpsipi",4.1E-5/1.026E-3,0.,0.1); //BF(jpsi_pi) = (4.1+-0.4)*10^-5 / BF(jpsi K) = (1.026+-0.031)*10^-3
  RooProduct n_jpsipi("n_jpsipi","n_jpsipi",RooArgList(n_signal,f_jpsipi));
  
  RooRealVar f_nonprompt("f_nonprompt","f_nonprompt",2.50259e-01,0,1);
  RooProduct n_nonprompt("n_nonprompt","n_nonprompt",RooArgList(n_signal,f_nonprompt));

  f_swap.setConstant(kTRUE);
  f_jpsipi.setConstant(kTRUE);
  //f_nonprompt.setConstant(kTRUE);
    
  RooAddPdf* model;
  RooAddPdf* pdf_m_signal_copy;

  switch(channel)
    {
    default:
    case 1:// B+ -> J/psi K+
      model = new RooAddPdf("model","model", RooArgList(*pdf_m_signal, *pdf_m_combinatorial, pdf_m_nonprompt_erf, pdf_m_jpsipi),RooArgList(n_signal, n_combinatorial, n_nonprompt, n_jpsipi));
      break;
    case 2:// B0 -> J/psi K*
      pdf_m_signal_copy = new RooAddPdf(*pdf_m_signal, "pdf_m_signal_copy");
      pdf_m_signal = new RooAddPdf("pdf_m_signal","pdf_m_signal",RooArgList(k_pi_swap,*pdf_m_signal_copy),RooArgList(f_swap));
      model = new RooAddPdf("model","model", RooArgList(*pdf_m_signal, *pdf_m_combinatorial), RooArgList(n_signal, n_combinatorial));
      break;
    case 3://B0 -> J/psi Ks
      model = new RooAddPdf("model","model", RooArgList(*pdf_m_signal, *pdf_m_combinatorial, pdf_m_nonprompt_erf),RooArgList(n_signal, n_combinatorial, n_nonprompt));
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

void setup_bins(TString measure, int channel, TString bins, TString* var1_name, int* n_var1_bins, TString* var2_name, int* n_var2_bins, double** var1_bins, double** var2_bins)
{
  int n_pt_bins=1;
  double* pt_bins = total_pt_bin_edges;
  
  int n_y_bins=1;
  double* y_bins = total_y_bin_edges;
  
  switch (channel)
    {
    default:
    case 1:
      pt_bins = ntkp_pt_bins;
      n_pt_bins = (sizeof(ntkp_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntkp_y_bins;
      n_y_bins = (sizeof(ntkp_y_bins) / sizeof(double)) - 1 ;
      break;

    case 2:
      pt_bins = ntkstar_pt_bins;
      n_pt_bins = (sizeof(ntkstar_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntkstar_y_bins;
      n_y_bins = (sizeof(ntkstar_y_bins) / sizeof(double)) - 1 ;
      break;

    case 3:
      pt_bins = ntks_pt_bins;
      n_pt_bins = (sizeof(ntks_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntks_y_bins;
      n_y_bins = (sizeof(ntks_y_bins) / sizeof(double)) - 1 ;
      break;

    case 4:
      pt_bins = ntphi_pt_bins;
      n_pt_bins = (sizeof(ntphi_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntphi_y_bins;
      n_y_bins = (sizeof(ntphi_y_bins) / sizeof(double)) - 1 ;
      break;

    case 5:
      pt_bins = ntmix_pt_bins;
      n_pt_bins = (sizeof(ntmix_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntmix_y_bins;
      n_y_bins = (sizeof(ntmix_y_bins) / sizeof(double)) - 1 ;
      break;

    case 6:
      pt_bins = ntlambda_pt_bins;
      n_pt_bins = (sizeof(ntlambda_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ntlambda_y_bins;
      n_y_bins = (sizeof(ntlambda_y_bins) / sizeof(double)) - 1 ;
      break;
    }
  
  if(measure == "ratio")
    {
      pt_bins = ratio_pt_bins;
      n_pt_bins = (sizeof(ratio_pt_bins) / sizeof(double)) - 1 ;
      y_bins = ratio_y_bins;
      n_y_bins = (sizeof(ratio_y_bins) / sizeof(double)) - 1 ;
    }
  
  if(bins == "full")
    {
      *n_var1_bins= 1;
      *var1_name = "pt";
      *var1_bins = total_pt_bin_edges;

      *n_var2_bins= 1;
      *var2_name = "y";
      *var2_bins = total_y_bin_edges;
    }
  else
    if(bins == "pt")
      {
        *n_var1_bins = n_pt_bins;
        *var1_name = "pt";
        *var1_bins = pt_bins;
	
        *n_var2_bins = 1;
        *var2_name = "y";
        *var2_bins = total_y_bin_edges;
      }
    else
      if(bins == "pt_y")
        {
          *n_var1_bins = n_pt_bins;
          *var1_name = "pt";
          *var1_bins = pt_bins;

          *n_var2_bins = n_y_bins;
          *var2_name = "y";
          *var2_bins = y_bins;
        }
      else
        if(bins == "y")
          {
            *n_var1_bins = n_y_bins;
            *var1_name = "y";
            *var1_bins = y_bins;

            *n_var2_bins = 1;
            *var2_name = "pt";
            *var2_bins = total_pt_bin_edges;
          }
        else
          if(bins == "y_pt")
            {
              *n_var1_bins = n_y_bins;
              *var1_name = "y";
              *var1_bins = y_bins;

              *n_var2_bins = n_pt_bins;
              *var2_name = "pt";
              *var2_bins = pt_bins;
            }
          else
            {
              printf("ERROR: The variables you asked for are not deffined. Only full, pt, y, pt_y, y_pt are deffined. Please check in the code.\n");
              return;
            }
}

double var_mean_value(RooWorkspace& w, TString var1_name, double var1_min, double var1_max, TString var2_name, double var2_min, double var2_max)
{
  RooRealVar var1 = *(w.var(var1_name));
  RooRealVar var1_low("var1_low","var1_low",var1_min);
  RooRealVar var1_high("var1_high","var1_high",var1_max);
  RooRealVar var2 = *(w.var(var2_name));
  RooRealVar var2_low("var2_low","var2_low",var2_min);
  RooRealVar var2_high("var2_high","var2_high",var2_max);
  
  RooAbsData* data_original;
  RooAbsData* data_cut;
  double mean_value;

  data_original = w.data("data");

  TString cut_formula = var1_name + ">var1_low && " + var1_name + "<var1_high && " + var2_name + ">var2_low && " + var2_name + "<var2_high";
  
  RooFormulaVar var_cut("var_cut", cut_formula, RooArgList(var1,var1_low,var1_high,var2,var2_low,var2_high));
  
  data_cut = data_original->reduce(var_cut);

  mean_value = (double) data_cut->meanVar(var1)->getVal();
  
  return mean_value;
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
  
  //model->paramOn(frame_m,Layout(0.60,1.00,0.57));
  
  model->plotOn(frame_m,Name("signal"),Precision(2E-4),Components("pdf_m_signal"),LineColor(8),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(8), VLines(), DrawOption("F"));
  
  if(channel==1 || channel==2 || channel==3 || channel==4 || channel==6)
    model->plotOn(frame_m,Name("combinatorial"),Precision(2E-4),Components("pdf_m_combinatorial"),LineColor(9),LineWidth(2),LineStyle(2));
  else
    model->plotOn(frame_m,Name("combinatorial"),Precision(2E-4),Components("pdf_m_combinatorial_bern"),LineColor(kCyan+1),LineWidth(2),LineStyle(2));

  if(channel==2) // k pi swap component
    model->plotOn(frame_m,Name("kpiswap"),Precision(2E-4),Components("k_pi_swap"),LineColor(kOrange),LineWidth(2),LineStyle(7),FillStyle(3008),FillColor(kOrange), VLines(), DrawOption("F"));

  if(channel==1 || channel==3)
    model->plotOn(frame_m,Name("nonprompt"),Precision(2E-4),Components("pdf_m_nonprompt_erf"),LineColor(kViolet),LineWidth(2),LineStyle(2));
  
  if(channel==1)
    model->plotOn(frame_m,Name("jpsipi"),Precision(2E-4),Components("pdf_m_jpsipi"),LineColor(kRed),LineWidth(2),LineStyle(7),FillStyle(3008),FillColor(kRed), VLines(), DrawOption("F"));
  
  if(channel==5)
    model->plotOn(frame_m,Precision(2E-4),Components("pdf_m_x3872"),LineColor(kOrange),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kOrange), VLines(), DrawOption("F"));
  
  frame_m->SetTitle("");
  frame_m->GetXaxis()->SetLabelSize(0.0);
  frame_m->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass.getMax()-mass.getMin())*1000./channel_to_nbins(channel)));
  frame_m->GetYaxis()->SetTitleFont(42);
  frame_m->GetYaxis()->SetTitleSize(0.05);
  frame_m->GetYaxis()->SetTitleOffset(0.9);
  frame_m->GetYaxis()->SetLabelFont(42);
  frame_m->GetYaxis()->SetLabelSize(0.045);
  frame_m->GetYaxis()->SetLabelOffset(0.01);
    
  RooHist* pull_hist = frame_m->pullHist("theData","thePdf");
  
  RooPlot* pull_plot = mass.frame();
  
  RooGenericPdf* line_ref = new RooGenericPdf("ref_0", "ref_0", RooConst(0.));
  line_ref->plotOn(pull_plot, LineStyle(7), LineColor(13), LineWidth(2));

  TLine *line = new TLine(mass.getMin(),0,mass.getMax(),0);
  line->SetLineColor(kBlack);
  line->Draw();

  pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
  pull_plot->SetTitle("");
  pull_plot->GetXaxis()->SetTitle(channel_to_xaxis_title(channel));
  pull_plot->GetXaxis()->SetTitleFont(42);
  pull_plot->GetXaxis()->SetTitleSize(0.17);
  pull_plot->GetXaxis()->SetTitleOffset(1.09);

  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelSize(0.17);
  pull_plot->GetXaxis()->SetLabelOffset(0.01);
  pull_plot->GetXaxis()->SetTickLength(0.15);
  
  pull_plot->GetYaxis()->SetTitle("Pull hist");
  pull_plot->GetYaxis()->SetTitleFont(42);  
  pull_plot->GetYaxis()->SetTitleSize(0.20);
  pull_plot->GetYaxis()->SetTitleOffset(0.15);

  pull_plot->GetYaxis()->SetLabelFont(42);
  pull_plot->GetYaxis()->SetLabelSize(0.14);
  pull_plot->GetYaxis()->SetLabelOffset(0.005);
  
  pull_plot->GetYaxis()->SetNdivisions(305);
  
  TCanvas *c1 = canvasDressing("c1"); c1->cd();
  
  TPad *p1 = new TPad("p1","p1",0.0,0.27,0.99,0.99);
  //p1->SetLogy();
  p1->SetBorderMode(0); 
  p1->SetFrameBorderMode(0); 
  p1->SetBorderSize(2);
  p1->SetBottomMargin(0.0);
  p1->Draw(); 
     
  TPad *p2 = new TPad("p2","p2",0.0,0.065,0.99,0.24);
  p2->SetTopMargin(0.);    
  p2->SetBorderMode(0);
  p2->SetBorderSize(2); 
  p2->SetFrameBorderMode(0); 
  p2->Draw();
  
  double chi_square = frame_m->chiSquare("thePdf","theData");
 
  TLatex* tex8 = new TLatex(0.17, 0.5, Form("#chi^{2} = %.3lf", chi_square));
  tex8->SetNDC(kTRUE);
  tex8->SetTextFont(42);
  tex8->SetTextSize(0.035);
  tex8->Draw();

  ///////////////////
  //legend position//
  ///////////////////
  double x_1 = 0.60;
  double x_2 = 0.95;
  double y_1;
  double y_2= 0.88;
  double y_space = 0.05;

  int nitems = 6;
  y_1 = y_2 - y_space*nitems;

  p1->cd();

  model->paramOn(frame_m,Layout(x_1,x_2,y_1));
  frame_m->Draw();

  histo_data->Draw("Esame");
  Legend(channel, pt_low, pt_high, y_low, y_high, 1); //to show the bin size and other common elements

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
      leg->AddEntry("kpiswap","Swapped K^{#pm} #pi^{#mp} signal", "F");
      break;
    case 3:
      break;
    case 4:
      leg->AddEntry("signal","B_{s}^{0} #rightarrow J/#psi #Phi Signal", "F");
      leg->AddEntry("combinatorial","Combinatorial background", "L");
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

RooRealVar* bin_mass_fit(RooWorkspace& w, int channel, double pt_min, double pt_max, double y_min, double y_max, std::string choice, std::string choice2, double mass_min, double mass_max, Bool_t verb)
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
  
  //if mass_min and mass_max are provided as input, it sets that mass window in the workspace ws_cut.
  //if mass_min or mass_max are not provided it uses mass window from workspace w.
  if(mass_min!=0.0 && mass_max!=0.0)
    set_up_workspace_variables(ws_cut,channel,mass_min,mass_max);

  set_up_workspace_variables(ws_cut,channel);
  
  RooFormulaVar cut("cut","pt>pt_low && pt<pt_high && ((y>y_low && y<y_high) || (y>-y_high && y<-y_low))", RooArgList(pt,pt_low,pt_high,y,y_low,y_high));
  
  data_cut = data_original->reduce(cut);
  read_data_cut(ws_cut,data_cut);

  build_pdf(ws_cut,channel, choice, choice2);
  
  model_cut = ws_cut.pdf("model");
   
  model_cut->fitTo(*data_cut,Verbose(verb),Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));
  
  TString mass_info = "";
  if(mass_min!=0.0 && mass_max!=0.0)
    mass_info = TString::Format("_mass_from_%.2f_to_%.2f",mass_min,mass_max);
  
  TString base_dir = TString::Format(VERSION) + "/mass_fits/"; 

  TString syst_info = "";
  
if(choice != "" && choice2 != "")
  {
    base_dir += "syst/"; 
    syst_info = "_syst_" + choice + "_" + choice2;
  }
  TString dir = "";
  
  dir = base_dir + channel_to_ntuple_name(channel) + "/" + channel_to_ntuple_name(channel) + syst_info + "_mass_fit_" + TString::Format("pt_from_%d_to_%d_y_from_%.2f_to_%.2f",(int)pt_min,(int)pt_max,y_min,y_max) + mass_info;
  
  plot_mass_fit(ws_cut, channel, dir, (int) pt_max, (int) pt_min, y_min, y_max);

  signal_res = ws_cut.var("n_signal");
  
  return signal_res;
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
  double passed_pt = 0;
  double passed_y = 0;
  double passed_mu1 = 0;
  double passed_mu2 = 0;

  for (int evt=0;evt<tin->GetEntries();evt++)
    {
      tin->GetEntry(evt);
      
      if (fabs(eta_b) > 2.4) continue; //B mesons inside the detector region eta < 2.4

      if (fabs(y_b)<y_min || fabs(y_b)>y_max) continue; // within the y binning
      passed_y ++;

      if (pt_b<pt_min || pt_b>pt_max) continue; //within the pt bin
      passed_pt ++;

      hist_tot->Fill(pt_b);
      total ++;

      bool muon1Filter = fabs(eta_mu1)<2.4 && pt_mu1>2.8;
      bool muon2Filter = fabs(eta_mu2)<2.4 && pt_mu2>2.8;
      
      if (muon1Filter)
	{
	passed_mu1 ++;
	if(muon2Filter)
	  {
	    hist_passed->Fill(pt_b);//count only the events with the muon selection above
	    passed_mu2 ++;
	    passed ++;
	  }
	}
    }
  
  //debug
  /*
  std::cout << "debug: evt number: " << tin->GetEntries() << std::endl;
  std::cout << "debug: passed_y number: " << passed_y << std::endl;
  std::cout << "debug: passed_pt number: " << passed_pt << std::endl;
  std::cout << "debug: total number: " << total << std::endl;
  std::cout << "debug: passed_mu1 number: " << passed_mu1 << std::endl;
  std::cout << "debug: passed_mu2 number: " << passed_mu2 << std::endl;
  std::cout << "debug: passed number: " << passed << std::endl;
  std::cout << "debug: passed: " << hist_passed->GetBinContent(1) << std::endl;
  std::cout << "debug: total: " << hist_tot->GetBinContent(1) << std::endl;
  */
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
  TTree *tin_no_cuts = (TTree*)fin_no_cuts->Get(channel_to_ntuple_name(channel));
  
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
      
      if (fabs(eta_b) > 2.4) continue; //B mesons inside the detector region eta < 2.4
      if (fabs(y_b)<y_min || fabs(y_b)>y_max) continue; // within the y binning
      if (pt_b<pt_min || pt_b>pt_max) continue; //within the pt bin

      bool muon1Filter = fabs(eta_mu1) < 2.4 && pt_mu1 > 2.8;
      bool muon2Filter = fabs(eta_mu2) < 2.4 && pt_mu2 > 2.8;
 
      if (muon1Filter && muon2Filter) hist_tot->Fill(pt_b); //count only the events with the muon selection above
    }
      
    //--------------------------------read monte carlo with cuts------------------------
    TString mc_input_with_cuts = TString::Format(BASE_DIR) + "reduced_selected_myloop_new_mc_truth_" + channel_to_ntuple_name(channel) + "_with_cuts.root";
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
	
	if (fabs(eta_b) > 2.4) continue; //B mesons inside the detector region eta < 2.4
	if (fabs(y_b)<y_min || fabs(y_b)>y_max) continue; // within the y binning
	if (pt_b<pt_min || pt_b>pt_max) continue; //within the pt bin
		
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
  
  RooRealVar* kstar_to_k_pi = new RooRealVar("kstar","kstar",0.66503);
  kstar_to_k_pi->setError(1.4e-4);

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
/*
void read_vector(TString measure, int channel, TString vector, TString var1_name , TString var2_name, int n_var1_bins, int n_var2_bins,  double* var1_bins, double* var2_bins, double* array, double* err_lo, double* err_hi)
{
  //input_file_name
  TString bins_str = "";
  TString in_file_name = "";
  TString dir = TString::Format(VERSION) + "/";
  TString ntuple_name = "";

  if(vector == "yield")
    {
      dir += "signal_yield_root/" + channel_to_ntuple_name(channel) + "/";
      measure += "_";
    }
  else
    if(vector == "combined_syst" || vector == "signal_pdf" || vector == "cb_pdf" || vector == "mass_window")
    {
      dir += "signal_yield_root/syst/" + channel_to_ntuple_name(channel) + "/";
      measure += "_";
    }
    else
      if(vector == "preeff" || vector == "recoeff" || vector == "totaleff")
	{
	  dir += "efficiencies_root/" + channel_to_ntuple_name(channel) + "/";
	  measure += "_";
	}
      else
	if(vector == "x_sec")
	  {
	    dir += "x_sec/";
	    measure = "";
	  }
	else
	  if(vector == "BsBu" || vector == "fsfu" || vector == "BsBd" || vector == "fsfd" || vector == "BdBu" || vector == "fdfu")
	    {
	      dir += "ratio/";
	      measure = "";
	    }
  
  if(vector == "BsBu" || vector == "fsfu" || vector == "BsBd" || vector == "fsfd" || vector == "BdBu" || vector == "fdfu")
    ntuple_name = "";
  else
    ntuple_name = channel_to_ntuple_name(channel) + "_";
  
  for(int j=0; j<n_var2_bins; j++)
    {
      //file_name
      if(var1_name == "pt")
	bins_str = TString::Format("%.2f_to_%.2f", var2_bins[j], var2_bins[j+1]);
      else
	bins_str = TString::Format("%d_to_%d", (int)var2_bins[j], (int)var2_bins[j+1]);
      
      in_file_name = dir + vector + "_vector_" + measure + ntuple_name + var1_name + "_bins_" + var2_name + "_from_" + bins_str + ".root";
      
      if(n_var1_bins == 1 && n_var2_bins == 1)
	in_file_name = dir + vector + "_vector_" + measure + ntuple_name + "full_bins.root";

      //open file
      TFile* fin = new TFile(in_file_name);
      TVectorD *in_val = (TVectorD*)fin->Get("val");
      
      for(int i=0; i<n_var1_bins; i++)
	{
	  //copy to output_array
	  *(array + j*n_var1_bins + i) = in_val[0][i];
	}

      if(err_lo != NULL && err_hi != NULL)
	{
	  TVectorD *in_err_lo = (TVectorD*)fin->Get("err_lo");
	  TVectorD *in_err_hi = (TVectorD*)fin->Get("err_hi");
	  
	  for(int i=0; i<n_var1_bins; i++)
	    {
	      //copy to output_array
	      *(err_lo + j*n_var1_bins + i) = in_err_lo[0][i];
	      *(err_hi + j*n_var1_bins + i) = in_err_hi[0][i];
	    }
	}
      delete fin;
    }
}
*/

void read_vector(int channel, TString vector, TString var1_name , TString var2_name, int n_var1_bins, int n_var2_bins,  double* var1_bins, double* var2_bins, double* array, double* err_lo, double* err_hi)
{
  //input_file_name
  TString bins_str = "";
  TString in_file_name = "";
  TString dir = TString::Format(VERSION) + "/";

  if(vector == "yield")
      dir += "signal_yield_root/";
  else
    if(vector == "combined_syst" || vector == "signal_pdf" || vector == "cb_pdf" || vector == "mass_window")
      dir += "signal_yield_root/syst/";
    else
      if(vector == "preeff" || vector == "recoeff" || vector == "totaleff")
	dir += "efficiencies_root/";
  
  dir += channel_to_ntuple_name(channel) + "/";
  
  for(int j=0; j<n_var2_bins; j++)
    {
      for(int i=0; i<n_var1_bins; i++)
	{
	  //file_name
	  if(var1_name == "pt")
	    bins_str = TString::Format("pt_from_%d_to_%d_y_from_%.2f_to_%.2f", (int)var1_bins[i], (int)var1_bins[i+1], var2_bins[j], var2_bins[j+1]);
	  else
	    bins_str = TString::Format("pt_from_%d_to_%d_y_from_%.2f_to_%.2f", (int)var2_bins[j], (int)var2_bins[j+1], var1_bins[i], var1_bins[i+1]);
	  
	  in_file_name = dir + vector + "_" + channel_to_ntuple_name(channel) + "_" + bins_str + ".root";
	  
	  //debug
	  std::cout << "Reading: " << in_file_name << std::endl;
	  
	  //open file
	  TFile* fin = new TFile(in_file_name);
	  
	  if(fin->IsZombie())
	    {
	      std::cout << "The file " << in_file_name << " was not found." << std::endl;
	      std::cout << "No problem, the " << vector << " for the bin with " << bins_str << " will be calculated" << std::endl;
	      
	      ///////////////////////////////
	      //calculate yield or eff syst//
	      ///////////////////////////////

	      TString line = "";
	      TString command = "calculate_bin_";
	      TString opt = "";

	      TString pt_min;
	      TString pt_max;
	      TString y_min;
	      TString y_max;
	      
	      if(var1_name == "pt")
                {
		  pt_min = TString::Format("%d", (int)var1_bins[i]);
		  pt_max = TString::Format("%d", (int)var1_bins[i+1]);
		  y_min = TString::Format("%.2f", var2_bins[j]);
		  y_max = TString::Format("%.2f", var2_bins[j+1]);
                }
              else
                {
                  pt_min = TString::Format("%d", (int)var2_bins[j]);
                  pt_max = TString::Format("%d", (int)var2_bins[j+1]);
                  y_min = TString::Format("%.2f", var1_bins[i]);
                  y_max = TString::Format("%.2f", var1_bins[i+1]);
                }

              opt = " --channel " + TString::Format("%d", channel) + " --ptmin " + pt_min + " --ptmax " + pt_max + " --ymin " +  y_min + " --ymax " + y_max;

              if(vector == "yield")
                {
                  command += "yield";
                  line = command + opt;
                }
              else
                if(vector == "preeff" || vector == "recoeff" || vector == "totaleff")
                  {
                    command += "eff";
                    line = command + opt + " --eff " + vector;
                  }
                else
                  if(vector == "combined_syst" || vector == "signal_pdf" || vector == "cb_pdf" || vector == "mass_window")
                    {
                      command += "syst";
                      line = command + opt + " --syst " + vector;
                    }
	      
              gSystem->Exec(line);
	      //////////////////////////////////////////////////////////////////////////////////	      

	      fin = new TFile(in_file_name);
	    }
	    
	    TVectorD *in_val = (TVectorD*)fin->Get("val");
	    
	    //copy to output_array
	    *(array + j*n_var1_bins + i) = in_val[0][0];
	    
	    if(err_lo != NULL && err_hi != NULL)
	      {
		TVectorD *in_err_lo = (TVectorD*)fin->Get("err_lo");
		TVectorD *in_err_hi = (TVectorD*)fin->Get("err_hi");

		//copy to output_array
		*(err_lo + j*n_var1_bins + i) = in_err_lo[0][0];
		*(err_hi + j*n_var1_bins + i) = in_err_hi[0][0];
	      }
	    	    
	    delete fin;
	    delete in_val;
	}
    }
}
	    

  void plot_eff(TString measure, TString eff_name, int channel, int n_var1_bins, TString var2_name, double var2_min, double var2_max, TString x_axis_name, TString b_title, double* var1_bin_centre, double* var1_bin_centre_lo, double* var1_bin_centre_hi, double* eff_array, double* eff_err_lo_array, double* eff_err_hi_array)
  {  
  TCanvas ce;
  TGraphAsymmErrors* graph_pre_eff = new TGraphAsymmErrors(n_var1_bins, var1_bin_centre, eff_array, var1_bin_centre_lo, var1_bin_centre_hi, eff_err_lo_array, eff_err_hi_array);
  
  TString eff_title = "";
  TString ntuple = "";

  if(eff_name != "ratioeff")
    ntuple = channel_to_ntuple_name(channel) + "_";
  
  if(eff_name == "preeff")
    eff_title = b_title + " Pre-filter efficiency";
  if(eff_name == "recoeff")
    eff_title = b_title + " Reconstruction efficiency";
  if(eff_name == "totaleff")
    eff_title = b_title + " Overall efficiency";
  if(eff_name == "ratioeff")
    eff_title = b_title + " Efficiency ratio";

  graph_pre_eff->SetTitle(eff_title);
  graph_pre_eff->GetXaxis()->SetTitle(x_axis_name);
  graph_pre_eff->Draw("AP");

  if(measure != "x_sec" && measure != "ratio")
    return;
  
  TString save_eff = TString::Format(VERSION) + "/efficiencies/" + measure + "_" + eff_name + "_" + ntuple + var2_name + TString::Format("_from_%.2f_to_%.2f", var2_min, var2_max) + ".png";
  
  if(n_var1_bins == 1)
    save_eff = TString::Format(VERSION) + "/efficiencies/" + measure + "_" + eff_name + "_" + ntuple + "full_bins.png";
  
  ce.SaveAs(save_eff);
} 

void print_table(TString title, int n_var1_bins, int n_var2_bins, TString var1_name, TString var2_name, double* var1_bin_edges, double* var2_bin_edges, double* array, double* stat_err_lo, double* stat_err_hi, double* syst_err_lo, double* syst_err_hi)
{
  std::cout << title << " : " << std::endl;
  for(int j=0; j<n_var2_bins; j++)
    {
      std::cout << "BIN: " << var2_name << " " << var2_bin_edges[j] << " to " << var2_bin_edges[j+1] << " : " << std::endl;

      for(int i=0; i<n_var1_bins; i++)
        {
	  if(syst_err_lo == NULL || syst_err_hi == NULL)
	    {
	      std::cout << "BIN: " << var1_name << " " << var1_bin_edges[i] << " to " << var1_bin_edges[i+1] << " : " << *(array + j*n_var1_bins + i) << " +" << *(stat_err_hi + j*n_var1_bins + i) << " -" << *(stat_err_lo + j*n_var1_bins + i) << " (stat)" << std::endl;
	    }
	  else
	    {
	      std::cout << "BIN: " << var1_name << " " << var1_bin_edges[i] << " to " << var1_bin_edges[i+1] << " : " << *(array + j*n_var1_bins + i) << " +" << *(stat_err_hi + j*n_var1_bins + i) << " -" << *(stat_err_lo + j*n_var1_bins + i) << " (stat) +" << *(syst_err_hi + j*n_var1_bins + i) << " -" << *(syst_err_lo + j*n_var1_bins + i) << " (syst)" << std::endl;
	    }
        }
      std::cout << std::endl;
    }
}

void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, std::vector<std::vector<double> > numbers, std::string caption)
{
  std::ofstream file;

  file.open(filename + ".tex");
  
  // Create table                                                                                                                                
  file << "\\begin{table}" << std::endl;

  //setup table size                                                                                                                             
  std::string col="c";

  for(int i=1; i<n_col; i++)
    col+="|c";

  file << "\\begin{center}" << std::endl;
  file << "\\begin{tabular}{"+col+"}" << std::endl;
  file << "\\toprule" << std::endl;

  for(int c=0; c<n_col-1; c++)
    file << col_name[c] << " & ";

  file << col_name[n_col-1] << " \\\\" << std::endl;
  file << "\\midrule" << std::endl;

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
  file << "\\end{center}" << std::endl;
  file << "\\end{table}" << std::endl;

  std::string line;

  std::ifstream infile;
  infile.open(filename + ".tex");

  if(infile.is_open())
    {
      while(getline (infile,line))
	{
	  std::cout << line << std::endl;
	}
      infile.close();
    }

  infile.close();

  //system(("pdflatex " + filename + ".tex").c_str());
  //system(("gnome-open " + filename + ".pdf").c_str());
}

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
  
  TString dir_mc = TString::Format(VERSION) + "/";
  dir_mc += "mc_study/" + channel_to_ntuple_name(channel) + "/" + channel_to_ntuple_name(channel) + "_mc_study_" + TString::Format("pt_from_%d_to_%d_y_from_%.2f_to_%.2f",(int)pt_min,(int)pt_max,y_min,y_max);
  
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
