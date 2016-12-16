//this sript is used to estimate the ammout of "signal" in the B0->jpsi K*0 channel, in which the k and pi are swapped.
//for this we look at an MC sample processed using myloop_new .cc where we separate between true signal and swapped signal.
//here we fit each of these two categories separatly and extract the relative magnitude and width of each gaussian.
//This is then introduced in the pdf to extract the signal of B0->jpsi K*0, in other sripts in the cross sections studies.

#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

using namespace RooFit;

int main(int argc, char** argv)
{
  std::string output ="";
  int log = 0;

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--output")
        {
          convert << argv[++i];
          convert >> output;
        }
      if(argument == "--log")
        {
          convert << argv[++i];
          convert >> log;
        }
    }

  double mass_window = 0.5;
  double mass_min = B0_MASS - mass_window;
  double mass_max = B0_MASS + mass_window;
  double pt_min = 0;
  double pt_max = 300;
  double y_min = -3;
  double y_max = 3;

  RooRealVar mass("mass","mass",mass_min,mass_max);
  RooRealVar pt("pt","pt",pt_min,pt_max);
  RooRealVar y("y","y",y_min,y_max);
  
  TString input_file = TString::Format(BASE_DIR) + "selected_mc_ntkstar_with_cuts.root";

  TFile *fin = new TFile(input_file);
  TString directory = "";
  
  //set the ntuple name to extract the full signal and put it in RooDataSet data_full
  TString nt = "ntkstar";
  TTree *tr = (TTree*)fin->Get(nt);
  RooDataSet* data_full = new RooDataSet("data_full","data_full", tr, RooArgSet(mass,pt,y) );

  //set the ntuple name to extract the true signal and put it in RooDataset data_signal
  nt = "ntkstar_true";
  tr = (TTree*)fin->Get(nt);
  RooDataSet* data_signal = new RooDataSet("data_signal","data_signal", tr, RooArgSet(mass,pt,y) );

  //set the ntuple name to extract the swapped signal and put it in RooDataset data_swapped
  nt = "ntkstar_swap";
  tr = (TTree*)fin->Get(nt);
  RooDataSet* data_swapped = new RooDataSet("data_swapped","data_swapped", tr, RooArgSet(mass,pt,y) );
  
  //create histograms
  TH1D* histo_full = (TH1D*)data_full->createHistogram("histo_full", mass, Binning(channel_to_nbins(2), mass.getMin(), mass.getMax() ));
  TH1D* histo_signal = (TH1D*)data_signal->createHistogram("histo_signal", mass, Binning(channel_to_nbins(2), mass.getMin(), mass.getMax() ));
  TH1D* histo_swapped = (TH1D*)data_swapped->createHistogram("histo_swapped", mass, Binning(channel_to_nbins(2), mass.getMin(), mass.getMax() ));
  
  for (int i=1; i<=channel_to_nbins(2); i++)
    {
      if (histo_full->GetBinContent(i)==0) histo_full->SetBinError(i,0.);
      if (histo_signal->GetBinContent(i)==0) histo_signal->SetBinError(i,0.);
      if (histo_swapped->GetBinContent(i)==0) histo_swapped->SetBinError(i,0.);
    }
  
  //========make fits of each category===========

  RooRealVar mean("mean","mean", B0_MASS);
  mean.setConstant(kTRUE);
  
  RooRealVar m_exp("m_exp","m_exp",0 ,-4., 0.);
  m_exp.setConstant(kTRUE);
  RooExponential combinatorial_exp("combinatorial_exp","combinatorial_exp", mass, m_exp);

  RooRealVar r1("r1","r1",0.1,0.0,1);

  //true signal fit
  RooRealVar sigma_signal1("sigma_signal1","sigma_signal1",0.015, 0.001, 0.5);
  RooRealVar sigma_signal2("sigma_signal2","sigma_signal2",0.015, 0.001, 0.5);
  RooGaussian gauss_signal1("gauss_signal1","gauss_signal1", mass, mean, sigma_signal1);
  RooGaussian gauss_signal2("gauss_signal2","gauss_signal2", mass, mean, sigma_signal2);
  RooRealVar r_signal("r_signal","r_signal", 0.5, 0, 1);
  RooAddPdf signal("signal","signal",RooArgSet(gauss_signal1,gauss_signal2), r_signal);

  //RooAddPdf model_signal("model_signal","model_signal",RooArgSet(combinatorial_exp, signal), r1);

  RooPlot* frame1 = mass.frame(Title("True signal fit"));
  data_signal->plotOn(frame1,Name("theSignal"),Binning(channel_to_nbins(2)));
    
  signal.fitTo(*data_signal);
  signal.paramOn(frame1);
  signal.plotOn(frame1, LineColor(7), LineWidth(4), LineStyle(1));
  //model_signal.plotOn(frame1,Components("combinatorial_exp"),LineColor(9),LineWidth(4),LineStyle(2));
  signal.plotOn(frame1,Components("gauss_signal1"),LineColor(8),LineWidth(4),LineStyle(2));
  signal.plotOn(frame1,Components("gauss_signal2"),LineColor(5),LineWidth(4),LineStyle(2));
  
  TCanvas c1;
  if(log)
    c1.SetLogy();
  c1.cd();
  frame1->Draw();
  directory = output + "mass_signal.png";
  c1.SaveAs(directory);

  //swapped signal fit
  RooRealVar sigma_swapped1("sigma_swapped1","sigma_swapped1",0.015, 0.001, 0.5);
  RooRealVar sigma_swapped2("sigma_swapped2","sigma_swapped2",0.015, 0.001, 0.5);
  RooGaussian gauss_swapped1("gauss_swapped1","gauss_swapped1", mass, mean, sigma_swapped1);
  RooGaussian gauss_swapped2("gauss_swapped2","gauss_swapped2", mass, mean, sigma_swapped2);
  RooRealVar r_swapped("r_swapped","r_swapped", 0.5, 0, 1);
  RooAddPdf swapped("swapped","swapped",RooArgSet(gauss_swapped1,gauss_swapped2), r_swapped);

  //RooAddPdf model_swapped("model_swapped","model_swapped",RooArgSet(combinatorial_exp, swapped), r1);

  RooPlot* frame2 = mass.frame(Title("Swapped signal fit"));
  data_swapped->plotOn(frame2,Name("theSwapped"),Binning(channel_to_nbins(2)));
 
  swapped.fitTo(*data_swapped);
  swapped.paramOn(frame2);
  swapped.plotOn(frame2, LineColor(7), LineWidth(4), LineStyle(1));
  //model_swapped.plotOn(frame2,Components("combinatorial_exp"),LineColor(9),LineWidth(4),LineStyle(2));
  swapped.plotOn(frame2,Components("gauss_swapped1"),LineColor(8),LineWidth(4),LineStyle(2));
  swapped.plotOn(frame2,Components("gauss_swapped2"),LineColor(5),LineWidth(4),LineStyle(2));
 
  TCanvas c2;
  if(log)
    c2.SetLogy();
  c2.cd();
  frame2->Draw();
  directory = output + "mass_swapped.png";
  c2.SaveAs(directory);

  //set signal and swapped fits constant
  sigma_signal1.setConstant(kTRUE);
  sigma_signal2.setConstant(kTRUE);
  r_signal.setConstant(kTRUE);

  sigma_swapped1.setConstant(kTRUE);
  sigma_swapped2.setConstant(kTRUE);
  r_swapped.setConstant(kTRUE);

  //full fit
  RooRealVar r2("r2","r2",0.1,0.0,1);
  
  RooAddPdf model_full("model_full","model_full",RooArgSet(swapped, signal),r1);

  RooPlot* frame3 = mass.frame(Title("Full signal fit"));
  data_full->plotOn(frame3,Name("theFull"),Binning(channel_to_nbins(2)));
  
  model_full.fitTo(*data_full);
  model_full.paramOn(frame3);
  model_full.plotOn(frame3, LineColor(7), LineWidth(4), LineStyle(1));
  //  model_full.plotOn(frame3,Components("combinatorial_exp"),LineColor(9),LineWidth(4),LineStyle(2));
  model_full.plotOn(frame3,Components("swapped"),LineColor(8),LineWidth(4),LineStyle(2));
  model_full.plotOn(frame3,Components("signal"),LineColor(5),LineWidth(4),LineStyle(2));

  TCanvas c3;
  if(log)
    c3.SetLogy();
  c3.cd();
  frame3->Draw();
  directory = output + "mass_full.png";
  c3.SaveAs(directory);

  double number_signal = (double)data_signal->sumEntries();
  double number_swap   = (double)data_swapped->sumEntries();

  std::cout << "Numbers from counting the entries:" << std::endl;

  std::cout << "number signal: " << number_signal << std::endl;
  std::cout << "number swap: " << number_swap << std::endl;
  std::cout << "number ratio: " << number_swap/number_signal << std::endl;

  std::cout << " " << std::endl;
  std::cout << "Numbers from the fits:" << std::endl;

  //std::cout << "ratios:" << " r1" << r1.getVal() << " r2" << r2.getVal() << std::endl;
  std::cout << "ratio: " << r1.getVal() << std::endl;
  
  delete data_signal;
  delete data_swapped;
  delete data_full;
  delete fin;
}
