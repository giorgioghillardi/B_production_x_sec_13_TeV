#include <sstream>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TSystem.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include "UserCode/B_production_x_sec_13_TeV/interface/myloop.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/plotDressing.h"
#include "TMath.h"
using namespace RooFit;

#define MASS_MIN_2 5.0
#define MASS_MAX_2 5.6


//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda

void plot_pt_dist(RooWorkspace& w, int channel, TString directory);
void plot_mass_dist(RooWorkspace& w, int channel, TString directory);

void read_data(RooWorkspace& w, TString filename,int channel);
void read_data_cut(RooWorkspace& w, RooDataSet* data);
void set_up_workspace_variables(RooWorkspace& w, int channel);
TString channel_to_ntuple_name(int channel);

void sideband_sub(RooWorkspace& w, double left, double right);

void data_selection(TString fin1,TString data_selection_output_file,int channel);

//input example: data_selection --channel 1 --input /some/place/ --sub 1 --showdist 1
int main(int argc, char** argv)
{
  int channel = 0;
  std::string input_file = "/lstore/cms/brunogal/input_for_B_production_x_sec_13_TeV/myloop_data.root";
  bool side_sub = 0, show_dist = 0;

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--channel")
        {
          convert << argv[++i];
          convert >> channel;
        }

      if(argument == "--input")
        {
          convert << argv[++i];
          convert >> input_file;
        }

      if(argument == "--sub")
	{
	  convert << argv[++i];
	  convert >> side_sub;
	}

      if(argument == "--showdist")
	{
	  convert << argv[++i];
	  convert >> show_dist;
	}
    }

  if(channel==0)
    {
      std::cout << "No channel was provided as input. Please use --channel. Example: data_selection --channel 1" << std::endl;
      return 0;
    }

  TString data_selection_output_file="";
  data_selection_output_file= "selected_data_" + channel_to_ntuple_name(channel) + ".root";
  
  data_selection(input_file,data_selection_output_file,channel);

  if(show_dist)
    { 
      RooWorkspace* ws = new RooWorkspace("ws","Bmass");
      TString pt_dist_directory="";
      TString mass_dist_directory="";
 
      gSystem->Exec("mkdir -p full_dataset_mass_pt_dist/");

      //set up mass and pt variables inside ws  
      set_up_workspace_variables(*ws,channel);
      
      //read data from the selected data file, and import it as a dataset into the workspace.
      read_data(*ws, data_selection_output_file,channel);
      
      RooAbsData* data = ws->data("data");
      
      pt_dist_directory = "full_dataset_mass_pt_dist/" + channel_to_ntuple_name(channel) + "_pt";
      plot_pt_dist(*ws,channel,pt_dist_directory);
      
      mass_dist_directory = "full_dataset_mass_pt_dist/" + channel_to_ntuple_name(channel) + "_mass";
      plot_mass_dist(*ws,channel,mass_dist_directory);
    }

  if(side_sub)
    {
      RooWorkspace* ws = new RooWorkspace("ws","Bmass");
      //set up mass and pt variables inside ws  
      set_up_workspace_variables(*ws,channel);

      //read data from the selected data file, and import it as a dataset into the workspace.
      read_data(*ws, data_selection_output_file,channel);

      switch(channel)
	{
	case 2:
	  sideband_sub(*ws, 5.1, 5.4);
	  break;
	case 4:
	  sideband_sub(*ws, 5.25, 5.45);
	  break;
	default:
	  std::cout << "WARNING! UNDEFINED LIMITS FOR PEAK REGION" << std::endl;
	}
    }
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
  c2.SaveAs(directory + ".root");
  c2.SaveAs(directory + ".png");
}

void plot_mass_dist(RooWorkspace& w, int channel, TString directory)
{
  //full dataset mass distribution
  RooRealVar mass = *(w.var("mass"));
  RooAbsData* data = w.data("data");

  TCanvas c2;
  TH1D* mass_dist = (TH1D*)data->createHistogram("mass_dist",mass);
  mass_dist->Draw();
  c2.SaveAs(directory + ".root");
  c2.SaveAs(directory + ".png");
}

void read_data(RooWorkspace& w, TString filename,int channel)
{
  std::cout<<std::endl<<"READ_DATA!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;

  TFile* f = new TFile(filename);
  TNtupleD* _nt = (TNtupleD*)f->Get(channel_to_ntuple_name(channel));
 
  RooArgSet arg_list(*(w.var("mass")) , *(w.var("pt")) , *(w.var("y")) , *(w.var("mu1pt")) , *(w.var("mu2pt")) , *(w.var("mu1eta")) , *(w.var("mu2eta")) , *(w.var("lxy")) , *(w.var("errxy")) );

  arg_list.add(*(w.var("vtxprob")));
  arg_list.add(*(w.var("lerrxy")));

  RooDataSet* data = new RooDataSet("data","data",_nt,arg_list);
  
  w.import(*data);
}

void read_data_cut(RooWorkspace& w, RooDataSet* data)
{
  w.import(*data);
}

void set_up_workspace_variables(RooWorkspace& w, int channel)
{
  std::cout<<std::endl<<"SET_UP_WORKSPACE_VARIABLES!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
  double mass_min, mass_max;
  double pt_min, pt_max;
  double y_min, y_max;
  double mu1pt_min, mu1pt_max;
  double mu2pt_min, mu2pt_max;
  double mu1eta_min, mu1eta_max;
  double mu2eta_min, mu2eta_max;
  double lxy_min, lxy_max;
  double errxy_min, errxy_max;
  double vtxprob_min, vtxprob_max;
  double lerrxy_min, lerrxy_max;

  pt_min=0;
  pt_max=400;
  
  y_min=-3;
  y_max=3;

  mu1pt_min=0;
  mu1pt_max=80;

  mu2pt_min=0;
  mu2pt_max=90;

  mu1eta_min=-3;
  mu1eta_max=3;

  mu2eta_min=-3;
  mu2eta_max=3;

  lxy_min=0;
  lxy_max=3.5;

  errxy_min=0;
  errxy_max=0.05;

  vtxprob_min=0;
  vtxprob_max=1;

  lerrxy_min=0;
  lerrxy_max=42;

  switch (channel) {
  default:
  case 1:
    mass_min = 5.0; mass_max = 6.0;
    break;
  case 2:
    mass_min = MASS_MIN_2; mass_max = MASS_MAX_2;
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
  RooRealVar y("y","y",y_min,y_max);
  RooRealVar mu1pt("mu1pt","mu1pt",mu1pt_min,mu1pt_max);
  RooRealVar mu2pt("mu2pt","mu2pt",mu2pt_min,mu2pt_max);
  RooRealVar mu1eta("mu1eta","mu1eta",mu1eta_min,mu1eta_max);
  RooRealVar mu2eta("mu2eta","mu2eta",mu2eta_min,mu2eta_max);
  RooRealVar lxy("lxy","lxy",lxy_min,lxy_max);
  RooRealVar errxy("errxy","errxy",errxy_min,errxy_max);
  RooRealVar vtxprob("vtxprob","vtxprob",vtxprob_min,vtxprob_max);
  RooRealVar lerrxy("lerrxy","lerrxy",lerrxy_min,lerrxy_max);

  w.import(mass);
  w.import(pt);
  w.import(y);
  w.import(mu1pt);
  w.import(mu2pt);
  w.import(mu1eta);
  w.import(mu2eta);
  w.import(lxy);
  w.import(errxy);
  w.import(vtxprob);
  w.import(lerrxy);
}

void data_selection(TString fin1, TString data_selection_output_file,int channel)
{
  std::cout<<std::endl<<"DATA_SELECTION!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
  TFile *fout = new TFile(data_selection_output_file,"recreate");

  TNtupleD *_nt1 = new TNtupleD("ntkp","ntkp","mass:pt:eta:y:mu1pt:mu2pt:mu1eta:mu2eta:lxy:errxy:vtxprob:lerrxy");
  TNtupleD *_nt2 = new TNtupleD("ntkstar","ntkstar","mass:pt:eta:y:mu1pt:mu2pt:mu1eta:mu2eta:lxy:errxy:vtxprob:lerrxy");
  TNtupleD *_nt3 = new TNtupleD("ntks","ntks","mass:pt:eta:y:mu1pt:mu2pt:mu1eta:mu2eta:lxy:errxy:vtxprob:lerrxy");
  TNtupleD *_nt4 = new TNtupleD("ntphi","ntphi","mass:pt:eta:y:mu1pt:mu2pt:mu1eta:mu2eta:lxy:errxy:vtxprob:lerrxy");
  TNtupleD *_nt5 = new TNtupleD("ntmix","ntmix","mass:pt:eta:y:mu1pt:mu2pt:mu1eta:mu2eta:lxy:errxy:vtxprob:lerrxy");
  TNtupleD *_nt6 = new TNtupleD("ntlambda","ntlambda","mass:pt:eta:y:mu1pt:mu2pt:mu1eta:mu2eta:lxy:errxy:vtxprob:lerrxy");

  /*
    switch (channel) {
    case 1:
    default:
    _nt1 = new TNtupleD("ntkp","ntkp","mass:pt:eta");
    break;
    case 2:
    _nt2 = new TNtupleD("ntkstar","ntkstar","mass:pt:eta");
    break;
    case 3:
    _nt3 = new TNtupleD("ntks","ntks","mass:pt:eta");
    break;
    case 4:
    _nt4 = new TNtupleD("ntphi","ntphi","mass:pt:eta");
    break;
    case 5:
    _nt5 = new TNtupleD("ntmix","ntmix","mass:pt:eta");
    break;
    case 6:
    _nt6 = new TNtupleD("ntlambda","ntlambda","mass:pt:eta");
    break;
    }
  */

  ReducedBranches br;
  TChain* tin;
  
  int n_br_queued = 0;
  ReducedBranches br_queue[32];
  TLorentzVector v4_tk1, v4_tk2;
  
  std::cout << "selecting data from channel " << channel << std::endl;

  tin = new TChain(channel_to_ntuple_name(channel));

  tin->Add(fin1);
  br.setbranchadd(tin);

  for (int evt=0;evt<tin->GetEntries();evt++) {
    tin->GetEntry(evt);
        
    if (channel==1) { // cuts for B+ -> J/psi K+
      if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
      if (br.vtxprob<=0.1) continue;//original cut 0.1
      if (br.tk1pt<=1.6) continue;//original cut 1.6
      if (br.lxy/br.errxy<=3.0) continue;//original cut 3.0
      if (br.cosalpha2d<=0.99) continue;//original cut 0.99
            
      _nt1->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu2pt,br.mu1eta,br.mu2eta,br.lxy,br.errxy,br.vtxprob, br.lxy/br.errxy);
	    
    }else
      if (channel==2) { // cuts for B0 -> J/psi K*
	if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
	if (br.vtxprob<=0.2) continue;//original cut 0.1
	if (br.lxy/br.errxy<=4.5) continue;//original cut 3.0
	if (br.cosalpha2d<=0.996) continue;//original cut 0.99
	if (fabs(br.tktkmass-KSTAR_MASS)>=0.05) continue;//original cut 0.05
            
	v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
	v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
	if (fabs((v4_tk1+v4_tk2).Mag()-PHI_MASS)<=0.01) continue;//original cut 0.01
            
	if (n_br_queued==0) {
	  memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
	  n_br_queued++;
	}else
	  if (br.run == br_queue[n_br_queued-1].run && br.event == br_queue[n_br_queued-1].event) { // same event
	    memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
	    n_br_queued++;
	    if (n_br_queued>=32) printf("Warning: maximum queued branches reached.\n");
	  }            
	if (br.run != br_queue[n_br_queued-1].run || br.event != br_queue[n_br_queued-1].event || evt==tin->GetEntries()-1) {
	  for (int i=0; i<n_br_queued; i++) {
                    
	    bool isBestKstarMass = true;
	    for (int j=0; j<n_br_queued; j++) {
	      if (j==i) continue;
	      if (br_queue[i].mu1idx==br_queue[j].mu1idx &&
		  br_queue[i].mu2idx==br_queue[j].mu2idx &&
		  br_queue[i].tk1idx==br_queue[j].tk1idx &&
		  br_queue[i].tk2idx==br_queue[j].tk2idx) {
                        
		if (fabs(br_queue[j].tktkmass-KSTAR_MASS)<fabs(br_queue[i].tktkmass-KSTAR_MASS)) {
		  isBestKstarMass = false;
		  continue;
		}
	      }
	    }
                                 
	    if (isBestKstarMass){
	      _nt2->Fill(br_queue[i].mass,br_queue[i].pt,br_queue[i].eta,br_queue[i].y,
			 br_queue[i].mu1pt,br_queue[i].mu2pt,br_queue[i].mu1eta,br_queue[i].mu2eta,
			 br_queue[i].lxy,br_queue[i].errxy,br_queue[i].vtxprob, br_queue[i].lxy/br_queue[i].errxy);

	    }

	  }                
	  n_br_queued = 0;
	  memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
	  n_br_queued++;
	}
      }else
	if (channel==3) { // cuts for B0 -> J/psi Ks
	  if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
	  if (br.vtxprob<=0.1) continue;//original cut 0.1
	  if (br.lxy/br.errxy<=3.0) continue;//original cut 3.0
	  if (br.tktkblxy/br.tktkberrxy<=3.0) continue;//original cut 3.0
	  if (br.cosalpha2d<=0.99) continue;//original cut 0.99
	  if (fabs(br.tktkmass-KSHORT_MASS)>=0.015) continue;//original cut 0.015
                
	  _nt3->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu2pt,br.mu1eta,br.mu2eta,br.lxy,br.errxy,br.vtxprob,br.lxy/br.errxy);

	}else
	  if (channel==4) { // cuts for Bs -> J/psi phi
	    if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
	    if (br.vtxprob<=0.2) continue;//original cut 0.1
	    if (br.lxy/br.errxy<=4.5) continue;//original cut 3.0
	    if (br.cosalpha2d<=0.996) continue;//original cut 0.99
	    if (fabs(br.tktkmass-PHI_MASS)>=0.010) continue;//original cut 0.010
            
	    v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
	    v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
	    if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;//original cut 0.05
	    v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
	    v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
	    if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
                
	    _nt4->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu2pt,br.mu1eta,br.mu2eta,br.lxy,br.errxy,br.vtxprob, br.lxy/br.errxy);

	  }else
	    if (channel==5) { // cuts for psi(2S)/X(3872) -> J/psi pipi
	      if (br.vtxprob<=0.2) continue;//original cut 0.2
	      if (fabs(br.tk1eta)>=1.6) continue;//original cut 1.6
	      if (fabs(br.tk2eta)>=1.6) continue;//original cut 1.6
            
	      _nt5->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu2pt,br.mu1eta,br.mu2eta,br.lxy,br.errxy,br.vtxprob, br.lxy/br.errxy);

	    }else
	      if (channel==6) {//cuts for lambda
		if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
		if (br.vtxprob<=0.1) continue;//original cut 0.1
		if (br.lxy/br.errxy<=3.0) continue;//original cut 3.0
		if (br.tktkblxy/br.tktkberrxy<=3.0) continue;//original cut 3.0
		if (br.cosalpha2d<=0.99) continue;//original cut 0.99
		if (fabs(br.tktkmass-LAMBDA_MASS)>=0.015) continue;//original cut 0.015
            
		v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
		v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
		if (fabs((v4_tk1+v4_tk2).Mag()-KSHORT_MASS)<=0.015) continue;//original cut 0.015
            
		_nt6->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu2pt,br.mu1eta,br.mu2eta,br.lxy,br.errxy,br.vtxprob, br.lxy/br.errxy);

	      }
  }//end of the for for the events
  fout->Write();
  fout->Close();
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

void sideband_sub(RooWorkspace& w, double left, double right)
{
  std::cout<<std::endl<<"SIDEBAND_SUB!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
  //Create appropriate variables and data sets (the pt isn't imported from the RooWorkspace because its range will change)                      
  RooRealVar pt = *(w.var("pt"));                                                                                                            
  RooRealVar y = *(w.var("y"));                                                                                                            
  RooRealVar mu1pt = *(w.var("mu1pt"));
  RooRealVar mu2pt = *(w.var("mu2pt"));
  RooRealVar mu1eta = *(w.var("mu1eta"));
  RooRealVar mu2eta = *(w.var("mu2eta"));
  RooRealVar lxy = *(w.var("lxy"));
  RooRealVar errlxy = *(w.var("errxy"));
  RooRealVar vtxprob = *(w.var("vtxprob"));
  RooRealVar lerrxy = *(w.var("lerrxy"));  
  //RooRealVar pt("pt","pt",0.,150.);
  RooRealVar mass = *(w.var("mass"));
  RooDataSet* data =(RooDataSet*) w.data("data");
  RooDataSet* reduceddata_side;  RooDataSet* reduceddata_aux;
  RooDataSet* reduceddata_central;

  //Make selection for the different bands using mass as the selection variable                                                                  

  reduceddata_side = (RooDataSet*) data->reduce(Form("mass<%lf", left));
  reduceddata_aux = (RooDataSet*) data->reduce(Form("mass>%lf",right));

  reduceddata_side->append(*reduceddata_aux);

  reduceddata_central = (RooDataSet*) data->reduce(Form("mass>%lf",left));
  reduceddata_central = (RooDataSet*) reduceddata_central->reduce(Form("mass<%lf",right));
  
  /*RooRealVar mean("mean", "mean", 0., 5.);                                                                                                    
    RooRealVar sigma("sigma", "sigma", -1000., 1000.);*/

  RooRealVar lambda("lambda", "lambda",-5., -20., 0.);

  RooExponential fit_side("fit_side", "fit_side_exp", mass, lambda);

  mass.setRange("all", mass.getMin(),mass.getMax());
  mass.setRange("right",right,mass.getMax());
  mass.setRange("left",mass.getMin(),left);
  mass.setRange("peak",left,right);
  
  std::cout<<"mass minimum: "<<mass.getMin()<<std::endl;
  std::cout<<"mass maximum: "<<mass.getMax()<<std::endl;
  

  fit_side.fitTo(*reduceddata_side,Range("left,right"));
//  RooRealVar* nll = (RooRealVar*) fit_side.createNLL(*reduceddata_side, Range("left,right"));

  RooPlot* massframe = mass.frame();
  reduceddata_side->plotOn(massframe);
  fit_side.plotOn(massframe, Range("all"));

  TCanvas d;
  massframe->Draw();
  d.SaveAs("fit_side.png");

  std::cout << std::endl << "chisquare: " << massframe->chiSquare() << std::endl;
//  std::cout << "LogLikelihood: " << nll->getVal() << std::endl;

  //Integrating the background distribution                            

  RooAbsReal* int_fit_side_left = fit_side.createIntegral(mass, "left");
  RooAbsReal* int_fit_side_right = fit_side.createIntegral(mass, "right");
  RooAbsReal* int_fit_peak = fit_side.createIntegral(mass, "peak");

  std::cout<< std::endl << "Integral left band: " << int_fit_side_left->getVal() << std::endl;
  std::cout<< std::endl << "Integral right band: " << int_fit_side_right->getVal() << std::endl;

  double factor = (int_fit_peak->getVal())/(int_fit_side_left->getVal()+int_fit_side_right->getVal());

  std::cout << std::endl << "Factor: " << factor << std::endl;

  //Build and draw signal and background distributionsx                

  TH1D* pt_dist_side = (TH1D*) reduceddata_side->createHistogram("pt_dist_side",pt);
  pt_dist_side->SetMarkerColor(kBlue);
  pt_dist_side->SetLineColor(kBlue);
  pt_dist_side->SetNameTitle("pt_dist_side", "Signal and Background Distributions - p_{T} (B) ");

  TH1D* pt_dist_peak = (TH1D*) reduceddata_central->createHistogram("pt_dist_peak", pt);
  pt_dist_peak->SetMarkerColor(kRed);
  pt_dist_peak->SetLineColor(kRed);
  pt_dist_peak->SetNameTitle("pt_dist_peak", "Signal and Background Distributions - p_{T} (B)");

  TH1D* pt_dist_total = (TH1D*) data->createHistogram("pt_dist_total",pt);
  pt_dist_total->SetMarkerColor(kBlack);
  pt_dist_total->SetLineColor(kBlack);

  pt_dist_peak->Add(pt_dist_side, -factor);
  pt_dist_side->Add(pt_dist_side, factor);

  TCanvas c;

  pt_dist_total->Draw();
  pt_dist_side->Draw("same");
  pt_dist_peak->Draw("same");
  pt_dist_peak->SetXTitle("p_{T} [GeV]");
  pt_dist_side->SetXTitle("p_{T} [GeV]");
  pt_dist_total->SetXTitle("p_{T} [GeV]");

  TLegend *leg = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg->AddEntry("pt_dist_total", "Total", "l");
  leg->AddEntry("pt_dist_peak", "Signal", "l");
  leg->AddEntry("pt_dist_side", "Background", "l");
  leg->Draw("same");

  c.SetLogy();
  c.SaveAs("pt_sideband_sub.png");


  TH1D* mu1pt_dist_side = (TH1D*) reduceddata_side->createHistogram("mu1pt_dist_side",mu1pt);
  mu1pt_dist_side->SetMarkerColor(kBlue);
  mu1pt_dist_side->SetLineColor(kBlue);
  mu1pt_dist_side->SetNameTitle("mu1pt_dist_side", "Signal and Background Distributions - p_{T} (#mu_{1}) ");

  TH1D* mu1pt_dist_peak = (TH1D*) reduceddata_central->createHistogram("mu1pt_dist_peak", mu1pt);
  mu1pt_dist_peak->SetMarkerColor(kRed);
  mu1pt_dist_peak->SetLineColor(kRed);
  mu1pt_dist_peak->SetNameTitle("mu1pt_dist_peak", "Signal and Background Distributions - p_{T} (#mu_{1})");

  TH1D* mu1pt_dist_total = (TH1D*) data->createHistogram("mu1pt_dist_total",mu1pt);
  mu1pt_dist_total->SetMarkerColor(kBlack);
  mu1pt_dist_total->SetLineColor(kBlack);

  mu1pt_dist_peak->Add(mu1pt_dist_side, -factor);
  mu1pt_dist_side->Add(mu1pt_dist_side, factor);

  TCanvas c1;

  mu1pt_dist_total->Draw();
  mu1pt_dist_side->Draw("same");
  mu1pt_dist_peak->Draw("same");
  mu1pt_dist_peak->SetXTitle("p_{T} [GeV]");
  mu1pt_dist_side->SetXTitle("p_{T} [GeV]");
  mu1pt_dist_total->SetXTitle("p_{T} [GeV]");

  TLegend *leg1 = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg1->AddEntry("mu1pt_dist_total", "Total", "l");
  leg1->AddEntry("mu1pt_dist_peak", "Signal", "l");
  leg1->AddEntry("mu1pt_dist_side", "Background", "l");
  leg1->Draw("same");

  c1.SetLogy();
  c1.SaveAs("mu1pt_sideband_sub.png");

  TH1D* mu2pt_dist_side = (TH1D*) reduceddata_side->createHistogram("mu2pt_dist_side",mu2pt);
  mu2pt_dist_side->SetMarkerColor(kBlue);
  mu2pt_dist_side->SetLineColor(kBlue);
  mu2pt_dist_side->SetNameTitle("mu2pt_dist_side", "Signal and Background Distributions - p_{T} (#mu_{2}) ");

  TH1D* mu2pt_dist_peak = (TH1D*) reduceddata_central->createHistogram("mu2pt_dist_peak", mu2pt);
  mu2pt_dist_peak->SetMarkerColor(kRed);
  mu2pt_dist_peak->SetLineColor(kRed);
  mu2pt_dist_peak->SetNameTitle("mu2pt_dist_peak", "Signal and Background Distributions - p_{T} (#mu_{2})");

  TH1D* mu2pt_dist_total = (TH1D*) data->createHistogram("mu2pt_dist_total",mu2pt);
  mu2pt_dist_total->SetMarkerColor(kBlack);
  mu2pt_dist_total->SetLineColor(kBlack);

  mu2pt_dist_peak->Add(mu2pt_dist_side, -factor);
  mu2pt_dist_side->Add(mu2pt_dist_side, factor);

  TCanvas c2;

  mu2pt_dist_total->Draw();
  mu2pt_dist_side->Draw("same");
  mu2pt_dist_peak->Draw("same");
  mu2pt_dist_peak->SetXTitle("p_{T} [GeV]");
  mu2pt_dist_side->SetXTitle("p_{T} [GeV]");
  mu2pt_dist_total->SetXTitle("p_{T} [GeV]");

  TLegend *leg2 = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg2->AddEntry("mu2pt_dist_total", "Total", "l");
  leg2->AddEntry("mu2pt_dist_peak", "Signal", "l");
  leg2->AddEntry("mu2pt_dist_side", "Background", "l");
  leg2->Draw("same");

  c2.SetLogy();
  c2.SaveAs("mu2pt_sideband_sub.png");

  TH1D* mu1eta_dist_side = (TH1D*) reduceddata_side->createHistogram("mu1eta_dist_side",mu1eta);
  mu1eta_dist_side->SetMarkerColor(kBlue);
  mu1eta_dist_side->SetLineColor(kBlue);
  mu1eta_dist_side->SetNameTitle("mu1eta_dist_side", "Signal and Background Distributions - #eta (#mu_{1}) ");

  TH1D* mu1eta_dist_peak = (TH1D*) reduceddata_central->createHistogram("mu1eta_dist_peak", mu1eta);
  mu1eta_dist_peak->SetMarkerColor(kRed);
  mu1eta_dist_peak->SetLineColor(kRed);
  mu1eta_dist_peak->SetNameTitle("mu1eta_dist_peak", "Signal and Background Distributions - #eta (#mu_{1})");

  TH1D* mu1eta_dist_total = (TH1D*) data->createHistogram("mu1eta_dist_total",mu1eta);
  mu1eta_dist_total->SetMarkerColor(kBlack);
  mu1eta_dist_total->SetLineColor(kBlack);

  mu1eta_dist_peak->Add(mu1eta_dist_side, -factor);
  mu1eta_dist_side->Add(mu1eta_dist_side, factor);

  TCanvas c3;

  mu1eta_dist_total->Draw();
  mu1eta_dist_side->Draw("same");
  mu1eta_dist_peak->Draw("same");
  mu1eta_dist_peak->SetXTitle("#eta");
  mu1eta_dist_side->SetXTitle("#eta");
  mu1eta_dist_total->SetXTitle("#eta");

  TLegend *leg3 = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg3->AddEntry("mu1eta_dist_total", "Total", "l");
  leg3->AddEntry("mu1eta_dist_peak", "Signal", "l");
  leg3->AddEntry("mu1eta_dist_side", "Background", "l");
  leg3->Draw("same");

  c3.SetLogy();
  c3.SaveAs("mu1eta_sideband_sub.png");

  TH1D* mu2eta_dist_side = (TH1D*) reduceddata_side->createHistogram("mu2eta_dist_side",mu2eta);
  mu2eta_dist_side->SetMarkerColor(kBlue);
  mu2eta_dist_side->SetLineColor(kBlue);
  mu2eta_dist_side->SetNameTitle("mu2eta_dist_side", "Signal and Background Distributions - #eta (#mu_{2}) ");

  TH1D* mu2eta_dist_peak = (TH1D*) reduceddata_central->createHistogram("mu2eta_dist_peak", mu2eta);
  mu2eta_dist_peak->SetMarkerColor(kRed);
  mu2eta_dist_peak->SetLineColor(kRed);
  mu2eta_dist_peak->SetNameTitle("mu2eta_dist_peak", "Signal and Background Distributions - #eta (#mu_{2})");

  TH1D* mu2eta_dist_total = (TH1D*) data->createHistogram("mu1eta_dist_total",mu2eta);
  mu2eta_dist_total->SetMarkerColor(kBlack);
  mu2eta_dist_total->SetLineColor(kBlack);

  mu2eta_dist_peak->Add(mu2eta_dist_side, -factor);
  mu2eta_dist_side->Add(mu2eta_dist_side, factor);

  TCanvas c4;

  mu2eta_dist_total->Draw();
  mu2eta_dist_side->Draw("same");
  mu2eta_dist_peak->Draw("same");
  mu2eta_dist_peak->SetXTitle("#eta");
  mu2eta_dist_side->SetXTitle("#eta");
  mu2eta_dist_total->SetXTitle("#eta");

  TLegend *leg4 = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg4->AddEntry("mu2eta_dist_total", "Total", "l");
  leg4->AddEntry("mu2eta_dist_peak", "Signal", "l");
  leg4->AddEntry("mu2eta_dist_side", "Background", "l");
  leg4->Draw("same");

  c4.SetLogy();
  c4.SaveAs("mu2eta_sideband_sub.png");

  TH1D* y_dist_side = (TH1D*) reduceddata_side->createHistogram("y_dist_side",y);
  y_dist_side->SetMarkerColor(kBlue);
  y_dist_side->SetLineColor(kBlue);
  y_dist_side->SetNameTitle("y_dist_side", "Signal and Background Distributions - y (B) ");

  TH1D* y_dist_peak = (TH1D*) reduceddata_central->createHistogram("y_dist_peak", y);
  y_dist_peak->SetMarkerColor(kRed);
  y_dist_peak->SetLineColor(kRed);
  y_dist_peak->SetNameTitle("y_dist_peak", "Signal and Background Distributions - y (B)");

  TH1D* y_dist_total = (TH1D*) data->createHistogram("y_dist_total",y);
  y_dist_total->SetMarkerColor(kBlack);
  y_dist_total->SetLineColor(kBlack);

  y_dist_peak->Add(y_dist_side, -factor);
  y_dist_side->Add(y_dist_side, factor);

  TCanvas c5;

  y_dist_total->Draw();
  y_dist_side->Draw("same");
  y_dist_peak->Draw("same");
  y_dist_peak->SetXTitle("y");
  y_dist_side->SetXTitle("y");
  y_dist_total->SetXTitle("y");

  TLegend *leg5 = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg5->AddEntry("y_dist_total", "Total", "l");
  leg5->AddEntry("y_dist_peak", "Signal", "l");
  leg5->AddEntry("y_dist_side", "Background", "l");
  leg5->Draw("same");

  c5.SetLogy();
  c5.SaveAs("y_sideband_sub.png");

  TH1D* vtxprob_dist_side = (TH1D*) reduceddata_side->createHistogram("vtxprob_dist_side",vtxprob);
  vtxprob_dist_side->SetMarkerColor(kBlue);
  vtxprob_dist_side->SetLineColor(kBlue);
  vtxprob_dist_side->SetNameTitle("vtxprob_dist_side", "Signal and Background Distributions -  #chi^{2} prob ");

  TH1D* vtxprob_dist_peak = (TH1D*) reduceddata_central->createHistogram("vtxprob_dist_peak", vtxprob);
  vtxprob_dist_peak->SetMarkerColor(kRed);
  vtxprob_dist_peak->SetLineColor(kRed);
  vtxprob_dist_peak->SetNameTitle("vtxprob_dist_peak", "Signal and Background Distributions - #chi^{2} prob");

  TH1D* vtxprob_dist_total = (TH1D*) data->createHistogram("vtxprob_dist_total",vtxprob);
  vtxprob_dist_total->SetMarkerColor(kBlack);
  vtxprob_dist_total->SetLineColor(kBlack);

  vtxprob_dist_peak->Add(vtxprob_dist_side, -factor);
  vtxprob_dist_side->Add(vtxprob_dist_side, factor);

  TCanvas c6;

  vtxprob_dist_total->Draw();
  vtxprob_dist_side->Draw("same");
  vtxprob_dist_peak->Draw("same");
  vtxprob_dist_peak->SetXTitle("#chi^{2} prob");
  vtxprob_dist_side->SetXTitle("#chi^{2} prob");
  vtxprob_dist_total->SetXTitle("#chi^{2} prob");

  TLegend *leg6 = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg6->AddEntry("vtxprob_dist_total", "Total", "l");
  leg6->AddEntry("vtxprob_dist_peak", "Signal", "l");
  leg6->AddEntry("vtxprob_dist_side", "Background", "l");
  leg6->Draw("same");

  c6.SetLogy();
  c6.SaveAs("vtxprob_sideband_sub.png");
   
  TH1D* lxy_dist_side = (TH1D*) reduceddata_side->createHistogram("lxy_dist_side",lxy);
  lxy_dist_side->SetMarkerColor(kBlue);
  lxy_dist_side->SetLineColor(kBlue);
  lxy_dist_side->SetNameTitle("lxy_dist_side", "Signal and Background Distributions -  l_{xy} ");

  TH1D* lxy_dist_peak = (TH1D*) reduceddata_central->createHistogram("lxy_dist_peak", lxy);
  lxy_dist_peak->SetMarkerColor(kRed);
  lxy_dist_peak->SetLineColor(kRed);
  lxy_dist_peak->SetNameTitle("lxy_dist_peak", "Signal and Background Distributions - l_{xy} ");

  TH1D* lxy_dist_total = (TH1D*) data->createHistogram("lxy_dist_total",lxy);
  lxy_dist_total->SetMarkerColor(kBlack);
  lxy_dist_total->SetLineColor(kBlack);

  lxy_dist_peak->Add(lxy_dist_side, -factor);
  lxy_dist_side->Add(lxy_dist_side, factor);

  TCanvas c7;

  lxy_dist_total->Draw();
  lxy_dist_side->Draw("same");
  lxy_dist_peak->Draw("same");
  lxy_dist_peak->SetXTitle("l_{xy} [cm]");
  lxy_dist_side->SetXTitle("l_{xy} [cm]");
  lxy_dist_total->SetXTitle("l_{xy} [cm]");

  TLegend *leg7 = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg7->AddEntry("lxy_dist_total", "Total", "l");
  leg7->AddEntry("lxy_dist_peak", "Signal", "l");
  leg7->AddEntry("lxy_dist_side", "Background", "l");
  leg7->Draw("same");

  c7.SetLogy();
  c7.SaveAs("lxy_sideband_sub.png");

  TH1D* errlxy_dist_side = (TH1D*) reduceddata_side->createHistogram("errlxy_dist_side",errlxy);
  errlxy_dist_side->SetMarkerColor(kBlue);
  errlxy_dist_side->SetLineColor(kBlue);
  errlxy_dist_side->SetNameTitle("errlxy_dist_side", "Signal and Background Distributions - #sigma l_{xy} ");

  TH1D* errlxy_dist_peak = (TH1D*) reduceddata_central->createHistogram("errlxy_dist_peak", errlxy);
  errlxy_dist_peak->SetMarkerColor(kRed);
  errlxy_dist_peak->SetLineColor(kRed);
  errlxy_dist_peak->SetNameTitle("errlxy_dist_peak", "Signal and Background Distributions - #sigma l_{xy} ");

  TH1D* errlxy_dist_total = (TH1D*) data->createHistogram("errlxy_dist_total",errlxy);
  errlxy_dist_total->SetMarkerColor(kBlack);
  errlxy_dist_total->SetLineColor(kBlack);

  errlxy_dist_peak->Add(errlxy_dist_side, -factor);
  errlxy_dist_side->Add(errlxy_dist_side, factor);

  TCanvas c8;

  errlxy_dist_total->Draw();
  errlxy_dist_side->Draw("same");
  errlxy_dist_peak->Draw("same");
  errlxy_dist_peak->SetXTitle("#sigma l_{xy} [cm]");
  errlxy_dist_side->SetXTitle("#sigma l_{xy} [cm]");
  errlxy_dist_total->SetXTitle("#sigma l_{xy} [cm]");

  TLegend *leg8 = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg8->AddEntry("errlxy_dist_total", "Total", "l");
  leg8->AddEntry("errlxy_dist_peak", "Signal", "l");
  leg8->AddEntry("errlxy_dist_side", "Background", "l");
  leg8->Draw("same");

  c8.SetLogy();
  c8.SaveAs("errlxy_sideband_sub.png");


  TH1D* lerrxy_dist_side = (TH1D*) reduceddata_side->createHistogram("lerrxy_dist_side",lerrxy);
  lerrxy_dist_side->SetMarkerColor(kBlue);
  lerrxy_dist_side->SetLineColor(kBlue);
  lerrxy_dist_side->SetNameTitle("lerrxy_dist_side", "Signal and Background Distributions - l_{xy}/#sigma l_{xy} ");

  TH1D* lerrxy_dist_peak = (TH1D*) reduceddata_central->createHistogram("lerrxy_dist_peak", lerrxy);
  lerrxy_dist_peak->SetMarkerColor(kRed);
  lerrxy_dist_peak->SetLineColor(kRed);
  lerrxy_dist_peak->SetNameTitle("lerrxy_dist_peak", "Signal and Background Distributions - l_{xy}/#sigma l_{xy} ");

  TH1D* lerrxy_dist_total = (TH1D*) data->createHistogram("lerrxy_dist_total",lerrxy);
  lerrxy_dist_total->SetMarkerColor(kBlack);
  lerrxy_dist_total->SetLineColor(kBlack);

  lerrxy_dist_peak->Add(lerrxy_dist_side, -factor);
  lerrxy_dist_side->Add(lerrxy_dist_side, factor);

  TCanvas c9;

  lerrxy_dist_total->Draw();
  lerrxy_dist_side->Draw("same");
  lerrxy_dist_peak->Draw("same");
  lerrxy_dist_peak->SetXTitle("l_{xy}/#sigma l_{xy} ");
  lerrxy_dist_side->SetXTitle("l_{xy}/#sigma l_{xy} ");
  lerrxy_dist_total->SetXTitle("l_{xy}/#sigma l_{xy} ");

  TLegend *leg9 = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg9->AddEntry("lerrxy_dist_total", "Total", "l");
  leg9->AddEntry("lerrxy_dist_peak", "Signal", "l");
  leg9->AddEntry("lerrxy_dist_side", "Background", "l");
  leg9->Draw("same");

  c9.SetLogy();
  c9.SaveAs("lerrxy_sideband_sub.png");


}
