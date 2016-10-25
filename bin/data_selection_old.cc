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
#include "UserCode/B_production_x_sec_13_TeV/interface/channel.h"
#include "TMath.h"
using namespace RooFit;

#define BASE_DIR "/lstore/cms/brunogal/input_for_B_production_x_sec_13_TeV/"

#define SHOW_DIST 0
//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda

/*
void plot_pt_dist(RooWorkspace& w, int channel, TString directory);
void plot_mass_dist(RooWorkspace& w, int channel, TString directory);
void read_data(RooWorkspace& w, TString filename,int channel);
void read_data_cut(RooWorkspace& w, RooDataSet* data);
void set_up_workspace_variables(RooWorkspace& w, int channel);
*/
void data_selection(TString fin1, int run_on_mc ,TString data_selection_output_file,int channel);

//input example: data_selection_old --channel 1 --mc 1 --input /some/place/
int main(int argc, char** argv)
{
  int channel = 0;
  TString input_file = TString::Format(BASE_DIR) + "myloop_data.root";
  int run_on_mc =0;

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--channel")
        {
          convert << argv[++i];
          convert >> channel;
        }

      if(argument == "--mc")
        {
          convert << argv[++i];
          convert >> run_on_mc;
        }

      if(argument == "--input")
        {
          convert << argv[++i];
          convert >> input_file;
        }
    }

  if(channel==0)
    {
      std::cout << "No channel was provided as input. Please use --channel. Example: data_selection --channel 1" << std::endl;
      return 0;
    }

  TString data_selection_output_file="";
  if(run_on_mc)
    {
      data_selection_output_file= TString::Format(BASE_DIR) + "selected_mc_" + channel_to_ntuple_name(channel) + "_bmuonfilter_with_cuts.root";
    }
  else
    {
      data_selection_output_file= TString::Format(BASE_DIR) + "selected_data_" + channel_to_ntuple_name(channel) + ".root";
    }
  
  data_selection(input_file, run_on_mc ,data_selection_output_file,channel);
  /*    
  if(SHOW_DIST)
    { 
      RooWorkspace* ws = new RooWorkspace("ws","Bmass");
      TString pt_dist_directory="";
      TString mass_dist_directory="";
      
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
  */
}

void data_selection(TString fin1, int run_on_mc, TString data_selection_output_file,int channel){

    TFile *fout = new TFile(data_selection_output_file,"recreate");

    TNtupleD *_nt1 = new TNtupleD("ntkp","ntkp","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
    TNtupleD *_nt2 = new TNtupleD("ntkstar","ntkstar","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
    TNtupleD *_nt3 = new TNtupleD("ntks","ntks","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
    TNtupleD *_nt4 = new TNtupleD("ntphi","ntphi","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
    TNtupleD *_nt5 = new TNtupleD("ntmix","ntmix","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
    TNtupleD *_nt6 = new TNtupleD("ntlambda","ntlambda","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");

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
	if(run_on_mc)
	  {
	     if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1) continue;
	  }
	else
	  {
	    if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1 && br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
	  }
	
	if (br.tk1pt<=1.6) continue;
	if (br.vtxprob<=0.2) continue; //original cut 0.1
	if (br.lxy/br.errxy<=4.5) continue; //original cut 3.0
	if (br.cosalpha2d<=0.996) continue; //original cut 0.99
	            
	_nt1->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	    
      }else
        if (channel==2) { // cuts for B0 -> J/psi K*
	  if(run_on_mc)
	    {
	      if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1) continue;
	    }
	  else
	  {
	    if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1 && br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
	  }
	  if (br.vtxprob<=0.2) continue; //original cut 0.1
	  if (br.lxy/br.errxy<=4.5) continue; //original cut 3.0
	  if (br.cosalpha2d<=0.996) continue; //original cut 0.99
	  if (fabs(br.tktkmass-KSTAR_MASS)>=0.05) continue;
            
	  v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
	  v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
	  if (fabs((v4_tk1+v4_tk2).Mag()-PHI_MASS)<=0.01) continue;
            
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
		_nt2->Fill(br_queue[i].mass,br_queue[i].pt,br_queue[i].eta,br_queue[i].y,br_queue[i].mu1pt,br_queue[i].mu1eta,br_queue[i].mu2pt,br_queue[i].mu2eta);

	      }

	    }                
	    n_br_queued = 0;
	    memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
	    n_br_queued++;
	  }
	}else
	  if (channel==3) { // cuts for B0 -> J/psi Ks
	    if(run_on_mc)
	      {
		if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1) continue;
	      }
	    else
	  {
	    if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1 && br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
	  }
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.tktkblxy/br.tktkberrxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-KSHORT_MASS)>=0.015) continue;
                
            _nt3->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);

	  }else
	    if (channel==4) { // cuts for Bs -> J/psi phi
	      if(run_on_mc)
		{
		  if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1) continue;
		}
	      else
		{
		  if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1 && br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
		}
	      if (br.vtxprob<=0.2) continue; //original cut 0.1
	      if (br.lxy/br.errxy<=4.5) continue; //original cut 3.0
	      if (br.cosalpha2d<=0.996) continue; //original cut 0.99
	      if (fabs(br.tktkmass-PHI_MASS)>=0.010) continue;
            
	      v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
	      v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
	      if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
	      v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
	      v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
	      if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
                
	      _nt4->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);

	    }else
	      if (channel==5) { // cuts for psi(2S)/X(3872) -> J/psi pipi
		if (br.vtxprob<=0.2) continue;
		if (fabs(br.tk1eta)>=1.6) continue;
		if (fabs(br.tk2eta)>=1.6) continue;
            
		_nt5->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);

	      }else
		if (channel==6) {//cuts for lambda
		  if(run_on_mc)
		    {
		      if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1) continue;
		    }
		  else
		    {
		      if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1 && br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
		    }
		  if (br.vtxprob<=0.1) continue;
		  if (br.lxy/br.errxy<=3.0) continue;
		  if (br.tktkblxy/br.tktkberrxy<=3.0) continue;
		  if (br.cosalpha2d<=0.99) continue;
		  if (fabs(br.tktkmass-LAMBDA_MASS)>=0.015) continue;
            
		  v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
		  v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
		  if (fabs((v4_tk1+v4_tk2).Mag()-KSHORT_MASS)<=0.015) continue;
            
		  _nt6->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);

		}
    }//end of the for for the events
    fout->Write();
    fout->Close();
}

/*
void plot_pt_dist(RooWorkspace& w, int channel, TString directory)
{
  //full dataset pt distribution
  RooRealVar pt = *(w.var("pt"));
  RooAbsData* data = w.data("data");
  
  TCanvas c2;
  TH1D* pt_dist = (TH1D*)data->createHistogram("pt_dist",pt);
  pt_dist->Draw();
  c2.SetLogy();
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
  c2.SaveAs(directory + ".png");
}

void read_data(RooWorkspace& w, TString filename,int channel)
{
  TFile* f = new TFile(filename);
  TNtupleD* _nt = (TNtupleD*)f->Get(channel_to_ntuple_name(channel));
 
  RooDataSet* data = new RooDataSet("data","data",_nt,RooArgSet( *(w.var("mass")) , *(w.var("pt")) ));
  
  w.import(*data);
}

void read_data_cut(RooWorkspace& w, RooDataSet* data)
{
  w.import(*data);
}

void set_up_workspace_variables(RooWorkspace& w, int channel)
{
  double mass_min, mass_max;
  double pt_min, pt_max;

  pt_min=0;
  pt_max=400;
  
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

  w.import(mass);
  w.import(pt);
}
*/
