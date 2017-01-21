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
#include "UserCode/B_production_x_sec_13_TeV/interface/channel.h"
#include "TMath.h"
using namespace RooFit;

#define BASE_DIR "/lstore/cms/brunogal/input_for_B_production_x_sec_13_TeV/"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda

//input example: reduce_mc --channel 1 --gen 1 --input /some/place/
int main(int argc, char** argv)
{
  int channel = 1;
  int gen = 1;
  TString input_file = "";

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--channel")
        {
          convert << argv[++i];
          convert >> channel;
        }
      
      if(argument == "--gen")
	{
	  convert << argv[++i];
          convert >> gen;
	}
      
      if(argument == "--input")
        {
          convert << argv[++i];
          convert >> input_file;
        }
    }
  
  if(input_file == "")
    {
      printf("ERROR: no input file was indicated. Please use --input /some/file \n");
    }
    
  TString data_selection_output_file = "reduced_" + input_file;
  
  TFile *fout = new TFile(data_selection_output_file,"recreate");
  
  TNtupleD *_nt1 = new TNtupleD("ntkp","ntkp","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  TNtupleD *_nt2 = new TNtupleD("ntkstar","ntkstar","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  TNtupleD *_nt3 = new TNtupleD("ntks","ntks","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  TNtupleD *_nt4 = new TNtupleD("ntphi","ntphi","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  TNtupleD *_nt5 = new TNtupleD("ntmix","ntmix","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  TNtupleD *_nt6 = new TNtupleD("ntlambda","ntlambda","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  
  ReducedGenBranches br;
  TChain* tin;

  std::cout << "Reducing the size of channel " << channel << std::endl;

  TString ntuple_name = channel_to_ntuple_name(channel);
  
  if(gen)
    ntuple_name += "_gen";

  tin = new TChain(ntuple_name);

  tin->Add(input_file);

  br.setbranchadd(tin);

  int n_entries = tin->GetEntries();

  for (int evt=0; evt<n_entries; evt++)
    {
      if (evt%1000==0 || evt==n_entries-1) printf("processing %d/%d (%.2f%%).\n",evt,n_entries-1,(double)evt/(double)(n_entries-1)*100.);

      tin->GetEntry(evt);
    
      switch(channel)
	{
	default:
	case 1:
	  _nt1->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	  break;
	case 2:
	  _nt2->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	  break;
	case 3:
	  _nt3->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	  break;
	case 4:
	  _nt4->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	  break;
	case 5:
	  _nt5->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	  break;
	case 6:
	  _nt6->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	  break;
	}
    }

  fout->Write();
  fout->Close();
}
