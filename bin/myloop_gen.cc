#include <TChain.h>
#include <iostream>
#include <sstream>
#include <TTree.h>
#include <TH1D.h>
#include <TNtupleD.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "UserCode/B_production_x_sec_13_TeV/interface/format.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/myloop.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/channel.h"

//myloop_gen --channel 2 --bfilter 1 --output some/place
int main(int argc, char** argv)
{
  int channel = 1;
  int bfilter = 1;
  std::string output_dir ="";

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--channel")
        {
          convert << argv[++i];
          convert >> channel;
        }
      if(argument == "--bfilter")
        {
          convert << argv[++i];
          convert >> bfilter;
        }
      if(argument == "--output")
        {
          convert << argv[++i];
          convert >> output_dir;
        }
    }
  
  if(channel==0)
    {
      std::cout << "The option --channel can be used to choose a channel. Example myloop_gen --channel 2" << std::endl;
    }
  
  TChain *root = new TChain("analysis/root");
  
  switch(channel)
    {
    default:
    case 1:
      if(bfilter)
	{
	  //for Bfilter processed with Bfinder_mc.cc
	  root->Add("/gstore/t3cms/store/user/bfontana/Bfinder_Bp_MC2016/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/Pincopallino_Bp/170920_100321/0000/bfinder_*.root");
	}
      else
	{
	  //for BMuonfilter processed with Bfinder_mc.cc
	  root->Add("/gstore/t3cms/store/user/bfontana/Bfinder_Bp_MC2016_nocuts/BuToJpsiK_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/Pincopallino/170921_152917/0000/bfinder_*.root");
	}
      break;
    case 2:
      if(bfilter)
	{
	  //for Bfilter processed with Bfinder_mc.cc
	  //root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bd_Bfilter_v1/BdToJpsiKstarV2_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bd_Bfilter_v1/160812_115744/0000/Bfinder_mc_*.root");
	  
	  root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bd_Bfilter_ext_v1/BdToJpsiKstarV2_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bd_Bfilter_ext_v1/170121_121833/0000/Bfinder_mc_*.root");
	}
      else
	{
	  //for BMuonfilter processed with Bfinder_mc.cc
	  root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bd_muonfilter_v1/BdToJpsiKstarV2_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bd_muonfilter_v1/160812_135133/0000/Bfinder_mc_*.root");
	}
      break;
    case 3:
      break;
    case 4:
      if(bfilter)
	{
	  //for Bfilter processed with Bfinder_mc.cc
	  //root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bs_Bfilter_v3/BsToJpsiPhi_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bs_Bfilter_v3/160812_235643/0000/Bfinder_mc_*.root");
	  
	  root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bs_Bfilter_ext_v1/BsToJpsiPhiV2_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bs_Bfilter_ext_v1/170120_150004/0000/Bfinder_mc_*.root");
	}
      else
	{
	  //for BMuonfilter processed with Bfinder_mc.cc
	  root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bs_muonfilter_v1/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bs_muonfilter_v1/160812_151233/0000/Bfinder_mc_*.root");
	}
      break;
    case 5:
      break;
    case 6:
      break;
    }
  
  int n_entries = root->GetEntries();
  printf("Going to process %d entries.\n",n_entries);
  
  //-----------------------------------------------------------------
  // setting memory addresses to the branches
  // create new output trees
  
  GenInfoBranches *GenInfo = new GenInfoBranches;
  
  GenInfo->setbranchadd(root);
  
  TString bf = "";
  
  if(bfilter)
    bf = "bfilter";
  else
    bf = "bmuonfilter";
  
  TString file_name = output_dir + "myloop_gen_" + channel_to_ntuple_name(channel) + "_" + bf + ".root";
  TFile *fout = new TFile(file_name,"recreate");
  
  ReducedGenBranches *br = new ReducedGenBranches();
  TString tree_name = channel_to_ntuple_name(channel) + "_gen";

  TTree *nt = new TTree(tree_name,tree_name);
  br->regTree(nt);
  
  int idx_b    = -1;
  int idx_jpsi = -1;
  int idx_tktk = -1;
  int idx_tk1  = -1;
  int idx_tk2  = -1;
  int idx_mu1  = -1;
  int idx_mu2  = -1;
  
  int b_counter = 0;

  for (int evt=0; evt<n_entries; evt++) 
    {
      if (evt%1000==0 || evt==n_entries-1) printf("processing %d/%d (%.2f%%).\n",evt,n_entries-1,(double)evt/(double)(n_entries-1)*100.);
      
      root->GetEntry(evt);

      // Look for indices of the whole decay tree
      for (int idx = 0; idx < GenInfo->size; idx++) 
	{
	  switch(channel)
	    {
	    default:
	    case 1:
	      idx_b    = idx;
	      idx_jpsi = GenInfo->da1[idx_b];
	      idx_tk1  = GenInfo->da2[idx_b];
	      idx_mu1  = GenInfo->da1[idx_jpsi];
	      idx_mu2  = GenInfo->da2[idx_jpsi];
	      
	      if (abs(GenInfo->pdgId[idx_b])!=521) continue; // not B+
	      if (abs(GenInfo->pdgId[idx_jpsi])!=443) continue; // not J/psi
	      if (abs(GenInfo->pdgId[idx_tk1])!=321) continue; // not K+-
	      if (abs(GenInfo->pdgId[idx_mu1])!=13) continue; // not mu+-
	      if (abs(GenInfo->pdgId[idx_mu2])!=13) continue; // not mu+-
	      break;
	    
	    case 2:
	      idx_b    = idx;
	      idx_jpsi = GenInfo->da1[idx_b];
	      idx_tktk = GenInfo->da2[idx_b];
	      idx_tk1  = GenInfo->da1[idx_tktk];
	      idx_tk2  = GenInfo->da2[idx_tktk];
	      idx_mu1  = GenInfo->da1[idx_jpsi];
	      idx_mu2  = GenInfo->da2[idx_jpsi];
	     	           
	      if (abs(GenInfo->pdgId[idx_b])!=511) continue; // not B0
	      if (abs(GenInfo->pdgId[idx_jpsi])!=443) continue; // not J/psi
	      if (abs(GenInfo->pdgId[idx_tktk])!=313) continue; // not K*0
	      if ((GenInfo->pdgId[idx_tk1]!=321 || GenInfo->pdgId[idx_tk2]!=-211) && //not k+pi-
		  (GenInfo->pdgId[idx_tk1]!=-321 || GenInfo->pdgId[idx_tk2]!=211) && //not k-pi+
		  (GenInfo->pdgId[idx_tk1]!=211 || GenInfo->pdgId[idx_tk2]!=-321) && //not pi+k-
		  (GenInfo->pdgId[idx_tk1]!=-211 || GenInfo->pdgId[idx_tk2]!=321)) continue; //not k+pi-
     	      if (abs(GenInfo->pdgId[idx_mu1])!=13) continue; // not mu+-
	      if (abs(GenInfo->pdgId[idx_mu2])!=13) continue; // not mu+-
	      break;
	      
	    case 3:
	      break;
	    case 4:
	      idx_b    = idx;
	      idx_jpsi = GenInfo->da1[idx_b];
	      idx_tktk = GenInfo->da2[idx_b];
	      idx_tk1  = GenInfo->da1[idx_tktk];
	      idx_tk2  = GenInfo->da2[idx_tktk];
	      idx_mu1  = GenInfo->da1[idx_jpsi];
	      idx_mu2  = GenInfo->da2[idx_jpsi];
	      
	      if(abs(GenInfo->pdgId[idx_b])!=531) continue; // not Bs
	      if(abs(GenInfo->pdgId[idx_jpsi])!=443) continue; // not J/psi
	      if(abs(GenInfo->pdgId[idx_tktk])!=333) continue; // not phi
	      if((GenInfo->pdgId[idx_tk1]!=321 || GenInfo->pdgId[idx_tk2]!=-321) && (GenInfo->pdgId[idx_tk1]!=-321 || GenInfo->pdgId[idx_tk2]!=321)) continue; //not k+k- and k-k+
	      if(abs(GenInfo->pdgId[idx_mu1])!=13) continue; // not mu+-
	      if(abs(GenInfo->pdgId[idx_mu2])!=13) continue; // not mu+
	      break;

	    case 5:
	      break;
	    case 6:
	      break;
	    }
	  
	  b_counter++;

	  TLorentzVector v4_b, v4_uj;
	  v4_b.SetPtEtaPhiM(GenInfo->pt[idx_b],GenInfo->eta[idx_b],GenInfo->phi[idx_b],GenInfo->mass[idx_b]);
	  v4_uj.SetPtEtaPhiM(GenInfo->pt[idx_jpsi],GenInfo->eta[idx_jpsi],GenInfo->phi[idx_jpsi],GenInfo->mass[idx_jpsi]);
	  
	  br->mass    = GenInfo->mass[idx_b];
	  br->pt      = GenInfo->pt[idx_b];
	  br->eta     = GenInfo->eta[idx_b];
	  br->phi     = GenInfo->phi[idx_b];
	  br->y       = v4_b.Rapidity();
	  br->vx      = GenInfo->vx[idx_b];
	  br->vy      = GenInfo->vy[idx_b];
	  br->vz      = GenInfo->vz[idx_b];
	  
	  br->ujmass  = GenInfo->mass[idx_jpsi];
	  br->ujpt    = GenInfo->pt[idx_jpsi];
	  br->ujeta   = GenInfo->eta[idx_jpsi];
	  br->ujphi   = GenInfo->phi[idx_jpsi];
	  br->ujy     = v4_uj.Rapidity();
	  br->ujvx    = GenInfo->vx[idx_jpsi];
	  br->ujvy    = GenInfo->vy[idx_jpsi];
	  br->ujvz    = GenInfo->vz[idx_jpsi];
	  
	  br->mu1pt   = GenInfo->pt[idx_mu1];
	  br->mu1eta  = GenInfo->eta[idx_mu1];
	  br->mu1phi  = GenInfo->phi[idx_mu1];

	  br->mu2pt   = GenInfo->pt[idx_mu2];
	  br->mu2eta  = GenInfo->eta[idx_mu2];
	  br->mu2phi  = GenInfo->phi[idx_mu2];
	  
	  br->tk1pt   = GenInfo->pt[idx_tk1];
	  br->tk1eta  = GenInfo->eta[idx_tk1];
	  br->tk1phi  = GenInfo->phi[idx_tk1];
	  br->tk1charge  = GenInfo->pdgId[idx_tk1]/abs(GenInfo->pdgId[idx_tk1]);

	  /*
	  br->tk2pt   = GenInfo->pt[idx_tk2];
	  br->tk2eta  = GenInfo->eta[idx_tk2];
	  br->tk2phi  = GenInfo->phi[idx_tk2];
	  br->tk2charge  = GenInfo->pdgId[idx_tk2]/abs(GenInfo->pdgId[idx_tk2]);
	  */

	  nt->Fill();
	} //end of GenInfo loop
    } // end of evt loop
  
  //debug: print the number of signal.
  std::cout << "Channel : " << channel << " Number of signal : " << b_counter << std::endl;
    
    
  fout->Write();
  fout->Close();
  
  delete GenInfo;
}
