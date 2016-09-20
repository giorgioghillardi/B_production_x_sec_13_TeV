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

int main(int argc, char** argv)
{
  int channel = 0;

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--channel")
        {
          convert << argv[++i];
          convert >> channel;
        }
    }
  
  if(channel==0)
    {
      std::cout << "The option --channel can be used to choose a channel. Example myloop_gen --channel 2" << std::endl;
    }
  
  TChain *root = new TChain("demo/root");
  
  switch(channel)
    {
    default:
    case 1:
      //for Bfilter processed with DumpGenInfo.py
      //root->Add("/gstore/t3cms/store/user/martinsg/no_pre_filter_MC_v1/BuToJpsiKV2_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_no_pre_filter_v1/160719_165646/0000/no_pre_filter_MC_*.root");

      //for Bfilter processed with Bfinder_mc.cc
      root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bu_Bfilter_v1/BuToJpsiKV2_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bu_Bfilter_v1/160812_095709/0000/Bfinder_mc_*.root");
      break;
    case 2:
      //for Bfilter processed with Bfinder_mc.cc
      root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bd_Bfilter_v1/BdToJpsiKstarV2_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bd_Bfilter_v1/160812_115744/0000/Bfinder_mc_*.root");
      break;
    case 3:
      break;
    case 4:
      //for Bfilter processed with Bfinder_mc.cc
      root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bs_Bfilter_v3/BsToJpsiPhi_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bs_Bfilter_v3/160812_235643/0000/Bfinder_mc_*.root");
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
  
  TString file_name = "myloop_gen_" + channel_to_ntuple_name(channel) + "_bfilter.root";
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
	      if (GenInfo->pdgId[idx_jpsi]!=443) continue; // not J/psi
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
	      if (GenInfo->pdgId[idx_jpsi]!=443) continue; // not J/psi
	      if (abs(GenInfo->pdgId[idx_tktk])!=313) continue; // not K*0
	      if ((GenInfo->pdgId[idx_tk1]!=321 || GenInfo->pdgId[idx_tk2]!=-211) && (GenInfo->pdgId[idx_tk1]!=-321 || GenInfo->pdgId[idx_tk2]!=211)) continue; //not k+pi- and k-pi+
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
	     	           
	      if (abs(GenInfo->pdgId[idx_b])!=531) continue; // not Bs
	      if (GenInfo->pdgId[idx_jpsi]!=443) continue; // not J/psi
	      if (abs(GenInfo->pdgId[idx_tktk])!=333) continue; // not phi
	      if ((GenInfo->pdgId[idx_tk1]!=321 || GenInfo->pdgId[idx_tk2]!=-321) && (GenInfo->pdgId[idx_tk1]!=-321 || GenInfo->pdgId[idx_tk2]!=321)) continue; //not k+k- and k-k+
     	      if (abs(GenInfo->pdgId[idx_mu1])!=13) continue; // not mu+-
	      if (abs(GenInfo->pdgId[idx_mu2])!=13) continue; // not mu+
	      break;

	    case 5:
	      break;
	    case 6:
	      break;
	    }

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
  
  fout->Write();
  fout->Close();
  
  delete GenInfo;
}

/*
for (int evt=0; evt<n_entries; evt++) 
    {
      if (evt%1000==0 || evt==n_entries-1) printf("processing %d/%d (%.2f%%).\n",evt,n_entries-1,(double)evt/(double)(n_entries-1)*100.);
      
      root->GetEntry(evt);
      
      // Look for indices of the whole decay tree
      for (int idx = 0; idx < GenInfo->size; idx++) 
	{
	  if (abs(GenInfo->pdgId[idx])==521) // B+ find
	    {
	      int idx_bp   = idx;
	      int idx_jpsi = GenInfo->da1[idx_bp];
	      int idx_kp   = GenInfo->da2[idx_bp];
	      int idx_mu1  = GenInfo->da1[idx_jpsi];
	      int idx_mu2  = GenInfo->da2[idx_jpsi];
	      
	    if (GenInfo->pdgId[idx_jpsi]!=443) continue; // not J/psi
	    if (abs(GenInfo->pdgId[idx_kp])!=321) continue; // not K+-
	    if (abs(GenInfo->pdgId[idx_mu1])!=13) continue; // not mu+-
	    if (abs(GenInfo->pdgId[idx_mu2])!=13) continue; // not mu+-
	    
	    TLorentzVector v4_bp, v4_uj;
	    v4_bp.SetPtEtaPhiM(GenInfo->pt[idx_bp],GenInfo->eta[idx_bp],GenInfo->phi[idx_bp],GenInfo->mass[idx_bp]);
	    v4_uj.SetPtEtaPhiM(GenInfo->pt[idx_jpsi],GenInfo->eta[idx_jpsi],GenInfo->phi[idx_jpsi],GenInfo->mass[idx_jpsi]);
	    
	    br->mass    = GenInfo->mass[idx_bp];
	    br->pt      = GenInfo->pt[idx_bp];
	    br->eta     = GenInfo->eta[idx_bp];
	    br->phi     = GenInfo->phi[idx_bp];
	    br->y       = v4_bp.Rapidity();
	    br->vx      = GenInfo->vx[idx_bp];
	    br->vy      = GenInfo->vy[idx_bp];
	    br->vz      = GenInfo->vz[idx_bp];
	    
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
	    
	    br->tk1pt   = GenInfo->pt[idx_kp];
	    br->tk1eta  = GenInfo->eta[idx_kp];
	    br->tk1phi  = GenInfo->phi[idx_kp];
	    br->tk1charge  = GenInfo->pdgId[idx_kp]/abs(GenInfo->pdgId[idx_kp]);
	    
	    nt->Fill();
	  }
      } //end of GenInfo loop
  } // end of evt loop
 */
