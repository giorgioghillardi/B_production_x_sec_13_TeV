#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

using namespace RooFit;

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda
// channel = 7: Bc -> J/psi Pi

//input example: data_selection --channel 1 --input /some/place/file.root --output /some/place
int main(int argc, char** argv)
{
  int channel = 1;

  TString input_file = "";
  TString output_dir = "";
  
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
      if(argument == "--output")
        {
          convert << argv[++i];





        convert >> output_dir;
        }
    }

  if(input_file=="")
    {
      std::cout << "No input file provided. Please choose one using --input." << std::endl;
      return 0;
    }

  //output file name
  TString data_selection_output_file = output_dir + "selected_" + input_file;

  //output ntuple
  //TNtupleD *nt = new TNtupleD(channel_to_ntuple_name(channel),channel_to_ntuple_name(channel),"mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");

  TLorentzVector v4_tk1, v4_tk2;

  TChain* tin;
  TString tree_name = channel_to_ntuple_name(channel);
  
  tin = new TChain(tree_name);
  tin->Add(input_file);
  
  ReducedBranches br;
  br.setbranchadd(tin);

  //output file structure
  TFile *fout = new TFile(data_selection_output_file,"recreate");
  
  TTree *tree_out = new TTree(tree_name,tree_name);
  br.regTree(tree_out);
  
  std::cout << "selecting data from channel " << channel << std::endl;
 
  int n_entries = (int) tin->GetEntries();
  int percent = (int)(0.01*n_entries);

  std::vector<std::string> particle_flow_string = {"signal","tktk_mass_window","tktk_veto"};
  std::vector<int> particle_flow_number;

  for(int i=0; i < (int)particle_flow_string.size(); i++)
    particle_flow_number.push_back(0);
  
  TH1D* ditrack_mass; 
  TH1D* ditrack_veto;
  TCanvas* cv1;

  double phi_window = 0.010;
  double phi_veto = 0.010;
  
  double k_star_window = 0.050;
  double k_star_veto = 0.050;

  printf("Filling histograms without mass or veto cuts \n");

  cv1 = new TCanvas("cv1", "First canvas", 1920, 1080);
  ditrack_mass = new TH1D("ditrack_mass","tktk_mass_" + channel_to_ntuple_name(channel), 100, 0.7, 1.2);
  ditrack_veto = new TH1D("ditrack_veto","tktk_veto_" + channel_to_ntuple_name(channel), 100, 0.7, 1.2);

  for (int evt=0; evt<n_entries; evt++) 
    {
      //show the percentage of processed events
      if (evt%percent==0 || evt==n_entries-1) printf("processing %d/%d (%.2f%%).\n",evt,n_entries-1,(double)evt/(double)(n_entries-1)*100.);

      tin->GetEntry(evt);

      switch(channel)
	{
	case 2:
	  ditrack_mass->Fill(br.tktkmass);

	  v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
	  v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
	  
	  ditrack_veto->Fill(fabs((v4_tk1+v4_tk2).Mag()));
	  break;

	case 4:
	  ditrack_mass->Fill(br.tktkmass);
	  
	  v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
	  v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
	  ditrack_veto->Fill(fabs((v4_tk1+v4_tk2).Mag()));
 	  
      	  v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
	  v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
	  ditrack_veto->Fill(fabs((v4_tk1+v4_tk2).Mag()));
	  break;
	}
    }//end of event cicle
  
  cv1->Divide(2,1);
  cv1->cd(1);
  ditrack_mass->Draw();
  cv1->cd(2);
  ditrack_veto->Draw();
  cv1->SaveAs("tktk_mass_and_veto_" + channel_to_ntuple_name(channel) + "_no_window_no_veto.png");

  delete cv1;
  delete ditrack_mass;
  delete ditrack_veto;
  
  printf("Filling histograms with mass or veto cuts \n");

  cv1 = new TCanvas("cv1", "First canvas", 1920, 1080);
  ditrack_mass = new TH1D("ditrack_mass","tktk_mass_" + channel_to_ntuple_name(channel), 100, 0.7, 1.2);
  ditrack_veto = new TH1D("ditrack_veto","tktk_veto_" + channel_to_ntuple_name(channel), 100, 0.7, 1.2);

  for (int evt=0; evt<n_entries; evt++) 
    {
      //show the percentage of processed events
      if (evt%percent==0 || evt==n_entries-1) printf("processing %d/%d (%.2f%%).\n",evt,n_entries-1,(double)evt/(double)(n_entries-1)*100.);

      tin->GetEntry(evt);

      particle_flow_number[0]++;

      switch(channel)
	{
	case 2:
	  if (fabs(br.tktkmass-KSTAR_MASS)>=k_star_window) continue;
	  particle_flow_number[1]++;
	  
     	  v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
	  v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
	  if (fabs((v4_tk1+v4_tk2).Mag()-PHI_MASS)<=phi_veto) continue;

	  particle_flow_number[2]++;
	  break;

	case 4:
	  if (fabs(br.tktkmass-PHI_MASS)>=phi_window) continue;
	  particle_flow_number[1]++;
	  
	  v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
	  v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
	  if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=k_star_veto) continue;
 	  
      	  v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
	  v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
	  if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=k_star_veto) continue;

	  particle_flow_number[2]++;
	  break;
	}
      
      ditrack_mass->Fill(br.tktkmass);
      ditrack_veto->Fill(fabs((v4_tk1+v4_tk2).Mag()));
      
      tree_out->Fill();

    }//end of event cicle

  cv1->Divide(2,1);
  cv1->cd(1);
  ditrack_mass->Draw();
  cv1->cd(2);
  ditrack_veto->Draw();
  cv1->SaveAs("tktk_mass_and_veto_" + channel_to_ntuple_name(channel) + ".png");
  
  delete cv1;
  delete ditrack_mass;
  delete ditrack_veto;
  
  std::cout << "cut   :   number   :   eff" << std::endl;
  
  for(int i=1; i < (int)particle_flow_number.size(); i++)
    {
      std::cout << particle_flow_string[i] << " : " << particle_flow_number[i] << " : " << (double)particle_flow_number[i]/(double)particle_flow_number[i-1] << std::endl;
    }

  fout->Write();
  fout->Close();
}
