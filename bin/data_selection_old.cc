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

//input example: data_selection_old --channel 2 --mc 1 --input /some/place/ --output /some/place
int main(int argc, char** argv)
{
  int channel = 0;
  TString input_file = TString::Format(BASE_DIR) + "myloop_data.root";
  TString output_dir = "";
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
      if(argument == "--output")
        {
          convert << argv[++i];
          convert >> output_dir;
        }
    }

  if(channel==0)
    {
      std::cout << "No channel was provided as input. Please use --channel. Example: data_selection --channel 1" << std::endl;
      return 0;
    }

  //output file name
  TString data_selection_output_file="";

  if(output_dir=="")
    {
      if(run_on_mc)
	data_selection_output_file= TString::Format(BASE_DIR) + "selected_mc_" + channel_to_ntuple_name(channel) + "_bmuonfilter_with_cuts.root";
      else
	data_selection_output_file= TString::Format(BASE_DIR) + "selected_data_" + channel_to_ntuple_name(channel) + ".root";
    }
  else
    {
      if(run_on_mc)
	data_selection_output_file= output_dir + "selected_mc_" + channel_to_ntuple_name(channel) + "_bmuonfilter_with_cuts.root";
      else
	data_selection_output_file= output_dir + "selected_data_" + channel_to_ntuple_name(channel) + ".root";
    }

  /*
  TNtupleD *_nt1 = new TNtupleD("ntkp","ntkp","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  TNtupleD *_nt2 = new TNtupleD("ntkstar","ntkstar","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  TNtupleD *_nt3 = new TNtupleD("ntks","ntks","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  TNtupleD *_nt4 = new TNtupleD("ntphi","ntphi","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  TNtupleD *_nt5 = new TNtupleD("ntmix","ntmix","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  TNtupleD *_nt6 = new TNtupleD("ntlambda","ntlambda","mass:pt:eta:y:mu1pt:mu1eta:mu2pt:mu2eta");
  */
   
  int n_br_queued = 0;
  ReducedBranches br_queue[32];
  TLorentzVector v4_tk1, v4_tk2;

  TChain* tin;
  tin = new TChain(channel_to_ntuple_name(channel));
  tin->Add(input_file);  
  
  ReducedBranches br;
  br.setbranchadd(tin);

  //output file structure
  TFile *fout = new TFile(data_selection_output_file,"recreate");
  
  TTree *tree_out = new TTree(channel_to_ntuple_name(channel),channel_to_ntuple_name(channel));
  br.regTree(tree_out);
  
  std::cout << "selecting data from channel " << channel << std::endl;
  std::cout << "start of the evt cicle, with: " << tin->GetEntries() << " events." << std::endl;
 
  for (int evt=0;evt<tin->GetEntries();evt++) 
    {
      //show the percentage of processed events
      if(evt%1000 == 0)
	std::cout << "processed: " << evt << " / " << tin->GetEntries() << std::endl;
	 
      tin->GetEntry(evt);
        
      if (channel==1) // cuts for B+ -> J/psi K+
	{
	  if(run_on_mc)
	    {if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1) continue;}
	  else
	    {if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1 && br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;}
	  
	  if (br.tk1pt<=1.6) continue;
	  if (br.vtxprob<=0.2) continue; //original cut 0.1
	  if (br.lxy/br.errxy<=4.5) continue; //original cut 3.0
	  if (br.cosalpha2d<=0.996) continue; //original cut 0.99
	    
	  //_nt1->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);	    
	}
	
      if (channel==2)  // cuts for B0 -> J/psi K*
	{
	  if(run_on_mc)
	    {if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1) continue;}
	  else
	    {if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v1]!=1 && br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;}
	  
	  if (br.vtxprob<=0.2) continue; //original cut 0.1
	  if (br.lxy/br.errxy<=4.5) continue; //original cut 3.0
	  if (br.cosalpha2d<=0.996) continue; //original cut 0.99
	  
	  if (fabs(br.tktkmass-KSTAR_MASS)>=0.05) continue;
 
	  v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
	  v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
	  if (fabs((v4_tk1+v4_tk2).Mag()-PHI_MASS)<=0.01) continue;
 
	  /*
	  if (n_br_queued==0)
	    {
	      memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
	      n_br_queued++;
	    }
	  else
	    if (br.run == br_queue[n_br_queued-1].run && br.event == br_queue[n_br_queued-1].event)  // same event
	      {
		memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
		n_br_queued++;

		if (n_br_queued>=32) printf("Warning: maximum queued branches reached.\n");
	      }
	  
	  if (br.run != br_queue[n_br_queued-1].run || br.event != br_queue[n_br_queued-1].event || evt==tin->GetEntries()-1) 
	    {
	      for (int i=0; i<n_br_queued; i++)
		{
		  bool isBestKstarMass = true;
		    
		  for (int j=0; j<n_br_queued; j++) 
		    {
		      if (j==i) continue;
		      if (br_queue[i].mu1idx==br_queue[j].mu1idx && br_queue[i].mu2idx==br_queue[j].mu2idx && br_queue[i].tk1idx==br_queue[j].tk1idx && br_queue[i].tk2idx==br_queue[j].tk2idx)
			{	      
			  if (fabs(br_queue[j].tktkmass-KSTAR_MASS)<fabs(br_queue[i].tktkmass-KSTAR_MASS))
			    {
			      isBestKstarMass = false;
			      continue;
			    }
			}
		    }
		  
		  if (isBestKstarMass)
		    {
		      _nt2->Fill(br_queue[i].mass,br_queue[i].pt,br_queue[i].eta,br_queue[i].y,br_queue[i].mu1pt,br_queue[i].mu1eta,br_queue[i].mu2pt,br_queue[i].mu2eta);
		    }
		}	
	      n_br_queued = 0;
	      memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
	      n_br_queued++;
	    }
	  */

	  //just to test.
	  //_nt2->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	}
      
      if (channel==3)  // cuts for B0 -> J/psi Ks 
	{
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
	  if (br.cosalpha2d<=0.99) continue;
	  
	  if (br.tktkblxy/br.tktkberrxy<=3.0) continue;
	  if (fabs(br.tktkmass-KSHORT_MASS)>=0.015) continue;
	    
	  //_nt3->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	}
	
      if (channel==4)  // cuts for Bs -> J/psi phi
	{
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
	  
	  //_nt4->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta); 
	}
	
      if (channel==5)  // cuts for psi(2S)/X(3872) -> J/psi pipi
	{
	  if (br.vtxprob<=0.2) continue;
	  if (fabs(br.tk1eta)>=1.6) continue;
	  if (fabs(br.tk2eta)>=1.6) continue;
	  
	  //_nt5->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	}
	
      if (channel==6) //cuts for lambda
	{
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
	  
	  //_nt6->Fill(br.mass,br.pt,br.eta,br.y,br.mu1pt,br.mu1eta,br.mu2pt,br.mu2eta);
	}
      tree_out->Fill();
    }//end of the for for the events
  
  fout->Write();
  fout->Close();
}
