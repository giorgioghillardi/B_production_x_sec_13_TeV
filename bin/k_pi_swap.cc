//this sript is used to estimate the ammout of "signal" in the B0->jpsi K*0 channel, in which the k and pi are swapped.
//for this we look at an MC sample processed using myloop_new .cc where we separate between true signal and swapped signal.
//here we fit each of these two categories separatly and extract the relative magnitude and width of each gaussian.
//This is then introduced in the pdf to extract the signal of B0->jpsi K*0, in other sripts in the cross sections studies.

#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

using namespace RooFit;

int main(int argc, char** argv)
{
  //int channel = 2;
  std::string dir ="";

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--output")
        {
          convert << argv[++i];
          convert >> dir;
        }
    }

  double mass_min = 4.0;
  double mass_max = 6.5;
  double pt_min = 0;
  double pt_max = 300;
  double y_min = -3;
  double y_max = 3;

  RooRealVar mass("mass","mass",mass_min,mass_max);
  RooRealVar pt("pt","pt",pt_min,pt_max);
  RooRealVar y("y","y",y_min,y_max);
  
  TString input_file = "selected_mc_ntkstar_with_cuts.root";

  TFile *fin = new TFile(input_file);

  TString ntuple_signal = "ntkstar_true";
  TTree *tree_signal = (TTree*)fin->Get(ntuple_signal);
  RooDataSet* data_signal = new RooDataSet("data_signal","data_signal", tree_signal, RooArgSet(mass,pt,y) );

  TCanvas c1;
  TH1D* mass_signal = (TH1D*)data_signal->createHistogram("mass_signal", mass);
  c1.SetLogy();
  mass_signal->Draw();
  
  TString directory = dir + "mass_signal.png";
  c1.SaveAs(directory);

  //TString ntuple_signal = "ntkstar_swap";

  delete data_signal;

  //ReducedBranches br;  
  //br.setbranchadd(tin);  
}  
  /*
    GenInfoBranches *GenInfo = new GenInfoBranches;
  EvtInfoBranches *EvtInfo = new EvtInfoBranches;
  VtxInfoBranches *VtxInfo = new VtxInfoBranches;
  MuonInfoBranches *MuonInfo = new MuonInfoBranches;
  TrackInfoBranches *TrackInfo = new TrackInfoBranches;
  BInfoBranches *BInfo = new BInfoBranches;

  GenInfo->setbranchadd(tin);
  EvtInfo->setbranchadd(tin);
  VtxInfo->setbranchadd(tin);
  MuonInfo->setbranchadd(tin);
  TrackInfo->setbranchadd(tin);
  BInfo->setbranchadd(tin);

  TString directory = "";
  directory = "k_pi_true_and_swapped_signal.root";

  if(dir != "")
    directory = dir + directory;

  TFile *fout = new TFile(directory,"recreate");  

  TTree *nt_true = new TTree("nttrue","nttrue");
  TTree *nt_swapped = new TTree("ntswapped","ntswapped");

  br.regTree(nt_true);
  br.regTree(nt_swapped);
  
  
  for (int evt=0; evt<tin->GetEntries(); evt++)
    {
      //processing percentage
      if (evt%1000==0 || evt==tin->GetEntries()-1) printf("processing %d/%d (%.2f%%).\n",evt,(int)tin->GetEntries()-1,(double)evt/(double)(tin->GetEntries()-1)*100.);

      tin->GetEntry(evt);

      int ujidx  = -1;
      int tk1idx = -1;
      int tk2idx = -1;
      int mu1idx = -1;
      int mu2idx = -1;

      // Start of BInfo loop
      for (int bidx = 0; bidx < BInfo->size; bidx++)
	{
	  int b_type = BInfo->type[bidx];
	  
	  //the indices to run over the Binfo. These are used to identify the true signal and the swapped signal, when running on MC.
	  ujidx = BInfo->rfuj_index[bidx];
	  tk1idx = BInfo->rftk1_index[bidx];
	  tk2idx = BInfo->rftk2_index[bidx];
	  mu1idx = BInfo->uj_rfmu1_index[ujidx];
	  mu2idx = BInfo->uj_rfmu2_index[ujidx];
	  
	  if (b_type != 4) continue; // skip any non Kstar
	  if (abs(GenInfo->pdgId[MuonInfo->geninfo_index[mu1idx]]) != 13) continue; //skip any mu that was not generated as a mu+-
	  if (abs(GenInfo->pdgId[MuonInfo->geninfo_index[mu2idx]]) != 13) continue; //skip any mu that was not generated as a mu+-
	  if (GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]] != GenInfo->mo1[MuonInfo->geninfo_index[mu2idx]]) continue; //skip if the two muons don't have the same mother particle index
	  if (GenInfo->pdgId[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]] != 443) continue; //skip if the mother of the muons is not jpsi, this is redundant, in principle all come from jpsi
	  if ((GenInfo->pdgId[TrackInfo->geninfo_index[tk1idx]]!=321 || GenInfo->pdgId[TrackInfo->geninfo_index[tk2idx]]!=-211) &&
	      (GenInfo->pdgId[TrackInfo->geninfo_index[tk1idx]]!=-321 || GenInfo->pdgId[TrackInfo->geninfo_index[tk2idx]]!=211)) continue;//skip anything that is not k+pi- or k-pi+            
	  if (GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]] != GenInfo->mo1[TrackInfo->geninfo_index[tk2idx]]) continue; //skip if the two tracks don't have the same mother particle index
	  if (abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]]]) != 313) continue; //skip if the mother of the tracks is not K*0
	  if (GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]] != GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]]]) continue; //skip if the index of the mother of the tracks is not the same as mother of the jpsi
	  if (abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]]) != 511) continue; //skip anything that is not a B0. probably redundant at this point in the decay chain. but it is reasonable to keep it.
	  
	  nt_true->Fill();
	  
	  //just as a test
	  if(br.mu1pt < 25) continue;
	  nt_swapped->Fill();

	}//end of BInfo loop
    }//end of evt loop

  fout->Write();
  fout->Close();

  delete GenInfo;
  delete EvtInfo;
  delete VtxInfo;
  delete MuonInfo;
  delete TrackInfo;
  delete BInfo;
}
  */
