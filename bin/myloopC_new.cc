#include <iostream>
#include <TChain.h>
#include <sstream>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "UserCode/B_production_x_sec_13_TeV/interface/format.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/myloop.h"
#include "UserCode/B_production_x_sec_13_TeV/interface/channel.h"

//myloop_new --channel 7 --mc 0 --truth 0 --cuts 1 --debug 0 --output prova
int main(int argc, char** argv)
{
  int channel = 1;
  int run_on_mc= 0;
  int mc_truth=0;
  int cuts = 1;
  int debug = 0;
  std::string dir ="";

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

      if(argument == "--truth")
        {
          convert << argv[++i];
          convert >> mc_truth;
        }

      if(argument == "--cuts")
	{
          convert << argv[++i];
          convert >> cuts;
        }
      
      if(argument == "--debug")
	{
          convert << argv[++i];
          convert >> debug;
        }
      
      if(argument == "--output")
        {
          convert << argv[++i];
          convert >> dir;
        }
    }
  
 
  
  TChain *root = new TChain("analysis/root");
  TChain *HltTree = new TChain("hltanalysis/HltTree");
 

    if(run_on_mc)
      {
	switch(channel)
	  {
	  default:
	  case 1:
	    //for BMuonFilter processed with Bfinder_mc	    
	    root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_test_v2/BuToJpsiKV2_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_test_v2/160811_210322/0000/Bfinder_mc_*.root");
	    HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_test_v2/BuToJpsiKV2_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_test_v2/160811_210322/0000/Bfinder_mc_*.root");
	    break;

	  case 2:
	    //for BMuonFilter processed with Bfinder_mc
	    root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bd_muonfilter_v1/BdToJpsiKstarV2_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bd_muonfilter_v1/160812_135133/0000/Bfinder_mc_*.root");
	    HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bd_muonfilter_v1/BdToJpsiKstarV2_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bd_muonfilter_v1/160812_135133/0000/Bfinder_mc_*.root");
	    break;

	  case 3:
	    break;

	  case 4:
	    //for BMuonFilter processed with Bfinder_mc.py
	    root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bs_muonfilter_v1/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bs_muonfilter_v1/160812_151233/0000/Bfinder_mc_*.root");
	    HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bs_muonfilter_v1/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bs_muonfilter_v1/160812_151233/0000/Bfinder_mc_*.root");
	    break;

	  case 5:
	    break;
	  case 6:
	    break;
	  case 7:
	    root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_BcMC2016/BcToJPsiBcPt8Y2p5_MuNoCut_13TeV-bcvegpy2-pythia8/Bfinder_Run2016D-23Sep2016-v1/170818_185252/0000/bfinder_*.root");
	      HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_BcMC2016/BcToJPsiBcPt8Y2p5_MuNoCut_13TeV-bcvegpy2-pythia8/Bfinder_Run2016D-23Sep2016-v1/170818_185252/0000/bfinder_*.root");

	  }
      }
    else
      { //the data contains all the B's from the different channels. 
	//the first is the D dataset and it's in another folder since we try as first with only  this one
	/*root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v1/Charmonium/Run2015D-Bfinder-promptreco-v1/160309_114238/0000/Bfinder_25ns_147.root"); //2015 data test
	//HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v1/Charmonium/Run2015D-Bfinder-promptreco-v1/160309_114238/0000/Bfinder_25ns_147.root");
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_25ns_promptreco_v1/Charmonium/Bfinder_Run2016D-23Sep2016-v1/170811_133948/0000/bfinder_147.root"); //2016 data test

	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_25ns_promptreco_v1/Charmonium/Bfinder_Run2016D-23Sep2016-v1/170811_133948/0000/bfinder_147.root");
       
	*/ 

 
	/* root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_25ns_promptreco_v1/Charmonium/Bfinder_Run2016D-23Sep2016-v1/170811_133948/0000/bfinder_*.root");
       root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_25ns_promptreco_v1/Charmonium/Bfinder_Run2016D-23Sep2016-v1/170811_133948/0001/bfinder_*.root");
       root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_25ns_promptreco_v1/Charmonium/Bfinder_Run2016D-23Sep2016-v1/170811_133948/0002/bfinder_*.root");
       
       
        root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0000/bfinder_*.root"); 
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0001/bfinder_*.root");
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0002/bfinder_*.root");
        root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0003/bfinder_*.root");
        root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0004/bfinder_*.root");
        root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0005/bfinder_*.root");
	
	*/
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/C/Charmonium/Bfinder_Run2016C-23Sep2016-v1/170811_175516/0000/bfinder_*.root");
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/C/Charmonium/Bfinder_Run2016C-23Sep2016-v1/170811_175516/0001/bfinder_*.root"); 
	
	/*
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/E/Charmonium/Bfinder_Run2016E-23Sep2016-v1/170811_175735/0000/bfinder_*.root");
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/E/Charmonium/Bfinder_Run2016E-23Sep2016-v1/170811_175735/0001/bfinder_*.root");
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/E/Charmonium/Bfinder_Run2016E-23Sep2016-v1/170811_175735/0002/bfinder_*.root");
	
	          
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/F/Charmonium/Bfinder_Run2016F-23Sep2016-v1/170811_180319/0000/bfinder_*.root");
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/F/Charmonium/Bfinder_Run2016F-23Sep2016-v1/170811_180319/0001/bfinder_*.root");


	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/G/Charmonium/Bfinder_Run2016G-23Sep2016-v1/170811_180434/0000/bfinder_*.root");
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/G/Charmonium/Bfinder_Run2016G-23Sep2016-v1/170811_180434/0001/bfinder_*.root");
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/G/Charmonium/Bfinder_Run2016G-23Sep2016-v1/170811_180434/0002/bfinder_*.root");
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/G/Charmonium/Bfinder_Run2016G-23Sep2016-v1/170811_180434/0003/bfinder_*.root");
	root->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/G/Charmonium/Bfinder_Run2016G-23Sep2016-v1/170811_180434/0004/bfinder_*.root");





	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_25ns_promptreco_v1/Charmonium/Bfinder_Run2016D-23Sep2016-v1/170811_133948/0000/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_25ns_promptreco_v1/Charmonium/Bfinder_Run2016D-23Sep2016-v1/170811_133948/0001/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_25ns_promptreco_v1/Charmonium/Bfinder_Run2016D-23Sep2016-v1/170811_133948/0002/bfinder_*.root");
	
       
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0000/bfinder_*.root"); 
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0001/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0002/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0003/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0004/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/B/Charmonium/Bfinder_Run2016B-23Sep2016-v3/170811_175216/0005/bfinder_*.root");
	



	*/

	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/C/Charmonium/Bfinder_Run2016C-23Sep2016-v1/170811_175516/0000/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/C/Charmonium/Bfinder_Run2016C-23Sep2016-v1/170811_175516/0001/bfinder_*.root"); 

	/*
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/E/Charmonium/Bfinder_Run2016E-23Sep2016-v1/170811_175735/0000/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/E/Charmonium/Bfinder_Run2016E-23Sep2016-v1/170811_175735/0001/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/E/Charmonium/Bfinder_Run2016E-23Sep2016-v1/170811_175735/0002/bfinder_*.root");
	
	          
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/F/Charmonium/Bfinder_Run2016F-23Sep2016-v1/170811_180319/0000/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/F/Charmonium/Bfinder_Run2016F-23Sep2016-v1/170811_180319/0001/bfinder_*.root");


	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/G/Charmonium/Bfinder_Run2016G-23Sep2016-v1/170811_180434/0000/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/G/Charmonium/Bfinder_Run2016G-23Sep2016-v1/170811_180434/0001/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/G/Charmonium/Bfinder_Run2016G-23Sep2016-v1/170811_180434/0002/bfinder_*.root");
	HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/G/Charmonium/Bfinder_Run2016G-23Sep2016-v1/170811_180434/0003/bfinder_*.root");
        HltTree->Add("/gstore/t3cms/store/user/gghillar/Bfinder_Run2016/G/Charmonium/Bfinder_Run2016G-23Sep2016-v1/170811_180434/0004/bfinder_*.root");  */
	

      }

       

    //-----------------------------------------------------------------
    // some basic integrity checks
    
    int n_entries = root->GetEntries();
    int percent = n_entries/100;

    if (n_entries!=HltTree->GetEntries())
      {
	printf("Error: # of entries are different in two main trees.\n");
	return 0;
      }
    printf("Going to process %d entries.\n",n_entries);
    
    //-----------------------------------------------------------------
    // setting memory addresses to the branches

    GenInfoBranches *GenInfo = new GenInfoBranches;    
    EvtInfoBranches *EvtInfo = new EvtInfoBranches;
    VtxInfoBranches *VtxInfo = new VtxInfoBranches;
    MuonInfoBranches *MuonInfo = new MuonInfoBranches;
    TrackInfoBranches *TrackInfo = new TrackInfoBranches;
    BInfoBranches *BInfo = new BInfoBranches;
    
    GenInfo->setbranchadd(root);
    EvtInfo->setbranchadd(root);
    VtxInfo->setbranchadd(root);
    MuonInfo->setbranchadd(root);
    TrackInfo->setbranchadd(root);
    BInfo->setbranchadd(root);
    
    ULong64_t HltTree_Event;
    int HltTree_Run;
    int hlt_size = N_HLT_BOOKINGS;

    if(run_on_mc)
      hlt_size =  N_HLT_BOOKINGS-1;// to change if there are different HLT
    
    int HLT_book[hlt_size];
    
    HltTree->SetBranchAddress("Event",&HltTree_Event);
    HltTree->SetBranchAddress("Run",&HltTree_Run);
    for (int i=0; i<hlt_size; i++)
        HltTree->SetBranchAddress(HLT_paths[i],&HLT_book[i]);
    
    TString directory = "";
    TString data = "";
    TString filter = "";
    
    if(cuts)
      filter = "with_cuts";
    else
      filter = "no_cuts";
    
    if(run_on_mc)
      {
	if(mc_truth)
	  data = "mc_truth";
	else
	  data = "mc";
      }
    else
      data = "data";
    
    directory = "myloop_new_" + data + "_" + channel_to_ntuple_name(channel) + "_" + filter + ".root";

    if(dir != "")
      directory = dir + directory;
    
    TFile *fout = new TFile(directory,"recreate");
    
    ReducedBranches brkp;
    ReducedBranches brpi;
    ReducedBranches brks;
    ReducedBranches brkstar;
    ReducedBranches brkstar_true;
    ReducedBranches brkstar_swap;
    ReducedBranches brphi;
    ReducedBranches brmix;
    ReducedBranches brlambda;
    
    TTree *ntkp = new TTree("ntkp","ntkp");
    TTree *ntpi = new TTree("ntpi","ntpi");
    TTree *ntks = new TTree("ntks","ntks");
    TTree *ntkstar = new TTree("ntkstar","ntkstar");
    TTree *ntkstar_true = new TTree("ntkstar_true","ntkstar_true");
    TTree *ntkstar_swap = new TTree("ntkstar_swap","ntkstar_swap");
    TTree *ntphi = new TTree("ntphi","ntphi");
    TTree *ntmix = new TTree("ntmix","ntmix");
    TTree *ntlambda = new TTree("ntlambda","ntlambda");
    
    brkp.regTree(ntkp);
    brpi.regTree(ntpi);
    brks.regTree(ntks);
    brkstar.regTree(ntkstar);
    brkstar_true.regTree(ntkstar_true);
    brkstar_swap.regTree(ntkstar_swap);
    brphi.regTree(ntphi);
    brmix.regTree(ntmix);
    brlambda.regTree(ntlambda);

    std::vector<ReducedBranches> selected_bees;

    //debug:
    std::vector<std::string> particle_flow_string = {"event","signal","mu1pt","mu2pt","mu1eta","mu2eta","mu1soft","mu2soft","jpsi_mass","jpsi_pt","tk1pt","tk2pt","tk1eta","tk2eta","tk1chi_ndf","tk2chi_ndf","tk1hits","tk2hits","tktk_mass","tktk_veto","HLT_selection","vtx_prob","lxy_errlxy","cos2D"};
    std::vector<int> particle_flow_number;

    for(int i=0; i < (int)particle_flow_string.size(); i++)
      particle_flow_number.push_back(0);
    
    particle_flow_number[0]= n_entries;
    for (int evt=0; evt< n_entries; evt++)   
    
      {
	if (evt%percent==0 || evt==n_entries-1) printf("processing %d/%d (%.2f%%).\n",evt,n_entries-1,(double)evt/(double)(n_entries-1)*100.);
        
        root->GetEntry(evt);
        HltTree->GetEntry(evt);
	
        // verify the Run # and Event #
        if (EvtInfo->EvtNo!=(int)HltTree_Event || EvtInfo->RunNo!=HltTree_Run) 
	  {
            printf("Error: mismatch of event # and run #.\n");
            return 0;
	  }
	
	int ujidx  = -1;
	int tk1idx = -1;
	int tk2idx = -1;
	int mu1idx = -1;
	int mu2idx = -1;

	//clear the vector that stores the selected bees in each event
	selected_bees.clear();
	
	// Start of BInfo loop
	for (int bidx = 0; bidx < BInfo->size; bidx++)//bidx is the index that will run over the Binfo
	  {
	    int b_type = BInfo->type[bidx];
	    
	    //the indices to run over the Binfo. These are used to identify the signal when running on MC.
	    ujidx = BInfo->rfuj_index[bidx]; 
	    tk1idx = BInfo->rftk1_index[bidx];
	    tk2idx = BInfo->rftk2_index[bidx];
	    mu1idx = BInfo->uj_rfmu1_index[ujidx];
	    mu2idx = BInfo->uj_rfmu2_index[ujidx];
		
	    switch(channel) //to translate the channel to the type. type is defined in Bfinder to destinguish the BParticles candidate
	      {
	      default:
	      case 1:
		if (b_type != 1) continue; // skip any non K+
		//to select the reconstructed Bees that we save. This way we only save signal.
		if (run_on_mc && mc_truth)
		  {
		    if (abs(GenInfo->pdgId[MuonInfo->geninfo_index[mu1idx]]) != 13) continue; //skip any mu that was not generated as a mu+-
		    if (abs(GenInfo->pdgId[MuonInfo->geninfo_index[mu2idx]]) != 13) continue; //skip any mu that was not generated as a mu+-
		    if (GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]] != GenInfo->mo1[MuonInfo->geninfo_index[mu2idx]]) continue; //skip if the two muons don't have the same index for the mother particle
		    if (abs(GenInfo->pdgId[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]) != 443) continue; //skip if the mother of the muons is not jpsi, this is redundant, in principle all come from jpsi
		    if (abs(GenInfo->pdgId[TrackInfo->geninfo_index[tk1idx]]) != 321) continue; //skip any tk that was not generated as a K+-
		    if (GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]] != GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]]) continue; //skip if the index of the mother of the track is not the same as  mother of the jpsi
		    if (abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]]) != 521) continue; //skip anything that is not a B+-. probably redundant at this point in the decay chain. but it is reasonable to keep it.
		  }
		break;

	      case 2:
		if (b_type != 4 && b_type != 5) continue; // skip any non Kstar
		//to select the reconstructed Bees that we save. This way we only save signal.
		if (run_on_mc && mc_truth)
		  {
		    if (abs(GenInfo->pdgId[MuonInfo->geninfo_index[mu1idx]]) != 13) continue; //skip any mu that was not generated as a mu+-
		    if (abs(GenInfo->pdgId[MuonInfo->geninfo_index[mu2idx]]) != 13) continue; //skip any mu that was not generated as a mu+-    
		    if (GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]] != GenInfo->mo1[MuonInfo->geninfo_index[mu2idx]]) continue; //skip if the two muons don't have the same mother particle index
		    if (abs(GenInfo->pdgId[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]) != 443) continue; //skip if the mother of the muons is not jpsi, this is redundant, in principle all come from jpsi	    
		    if ((GenInfo->pdgId[TrackInfo->geninfo_index[tk1idx]]!=321 || GenInfo->pdgId[TrackInfo->geninfo_index[tk2idx]]!=-211) && //not k+pi-
			(GenInfo->pdgId[TrackInfo->geninfo_index[tk1idx]]!=-321 || GenInfo->pdgId[TrackInfo->geninfo_index[tk2idx]]!=211) && //not k-pi+
			(GenInfo->pdgId[TrackInfo->geninfo_index[tk1idx]]!=211 || GenInfo->pdgId[TrackInfo->geninfo_index[tk2idx]]!=-321) && //not pi+k-
			(GenInfo->pdgId[TrackInfo->geninfo_index[tk1idx]]!=-211 || GenInfo->pdgId[TrackInfo->geninfo_index[tk2idx]]!=321)) continue;//not pi-k+
		    if (GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]] != GenInfo->mo1[TrackInfo->geninfo_index[tk2idx]]) continue; //skip if the two tracks don't have the same mother particle index
		    if (abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]]]) != 313) continue; //skip if the mother of the tracks is not K*0
		    if (GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]] != GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]]]) continue; //skip if the index of the mother of the tracks is not the same as mother of the jpsi
		    if (abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]]) != 511) continue; //skip anything that is not a B0. probably redundant at this point in the decay chain. but it is reasonable to keep it.
		  }
		break;
		    
	      case 3:
		if (b_type != 3) continue; // skip any non Kshort
		break;

	      case 4:
		if (b_type != 6) continue; // skip any non phi
		//to select the reconstructed Bees that we save. This way we only save signal.
		if (run_on_mc && mc_truth)
		  {	    
		    if (abs(GenInfo->pdgId[MuonInfo->geninfo_index[mu1idx]]) != 13) continue; //skip any mu that was not generated as a mu+-
		    if (abs(GenInfo->pdgId[MuonInfo->geninfo_index[mu2idx]]) != 13) continue; //skip any mu that was not generated as a mu+-    
		    if (GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]] != GenInfo->mo1[MuonInfo->geninfo_index[mu2idx]]) continue; //skip if the two muons don't have the same mother particle index
		    if (abs(GenInfo->pdgId[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]) != 443) continue; //skip if the mother of the muons is not jpsi, this is redundant, in principle all come from jpsi	    
		    if ((GenInfo->pdgId[TrackInfo->geninfo_index[tk1idx]]!=321 || GenInfo->pdgId[TrackInfo->geninfo_index[tk2idx]]!=-321) &&
			(GenInfo->pdgId[TrackInfo->geninfo_index[tk1idx]]!=-321 || GenInfo->pdgId[TrackInfo->geninfo_index[tk2idx]]!=321)) continue; //skip anything that is not k+k- or k-k+
		    if (GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]] != GenInfo->mo1[TrackInfo->geninfo_index[tk2idx]]) continue; //skip if the two tracks don't have the same mother particle index
		    if (abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]]]) != 333) continue; //skip if the mother of the tracks is not phi
		    if (GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]] != GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]]]) continue; //skip if the index of the mother of the tracks is not the same as mother of the jpsi
		    if (abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]]) != 531) continue; //skip anything that is not a Bs. probably redundant at this point in the decay chain. but it is reasonable to keep it.
		  }
		break;
	      case 5:
		if (b_type != 7) continue; // skip any non pipi
		break;
	      case 6:
		if (b_type != 8 && b_type!=9) continue; // skip any non lambda
		break;
	     
              case 7:
                if (b_type !=2) continue; // skip any non Bc to Jps Pi
		if (run_on_mc && mc_truth)
		  {
		    if (abs(GenInfo->pdgId[MuonInfo->geninfo_index[mu1idx]]) != 13) continue; //skip any mu that was not generated as a mu+-
		    if (abs(GenInfo->pdgId[MuonInfo->geninfo_index[mu2idx]]) != 13) continue; //skip any mu that was not generated as a mu+-
		    if (GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]] != GenInfo->mo1[MuonInfo->geninfo_index[mu2idx]]) continue; //skip if the two muons don't have the same index for the mother particle
		    if (abs(GenInfo->pdgId[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]) != 443) continue; //skip if the mother of the muons is not jpsi, this is redundant, in principle all come from jpsi
		    if (abs(GenInfo->pdgId[TrackInfo->geninfo_index[tk1idx]]) != 211) continue; //skip any tk that was not generated as a Pi+-
		    if (GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]] != GenInfo->mo1[TrackInfo->geninfo_index[tk1idx]]) continue; //skip if the index of the mother of the track is not the same as  mother of the jpsi
		    if (abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]]) != 541) continue; //skip anything that is not a Bc-. probably redundant at this point in the decay chain. but it is reasonable to keep it.
		  }
                break;
	      }
	    
	    if(run_on_mc)
	      particle_flow_number[1]++;
	    
	    // Find the target branching/ntuple to fill
	    ReducedBranches *br = NULL;
	    TTree *nt = NULL;
	    
	    switch (b_type)
	      {
	      case 1:
		br = &brkp; nt = ntkp; break; // K+
	      case 2:
		br = &brpi; nt = ntpi; break; // pi+
	      case 3:
		br = &brks; nt = ntks; break; // Kshort
	      case 4:
	      case 5:
		br = &brkstar; nt = ntkstar; break; // Kstar
	      case 6:
		br = &brphi; nt = ntphi; break; // phi
	      case 7:
		br  = &brmix; nt = ntmix; break; // low mass pipi
	      case 8:
	      case 9:
		br  = &brlambda; nt = ntlambda; break; // lambda
	      default:
		printf("Error: unknown BInfo->type.\n");
		return 0;
	      }
	    
	    // cleanup
	    memset(br,0x00,sizeof(ReducedBranches));
            
	    br->run = EvtInfo->RunNo;
	    br->event = EvtInfo->EvtNo;
	    br->type = b_type;
            
	    //-----------------------------------------------------------------
	    // save HLT paths
	    br->nhltbook = hlt_size;
	    for (int i=0;i<hlt_size;i++)
	      br->hltbook[i] = HLT_book[i];
              
	    //the user chooses to preforme the cuts or not. this is useful to calculate efficiencies. this affects both data and MC.
	    if(cuts)
	      {
		// Basic muon selections
		if (MuonInfo->pt[mu1idx]<=4) continue;
		if(run_on_mc)
		  particle_flow_number[2]++;

		if (MuonInfo->pt[mu2idx]<=4) continue;
		if(run_on_mc)
		  particle_flow_number[3]++;
		
		if (fabs(MuonInfo->eta[mu1idx])>=2.4) continue;
		if(run_on_mc)
		  particle_flow_number[4]++;
		
		if (fabs(MuonInfo->eta[mu2idx])>=2.4) continue;
		if(run_on_mc)
		  particle_flow_number[5]++;

		if (!MuonInfo->SoftMuID[mu1idx]) continue;
		if(run_on_mc)
		  particle_flow_number[6]++;

		if (!MuonInfo->SoftMuID[mu2idx]) continue;
		if(run_on_mc)
		  particle_flow_number[7]++;
		
		//-----------------------------------------------------------------
		// J/psi cut
		// KFC: May need to consider an y dependent cut?
		//add the jpsi vertex prob!?
		if (fabs(BInfo->uj_mass[ujidx]-JPSI_MASS)>=0.150) continue;
		if(run_on_mc)
		  particle_flow_number[8]++;
		
		if (BInfo->uj_pt[ujidx]<=10.0) continue; //was 8.0 before
		if(run_on_mc)
		  particle_flow_number[9]++;
		
		//-----------------------------------------------------------------
		// Basic track selections
		if (b_type==1 || b_type==2) // k, pi
		  {
		    //double tptcut = b_type==2?2.3:1.6; //CHECK!!!
		    //if (TrackInfo->pt[tk1idx]<= tptcut) continue; //was 0.8
		    if (TrackInfo->pt[tk1idx]<=2.3) continue; //was 0.8
		    if (fabs(TrackInfo->eta[tk1idx])>=2.5) continue;
		    if (TrackInfo->chi2[tk1idx]/TrackInfo->ndf[tk1idx]>=5.) continue;
		    if (TrackInfo->striphit[tk1idx]+TrackInfo->pixelhit[tk1idx]<5) continue;
		  }
		else
		  { // others (2 tracks)
		    if (TrackInfo->pt[tk1idx]<=0.7) continue;
		    if(run_on_mc)
		      particle_flow_number[10]++;
		    
		    if (TrackInfo->pt[tk2idx]<=0.7) continue;
		    if(run_on_mc)
		      particle_flow_number[11]++;

		    if (fabs(TrackInfo->eta[tk1idx])>=2.5) continue;
		    if(run_on_mc)
		      particle_flow_number[12]++;
		    
		    if (fabs(TrackInfo->eta[tk2idx])>=2.5) continue;
		    if(run_on_mc)
		      particle_flow_number[13]++;

		    if (TrackInfo->chi2[tk1idx]/TrackInfo->ndf[tk1idx]>=5.) continue;
		    if(run_on_mc)
		      particle_flow_number[14]++;
		    
		    if (TrackInfo->chi2[tk2idx]/TrackInfo->ndf[tk2idx]>=5.) continue;
		    if(run_on_mc)
		      particle_flow_number[15]++;
		    
		    if (TrackInfo->striphit[tk1idx]+TrackInfo->pixelhit[tk1idx]<5) continue;
		    if(run_on_mc)
		      particle_flow_number[16]++;

		    if (TrackInfo->striphit[tk2idx]+TrackInfo->pixelhit[tk2idx]<5) continue;
		    if(run_on_mc)
		      particle_flow_number[17]++;
		  }

		//---------------------------------------------------------------------
		// ditrack mass window selection
		double k_short_window = 0.015; //originally was 0.060
		double lambda_window = 0.015;  //originally was 0.010

		//double k_star_window = 0.200;  //try 0.060 //originally was 0.050
		//double k_star_veto = 0.050;

		//double phi_window = 0.015;    //try 0.010 //originally was 0.060
		//double phi_veto = 0.010;

		switch(b_type)
		  {
		  case 3: // Ks mode
		   if (fabs(BInfo->tktk_mass[bidx]-KSHORT_MASS)>=k_short_window) continue; 
		   break;

		  case 4: // Kstar mode
		  case 5: // Kstar mode
		    //if (fabs(BInfo->tktk_mass[bidx]-KSTAR_MASS)>=k_star_window) continue;
		    break;
		    
		  case 6: // phi mode
		    //if (fabs(BInfo->tktk_mass[bidx]-PHI_MASS)>=phi_window) continue;
		    break;
		    
		  case 8: // Lambda mode
		  case 9: // Lambda mode
		    if (fabs(BInfo->tktk_mass[bidx]-LAMBDA_MASS)>=lambda_window) continue; 
		    break;
		  }

		if(run_on_mc)
		  particle_flow_number[18]++;
		
		//------------------------------------------------------------------------------------
		// ditrack vetos
		TLorentzVector v4_tk1, v4_tk2;
		
		switch(b_type)
		  {
		  case 4: // Kstar mode
		  case 5: // Kstar mode
		    v4_tk1.SetPtEtaPhiM(TrackInfo->pt[tk1idx],TrackInfo->eta[tk1idx],TrackInfo->phi[tk1idx],KAON_MASS);
		    v4_tk2.SetPtEtaPhiM(TrackInfo->pt[tk2idx],TrackInfo->eta[tk2idx],TrackInfo->phi[tk2idx],KAON_MASS);
		    //if (fabs((v4_tk1+v4_tk2).Mag()-PHI_MASS)<=phi_veto) continue;
		    break;
		    
		  case 6: // phi mode
		    v4_tk1.SetPtEtaPhiM(TrackInfo->pt[tk1idx],TrackInfo->eta[tk1idx],TrackInfo->phi[tk1idx],KAON_MASS);
		    v4_tk2.SetPtEtaPhiM(TrackInfo->pt[tk2idx],TrackInfo->eta[tk2idx],TrackInfo->phi[tk2idx],PION_MASS);
		    //if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=k_star_veto) continue;
		    
		    v4_tk1.SetPtEtaPhiM(TrackInfo->pt[tk1idx],TrackInfo->eta[tk1idx],TrackInfo->phi[tk1idx],PION_MASS);
		    v4_tk2.SetPtEtaPhiM(TrackInfo->pt[tk2idx],TrackInfo->eta[tk2idx],TrackInfo->phi[tk2idx],KAON_MASS);
		    //if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=k_star_veto) continue;
		    break;
		    
		  case 8: // Lambda mode
		  case 9: // Lambda mode
		    v4_tk1.SetPtEtaPhiM(TrackInfo->pt[tk1idx],TrackInfo->eta[tk1idx],TrackInfo->phi[tk1idx],PION_MASS);
		    v4_tk2.SetPtEtaPhiM(TrackInfo->pt[tk2idx],TrackInfo->eta[tk2idx],TrackInfo->phi[tk2idx],PION_MASS);
		    if (fabs((v4_tk1+v4_tk2).Mag()-KSHORT_MASS)<=k_short_window) continue;
		    break;
		  }
		
		if(run_on_mc)
		  particle_flow_number[19]++;

	      } //end of cuts, there are more cuts later in the code
	    
            //--------------------------------------------------------------------------------------
            // Find the best pointing PV
            // KFC@20150713: keep the selecton code but PV is replaced with BS in the end.
            
	    TVector3 bvtx(BInfo->vtxX[bidx],BInfo->vtxY[bidx],BInfo->vtxZ[bidx]);
	    TVector3 bvtx_err(BInfo->vtxXE[bidx],BInfo->vtxYE[bidx],BInfo->vtxZE[bidx]);
	    TVector3 bmom(BInfo->px[bidx],BInfo->py[bidx],BInfo->pz[bidx]);
	    int vidx = -1;
	    double max_cosang = -1.;
            
	    for (int idx = 0; idx < VtxInfo->Size; idx++)
	      {
		TVector3 vtx(VtxInfo->x[idx],VtxInfo->y[idx],VtxInfo->z[idx]);    
		double cosang = bmom.Dot(bvtx-vtx)/(bmom.Mag()*(bvtx-vtx).Mag());
		
		if (cosang>max_cosang)
		  {
		    vidx = idx;
		    max_cosang = cosang;
		  }
	      }
	    
	    if (vidx==-1)
	      {
		printf("Error: no PV found. Run: %d, Event: %d.\n",EvtInfo->RunNo,EvtInfo->EvtNo);
		continue;
	      }
            
	    TVector3 BS(EvtInfo->PVx,EvtInfo->PVy,EvtInfo->PVz);
	    TVector3 BS_err(EvtInfo->PVxE,EvtInfo->PVyE,EvtInfo->PVzE);
	    TVector3 PV(VtxInfo->x[vidx],VtxInfo->y[vidx],VtxInfo->z[vidx]);
	    TVector3 PV_err(VtxInfo->xE[vidx],VtxInfo->yE[vidx],VtxInfo->zE[vidx]);
            
	    // KFC@20150713: Don't do PV finding for now. Always use the beamspot according to the agreement
	    PV = BS; PV_err = BS_err;
            
	    TLorentzVector v4_uj;
	    v4_uj.SetPtEtaPhiM(BInfo->uj_pt[ujidx],BInfo->uj_eta[ujidx],BInfo->uj_phi[ujidx],BInfo->uj_mass[ujidx]);
	    TLorentzVector v4_b;
	    v4_b.SetPtEtaPhiM(BInfo->pt[bidx],BInfo->eta[bidx],BInfo->phi[bidx],BInfo->mass[bidx]);
	    
	    //fill the br->gen info
	    if(run_on_mc && mc_truth)
	       {
		 TLorentzVector v4_b_gen;
		 
		 br->genmass = GenInfo->mass[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]];
		 br->genpt   = GenInfo->pt[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]];
		 br->geneta  = GenInfo->eta[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]];
		 br->genphi  = GenInfo->phi[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]];
		 
		 v4_b_gen.SetPtEtaPhiM(br->genpt,br->geneta,br->genphi,br->genmass);
		 br->geny    = v4_b_gen.Rapidity();
		 
		 br->genvx   = GenInfo->vx[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]];
		 br->genvy   = GenInfo->vy[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]];
		 br->genvz   = GenInfo->vz[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]]];
		 
		 br->genujvx = GenInfo->vx[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]];
		 br->genujvy = GenInfo->vy[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]];
		 br->genujvz = GenInfo->vz[GenInfo->mo1[MuonInfo->geninfo_index[mu1idx]]];
	       }
	    
	    TLorentzVector v4_tktk;
	    TVector3 tktkvtx, tktkvtx_err;
            
	    if (b_type!=1 && b_type!=2) // only for the ditrack channels
	      {  
		v4_tktk.SetPtEtaPhiM(BInfo->tktk_pt[bidx],BInfo->tktk_eta[bidx],BInfo->tktk_phi[bidx],BInfo->tktk_mass[bidx]);
		tktkvtx.SetXYZ(BInfo->tktk_vtxX[bidx],BInfo->tktk_vtxY[bidx],BInfo->tktk_vtxZ[bidx]);
		tktkvtx_err.SetXYZ(BInfo->tktk_vtxXE[bidx],BInfo->tktk_vtxYE[bidx],BInfo->tktk_vtxZE[bidx]);
	      }
	    
	    //-----------------------------------------------------------------
	    // Start to fill the B hadron information            
	    br->mass = BInfo->mass[bidx];
	    br->pt = BInfo->pt[bidx];
	    br->eta = BInfo->eta[bidx];
	    br->phi = BInfo->phi[bidx];
	    br->y = v4_b.Rapidity();
	    br->vx = BInfo->vtxX[bidx];
	    br->vy = BInfo->vtxY[bidx];
	    br->vz = BInfo->vtxZ[bidx];
	    br->PVx = PV.x();
	    br->PVy = PV.y();
	    br->PVz = PV.z();
	    br->lxy = (bvtx-PV).Perp();
	    br->lxyz = (bvtx-PV).Mag();
	    br->errxy = sqrt(bvtx_err.Perp2()+PV_err.Perp2());
	    br->errxyz = sqrt(bvtx_err.Mag2()+PV_err.Mag2());
	    br->vtxprob = TMath::Prob(BInfo->vtxchi2[bidx],BInfo->vtxdof[bidx]);
	    br->cosalpha2d = bmom.XYvector()*(bvtx-PV).XYvector()/(bmom.Perp()*(bvtx-PV).Perp());
	    br->cosalpha3d = bmom.Dot(bvtx-PV)/(bmom.Mag()*(bvtx-PV).Mag());
            
	    //-----------------------------------------------------------------
	    // calculate the proper decay time
	    TVector3 v_l = bvtx-PV, v_lerr2; // displace vector, error^2 vector for displacement
	    v_lerr2.SetX(bvtx_err.x()*bvtx_err.x()+PV_err.x()*PV_err.x());
	    v_lerr2.SetY(bvtx_err.y()*bvtx_err.y()+PV_err.y()*PV_err.y());
	    v_lerr2.SetZ(bvtx_err.z()*bvtx_err.z()+PV_err.z()*PV_err.z());
            
	    // B hadron mass for normalization of proper decay time
	    double default_bmass = BP_MASS; // B+ channels (type 1)
            if (b_type==2) default_bmass = BC_MASS; // Bc channels
	    if (b_type==3 || b_type==4 || b_type==5) default_bmass = B0_MASS; // B0 channels
	    if (b_type==6 || b_type==7) default_bmass = BS_MASS; // Bs channels
	    if (b_type==8 || b_type==9) default_bmass = LAMBDAB_MASS; // Lambdab channels
            
	    
            TVector3 v_p = v4_b.Vect();
	    TVector3 v_p2(v_p.x()*v_p.x(),v_p.y()*v_p.y(),v_p.z()*v_p.z());

	    br->ctau3d = v_l.Dot(v_p)*default_bmass/v_p.Mag2();
	    br->ctau3derr = sqrt(v_lerr2.Dot(v_p2))*default_bmass/v_p.Mag2();
	    br->ctau2d = v_l.XYvector()*v_p.XYvector()*default_bmass/v_p.Perp2();
	    br->ctau2derr = sqrt(v_lerr2.XYvector()*v_p2.XYvector())*default_bmass/v_p.Perp2();

	    //-----------------------------------------------------------------
	    // fill J/psi, tracks, muons, etc.
            
	    br->ujmass = BInfo->uj_mass[ujidx];
	    br->ujpt = BInfo->uj_pt[ujidx];
	    br->ujeta = BInfo->uj_eta[ujidx];
	    br->ujphi = BInfo->uj_phi[ujidx];
	    br->ujy = v4_uj.Rapidity();
	    br->ujvtxprob = TMath::Prob(BInfo->uj_vtxchi2[ujidx],BInfo->uj_vtxdof[ujidx]);
            
	    if (b_type!=1 && b_type!=2) // only for ditrack channels
	      {  
		br->tktkmass = BInfo->tktk_mass[bidx];
		br->tktkpt = BInfo->tktk_pt[bidx];
		br->tktketa = BInfo->tktk_eta[bidx];
		br->tktkphi = BInfo->tktk_phi[bidx];
		br->tktky = v4_tktk.Rapidity();
		br->tktkvtxprob = TMath::Prob(BInfo->tktk_vtxchi2[bidx],BInfo->tktk_vtxdof[bidx]);
		
		br->tktklxy = (tktkvtx-PV).Perp();
		br->tktklxyz = (tktkvtx-PV).Mag();
		br->tktkerrxy = sqrt(tktkvtx_err.Perp2()+PV_err.Perp2());
		br->tktkerrxyz = sqrt(tktkvtx_err.Mag2()+PV_err.Mag2());
		br->tktkblxy = (tktkvtx-bvtx).Perp();
		br->tktkblxyz = (tktkvtx-bvtx).Mag();
		br->tktkberrxy = sqrt(tktkvtx_err.Perp2()+bvtx_err.Perp2());
		br->tktkberrxyz = sqrt(tktkvtx_err.Mag2()+bvtx_err.Mag2());
	      }
	    
	    br->mu1idx = mu1idx;
	    br->mu1pt  = MuonInfo->pt[mu1idx];
	    br->mu1eta = MuonInfo->eta[mu1idx];
	    br->mu1phi = MuonInfo->phi[mu1idx];
	    br->mu2idx = mu2idx;
	    br->mu2pt  = MuonInfo->pt[mu2idx];
	    br->mu2eta = MuonInfo->eta[mu2idx];
	    br->mu2phi = MuonInfo->phi[mu2idx];
            
	    br->tk1idx = tk1idx;
	    br->tk1pt  = TrackInfo->pt[tk1idx];
	    br->tk1eta = TrackInfo->eta[tk1idx];
	    br->tk1phi = TrackInfo->phi[tk1idx];
	    br->tk1charge = TrackInfo->charge[tk1idx];
	    br->tk2idx = tk2idx;
	    br->tk2pt  = TrackInfo->pt[tk2idx];
	    br->tk2eta = TrackInfo->eta[tk2idx];
	    br->tk2phi = TrackInfo->phi[tk2idx];
	    br->tk2charge = TrackInfo->charge[tk2idx];
            
	    br->nhltmatch = N_HLT_MATCHINGS;

	    for (int i=0; i<N_HLT_MATCHINGS; i++)
	      {
		br->mu1hltpt[i]  = MuonInfo->MuTrgMatchTrgObjPt->at(mu1idx)[i];
		br->mu1hlteta[i] = MuonInfo->MuTrgMatchTrgObjEta->at(mu1idx)[i];
		br->mu1hltphi[i] = MuonInfo->MuTrgMatchTrgObjPhi->at(mu1idx)[i];
		br->mu2hltpt[i]  = MuonInfo->MuTrgMatchTrgObjPt->at(mu2idx)[i];
		br->mu2hlteta[i] = MuonInfo->MuTrgMatchTrgObjEta->at(mu2idx)[i];
		br->mu2hltphi[i] = MuonInfo->MuTrgMatchTrgObjPhi->at(mu2idx)[i];
	      }
            
	    if(cuts)
	      {
		// cuts that depend on complex variables defined above.
		
		//HLT selection
		//----------------------------------------------------------------
			if(b_type != 7) //not to use HLT filter in the jpsi pipi channel
		  {
		     if(run_on_mc)
		       {if (br->hltbook[ HLT_DoubleMu4_JpsiTrk_Displaced_v1]&&br->hltbook[ HLT_DoubleMu4_JpsiTrk_Displaced_v2]&&br->hltbook[ HLT_DoubleMu4_JpsiTrk_Displaced_v3]&&br->hltbook[ HLT_DoubleMu4_JpsiTrk_Displaced_v4]&&br->hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v5]&&br->hltbook[ HLT_DoubleMu4_JpsiTrk_Displaced_v6]) continue;
		      //is looking for all HLT version that we have on MC
		
                	
                        }
	             else
		       {if (br->hltbook[ HLT_DoubleMu4_JpsiTrk_Displaced_v1]&&br->hltbook[ HLT_DoubleMu4_JpsiTrk_Displaced_v2]&&br->hltbook[ HLT_DoubleMu4_JpsiTrk_Displaced_v3]&&br->hltbook[ HLT_DoubleMu4_JpsiTrk_Displaced_v4]&&br->hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v5]&&br->hltbook[ HLT_DoubleMu4_JpsiTrk_Displaced_v6]) continue;}
		     //is looking for all HLT version that we have on DATA
												

                	if(run_on_mc)
												
                 	particle_flow_number[20]++;
		      }
				
		if(b_type==1 || b_type==2 || b_type==4 || b_type==5 || b_type==6) //for K+, pi+, K*0, phi
		  {
		    if(br->vtxprob<=0.1) continue; //original cut 0.2
		    if(run_on_mc)
		      particle_flow_number[21]++;

		    if(br->lxy/br->errxy<=3.0) continue; //original cut 4.5
		    if(run_on_mc)
		      particle_flow_number[22]++;

		    if(br->cosalpha2d<=0.99) continue; //original cut 0.996
		    if(run_on_mc)
		      particle_flow_number[23]++;
		  }
		
		if(b_type==3 || b_type==8 || b_type==9) // Ks and lambda
		  {
		    if (br->vtxprob<=0.1) continue;
		    if (br->lxy/br->errxy<=3.0) continue;
		    if (br->cosalpha2d<=0.99) continue;
		    if (br->tktkblxy/br->tktkberrxy<=3.0) continue;
		  }
		
		if(b_type==7) // pipi
		  {
		    if (br->vtxprob<=0.2) continue;
		    if (fabs(br->tk1eta)>=1.6) continue;
		    if (fabs(br->tk2eta)>=1.6) continue;
		  }
		//-----------------------------------------------------------------
	      }//end of cuts
	    
	    if(b_type == 4 || b_type == 5)
	      {
		selected_bees.push_back(*br); //copy the selected bees of channel 2 into the selected_bees vector
	      }
	    else
	      nt->Fill();//fill the other channels the usual way
	    
	  } // end of BInfo loop
	
	//===========================================================================================================================
	//cicle over the selected_bees vector, and choose the best B0 candidate. This is for channel 2, to take care of the K pi swap
		
	//set the reducedbranches and ttree pointer to the right addresses
	ReducedBranches *br_cand = NULL;
	TTree *nt_cand = NULL;
	br_cand = &brkstar;
	nt_cand = ntkstar;
	
	//compare the i-th candidade with the j-th candidate: if the j is not closer to the KSTAR_MASS, keep the i as the best candidate.
	for(std::vector<int>::size_type i = 0; i != selected_bees.size(); i++)
	  {
	    if(selected_bees[i].type != 4 && selected_bees[i].type != 5)
	      printf("ERROR: the vector with the B0 candidates contains other channels!! \n");  
	    
	    bool isbestkstar = true; //start by assuming the i-th is the best mass.
	    
	    for(std::vector<int>::size_type j = 0; j != selected_bees.size(); j++)
	      {
		if(j==i) continue; //only compare with different candidates.
		if (selected_bees[i].mu1idx==selected_bees[j].mu1idx && selected_bees[i].mu2idx==selected_bees[j].mu2idx && selected_bees[i].tk1idx==selected_bees[j].tk1idx && selected_bees[i].tk2idx==selected_bees[j].tk2idx) //if a candidate is built with the same muons and tracks
		  {
		    if (fabs(selected_bees[j].tktkmass-KSTAR_MASS)<fabs(selected_bees[i].tktkmass-KSTAR_MASS)) //if the j-th is closer, then set isbestkstar to false, because it refers to the i-th candidate.
		      {
			isbestkstar = false;
			continue;
		      }
		  }
	      }
	    
	    if(isbestkstar) //isbestkstar is true if the i-th is closer than the j-th candidate, or if there is only one candidate.
	      {
		*br_cand = selected_bees[i]; //put the values of the selected B into br_cand
		nt_cand->Fill();
		
		if(run_on_mc && channel==2)//this runs for --truth 0 or --truth 1.
		  {
		    if(debug) printf("debug: evt %d \n", evt);
		    
		    if((selected_bees[i].type == 4 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]]==321 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]==-211) ||
		       (selected_bees[i].type == 4 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]]==-321 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]==211) ||
		       (selected_bees[i].type == 5 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]]==211 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]==-321) ||
		       (selected_bees[i].type == 5 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]]==-211 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]==321)) //true signal
		      {
			br_cand = &brkstar_true;
			nt_cand = ntkstar_true;
			*br_cand = selected_bees[i]; //put the values of the selected B into br_cand
			nt_cand->Fill();

			if(debug) printf("debug signal: type: %d B_pdgId: %d tk1_pdgId: %d tk2_pdgId: %d \n", selected_bees[i].type, GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[selected_bees[i].mu1idx]]]], GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]], GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]);
		      }
		    else 
		      if((selected_bees[i].type == 4 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]]==211 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]==-321) ||
			 (selected_bees[i].type == 4 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]]==-211 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]==321) ||
			 (selected_bees[i].type == 5 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]]==321 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]==-211) ||
			 (selected_bees[i].type == 5 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]]==-321 && GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]==211)) //swapped signal
			{
			  br_cand = &brkstar_swap;
			  nt_cand = ntkstar_swap;
			  *br_cand = selected_bees[i]; //put the values of the selected B into br_cand
			  nt_cand->Fill();
			  
			  if(debug) printf("debug swapped: type: %d B_pdgId: %d tk1_pdgId: %d tk2_pdgId: %d \n", selected_bees[i].type, GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[selected_bees[i].mu1idx]]]], GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]], GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]);
			}
		      else //some of the reconstructed muons and tracks are not associated with any gen particle.
			if(debug) printf("debug unidentified: type: %d B_pdgId: %d tk1_pdgId: %d tk2_pdgId: %d \n", selected_bees[i].type, GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[MuonInfo->geninfo_index[selected_bees[i].mu1idx]]]], GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk1idx]], GenInfo->pdgId[TrackInfo->geninfo_index[selected_bees[i].tk2idx]]);
		  }
	      }//end of is best kstar
	  }//end of the selected_bees loop	
      } // end of evt loop
    
    if(run_on_mc)
      {
	std::cout << "cut  :  number  :  efficiency" << std::endl;
	std::cout << particle_flow_string[0] << " : " << particle_flow_number[0] << " : " << ((double)particle_flow_number[0]/(double)particle_flow_number[0])*100 << std::endl;

	for(int i=1; i < (int)particle_flow_number.size(); i++)
	  {
	    std::cout << particle_flow_string[i] << " : " << particle_flow_number[i] << " : " << ((double)particle_flow_number[i]/(double)particle_flow_number[i-1])*100 << std::endl;
	  }
      }
    
    fout->Write();
    fout->Close();
    
    delete GenInfo;
    delete EvtInfo;
    delete VtxInfo;
    delete MuonInfo;
    delete TrackInfo;
    delete BInfo;
}
