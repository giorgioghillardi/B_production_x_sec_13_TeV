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
//#include "json_50nsMuonPhysV2_processed.h"

//myloop --channel 2 --mc 1 --cuts 1 --dir /some/place
int main(int argc, char** argv)
{
  int channel = 0;
  int run_on_mc= 0;
  int cuts = 1;
  int bfilter = 1;
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

      if(argument == "--cuts")
	{
          convert << argv[++i];
          convert >> cuts;
        }

      if(argument == "--bfilter")
	{
          convert << argv[++i];
          convert >> bfilter;
        }
      
      if(argument == "--dir")
        {
          convert << argv[++i];
          convert >> dir;
        }
    }

  /*
  if(dir == "")
    {
      std::cout << "No directory was specified. Saving in the current directory. To specify use --dir" << std::endl;
      return 0;
    }
  */
    TChain *root = new TChain("demo/root");
    TChain *HltTree = new TChain("hltanalysis/HltTree");

    if(run_on_mc)
      {
	switch(channel)
	  {
	  default:
	  case 1:
	    if(bfilter)
	      {
		//for Bfilter processed with Bfinder_mc.py	    
		root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bu_Bfilter_v1/BuToJpsiKV2_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bu_Bfilter_v1/160812_095709/0000/Bfinder_mc_*.root");
		HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bu_Bfilter_v1/BuToJpsiKV2_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bu_Bfilter_v1/160812_095709/0000/Bfinder_mc_*.root");
	      }
	    else
	      {
		//for BMuonFilter processed with Bfinder_mc	    
		root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_test_v2/BuToJpsiKV2_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_test_v2/160811_210322/0000/Bfinder_mc_*.root");
		HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_test_v2/BuToJpsiKV2_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_test_v2/160811_210322/0000/Bfinder_mc_*.root");
	      }
	    break;

	  case 2:
	    if(bfilter)
	      {
		//for Bfilter processed with Bfinder_mc.py
		root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bd_Bfilter_v1/BdToJpsiKstarV2_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bd_Bfilter_v1/160812_115744/0000/Bfinder_mc_*.root");
		HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bd_Bfilter_v1/BdToJpsiKstarV2_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bd_Bfilter_v1/160812_115744/0000/Bfinder_mc_*.root");
	      }
	    else
	      {
		//for BMuonFilter processed with Bfinder_mc
		root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bd_muonfilter_v1/BdToJpsiKstarV2_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bd_muonfilter_v1/160812_135133/0000/Bfinder_mc_*.root");
		HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bd_muonfilter_v1/BdToJpsiKstarV2_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bd_muonfilter_v1/160812_135133/0000/Bfinder_mc_*.root");
	      }
	    break;
	  case 3:
	    break;

	  case 4:
	    if(bfilter)
	      {
	    //for Bfilter processed with Bfinder_mc.py
	    root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bs_Bfilter_v3/BsToJpsiPhi_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bs_Bfilter_v3/160812_235643/0000/Bfinder_mc_*.root");
	    HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bs_Bfilter_v3/BsToJpsiPhi_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bs_Bfilter_v3/160812_235643/0000/Bfinder_mc_*.root");
	      }
	    else
	      {
		//for BMuonFilter processed with Bfinder_mc.py
		root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bs_muonfilter_v1/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bs_muonfilter_v1/160812_151233/0000/Bfinder_mc_*.root");
		HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_mc_Bs_muonfilter_v1/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/crab_Bfinder_mc_Bs_muonfilter_v1/160812_151233/0000/Bfinder_mc_*.root");
	      }
	    break;

	  case 5:
	    break;
	  case 6:
	    break;
	  }
      }
    else
      { //the data contains all the B's from the different channels. 
	root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v1/Charmonium/Run2015D-Bfinder-promptreco-v1/160309_114238/0000/Bfinder_25ns_*.root");
	root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v3/Charmonium/Run2015D-Bfinder-promptreco-v3/160308_233052/0001/Bfinder_25ns_*.root");                                                                                                                                 
	root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v3/Charmonium/Run2015D-Bfinder-promptreco-v3/160308_233052/0000/Bfinder_25ns_*.root");    
	root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v4/Charmonium/Run2015D-Bfinder-promptreco-v4/160315_105743/0000/Bfinder_25ns_*.root");
	root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v4/Charmonium/Run2015D-Bfinder-promptreco-v4/160315_105743/0001/Bfinder_25ns_*.root");
	root->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v4/Charmonium/Run2015D-Bfinder-promptreco-v4/160315_105743/0002/Bfinder_25ns_*.root");
	
	HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v1/Charmonium/Run2015D-Bfinder-promptreco-v1/160309_114238/0000/Bfinder_25ns_*.root");
	HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v3/Charmonium/Run2015D-Bfinder-promptreco-v3/160308_233052/0001/Bfinder_25ns_*.root");                                                                        
	HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v3/Charmonium/Run2015D-Bfinder-promptreco-v3/160308_233052/0000/Bfinder_25ns_*.root");    
	HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v4/Charmonium/Run2015D-Bfinder-promptreco-v4/160315_105743/0000/Bfinder_25ns_*.root");
	HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v4/Charmonium/Run2015D-Bfinder-promptreco-v4/160315_105743/0001/Bfinder_25ns_*.root");
	HltTree->Add("/gstore/t3cms/store/user/martinsg/Bfinder_25ns_promptreco_v4/Charmonium/Run2015D-Bfinder-promptreco-v4/160315_105743/0002/Bfinder_25ns_*.root");
      }
    
    //-----------------------------------------------------------------
    // some basic integrity checks
    
    int n_entries = root->GetEntries();
    if (n_entries!=HltTree->GetEntries())
      {
	printf("Error: # of entries are different in two main trees.\n");
	return 0;
      }
    printf("Going to process %d entries.\n",n_entries);
    
    //-----------------------------------------------------------------
    // setting memory addresses to the branches
    // create new output trees

    //GenInfoBranches *GenInfo = new GenInfoBranches;    
    EvtInfoBranches *EvtInfo = new EvtInfoBranches;
    VtxInfoBranches *VtxInfo = new VtxInfoBranches;
    MuonInfoBranches *MuonInfo = new MuonInfoBranches;
    TrackInfoBranches *TrackInfo = new TrackInfoBranches;
    BInfoBranches *BInfo = new BInfoBranches;
    
    //GenInfo->setbranchadd(root);
    EvtInfo->setbranchadd(root);
    VtxInfo->setbranchadd(root);
    MuonInfo->setbranchadd(root);
    TrackInfo->setbranchadd(root);
    BInfo->setbranchadd(root);
    
    ULong64_t HltTree_Event;
    int HltTree_Run;
    int HLT_book[N_HLT_BOOKINGS];
    
    HltTree->SetBranchAddress("Event",&HltTree_Event);
    HltTree->SetBranchAddress("Run",&HltTree_Run);
    for (int i=0;i<N_HLT_BOOKINGS;i++)
        HltTree->SetBranchAddress(HLT_paths[i],&HLT_book[i]);
    
    TString directory = "";
    TString filter = "";
    TString bf = "";
    
    if(cuts)
      filter = "with_cuts";
    else
      filter = "no_cuts";
    
    if(bfilter)
      bf = "bfilter";
    else
      bf = "bmuonfilter";
    
    directory = "myloop_new_" + channel_to_ntuple_name(channel) + "_" + bf + "_" + filter + ".root";

    if(dir != "")
      directory = dir + directory;
    
    TFile *fout = new TFile(directory,"recreate");
    
    ReducedBranches brkp;
    ReducedBranches brpi;
    ReducedBranches brks;
    ReducedBranches brkstar;
    ReducedBranches brphi;
    ReducedBranches brmix;
    ReducedBranches brlambda;
    
    TTree *ntkp = new TTree("ntkp","ntkp");
    TTree *ntpi = new TTree("ntpi","ntpi");
    TTree *ntks = new TTree("ntks","ntks");
    TTree *ntkstar = new TTree("ntkstar","ntkstar");
    TTree *ntphi = new TTree("ntphi","ntphi");
    TTree *ntmix = new TTree("ntmix","ntmix");
    TTree *ntlambda = new TTree("ntlambda","ntlambda");
    
    brkp.regTree(ntkp);
    brpi.regTree(ntpi);
    brks.regTree(ntks);
    brkstar.regTree(ntkstar);
    brphi.regTree(ntphi);
    brmix.regTree(ntmix);
    brlambda.regTree(ntlambda);

    for (int evt=0; evt<n_entries; evt++) 
      {
        if (evt%1000==0 || evt==n_entries-1) printf("processing %d/%d (%.2f%%).\n",evt,n_entries-1,(double)evt/(double)(n_entries-1)*100.);
        
        root->GetEntry(evt);
        HltTree->GetEntry(evt);

	/* 
        //-----------------------------------------------------------------
        // Do json file filtering
        bool be_skip = true;
        for (int i=0; i<N_LUMI_PROCESSED; i++) {
            if (EvtInfo->RunNo==lumi_processed[i][0] &&
                EvtInfo->LumiNo>=lumi_processed[i][1] &&
                EvtInfo->LumiNo<=lumi_processed[i][2]) {
                be_skip = false;
                break;
            }
        }
        if (be_skip) continue; 
        */
        //-----------------------------------------------------------------

        // verify the Run # and Event #

        if (EvtInfo->EvtNo!=(int)HltTree_Event || EvtInfo->RunNo!=HltTree_Run) 
	  {
            printf("Error: mismatch of event # and run #.\n");
            return 0;
	  }
	
	// Start of BInfo loop
	for (int bidx = 0; bidx < BInfo->size; bidx++)
	  {   
	    int b_type = BInfo->type[bidx];
	    
	    if(run_on_mc)
	      {		
		switch(channel)
		  {
		  default:
		  case 1:
		    if (b_type != 1) continue; // skip any non K+
		    break;
		  case 2:
		    if (b_type != 4) continue; // skip any non Kstar
		    break;
		  case 3:
		    if (b_type != 3) continue; // skip any non Kshort
		    break;
		  case 4:
		    if (b_type != 6) continue; // skip any non phi
		    break;
		  case 5:
		    if (b_type != 7) continue; // skip any non pipi
		    break;
		  case 6:
		    if (b_type != 8) continue; // skip any non lambda
		    break;
		  }
	      }//end of if(run_on_mc)
	   
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
		br = &brmix; nt = ntmix; break; // low mass pipi
	      case 8:
	      case 9:
		br = &brlambda; nt = ntlambda; break; // lambda
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
            br->nhltbook = N_HLT_BOOKINGS;
            for (int i=0;i<N_HLT_BOOKINGS;i++)
                br->hltbook[i] = HLT_book[i];
            
            int ujidx = BInfo->rfuj_index[bidx];
            int tk1idx = BInfo->rftk1_index[bidx];
            int tk2idx = BInfo->rftk2_index[bidx];
            int mu1idx = BInfo->uj_rfmu1_index[ujidx];
            int mu2idx = BInfo->uj_rfmu2_index[ujidx];
            
	    //the user chooses to preforme the cuts o not. this is useful to calculate efficiencies. this affects both data and MC.
	    if(cuts)
	      {
		//-----------------------------------------------------------------
		// Basic muon selections
		if (MuonInfo->pt[mu1idx]<=4.) continue;
		if (MuonInfo->pt[mu2idx]<=4.) continue;
		if (fabs(MuonInfo->eta[mu1idx])>=2.4) continue;
		if (fabs(MuonInfo->eta[mu2idx])>=2.4) continue;
		if (!MuonInfo->SoftMuID[mu1idx]) continue;
		if (!MuonInfo->SoftMuID[mu2idx]) continue;
		
		//-----------------------------------------------------------------
		// Basic track selections
		if (b_type==1 || b_type==2) { // k, pi
		  if (TrackInfo->pt[tk1idx]<=0.8) continue;
		  if (fabs(TrackInfo->eta[tk1idx])>=2.5) continue;
		  if (TrackInfo->chi2[tk1idx]/TrackInfo->ndf[tk1idx]>=5.) continue;
		  if (TrackInfo->striphit[tk1idx]+TrackInfo->pixelhit[tk1idx]<5) continue;
		}else { // others (2 tracks)
		  if (TrackInfo->pt[tk1idx]<=0.7) continue;
		  if (TrackInfo->pt[tk2idx]<=0.7) continue;
		  if (fabs(TrackInfo->eta[tk1idx])>=2.5) continue;
		  if (fabs(TrackInfo->eta[tk2idx])>=2.5) continue;
		  if (TrackInfo->chi2[tk1idx]/TrackInfo->ndf[tk1idx]>=5.) continue;
		  if (TrackInfo->chi2[tk2idx]/TrackInfo->ndf[tk2idx]>=5.) continue;
		  if (TrackInfo->striphit[tk1idx]+TrackInfo->pixelhit[tk1idx]<5) continue;
		  if (TrackInfo->striphit[tk2idx]+TrackInfo->pixelhit[tk2idx]<5) continue;
		}
		
		//-----------------------------------------------------------------
		// J/psi cut
		// KFC: May need to consider an y dependent cut?
		if (fabs(BInfo->uj_mass[ujidx]-JPSI_MASS)>=0.150) continue;
		if (BInfo->uj_pt[ujidx]<=8.0) continue;
		
		//-----------------------------------------------------------------
		// ditrack selections
		if (b_type==3) // Ks mode
		  {
		    if (fabs(BInfo->tktk_mass[bidx]-KSHORT_MASS)>=0.060) continue;
		  }
		if (b_type==4 || b_type==5) // Kstar mode
		  {
		    if (fabs(BInfo->tktk_mass[bidx]-KSTAR_MASS)>=0.100) continue;
		  }
		if (b_type==6) // phi mode
		  {
		    if (fabs(BInfo->tktk_mass[bidx]-PHI_MASS)>=0.060) continue;
		  }
					
		if (b_type==8 || b_type==9) // Lambda mode
		  {
		    if (fabs(BInfo->tktk_mass[bidx]-LAMBDA_MASS)>=0.060) continue;
		  }
	      } //end of cuts
	    
            //-----------------------------------------------------------------
            // Find the best pointing PV
            //
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
            
            TLorentzVector v4_tktk;
            TVector3 tktkvtx, tktkvtx_err;
            
            if (b_type!=1 && b_type!=2) // other then K+, pi
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
            double default_bmass = BP_MASS; // B+ channels (type 1 or 2)
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
            
            if (b_type!=1 && b_type!=2) // other then K+, pi
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
            
            nt->Fill();
            
	  } // end of BInfo loop
      } // end of evt loop
    
    fout->Write();
    fout->Close();
    
    //delete GenInfo;
    delete EvtInfo;
    delete VtxInfo;
    delete MuonInfo;
    delete TrackInfo;
    delete BInfo;
}
