#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda
// channel = 7: Bc -> J/psi Pi
//-----------------------------------------------------------------

//input example: calculate_bin_eff --channel 1 --eff preeff --ptmin 10 --ptmax 20 --ymin 0.00 --ymax 0.50
int main(int argc, char** argv)
{
  int channel = 1;
  TString eff_name = "";
  double pt_min = -1;
  double pt_max = -1;
  double y_min = -1;
  double y_max = -1;

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

      if(argument == "--channel")
	{
	  convert << argv[++i];
	  convert >> channel;
	}
      if(argument == "--eff")
	{
	  convert << argv[++i];
	  convert >> eff_name;
	}      
      if(argument == "--ptmin")
	{
	  convert << argv[++i];
	  convert >> pt_min;
	}
      if(argument == "--ptmax")
	{
	  convert << argv[++i];
	  convert >> pt_max;
	}
      if(argument == "--ymin")
	{
	  convert << argv[++i];
	  convert >> y_min;
	}
      if(argument == "--ymax")
	{
	  convert << argv[++i];
	  convert >> y_max;
	}
    }
  
  if(eff_name != "preeff" && eff_name != "recoeff" && eff_name != "totaleff")
    {
      std::cout << "Error: Please indicate the efficiency (preeff or recoeff or totaleff)" << std::endl;
      return 0;
    }
  
  if(pt_min == -1 || pt_max == -1 ||y_min == -1 ||y_max == -1)
    {
      std::cout << "Error: The bin was not well defined. Please enter pt and y bin." << std::endl;
      return 0;
    }
  
  //to create the directories to save the files
  std::vector<std::string> dir_list;  
  
  dir_list.push_back(static_cast<const char*>(TString::Format(VERSION) + "/efficiencies"));
  dir_list.push_back(static_cast<const char*>(TString::Format(VERSION) + "/efficiencies_root/" + channel_to_ntuple_name(channel) + "/"));
  create_dir(dir_list);

  RooRealVar* eff_res = new RooRealVar("eff_res","eff_res",1);
  eff_res->setError(0);
  
  ///////////////////////////////////////////////////
  //calculate the efficiency for a bin of pt and y.//
  ///////////////////////////////////////////////////
  std::cout << "processing subsample: " << pt_min << " < " << "pt" << " < " << pt_max << " and " << y_min << " < " << "|y|" << " < " << y_max << std::endl;
  
  if(eff_name == "preeff")
    eff_res = prefilter_efficiency(channel, pt_min, pt_max, y_min, y_max);
  else
    if(eff_name == "recoeff")
      eff_res = reco_efficiency(channel, pt_min, pt_max, y_min, y_max);
    else
      if(eff_name == "totaleff")
	{
	  //read pre-filter eff values
	  TString eff_dir = TString::Format(VERSION) + "/efficiencies_root/" + channel_to_ntuple_name(channel) + "/preeff_" + channel_to_ntuple_name(channel) + "_pt_from_" + TString::Format("%d_to_%d", (int)pt_min, (int)pt_max) + "_y_from_" + TString::Format("%.2f_to_%.2f", y_min, y_max) + ".root";
	  
	  TFile* f = new TFile(eff_dir);

	  ////////////////////////////////////////////////////
	  TString line = "";
	  TString eff_str = "";
	  TString opt = "--channel " + TString::Format("%d", channel) + " --ptmin " + TString::Format("%d",(int)pt_min) + " --ptmax " + TString::Format("%d", (int)pt_max) + " --ymin " +  TString::Format("%.2f", y_min) + " --ymax " + TString::Format("%.2f", y_max);

	  if(f->IsZombie())
	    {
	      line = "calculate_bin_eff --eff";
	      eff_str = "preeff";
	      line += " " + eff_str + " " + opt;
	      gSystem->Exec(line);
	      f = new TFile(eff_dir);
	    }

	  TVectorD *pre_eff_val = (TVectorD*)f->Get("val");
	  TVectorD *pre_stat_lo = (TVectorD*)f->Get("err_lo");
	  TVectorD *pre_stat_hi = (TVectorD*)f->Get("err_hi");

	  //read reco eff values
	  eff_dir = TString::Format(VERSION) + "/efficiencies_root/" + channel_to_ntuple_name(channel) + "/recoeff_" + channel_to_ntuple_name(channel) + "_pt_from_" + TString::Format("%d_to_%d", (int)pt_min, (int)pt_max) + "_y_from_" + TString::Format("%.2f_to_%.2f", y_min, y_max) + ".root";
      
	  TFile* f2 = new TFile(eff_dir);

	  /////////////////////////////////////////////////////////
	  if(f2->IsZombie())
	    {
	      line = "calculate_bin_eff --eff";
	      eff_str = "recoeff";
	      line += " " + eff_str + " " + opt;
	      gSystem->Exec(line);
	      f2 = new TFile(eff_dir);
	    }
	  
	  TVectorD *reco_eff_val = (TVectorD*)f2->Get("val");
	  TVectorD *reco_stat_lo = (TVectorD*)f2->Get("err_lo");
	  TVectorD *reco_stat_hi = (TVectorD*)f2->Get("err_hi");
      
	  //put the total efficiency values in eff_res
      
	  eff_res->setVal(pre_eff_val[0][0] * reco_eff_val[0][0]);
      
	  double err_lo = -eff_res->getVal() * sqrt( pow(pre_stat_lo[0][0]/pre_eff_val[0][0],2) + pow(reco_stat_lo[0][0]/reco_eff_val[0][0],2) );
	  double err_hi = eff_res->getVal() * sqrt( pow(pre_stat_hi[0][0]/pre_eff_val[0][0],2) + pow(reco_stat_hi[0][0]/reco_eff_val[0][0],2) );
 
	  eff_res->setAsymError(err_lo , err_hi);

	  delete f;
	  delete f2;
	}

  //write efficiency and statistical error to file
  TString eff_file_name = TString::Format(VERSION) + "/efficiencies_root/" + channel_to_ntuple_name(channel) + "/" + eff_name + "_" + channel_to_ntuple_name(channel) + "_pt_from_" + TString::Format("%d_to_%d", (int)pt_min, (int)pt_max) + "_y_from_" + TString::Format("%.2f_to_%.2f", y_min, y_max) + ".root";
  
  TFile* eff_file = new TFile(eff_file_name,"recreate");
  
  TVectorD efficiency(1);
  TVectorD stat_err_lo(1);
  TVectorD stat_err_hi(1);

  efficiency[0] = eff_res->getVal();
  stat_err_lo[0] = -eff_res->getAsymErrorLo();
  stat_err_hi[0] = eff_res->getAsymErrorHi();
  
  efficiency.Write("val");
  stat_err_lo.Write("err_lo");
  stat_err_hi.Write("err_hi");

  eff_file->Close();
  delete eff_file; 
  
  /////////////////
  //for debugging//
  /////////////////
  /*
  TFile* f3 = new TFile(eff_file_name);
  TVectorD *eff_val = (TVectorD*)f3->Get("val");
  TVectorD *stat_lo = (TVectorD*)f3->Get("err_lo");
  TVectorD *stat_hi = (TVectorD*)f3->Get("err_hi");
  
  eff_val->Print();
  stat_lo->Print();
  stat_hi->Print();
  delete f3;
  */
}//end of calculate_bin_eff
