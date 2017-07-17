#include "UserCode/B_production_x_sec_13_TeV/interface/functions.h"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda
//-----------------------------------------------------------------

//input example: combine_syst --measure x_sec --channel 1 --bins pt_y
int main(int argc, char** argv)
{
  std::vector<std::string> x_sec_syst = {"signal_pdf","cb_pdf","mass_window"};
  std::vector<std::string> ratio_syst = {"signal_pdf","cb_pdf","mass_window"};

  TString measure = "";
  int channel = 1;
  TString bins = "";

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;
      
      if(argument == "--measure")
	{
	  convert << argv[++i];
	  convert >> measure;
	}
      if(argument == "--channel")
	{
	  convert << argv[++i];
	  convert >> channel;
	}
      if(argument == "--bins")
	{
	  convert << argv[++i];
	  convert >> bins;
	}
    }

  if(bins == "")
    {
      std::cout << "Error: Enter the --bins option." << std::endl;
      return 0;
    }

  std::vector<std::string> syst;

  //copy the right syst names
  if(measure == "x_sec")
    {
      syst.reserve((int)x_sec_syst.size());
      
      for(int i=0; i<(int)x_sec_syst.size(); i++)
        syst.push_back(x_sec_syst[i]);
    }
  else
    if(measure == "ratio")
      {
        syst.reserve((int)ratio_syst.size());
	
        for(int i=0; i<(int)ratio_syst.size(); i++)
          syst.push_back(ratio_syst[i]);
      }

  TString var1_name = "";
  int n_var1_bins=1;
  TString var2_name = "";
  int n_var2_bins=1;
  
  double* var1_bins = NULL;
  double* var2_bins = NULL;
  
  setup_bins(measure,channel, bins, &var1_name, &n_var1_bins, &var2_name, &n_var2_bins, &var1_bins, &var2_bins);
  
  /////////////////////////////////////////////////////////////////////
  /////////////////////////////create vector///////////////////////////
  /////////////////////////////////////////////////////////////////////

  TString out_file_name = "";
  TString in_file_name = "";
  TString bins_str = "";
  TString dir = "";
  TString vector = "signal_pdf";
  
  dir = TString::Format(VERSION) + "/signal_yield_root/syst/" + channel_to_ntuple_name(channel) + "/"; 

  double sqrt_err_lo[n_var1_bins];
  double sqrt_err_hi[n_var1_bins];
  
  for(int i=0; i<n_var1_bins; i++)
    {
      sqrt_err_lo[i] = 0;
      sqrt_err_hi[i] = 0;
    }
  
  for(int j=0; j<n_var2_bins; j++)
    {
      TVectorD out_val(n_var1_bins);
      TVectorD out_err_lo(n_var1_bins);
      TVectorD out_err_hi(n_var1_bins);

      if(var1_name == "pt")
	bins_str = TString::Format("%.2f_to_%.2f", var2_bins[j], var2_bins[j+1]);
      else
	bins_str = TString::Format("%d_to_%d", (int)var2_bins[j], (int)var2_bins[j+1]);
            
      /////////////////////vector_name cicle
      for(int k=0; k<(int)syst.size(); k++)
	{
	  in_file_name = dir + syst[k] + "_vector_" + measure + "_" + channel_to_ntuple_name(channel) + "_" + var1_name + "_bins_" + var2_name + "_from_" + bins_str + ".root";
 
	  if(bins == "full")
	    in_file_name = dir + syst[k] + "_vector_" + measure + "_" + channel_to_ntuple_name(channel) + "_full_bins.root";

	  //debug
	  std::cout << "read :" << in_file_name << std::endl;
      	  
	  //open input file
	  TFile* fin = new TFile(in_file_name);
	  TVectorD *in_err_lo = (TVectorD*)fin->Get("err_lo");
	  TVectorD *in_err_hi = (TVectorD*)fin->Get("err_hi");
	  delete fin;
	  
	  for(int i=0; i<n_var1_bins; i++)
	    {
	      sqrt_err_lo[i] += pow(in_err_lo[0][i],2);
	      sqrt_err_hi[i] += pow(in_err_hi[0][i],2);
	    }
	}//end of vector_name cicle
      
      //copy the value to the output
      for(int i=0; i<n_var1_bins; i++)
	{
	  out_val[i] = 1;
	  out_err_lo[i] = sqrt(sqrt_err_lo[i]);
	  out_err_hi[i] = sqrt(sqrt_err_hi[i]);
	}
      
      out_file_name = dir + "combined_syst_vector_" + measure + "_" + channel_to_ntuple_name(channel) + "_" + var1_name + "_bins_" + var2_name + "_from_" + bins_str + ".root";
      
      if(bins == "full")
	out_file_name = dir + "combined_syst_vector_" + measure + "_" + channel_to_ntuple_name(channel) + "_full_bins.root";

      //create output file
      TFile* fout = new TFile(out_file_name,"recreate");
	  
      //write to the output file
      out_val.Write("val");
      out_err_lo.Write("err_lo");
      out_err_hi.Write("err_hi");
        
      fout->Close();
      delete fout;
      
      //debug
      std::cout << "write :" << out_file_name << std::endl;
      out_val.Print();
      out_err_lo.Print();
      out_err_hi.Print();
    }//end of var2 cicle

}//end
