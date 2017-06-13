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

//input example: create_yield_eff_syst_vector --measure x_sec --channel 1 --bins pt_y --vector yield --calculate 0
int main(int argc, char** argv)
{
  TString measure = "";
  int channel = 1;
  TString bins = "";
  TString vector = "";
  int calculate = 0;

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
      if(argument == "--vector")
	{
	  convert << argv[++i];
	  convert >> vector;
	}
      if(argument == "--calculate")
	{
	  convert << argv[++i];
	  convert >> calculate;
	}
    }

  if(bins == "")
    {
      std::cout << "Error: Enter the --bins option." << std::endl;
      return 0;
    }
  
  if(vector == "")
    {
      std::cout << "Error: Enter the --vector option." << std::endl;
      return 0;
    }

  TString var1_name = "";
  int n_var1_bins=1;
  TString var2_name = "";
  int n_var2_bins=1;
  
  double* var1_bins = NULL;
  double* var2_bins = NULL;
  
  setup_bins(measure,channel, bins, &var1_name, &n_var1_bins, &var2_name, &n_var2_bins, &var1_bins, &var2_bins);
  
  TString line = "";
  TString command = "";
  TString opt = "";
  
  TString pt_min;
  TString pt_max;
  TString y_min;
  TString y_max;
  
  //calculate the yield or eff or syst for all the bins, if --calculate 1
  if(calculate)
    {
      for(int j=0; j<n_var2_bins; j++)
	{
	  for(int i=0; i<n_var1_bins; i++)
	    {
	      if(var1_name == "pt")
		{
		 pt_min = TString::Format("%d", (int)var1_bins[i]);
		 pt_max = TString::Format("%d", (int)var1_bins[i+1]);
		 y_min = TString::Format("%.2f", var2_bins[j]);
		 y_max = TString::Format("%.2f", var2_bins[j+1]);
		}
	      else
		{
		  pt_min = TString::Format("%d", (int)var2_bins[j]);
		  pt_max = TString::Format("%d", (int)var2_bins[j+1]);
		  y_min = TString::Format("%.2f", var1_bins[i]);
		  y_max = TString::Format("%.2f", var1_bins[i+1]);
		}
	      
	      opt = " --channel " + TString::Format("%d", channel) + " --ptmin " + pt_min + " --ptmax " + pt_max + " --ymin " +  y_min + " --ymax " + y_max; 
	  
	      if(vector == "yield")
		{
		  command = "calculate_bin_yield";
		  line = command + opt;
		}
	      else
		if(vector == "preeff" || vector == "recoeff" || vector == "totaleff")
		  {
		    command = "calculate_bin_eff";
		    line = command + opt + " --eff " + vector;
		  }
		else
		  if(vector == "signal_pdf" || vector == "cb_pdf" || vector == "mass_window")
		    {
		      command = "calculate_bin_syst";
		      line = command + opt + " --syst " + vector;
		    }
	      
	      gSystem->Exec(line);
	      
	    }//end of var1 cicle
	}//end of var2 cicle
    }//end of if(calculate)
  
  /////////////////////////////////////////////////////////////////////
  /////////////////////////////create vector///////////////////////////
  /////////////////////////////////////////////////////////////////////

  TString out_file_name = "";
  TString in_file_name = "";
  TString bins_str = "";
  TString dir = "";

  if(vector == "yield")
    dir = "signal_yield_root/";
  else
    if(vector == "signal_pdf" || vector == "cb_pdf" || vector == "mass_window")
      dir = "signal_yield_root/syst/";
    else
      dir = "efficiencies_root/";
      
  dir += channel_to_ntuple_name(channel) + "_" + TString::Format(VERSION) + "/"; 

  for(int j=0; j<n_var2_bins; j++)
    {
      TVectorD out_val(n_var1_bins);
      TVectorD out_err_lo(n_var1_bins);
      TVectorD out_err_hi(n_var1_bins);
	  
      for(int i=0; i<n_var1_bins; i++)
	{
	  if(var1_name == "pt")
	    in_file_name = dir + vector + "_" + channel_to_ntuple_name(channel) + "_pt_from_" + TString::Format("%d_to_%d", (int)var1_bins[i], (int)var1_bins[i+1]) + "_y_from_" + TString::Format("%.2f_to_%.2f", var2_bins[j], var2_bins[j+1]) + "_" + TString::Format(VERSION) + ".root";
	  else
	    in_file_name = dir + vector + "_" + channel_to_ntuple_name(channel) + "_pt_from_" + TString::Format("%d_to_%d", (int)var2_bins[j], (int)var2_bins[j+1]) + "_y_from_" + TString::Format("%.2f_to_%.2f", var1_bins[i], var1_bins[i+1]) + "_" + TString::Format(VERSION) + ".root";
	      
	  //debug
	  std::cout << "read :" << in_file_name << std::endl;

	  //open input file
	  TFile* fin = new TFile(in_file_name);
	  TVectorD *in_val = (TVectorD*)fin->Get("val");
	  TVectorD *in_err_lo = (TVectorD*)fin->Get("err_lo");
	  TVectorD *in_err_hi = (TVectorD*)fin->Get("err_hi");
	  delete fin;

	  out_val[i] = in_val[0][0];
	  out_err_lo[i] = in_err_lo[0][0];
	  out_err_hi[i] = in_err_hi[0][0];
	}//end of var1 cicle

      //output file name
      if(var1_name == "pt")
	bins_str = TString::Format("%.2f_to_%.2f", var2_bins[j], var2_bins[j+1]);
      else
	bins_str = TString::Format("%d_to_%d", (int)var2_bins[j], (int)var2_bins[j+1]);
	  
      out_file_name = dir + vector + "_vector_" + measure + "_" + channel_to_ntuple_name(channel) + "_" + var1_name + "_bins_" + var2_name + "_from_" + bins_str + "_" + TString::Format(VERSION) + ".root";

      if(bins == "full")
	out_file_name = dir + vector + "_vector_" + measure + "_" + channel_to_ntuple_name(channel) + "_full_bins_" + TString::Format(VERSION) + ".root";

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
