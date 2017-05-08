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

//input example: create_table --measure x_sec --channel 1 --bins pt --vector yield
int main(int argc, char** argv)
{
  TString measure = "";
  int channel = 1;
  std::string bins = "pt";
  TString vector = "";

  for(int i=1 ; i<argc ; ++i)
    {
      std::string argument = argv[i];
      std::stringstream convert;

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
      if(argument == "--measure")
	{
	  convert << argv[++i];
	  convert >> measure;
	}
      if(argument == "--vector")
	{
	  convert << argv[++i];
	  convert >> vector;
	}
    }
  
  if(measure == "" || vector == "")
    {
      std::cout << "ERROR: No --measure or --vector input was provided." << std::endl;
      return 0;
    }
  
  //create dir to save the tables
  std::vector<std::string> dir_list;
  dir_list.push_back("tables/");
  create_dir(dir_list);
  
  //set up the vectors
  TString var1_name = "";
  TString var2_name = "";
  int n_var1_bins=1;
  int n_var2_bins=1;
  double* var1_bins = NULL;
  double* var2_bins = NULL;
  
  setup_bins(measure, channel, bins, &var1_name, &n_var1_bins, &var2_name, &n_var2_bins, &var1_bins, &var2_bins);
  
  //initialize arrays for yield, efficiencies, etc 
  double val_array[n_var2_bins][n_var1_bins];
  double val_err_lo[n_var2_bins][n_var1_bins];
  double val_err_hi[n_var2_bins][n_var1_bins];
  
  TString b_title= "";
  switch(channel)
    {
    default:
    case 1:
      b_title = "B+";
      break;
    case 2:
      b_title = "B0";
      break;
    case 4:
      b_title = "Bs";
      break;
    }
  
  //read yields or effs or syst
  read_vector(measure, channel, vector, var1_name , var2_name, n_var1_bins, n_var2_bins, var1_bins, var2_bins, val_array[0], val_err_lo[0], val_err_hi[0]);
  
  std::vector<std::string> col_name;

  std::string var1_str = var1_name.Data();
  var1_str.append(" bins");
  
  std::string table_str = vector.Data();

  col_name.push_back(var1_str);
  col_name.push_back(table_str);
  col_name.push_back("errlo");
  col_name.push_back("errhi");

  for(int j=0; j<n_var2_bins; j++)
    {
      //maybe set the precision of the doubles here
      std::vector<std::string> labels;
      std::string lab;

      for(int i=0; i<n_var1_bins; i++)
	{
	  if(var1_name == "pt")
	    lab = std::to_string((int)var1_bins[i]) + " to " + std::to_string((int)var1_bins[i+1]);
	  else
	    lab = std::to_string((double)var1_bins[i]) + " to " + std::to_string((double)var1_bins[i+1]);

	  labels.push_back(lab);
	}
  
      std::vector<std::vector<double> > numbers;
      std::vector<double> aux;

      for(int i=0; i<n_var1_bins; i++)
	{
	  aux.push_back(val_array[j][i]);
	}
      numbers.push_back(aux);
      aux.clear();

      for(int i=0; i<n_var1_bins; i++)
	{
	  aux.push_back(val_err_lo[j][i]);
	}
      numbers.push_back(aux);
      aux.clear();

      for(int i=0; i<n_var1_bins; i++)
	{
	  aux.push_back(val_err_hi[j][i]);
	}
      numbers.push_back(aux);
      aux.clear();
      
      TString file_name = "";
      TString dir = "tables/";
      TString bins_str ="";
      TString ntuple_name = "";

      if(vector == "Bs_Bu" || vector == "fs_fu" || vector == "Bs_Bd" || vector == "fs_fd" || vector == "Bd_Bu" || vector == "fd_fu")
	ntuple_name = "";
      else
	ntuple_name = channel_to_ntuple_name(channel) + "_";
      
      if(var1_name == "pt")
	bins_str = TString::Format("%.2f_to_%.2f", var2_bins[j], var2_bins[j+1]);
      else
	bins_str = TString::Format("%d_to_%d", (int)var2_bins[j], (int)var2_bins[j+1]);

      file_name = dir + vector + "_" + measure + "_" + ntuple_name + var1_name + "_bins_" + var2_name + "_from_" + bins_str + "_" + TString::Format(VERSION);
      
      if(bins == "full")
	file_name = dir + vector + "_" + measure + "_" + ntuple_name + "full_bins_" + TString::Format(VERSION);
      
      TString bins_cap ="";
      
      if(var1_name == "pt")
	bins_cap = TString::Format("$%.2f$ to $%.2f$", var2_bins[j], var2_bins[j+1]);
      else
	bins_cap = TString::Format("$%d$ to $%d$", (int)var2_bins[j], (int)var2_bins[j+1]);

      TString caption = b_title + " " + vector + " " + var2_name + " from " + bins_cap;
      
      latex_table(file_name.Data(), col_name.size(), n_var1_bins + 1, col_name, labels, numbers, caption.Data());
    }//end of var2 cicle
  
}//end of create_table
