#include <TString.h>

TString channel_to_ntuple_name(int channel)
{
  //returns a TString with the ntuple name corresponding to the channel. It can be used to find the data on each channel saved in a file. or to write the name of a directory                                                                                                                               
  TString ntuple_name = "";
  
  switch(channel){
  default:
  case 1:
    ntuple_name="ntkp";
    break;
  case 2:
    ntuple_name="ntkstar";
    break;
  case 3:
   ntuple_name="ntks";
   break;
  case 4:
    ntuple_name="ntphi";
    break;
  case 5:
   ntuple_name="ntmix";
   break;
  case 6:
    ntuple_name="ntlambda";
    break;
  case 7:
    ntuple_name="ntpi";
    break;
 }
  return ntuple_name;
}

TString channel_to_xaxis_title(int channel)
{
  TString xaxis_title = "";
  
  switch (channel) {
  default:
  case 1:
    xaxis_title = "M_{J/#psi K^{#pm}} [GeV]";
    break;
  case 2:
    xaxis_title = "M_{J/#psi K^{#pm}#pi^{#mp}} [GeV]";
    break;
  case 3:
    xaxis_title = "M_{J/#psi K^{0}_{S}} [GeV]";
    break;
  case 4:
    xaxis_title = "M_{J/#psi K^{#pm}K^{#mp}} [GeV]";
    break;
  case 5:
    xaxis_title = "M_{J/#psi #pi^{#pm}#pi^{#mp}} [GeV]";
    break;
  case 6:
    xaxis_title = "M_{J/#psi #Lambda} [GeV]";
    break;
  }
  return xaxis_title;
}

int channel_to_nbins(int channel)
{
  int nbins;

  switch (channel) {
  default:
  case 1:
    nbins = 50;
    break;
  case 2:
    nbins = 50;
    break;
  case 3:
    nbins = 50;
    break;
  case 4:
    nbins = 50;
    break;
  case 5:
    nbins = 80;
    break;
  case 6:
    nbins = 50;
    break;
  }
  return nbins;
}
