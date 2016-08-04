#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include "PhysicsTools/Utilities/interface/SideBandSubtraction.h"
#include "myloop.h"
#include "plotDressing.h"
#include "TMath.h"
using namespace RooFit;

// General fitting options
#define NUMBER_OF_CPU       1
//#define DO_MINOS            kTRUE
#define SELECT_DATA         0
#define YIELD_SUB_SAMPLES   0

#define SOURCE1             "/lstore/cms/brunogal/25ns_analysis_v1/myloop_data_run2015D_v4_v1.root"
#define SOURCE2             "/lstore/cms/brunogal/25ns_analysis_v1/myloop_data_run2015D_v3_v1.root"
#define SOURCE3             "/lstore/cms/brunogal/25ns_analysis_v1/myloop_data_run2015C_v1_v1.root"

#define VERSION             "v5"
//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda

void data_selection(TString fin1,TString fin2,TString fin3,TString data_selection_output_file);
void fityield(RooDataSet* data, double** yield, double** errHi, double** errLo, int channel, RooRealVar mass, double mass_min, double mass_max, double mass_peak, int nbins, TString xaxis_title, TString ntuple_name, double pt_bin_lo, double pt_bin_hi);

void signal_yield_bin_center(int channel=1)
{
    double mass_min, mass_max, mass_peak;
    double pt_min, pt_max;
    int nbins;
    TString ntuple_name = "", xaxis_title = "";
 
    double ntkp_pt_bin_edges[]={10,20,30,40,50,60,70,80,90,100,120};
    double ntkstar_pt_bin_edges[]={0};
    double ntks_pt_bin_edges[]={10,15,20,25,30,35,40,45,50,55,70};
    double ntphi_pt_bin_edges[]={10,15,20,25,30,35,40,45,50,55,65,80};
    double ntmix_pt_bin_edges[]={0};
    double ntlambda_pt_bin_edges[]={0};
    double* pt_bin_edges;
    
    pt_min=5;
    pt_max=400;

    int nptbins;
         
    switch (channel) {
        case 1:
	  pt_bin_edges = ntkp_pt_bin_edges;
	  nptbins = (sizeof(ntkp_pt_bin_edges) / sizeof(double)) - 1 ; //if pt_bin_edges is an empty array, then nptbins is equal to 0
	  mass_min = 5.0; mass_max = 6.0;
	  mass_peak = BP_MASS;
	  nbins = 50;
	  ntuple_name = "ntkp";
	  xaxis_title = "M_{J/#psi K^{#pm}} [GeV]";
	  break;
        case 2:
	  pt_bin_edges = ntkstar_pt_bin_edges;
	  nptbins = (sizeof(ntkstar_pt_bin_edges) / sizeof(double)) - 1 ;
	  mass_min = 5.0; mass_max = 6.0;
	  mass_peak = B0_MASS;
	  nbins = 50;
	  ntuple_name = "ntkstar";
	  xaxis_title = "M_{J/#psi K^{#pm}#pi^{#mp}} [GeV]";
	  break;
        case 3:
	  pt_bin_edges = ntks_pt_bin_edges;
	  nptbins = (sizeof(ntks_pt_bin_edges) / sizeof(double)) - 1 ;
	  mass_min = 5.0; mass_max = 6.0;
	  mass_peak = B0_MASS;
	  nbins = 50;
	  ntuple_name = "ntks";
	  xaxis_title = "M_{J/#psi K^{0}_{S}} [GeV]";
	  break;
        case 4:
	  pt_bin_edges = ntphi_pt_bin_edges;
	  nptbins = (sizeof(ntphi_pt_bin_edges) / sizeof(double)) - 1 ;
	  mass_min = 5.0; mass_max = 6.0;
	  mass_peak = BS_MASS;
	  nbins = 50;
	  ntuple_name = "ntphi";
	  xaxis_title = "M_{J/#psi K^{#pm}K^{#mp}} [GeV]";
	  break;
        case 5:
	  pt_bin_edges = ntmix_pt_bin_edges;
	  nptbins = (sizeof(ntmix_pt_bin_edges) / sizeof(double)) - 1 ;
	  mass_min = 3.6; mass_max = 4.0;
	  mass_peak = PSI2S_MASS;
	  nbins = 80;
	  ntuple_name = "ntmix";
	  xaxis_title = "M_{J/#psi #pi^{#pm}#pi^{#mp}} [GeV]";
	  break;
        case 6:
	  pt_bin_edges = ntlambda_pt_bin_edges;
	  nptbins = (sizeof(ntlambda_pt_bin_edges) / sizeof(double)) - 1 ;
	  mass_min = 5.3; mass_max = 6.3;
	  mass_peak = LAMBDAB_MASS;
	  nbins = 50;
	  ntuple_name = "ntlambda";
	  xaxis_title = "M_{J/#psi #Lambda} [GeV]";
	  break;
    }
   
    TString data_selection_output_file = "selected_data.root";
   
    if(SELECT_DATA)
      data_selection(SOURCE1,SOURCE2,SOURCE3,data_selection_output_file);
   
    TFile* f = new TFile(data_selection_output_file);
    TNtupleD* _nt = (TNtupleD*)f->Get(ntuple_name);
  
    RooRealVar mass("mass","mass",mass_min,mass_max);
    RooRealVar pt("pt","pt",pt_min,pt_max);

    RooDataSet* data_original = new RooDataSet("data_original","data_original",_nt,RooArgSet(mass,pt));

    TString mass_signal_cut = TString::Format("abs(mass-%f)<0.1", mass_peak);
    TString mass_background_cut = TString::Format("abs(mass-%f)>0.1", mass_peak);

    RooDataSet* data_signal = (RooDataSet*) data_original->reduce(Cut(mass_signal_cut));
    RooDataSet* data_background = (RooDataSet*) data_original->reduce(Cut(mass_background_cut));
  
    //#########To cut the dataset in bins of pt, preform the fit for each bin, and extract the yield##########################################

    if(YIELD_SUB_SAMPLES && nptbins!=0)
      {
    RooThresholdCategory ptRegion("ptRegion", "region of pt", pt);

    ptRegion.addThreshold(*(pt_bin_edges),"below 1st bin");

    for(int i=0; i<nptbins; i++)
      {
       TString reg = TString::Format("PtBin%d",i+1);    
       ptRegion.addThreshold(*(pt_bin_edges+i+1),reg);    
      }
    data_original->addColumn(ptRegion);

    Roo1DTable * tab = data_original->table(ptRegion);
    tab->Print("v");
    delete tab;

    //to produce and process each pt subsample
    RooDataSet *data;
    RooRealVar* pt_mean;
    double* yield;
    double* errHi;
    double* errLo;
    double pt_bin_centre[nptbins];
    double pt_bin_edges_Lo[nptbins];
    double pt_bin_edges_Hi[nptbins];
    double yield_array[nptbins];
    double errLo_array[nptbins];
    double errHi_array[nptbins];
        
    for(int i=0; i<nptbins; i++)
      {
	cout << "processing subsample pt: " << i+1 << std::endl;

	TString ptcut(TString::Format("(ptRegion==ptRegion::PtBin%d)", i+1));
	
	data = new RooDataSet("data", "data", *(data_original->get()),Import(*data_original), Cut(ptcut));
	pt_mean = data->meanVar(pt);
	
	//pt_bin_centre[i] = *(pt_bin_edges+i) + (*(pt_bin_edges+i+1)-*(pt_bin_edges+i))/2;
	pt_bin_centre[i] = (double) pt_mean->getVal();

	pt_bin_edges_Lo[i] = pt_bin_centre[i] - *(pt_bin_edges+i);
	pt_bin_edges_Hi[i] = *(pt_bin_edges+i+1) - pt_bin_centre[i];
	
    fityield(data, &yield, &errHi, &errLo, channel, mass, mass_min, mass_max, mass_peak, nbins, xaxis_title, ntuple_name, *(pt_bin_edges+i), *(pt_bin_edges+i+1));

	yield_array[i] = *yield;
	errLo_array[i] = -(*errLo);
	errHi_array[i] = *errHi;
     }

    for(int i=0; i<nptbins; i++)
      {
	std::cout << "BIN: "<< *(pt_bin_edges+i) << " to " << *(pt_bin_edges+i+1) << " : " <<  yield_array[i] << " +" << errHi_array[i] << " -" << errLo_array[i] << std::endl;
      }

    TCanvas cz;
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(nptbins, pt_bin_centre, yield_array, pt_bin_edges_Lo, pt_bin_edges_Hi, errLo_array, errHi_array);
    graph->SetTitle("Raw signal yield in Pt bins");
    graph->SetFillColor(2);
    graph->SetFillStyle(3001);
    graph->Draw("a2");
    graph->Draw("p");
    cz.SetLogy();
    cz.SaveAs("signal_yield/signal_yield_" + ntuple_name + "_" + TString::Format(VERSION) + ".root");
    cz.SaveAs("signal_yield/signal_yield_" + ntuple_name + "_" + TString::Format(VERSION) + ".png");
      }
else
  {
    RooDataSet* data = data_original;

    /*
    //for sideband background subtraction
    TH1F massHist("massHist","massHist",5,0,6);
    TH1F ptHist("ptHist","ptHist",0,0,350);

    vector<TH1F*> basehistos(0);
    basehistos.push_back(&massHist);
    basehistos.push_back(&ptHist);
    //---------------------------------------
    */

    double n_signal_initial = data->sumEntries(TString::Format("abs(mass-%g)<0.015",mass_peak))
      - data->sumEntries(TString::Format("abs(mass-%g)<0.030&&abs(mass-%g)>0.015",mass_peak,mass_peak));

    double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
      
    //-----------------------------------------------------------------
    // signal PDF 
    
    RooRealVar m_mean("m_mean","m_mean",mass_peak,mass_min,mass_max);
    RooRealVar m_sigma1("m_sigma1","m_sigma1",0.015,0.001,0.050);
    RooRealVar m_sigma2("m_sigma2","m_sigma2",0.030,0.001,0.100);
    RooRealVar m_fraction("m_fraction","m_fraction",0.5);
    
    RooGaussian m_gaussian1("m_gaussian1","m_gaussian1",mass,m_mean,m_sigma1);
    RooGaussian m_gaussian2("m_gaussian2","m_gaussian2",mass,m_mean,m_sigma2);
    
    RooAddPdf pdf_m_signal("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));
    
    // use single Gaussian for J/psi Ks and J/psi Lambda due to low statistics
    if (channel==3 || channel==6) {
      m_sigma2.setConstant(kTRUE);
      m_fraction.setVal(1.);
    }
    
    //-----------------------------------------------------------------
    // combinatorial background PDF (exponential or bernstean poly.)
    
    RooRealVar m_exp("m_exp","m_exp",-0.3,-4.,+4.);
    RooExponential pdf_m_combinatorial_exp("pdf_m_combinatorial_exp","pdf_m_combinatorial_exp",mass,m_exp);
    
    RooRealVar m_par1("m_par1","m_par2",1.,0,+10.);
    RooRealVar m_par2("m_par2","m_par3",1.,0,+10.);
    RooRealVar m_par3("m_par3","m_par3",1.,0,+10.);
    RooBernstein pdf_m_combinatorial_bern("pdf_m_combinatorial_bern","pdf_m_combinatorial_bern",mass,RooArgList(RooConst(1.),m_par1,m_par2,m_par3));
    
    //-----------------------------------------------------------------
    // B->J/psi+X background PDF, for B+ -> J/psi K+ and B0->J/Psi Ks
    
    //	RooRealVar m_jpsix_mean("m_jpsix_mean","m_jpsix_mean",5.07,5.0,5.2);
    //	RooRealVar m_jpsix_sigma("m_jpsix_sigma","m_jpsix_sigma",0.05,0.01,0.10);
    //	RooGaussian pdf_m_jpsix("pdf_m_jpsix","pdf_m_jpsix",mass,m_jpsix_mean,m_jpsix_sigma);
	
    //RooRealVar x("x","x",5.13,5.11,5.15);
    //RooRealVar width("width","width",0.055,0.02,0.1);
    RooFormulaVar pdf_m_jpsix("pdf_m_jsix","2.7*erfc((mass-5.14)/(0.5*0.08))",{mass});
	
    //-----------------------------------------------------------------
    // X(3872) PDF, only for J/psi pipi fit
    
    RooRealVar m_x3872_mean("m_x3872_mean","m_x3872_mean",3.872,3.7,3.9);
    RooRealVar m_x3872_sigma("m_x3872_sigma","m_x3872_sigma",0.01,0.001,0.010);
    RooGaussian pdf_m_x3872("pdf_m_x3872","pdf_m_x3872",mass,m_x3872_mean,m_x3872_sigma);
    
    //-----------------------------------------------------------------
    // full model
    
    RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
    RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
    RooRealVar n_x3872("n_x3872","n_x3872",200.,0.,data->sumEntries());

    RooRealVar n_jpsix("n_jpsix","n_jpsix",data->sumEntries(TString::Format("mass>4.9&&mass<5.14")),data->sumEntries(TString::Format("mass>4.9&&mass<5.14")),data->sumEntries());
    
    RooAddPdf* model;
    
    if (channel==1 || channel ==3) // B+ -> J/psi K+, B0 -> J/psi Ks
      model = new RooAddPdf("model","model",
			    RooArgList(pdf_m_signal, pdf_m_combinatorial_exp, pdf_m_jpsix),
			    RooArgList(n_signal, n_combinatorial, n_jpsix));
    else
      if (channel==2 || channel==4 || channel==6) // B0 -> J/psi K*; Bs -> J/psi phi; Lambda_b -> J/psi Lambda
	model = new RooAddPdf("model","model",
			      RooArgList(pdf_m_signal, pdf_m_combinatorial_exp),
			      RooArgList(n_signal, n_combinatorial));
      else
	if (channel==5) // J/psi pipi
	  model = new RooAddPdf("model","model",
				RooArgList(pdf_m_signal, pdf_m_combinatorial_bern, pdf_m_x3872),
				RooArgList(n_signal, n_combinatorial, n_x3872));

    /*
    //for sideband subtration need to be corrected for channel 5, the background is not exponential
    SideBandSubtract sbs(model,&pdf_m_combinatorial_exp,data,&mass,basehistos,true);
    sbs.addSignalRegion(5.2,5.4);
    sbs.addSideBandRegion(5.5,5.7);
    sbs.doGlobalFit();
    // sbs.printResults();
    sbs.saveResults("sideband_results.root");
    //-----------------------------------------
    */

    model->fitTo(*data,Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));
	
    std::cout <<"SIGNAL: "<< n_signal.getVal() << " -" << n_signal.getAsymErrorLo() << " +" << n_signal.getAsymErrorHi() << std::endl;
   
    //full dataset pt distribution
    TCanvas c2;
    TH1D* pt_dist = (TH1D*)data->createHistogram("pt_dist",pt);
    pt_dist->Draw();
    c2.SetLogy();
    c2.SaveAs("full_dataset_mass_pt_histo/" + ntuple_name + "_" + TString::Format(VERSION) + ".root");
    c2.SaveAs("full_dataset_mass_pt_histo/" + ntuple_name + "_" + TString::Format(VERSION) + ".png");
    
    /*
    //signal pt distribution
    TCanvas c3;
    TH1D* pt_signal = (TH1D*)data_signal->createHistogram("pt_signal",pt);
    pt_signal->Draw();
    c3.SetLogy();
    c3.SaveAs("full_dataset_mass_pt_histo/" + ntuple_name + TString::Format("_pt_signal_") + TString::Format(VERSION) + ".png");

    //signal mass distribution
    TCanvas c4;
    TH1D* mass_signal = (TH1D*)data_signal->createHistogram("mass_signal",mass);
    mass_signal->Draw();
    c4.SaveAs("full_dataset_mass_pt_histo/" + ntuple_name + TString::Format("_mass_signal_") + TString::Format(VERSION) + ".png");
    */
   
    RooPlot* frame_m = mass.frame();
    
    TH1D* histo_data = (TH1D*)data->createHistogram("histo_data", mass, Binning(nbins,mass_min,mass_max));
    histo_data->Sumw2(false);
    histo_data->SetBinErrorOption(TH1::kPoisson);
    histo_data->SetMarkerStyle(20);
    histo_data->SetMarkerSize(0.8);
    histo_data->SetLineColor(kBlack);
    
    for (int i=1; i<=nbins; i++)
      if (histo_data->GetBinContent(i)==0) histo_data->SetBinError(i,0.);
    
    data->plotOn(frame_m,Name("theData"),Binning(nbins),Invisible());
    
    model->plotOn(frame_m,Name("thePdf"),Precision(2E-4));
    
    //model->paramOn(frame_m); //show all the parameters of the fit in the plot.
        
    model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_signal),LineColor(kRed),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kRed), VLines(), DrawOption("F"));

    if (channel==1 || channel==2 || channel==3 || channel==4 || channel==6)
      model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_combinatorial_exp),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    else
      model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_combinatorial_bern),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    
    if (channel==1 || channel==3)
      model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_jpsix),LineColor(kViolet),LineWidth(2),LineStyle(7));
    if (channel==5)
      model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_x3872),LineColor(kOrange),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kOrange), VLines(), DrawOption("F"));
    
    frame_m->SetTitle("");
    frame_m->GetXaxis()->SetTitle(xaxis_title);
    frame_m->GetXaxis()->SetLabelFont(42);
    frame_m->GetXaxis()->SetLabelOffset(0.01);
    frame_m->GetXaxis()->SetTitleSize(0.06);
    frame_m->GetXaxis()->SetTitleOffset(1.09);
    frame_m->GetXaxis()->SetLabelFont(42);
    frame_m->GetXaxis()->SetLabelSize(0.055);
    frame_m->GetXaxis()->SetTitleFont(42);
    frame_m->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass_max-mass_min)*1000./nbins));
    frame_m->GetYaxis()->SetLabelFont(42);
    frame_m->GetYaxis()->SetLabelOffset(0.01);
    frame_m->GetYaxis()->SetTitleOffset(1.14);
    frame_m->GetYaxis()->SetTitleSize(0.06);
    frame_m->GetYaxis()->SetTitleFont(42);
    frame_m->GetYaxis()->SetLabelFont(42);
    frame_m->GetYaxis()->SetLabelSize(0.055);
    
    RooHist* pull_hist = frame_m->residHist("theData","thePdf",true);
    RooPlot* pull_plot = mass.frame();
    pull_plot->addPlotable(pull_hist,"P");

    TCanvas c1; c1.cd();// = canvasDressing("c1");

    TPad *p1 = new TPad("p1","p1",0,0,1,1);
    p1->Draw();

    /*
    TPad *p1 = new TPad("p1","p1",0.03,0.27,0.99,0.99);
    // p1->SetBorderMode(0); 
    p1->Draw(); 
   
    TPad *p2 = new TPad("p2","p2",0.03,0.01,0.99,0.24); 
    // p2->SetTopMargin(0.);    
    // p2->SetBorderMode(0); 
    p2->SetTicks(1,2); 
    p2->Draw();
    */

    p1->cd(); frame_m->Draw(); histo_data->Draw("Esame"); Legend(channel,0,0,0);
    // p2->cd(); pull_plot->Draw();
         
    c1.SaveAs("full_dataset_mass_fit/" + ntuple_name + "_" + TString::Format(VERSION) + ".root");
    c1.SaveAs("full_dataset_mass_fit/" + ntuple_name + "_" + TString::Format(VERSION) + ".png");      
  }
}

void data_selection(TString fin1, TString fin2, TString fin3, TString data_selection_output_file){

    TString ntuple_name = "";

    TFile *fout = new TFile(data_selection_output_file,"recreate");

    TNtupleD *_nt1 = new TNtupleD("ntkp","ntkp","mass:pt");
    TNtupleD *_nt2 = new TNtupleD("ntkstar","ntkstar","mass:pt");
    TNtupleD *_nt3 = new TNtupleD("ntks","ntks","mass:pt");
    TNtupleD *_nt4 = new TNtupleD("ntphi","ntphi","mass:pt");
    TNtupleD *_nt5 = new TNtupleD("ntmix","ntmix","mass:pt");
    TNtupleD *_nt6 = new TNtupleD("ntlambda","ntlambda","mass:pt");
    
    TChain* tin;
    ReducedBranches br;
        
    int n_br_queued = 0;
    ReducedBranches br_queue[32];
    TLorentzVector v4_tk1, v4_tk2;

    for(int channel=1; channel<=6; channel++)
      {
	std::cout << "selecting data channel " << channel << std::endl;

	switch(channel){
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
	}

	tin = new TChain(ntuple_name);

	tin->Add(fin1);
	tin->Add(fin2);
	tin->Add(fin3);

	br.setbranchadd(tin);

	for (int evt=0;evt<tin->GetEntries();evt++) {
        tin->GetEntry(evt);
        
        if (channel==1) { // cuts for B+ -> J/psi K+
	  if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
	  if (br.vtxprob<=0.1) continue;
	  if (br.tk1pt<=1.6) continue;
	  if (br.lxy/br.errxy<=3.0) continue;
	  if (br.cosalpha2d<=0.99) continue;
            
	    _nt1->Fill(br.mass,br.pt);

	    //	    if (fabs(br.mass-BP_MASS)<=0.1) _nt11->Fill(br.mass,br.pt); //signal
	    //  else _nt111->Fill(br.mass,br.pt); //background

        }else
        if (channel==2) { // cuts for B0 -> J/psi K*
	  if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-KSTAR_MASS)>=0.05) continue;
            
	    v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-PHI_MASS)<=0.01) continue;
            
            if (n_br_queued==0) {
                memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
                n_br_queued++;
            }else
            if (br.run == br_queue[n_br_queued-1].run && br.event == br_queue[n_br_queued-1].event) { // same event
                memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
                n_br_queued++;
                if (n_br_queued>=32) printf("Warning: maximum queued branches reached.\n");
            }            
            if (br.run != br_queue[n_br_queued-1].run || br.event != br_queue[n_br_queued-1].event || evt==tin->GetEntries()-1) {
                for (int i=0; i<n_br_queued; i++) {
                    
                    bool isBestKstarMass = true;
                    for (int j=0; j<n_br_queued; j++) {
                        if (j==i) continue;
                        if (br_queue[i].mu1idx==br_queue[j].mu1idx &&
                            br_queue[i].mu2idx==br_queue[j].mu2idx &&
                            br_queue[i].tk1idx==br_queue[j].tk1idx &&
                            br_queue[i].tk2idx==br_queue[j].tk2idx) {
                        
                            if (fabs(br_queue[j].tktkmass-KSTAR_MASS)<fabs(br_queue[i].tktkmass-KSTAR_MASS)) {
                                isBestKstarMass = false;
                                continue;
                            }
                        }
                    }
                                 
                    if (isBestKstarMass){
		      _nt2->Fill(br_queue[i].mass,br_queue[i].pt);

		      //	      if (fabs(br_queue[i].mass-B0_MASS)<=0.1) _nt21->Fill(br_queue[i].mass,br_queue[i].pt); //signal
		      //	      else _nt211->Fill(br_queue[i].mass,br_queue[i].pt); //background
		    }

                }                
                n_br_queued = 0;
                memcpy(&br_queue[n_br_queued],&br,sizeof(ReducedBranches));
                n_br_queued++;
            }
	    }else
	  if (channel==3) { // cuts for B0 -> J/psi Ks
	    if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.tktkblxy/br.tktkberrxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-KSHORT_MASS)>=0.015) continue;
                
            _nt3->Fill(br.mass,br.pt);

	    // if (fabs(br.mass-B0_MASS)<=0.1) _nt31->Fill(br.mass,br.pt); //signal
	    // else _nt311->Fill(br.mass,br.pt); //background

        }else
        if (channel==4) { // cuts for Bs -> J/psi phi
	    if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-PHI_MASS)>=0.010) continue;
            
	    v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
                
            _nt4->Fill(br.mass,br.pt);

	    // if (fabs(br.mass-BS_MASS)<=0.1) _nt41->Fill(br.mass,br.pt); //signal
	    // else _nt411->Fill(br.mass,br.pt); //background

        }else
        if (channel==5) { // cuts for psi(2S)/X(3872) -> J/psi pipi
	    if (br.vtxprob<=0.2) continue;
            if (fabs(br.tk1eta)>=1.6) continue;
            if (fabs(br.tk2eta)>=1.6) continue;
            
            _nt5->Fill(br.mass,br.pt);

	    // if (fabs(br.mass-PSI2S_MASS)<=0.1) _nt51->Fill(br.mass,br.pt);
	    // else _nt511->Fill(br.mass,br.pt);

        }else
	  if (channel==6) {//cuts for lambda
	    if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.tktkblxy/br.tktkberrxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-LAMBDA_MASS)>=0.015) continue;
            
	    v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSHORT_MASS)<=0.015) continue;
            
            _nt6->Fill(br.mass,br.pt);

	    // if (fabs(br.mass-LAMBDAB_MASS)<=0.1) _nt61->Fill(br.mass,br.pt);
	    // else _nt611->Fill(br.mass,br.pt);
	}
	}
      }
    fout->Write();
    fout->Close();
}

void fityield(RooDataSet* data, double** yield, double** errHi, double** errLo, int channel, RooRealVar mass, double mass_min,
	      double mass_max, double mass_peak,int nbins, TString xaxis_title, TString ntuple_name,
	      double pt_bin_lo, double pt_bin_hi)
{   
        double n_signal_initial = data->sumEntries(TString::Format("abs(mass-%g)<0.015",mass_peak))
                            - data->sumEntries(TString::Format("abs(mass-%g)<0.030&&abs(mass-%g)>0.015",mass_peak,mass_peak));

	double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
    
	//-----------------------------------------------------------------
	// signal PDF 
    
	RooRealVar m_mean("m_mean","m_mean",mass_peak,mass_min,mass_max);
	RooRealVar m_sigma1("m_sigma1","m_sigma1",0.015,0.001,0.050);
	RooRealVar m_sigma2("m_sigma2","m_sigma2",0.030,0.001,0.100);
	RooRealVar m_fraction("m_fraction","m_fraction",0.5);
    
	RooGaussian m_gaussian1("m_gaussian1","m_gaussian1",mass,m_mean,m_sigma1);
	RooGaussian m_gaussian2("m_gaussian2","m_gaussian2",mass,m_mean,m_sigma2);
    
	RooAddPdf pdf_m_signal("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));
    
	// use single Gaussian for J/psi Ks and J/psi Lambda due to low statistics
	if (channel==3 || channel==6) {
	  m_sigma2.setConstant(kTRUE);
	  m_fraction.setVal(1.);
        }
    
	//-----------------------------------------------------------------
	// combinatorial background PDF (exponential or bernstean poly.)
    
	RooRealVar m_exp("m_exp","m_exp",-0.3,-4.,+4.);
	RooExponential pdf_m_combinatorial_exp("pdf_m_combinatorial_exp","pdf_m_combinatorial_exp",mass,m_exp);
    
	RooRealVar m_par1("m_par1","m_par2",1.,0,+10.);
	RooRealVar m_par2("m_par2","m_par3",1.,0,+10.);
	RooRealVar m_par3("m_par3","m_par3",1.,0,+10.);
	RooBernstein pdf_m_combinatorial_bern("pdf_m_combinatorial_bern","pdf_m_combinatorial_bern",mass,RooArgList(RooConst(1.),m_par1,m_par2,m_par3));
    
	//-----------------------------------------------------------------
	// B->J/psi+X background PDF, for B+ -> J/psi K+ and B0->J/Psi Ks
    
	//	RooRealVar m_jpsix_mean("m_jpsix_mean","m_jpsix_mean",5.07,5.0,5.2);
	//	RooRealVar m_jpsix_sigma("m_jpsix_sigma","m_jpsix_sigma",0.05,0.01,0.10);
	//	RooGaussian pdf_m_jpsix("pdf_m_jpsix","pdf_m_jpsix",mass,m_jpsix_mean,m_jpsix_sigma);
	
	//RooRealVar x("x","x",5.13,5.11,5.15);
	//RooRealVar width("width","width",0.055,0.02,0.1);
        RooFormulaVar pdf_m_jpsix("pdf_m_jsix","2.7*erfc((mass-5.14)/(0.5*0.08))",{mass});
	
	//-----------------------------------------------------------------
	// X(3872) PDF, only for J/psi pipi fit
    
	RooRealVar m_x3872_mean("m_x3872_mean","m_x3872_mean",3.872,3.7,3.9);
	RooRealVar m_x3872_sigma("m_x3872_sigma","m_x3872_sigma",0.01,0.001,0.010);
	RooGaussian pdf_m_x3872("pdf_m_x3872","pdf_m_x3872",mass,m_x3872_mean,m_x3872_sigma);
    
	//-----------------------------------------------------------------
	// full model
    
	RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
	RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
	RooRealVar n_x3872("n_x3872","n_x3872",200.,0.,data->sumEntries());

	RooRealVar n_jpsix("n_jpsix","n_jpsix",data->sumEntries(TString::Format("mass>4.9&&mass<5.14")),data->sumEntries(TString::Format("mass>4.9&&mass<5.14")),data->sumEntries());
    
	RooAddPdf* model;
    
	if (channel==1 || channel ==3) // B+ -> J/psi K+, B0 -> J/psi Ks
	  model = new RooAddPdf("model","model",
				RooArgList(pdf_m_signal, pdf_m_combinatorial_exp, pdf_m_jpsix),
				RooArgList(n_signal, n_combinatorial, n_jpsix));
	else
	  if (channel==2 || channel==4 || channel==6) // B0 -> J/psi K*; Bs -> J/psi phi; Lambda_b -> J/psi Lambda
	    model = new RooAddPdf("model","model",
				  RooArgList(pdf_m_signal, pdf_m_combinatorial_exp),
				  RooArgList(n_signal, n_combinatorial));
	  else
	    if (channel==5) // J/psi pipi
	      model = new RooAddPdf("model","model",
				    RooArgList(pdf_m_signal, pdf_m_combinatorial_bern, pdf_m_x3872),
				    RooArgList(n_signal, n_combinatorial, n_x3872));


	model->fitTo(*data,Minos(kTRUE),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));
	
	double y;
	double y_err_hi;
	double y_err_lo;

	y = n_signal.getVal();
	y_err_hi = n_signal.getAsymErrorHi();
	y_err_lo = n_signal.getAsymErrorLo();

	*yield = &y;
	*errHi = &y_err_hi;
	*errLo = &y_err_lo;

    //display the mass fit of each subsample
	
    TCanvas *c1 = canvasDressing("c1");
        
    RooPlot* frame_m = mass.frame();
    
    TH1D* histo_data = (TH1D*)data->createHistogram("histo_data", mass, Binning(nbins,mass_min,mass_max));
    histo_data->Sumw2(false);
    histo_data->SetBinErrorOption(TH1::kPoisson);
    histo_data->SetMarkerStyle(20);
    histo_data->SetMarkerSize(0.8);
    histo_data->SetLineColor(kBlack);
    for (int i=1; i<=nbins; i++)
        if (histo_data->GetBinContent(i)==0) histo_data->SetBinError(i,0.);
    
    data->plotOn(frame_m,Binning(nbins),Invisible());
    
    model->plotOn(frame_m,Precision(2E-4));
    
    //model->paramOn(frame_m); //show all the parameters of the fit in the plot.
        
    model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_signal),LineColor(kRed),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kRed), VLines(), DrawOption("F"));

    if (channel==1 || channel==2 || channel==3 || channel==4 || channel==6)
        model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_combinatorial_exp),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    else
        model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_combinatorial_bern),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    
    if (channel==1 || channel==3)
        model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_jpsix),LineColor(kViolet),LineWidth(2),LineStyle(7));
    if (channel==5)
        model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_x3872),LineColor(kOrange),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kOrange), VLines(), DrawOption("F"));
    
    frame_m->SetTitle("");
    frame_m->GetXaxis()->SetTitle(xaxis_title);
    frame_m->GetXaxis()->SetLabelFont(42);
    frame_m->GetXaxis()->SetLabelOffset(0.01);
    frame_m->GetXaxis()->SetTitleSize(0.06);
    frame_m->GetXaxis()->SetTitleOffset(1.09);
    frame_m->GetXaxis()->SetLabelFont(42);
    frame_m->GetXaxis()->SetLabelSize(0.055);
    frame_m->GetXaxis()->SetTitleFont(42);
    frame_m->GetYaxis()->SetTitle(TString::Format("Events / %g MeV",(mass_max-mass_min)*1000./nbins));
    frame_m->GetYaxis()->SetLabelFont(42);
    frame_m->GetYaxis()->SetLabelOffset(0.01);
    frame_m->GetYaxis()->SetTitleOffset(1.14);
    frame_m->GetYaxis()->SetTitleSize(0.06);
    frame_m->GetYaxis()->SetTitleFont(42);
    frame_m->GetYaxis()->SetLabelFont(42);
    frame_m->GetYaxis()->SetLabelSize(0.055);
    
    frame_m->Draw();
    histo_data->Draw("Esame");
    /*    
    if (channel==1)
        LegendChannelOne();
    if (channel==2)
        LegendChannelTwo();
    if (channel==4)
        LegendChannelFour();
    if (channel==5)
        LegendChannelFive();
     */
    Legend(channel,(int)pt_bin_lo,(int)pt_bin_hi,1);

    c1->SaveAs("pt_bin_mass_fit/" + ntuple_name + "_" + TString::Format(VERSION) + "/" + "mass_fit_" + ntuple_name + TString::Format("_bin_%d_%d", (int)pt_bin_lo, (int)pt_bin_hi) + ".root");
    c1->SaveAs("pt_bin_mass_fit/" + ntuple_name + "_" + TString::Format(VERSION) + "/" + "mass_fit_" + ntuple_name + TString::Format("_bin_%d_%d", (int)pt_bin_lo, (int)pt_bin_hi) + ".png"); 
}
