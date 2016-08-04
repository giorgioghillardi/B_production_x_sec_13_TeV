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
#include "myloop.h"
#include "plotDressing.h"

using namespace RooFit;

// General fitting options
#define NUMBER_OF_CPU       1
#define DO_MINOS            kTRUE
// 0 - w/o DISPLAY
// 1 - w/  DISPLAY
#define DISPLAY             1

#define SOURCE              "myloop_data_test.root"

//-----------------------------------------------------------------
// Definition of channel #
// channel = 1: B+ -> J/psi K+
// channel = 2: B0 -> J/psi K*
// channel = 3: B0 -> J/psi Ks
// channel = 4: Bs -> J/psi phi
// channel = 5: Jpsi + pipi
// channel = 6: Lambda_b -> Jpsi + Lambda

void myfitter(int channel=6)
{
    double mass_min, mass_max, mass_peak;
    double pt_min, pt_max;
    int nbins;
    TString ntuple_name = "", xaxis_title = "";
    
    pt_min=5;
    pt_max=200;

    switch (channel) {
        case 1:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = BP_MASS;
            nbins = 50;
	    ntuple_name = "ntkp";
            xaxis_title = "M_{J/#psi K^{#pm}} [GeV]";
            break;
        case 2:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = B0_MASS;
            nbins = 50;
            ntuple_name = "ntkstar";
            xaxis_title = "M_{J/#psi K^{#pm}#pi^{#mp}} [GeV]";
            break;
        case 3:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = B0_MASS;
            nbins = 50;
            ntuple_name = "ntks";
            xaxis_title = "M_{J/#psi K^{0}_{S}} [GeV]";
            break;
        case 4:
            mass_min = 5.0; mass_max = 6.0;
            mass_peak = BS_MASS;
            nbins = 50;
            ntuple_name = "ntphi";
            xaxis_title = "M_{J/#psi K^{#pm}K^{#mp}} [GeV]";
            break;
        case 5:
            mass_min = 3.6; mass_max = 4.0;
            mass_peak = PSI2S_MASS;
            nbins = 80;
            ntuple_name = "ntmix";
            xaxis_title = "M_{J/#psi #pi^{#pm}#pi^{#mp}} [GeV]";
            break;
        case 6:
            mass_min = 5.3; mass_max = 6.3;
            mass_peak = LAMBDAB_MASS;
            nbins = 50;
            ntuple_name = "ntlambda";
            xaxis_title = "M_{J/#psi #Lambda} [GeV]";
            break;
    }
    
    RooRealVar mass("mass","mass",mass_min,mass_max);
    RooRealVar pt("pt","pt",pt_min,pt_max);    

    TFile *fout = new TFile("myfitter_"+ntuple_name+"_v1.root","recreate");
    TNtupleD *_nt = new TNtupleD("_nt","_nt","mass");
    TNtupleD *_nt2 = new TNtupleD("_nt2","_nt2","pt");
    
    TFile *fin = new TFile(SOURCE);
    TTree *tin = (TTree*)fin->Get(ntuple_name);
    
    ReducedBranches br;
    br.setbranchadd(tin);
    
    int n_br_queued = 0;
    ReducedBranches br_queue[32];
    
    for (int evt=0;evt<tin->GetEntries();evt++) {
        tin->GetEntry(evt);
        
        if (channel==1) { // cuts for B+ -> J/psi K+
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.tk1pt<=1.6) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            
            _nt->Fill(&br.mass);
	    _nt2->Fill(&br.pt);
        }else
        if (channel==2) { // cuts for B0 -> J/psi K*
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-KSTAR_MASS)>=0.05) continue;
            
            TLorentzVector v4_tk1, v4_tk2;
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
		      _nt->Fill(&br_queue[i].mass);
		      _nt2->Fill(&br_queue[i].pt);
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
                
	    _nt->Fill(&br.mass);
	    _nt2->Fill(&br.pt);
        }else
        if (channel==4) { // cuts for Bs -> J/psi phi
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-PHI_MASS)>=0.010) continue;
            
            TLorentzVector v4_tk1, v4_tk2;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,KAON_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,KAON_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSTAR_MASS)<=0.05) continue;
            
	    _nt->Fill(&br.mass);
	    _nt2->Fill(&br.pt);
        }else
        if (channel==5) { // cuts for psi(2S)/X(3872) -> J/psi pipi
            if (br.vtxprob<=0.2) continue;
            if (fabs(br.tk1eta)>=1.6) continue;
            if (fabs(br.tk2eta)>=1.6) continue;
            
	    _nt->Fill(&br.mass);
	    _nt2->Fill(&br.pt);
        }else
        if (channel==6) {
            if (br.hltbook[HLT_DoubleMu4_JpsiTrk_Displaced_v2]!=1) continue;
            if (br.vtxprob<=0.1) continue;
            if (br.lxy/br.errxy<=3.0) continue;
            if (br.tktkblxy/br.tktkberrxy<=3.0) continue;
            if (br.cosalpha2d<=0.99) continue;
            if (fabs(br.tktkmass-LAMBDA_MASS)>=0.015) continue;
            
            TLorentzVector v4_tk1, v4_tk2;
            v4_tk1.SetPtEtaPhiM(br.tk1pt,br.tk1eta,br.tk1phi,PION_MASS);
            v4_tk2.SetPtEtaPhiM(br.tk2pt,br.tk2eta,br.tk2phi,PION_MASS);
            if (fabs((v4_tk1+v4_tk2).Mag()-KSHORT_MASS)<=0.015) continue;
            
            
            _nt->Fill(&br.mass);
	    _nt2->Fill(&br.pt);
        }
    }
    fin->Close();
    
    RooDataSet *data = new RooDataSet("data","data",_nt,RooArgSet(mass));
    
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
    // B->J/psi+X background PDF, only for B+ -> J/psi K+
    
    RooRealVar m_jpsix_mean("m_jpsix_mean","m_jpsix_mean",5.07,5.0,5.2);
    RooRealVar m_jpsix_sigma("m_jpsix_sigma","m_jpsix_sigma",0.05,0.01,0.10);
    RooGaussian pdf_m_jpsix("pdf_m_jpsix","pdf_m_jpsix",mass,m_jpsix_mean,m_jpsix_sigma);

    //-----------------------------------------------------------------
    // X(3872) PDF, only for J/psi pipi fit
    
    RooRealVar m_x3872_mean("m_x3872_mean","m_x3872_mean",3.872,3.7,3.9);
    RooRealVar m_x3872_sigma("m_x3872_sigma","m_x3872_sigma",0.01,0.001,0.010);
    RooGaussian pdf_m_x3872("pdf_m_x3872","pdf_m_x3872",mass,m_x3872_mean,m_x3872_sigma);
    
    //-----------------------------------------------------------------
    // full model
    
    RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
    RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());
    RooRealVar n_jpsix("n_jpsix","n_jpsix",200.,0.,data->sumEntries());
    RooRealVar n_x3872("n_x3872","n_x3872",200.,0.,data->sumEntries());
    
    RooAddPdf *model;
    
    if (channel==1) // B+ -> J/psi K+
        model = new RooAddPdf("model","model",
                    RooArgList(pdf_m_signal, pdf_m_combinatorial_exp, pdf_m_jpsix),
                    RooArgList(n_signal, n_combinatorial, n_jpsix));
    else
    if (channel==2 || channel==3 || channel==4 || channel==6) // B0 -> J/psi K*, Ks; Bs -> J/psi phi; Lambda_b -> J/psi Lambda
        model = new RooAddPdf("model","model",
                    RooArgList(pdf_m_signal, pdf_m_combinatorial_exp),
                    RooArgList(n_signal, n_combinatorial));
    else
    if (channel==5) // J/psi pipi
        model = new RooAddPdf("model","model",
                    RooArgList(pdf_m_signal, pdf_m_combinatorial_bern, pdf_m_x3872),
                    RooArgList(n_signal, n_combinatorial, n_x3872));
    
    
    model->fitTo(*data,Minos(DO_MINOS),NumCPU(NUMBER_OF_CPU),Offset(kTRUE));

#if DISPLAY
    TCanvas *c1 = canvasDressing("c1");

    // Display
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
    model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_signal),LineColor(kRed),LineWidth(2),LineStyle(kSolid),FillStyle(3008),FillColor(kRed), VLines(), DrawOption("F"));

    if (channel==1 || channel==2 || channel==3 || channel==4 || channel==6)
        model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_combinatorial_exp),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    else
        model->plotOn(frame_m,Precision(2E-4),Components(pdf_m_combinatorial_bern),LineColor(kCyan+1),LineWidth(2),LineStyle(2));
    
    if (channel==1)
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
    
    if (channel==1)
        LegendChannelOne();
    if (channel==2)
        LegendChannelTwo();
    if (channel==4)
        LegendChannelFour();
    if (channel==5)
        LegendChannelFive();

    TLatex *   tex = new TLatex(5.4,1343,TString::Format("Signal: %f", n_signal.getVal()));
    tex->SetLineWidth(2);
    tex->Draw();

    Printf("Signal: %f  +%f %f",n_signal.getVal(),n_signal.getAsymErrorHi(),n_signal.getAsymErrorLo());
    
#endif
    // c1->Print(ntuple_name+"_v1.gif");
    c1->SaveAs("mass_"+ntuple_name+"_v1.root");
    //c1->SaveAs("pt_"+ntuple_name+"_v1.gif");
    fout->Write();
    fout->Close();
    
}
