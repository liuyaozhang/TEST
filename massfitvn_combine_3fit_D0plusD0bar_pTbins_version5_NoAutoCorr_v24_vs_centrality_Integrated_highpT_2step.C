#include <iostream>

#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TString.h"

#include <vector>

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

#include "TProfile.h"
#include "TH2.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"

//gROOT->Reset();

int iparmassfit_poly3bkg_floatwidth[13] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
int iparvnfit1_poly3bkg_floatwidth[16] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

struct GlobalChi2_poly3bkg_floatwidth {
    GlobalChi2_poly3bkg_floatwidth(ROOT::Math::IMultiGenFunction & f1,
                                   ROOT::Math::IMultiGenFunction & f2) :
    fChi2_1(&f1), fChi2_2(&f2) {}
    
    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)
    double operator() (const double *par) const {
        double p1[13];
        for(int i = 0; i < 13; ++i) p1[i] = par[iparmassfit_poly3bkg_floatwidth[i]];

        double p2[16];
        for(int i = 0; i < 16; ++i) p2[i] = par[iparvnfit1_poly3bkg_floatwidth[i]];
        
        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }
    
    const  ROOT::Math::IMultiGenFunction * fChi2_1;
    const  ROOT::Math::IMultiGenFunction * fChi2_2;
};

TString detector[2] = {"tracker", "HF"};
TString range_ybin[2]   = {"|y|<1", "1<|y|<2"};
TString labelN_ybin[2] = {"1T1","1T2"};
TString label_ybin[2]  = {"y-1.0to1.0", "y1.0to2.0"};
TString label_pTselection[2]  = {"pT2_8GeV", "pT8_15GeV"};
TString label_centbin[5] = {"cent0to10", "cent10to20", "cent20to30", "cent30to50", "cent50to70"}; 
TString range_centbin[5] = {"cent. 0-10%", "cent. 10-20%", "cent. 20-30%", "cent. 30-50%", "cent. 50-70%"};

double centralitybin[6]= {0., 10., 20., 30.0, 50, 70};
double centbin[5]={5., 15., 25., 40., 60.};
double centbinwidth[5]={5., 5., 5., 10., 10.};

double momentumbin[3]= {8.0, 10, 15.0};
const Int_t n_pTbin = sizeof(momentumbin)/sizeof(double) - 2;

TString label_pTbin[n_pTbin+1]={"pT8.0to10.0", "pT10.0to15.0"};
const Int_t n_centbin = sizeof(centralitybin)/sizeof(double)-1;; 

Int_t number_detector=1; //detector
Int_t number_pTselection=1;
Int_t number_y=0;//rapidity

Bool_t isV24tracker;
Bool_t isV24HF;

TString folderForPlots = "f_update_build17_ExactCumulant_"+detector[number_detector]+"_NoAutoCorr_update_Jun30/v2_vs_cent_2step";
TString label2_ybin = "abseta"+labelN_ybin[number_y]+"_sumw2_"+label_pTselection[number_pTselection];

void massfitvn_combine_3fit_D0plusD0bar_pTbins_version5_NoAutoCorr_v24_vs_centrality_Integrated_highpT_2step(Double_t vnYrangeMin = -0.2, Double_t vnYrangeMax = 0.6, Bool_t isForward=false)
{   


    cout << "==================" << n_centbin << endl; 
    cout << "==================" << n_pTbin << endl; 
    TH1::StatOverflows(kTRUE);
    if(number_detector==0){
    isV24tracker=true;
    isV24HF=false;
    }
  else{
    isV24tracker=false;
    isV24HF=true;
    }

    gSystem->Exec("mkdir "+folderForPlots);
    gSystem->Exec("mkdir "+folderForPlots+"/v24_"+label2_ybin);
    
    TString folder_rootfiles = "";
    TString folder_plots = "";
    TString plots_label1 = "";
    TString plots_label2 = "";
    TString plots_label3 = "";

    folder_rootfiles = "v24_"+label2_ybin+"/v24vspt_"+label2_ybin+".root";
    folder_plots = folderForPlots+"/v24_"+label2_ybin+"/plot_";
    plots_label1 = "v_{2}^{sig}{4} = %.3f #pm %.3f";
    plots_label2 = "v_{2}{4} ";
    plots_label3 = "v_{2}{4} (D^{0} + #bar{D^{0}}) ";
    
    TString ofName; 
    TString ofpara; 
    if(isV24tracker) ofName = "v24_"+label_ybin[number_y]+"_integrated_sumw2_tracker"+".csv";
    if(isV24HF)      ofName = "v24_"+label_ybin[number_y]+"_integrated_sumw2_HF"+".csv";
    ofstream ofout(ofName);  
    
    double fit_range_low = 1.74;
    double fit_range_high = 2.0;
    double D0_mass = 1.8648;
    
    double vn[n_centbin];
    double vne[n_centbin]={0.};
    
    TF1* fmasssig[n_pTbin];
    TF1* fmassswap[n_pTbin];
    TF1* fmassbkg[n_pTbin];
    TF1* fmasstotal[n_pTbin];
    TF1* fvn[n_pTbin];

    TFile ofile(folderForPlots+"/"+folder_rootfiles,"RECREATE"); 
    TCanvas* c1[n_centbin][n_pTbin];
    for(int icent=0;icent<n_centbin;icent++){
     for(int i=0;i<n_pTbin;i++){
    
        c1[icent][i] = new TCanvas(Form("c1_%d_%d",icent, i),Form("c1_%d_%d",icent, i),800,400);
        c1[icent][i]->Divide(2,1);
        }
       }
    for(int icent=0;icent<n_centbin;icent++){
     for(int i=0;i<n_pTbin;i++){
        c1[icent][i]->cd(1)->SetTopMargin(0.06);
        c1[icent][i]->cd(1)->SetLeftMargin(0.18);
        c1[icent][i]->cd(1)->SetRightMargin(0.043);
        c1[icent][i]->cd(1)->SetBottomMargin(0.145);
        c1[icent][i]->cd(2)->SetTopMargin(0.06);
        c1[icent][i]->cd(2)->SetLeftMargin(0.18);
        c1[icent][i]->cd(2)->SetRightMargin(0.043);
        c1[icent][i]->cd(2)->SetBottomMargin(0.145);
        }
       }
    
    TLatex* texCMS = new TLatex;
    texCMS->SetNDC();
    texCMS->SetTextFont(42);
    texCMS->SetTextSize(0.05);
    texCMS->SetTextAlign(12);
     
    for(Int_t icent=0; icent<n_centbin; icent++){//icent
    
    Double_t sumyvn_icent=0; 
    Double_t sumyvne_icent=0; 
    Double_t sumy_icent=0;      

    TFile* file0 = TFile::Open("/afs/cern.ch/user/z/zhangli/public_work/Scalar-product/MC/histograms_3Dcuts_MC_default_differential_pT_cent0_10_20_30_50_70_new_5pTbins_highpT15.root");
    TFile* file1 = TFile::Open("/afs/cern.ch/user/z/zhangli/public_work/Scalar-product/batch/TestOutPutDir/data_build17_AverageParticle_ExactCumulant_MVA_selection_MC_new_Centtable_Jun30/v2_vs_cent/differential_pT_cent0_10_20_30_50_70_sumw2_5pTbins_highpT15/PbPbMB.ExactCumulant_total.root");
    ofile.cd(); 
    
    TLatex* tex = new TLatex;
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.045);
    tex->SetLineWidth(2);
 
    
    Float_t massbinning[14]={1.74,1.78,1.80,1.82,1.84,1.85,1.86,1.865,1.87,1.88,1.90,1.92,1.96,2.00};
    
    TH1D* hist111_check=((TH1D*)(gROOT->FindObject("hist111")));
    delete hist111_check;     

    TH1D* hist111 = new TH1D("hist111","",13,massbinning);
    hist111->SetLineWidth(0);
    hist111->SetStats(kFALSE);
    hist111->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
    hist111->GetYaxis()->SetTitle(plots_label3);
    hist111->GetXaxis()->CenterTitle();
    hist111->GetYaxis()->CenterTitle();
    hist111->GetXaxis()->SetTitleOffset(1.2);
    hist111->GetYaxis()->SetTitleOffset(1.4);
    hist111->GetXaxis()->SetLabelOffset(0.007);
    hist111->GetYaxis()->SetLabelOffset(0.007);
    hist111->GetXaxis()->SetTitleSize(0.047);
    hist111->GetYaxis()->SetTitleSize(0.055);
    hist111->GetXaxis()->SetTitleFont(42);
    hist111->GetYaxis()->SetTitleFont(42);
    hist111->GetXaxis()->SetLabelFont(42);
    hist111->GetYaxis()->SetLabelFont(42);
    hist111->GetXaxis()->SetLabelSize(0.04);
    hist111->GetYaxis()->SetLabelSize(0.04);
    hist111->SetMinimum(vnYrangeMin);
    hist111->SetMaximum(vnYrangeMax);


    // c2->cd();
    // hist->Draw();
    
    for(int i=0;i<n_pTbin;i++)//pTbin
    {

        TH1D* h_mc_match_signal;
        TH1D* h_mc_match_signal1;
        TH1D* h_mc_match_all; 
        TH1D* h_mc_match_all1; 
        if(i>6){
          h_mc_match_signal = (TH1D*)file0->Get("hist_D0_mass_mc_signal_"+label_pTbin[6]+"_"+label_ybin[number_y]+"_"+label_centbin[icent]);
	  h_mc_match_all = (TH1D*)file0->Get("hist_D0_mass_mc_signalPlusSwap_"+label_pTbin[6]+"_"+label_ybin[number_y]+"_"+label_centbin[icent]); 
	 if(isForward){
            TH1D* aux_h_mc_match_signal = (TH1D*)file0->Get("hist_D0_mass_mc_signal_"+label_pTbin[6]+"_y-2.0to-1.0_"+label_centbin[icent]);
            TH1D* aux_h_mc_match_all = (TH1D*)file0->Get("hist_D0_mass_mc_signalPlusSwap_"+label_pTbin[6]+"_y-2.0to-1.0_"+label_centbin[icent]);
            h_mc_match_signal->Sumw2();
            h_mc_match_all->Sumw2();
            h_mc_match_signal->Add(aux_h_mc_match_signal);
            h_mc_match_all->Add(aux_h_mc_match_all);
          }
        }else{
          h_mc_match_signal = (TH1D*)file0->Get("hist_D0_mass_mc_signal_"+label_pTbin[i]+"_"+label_ybin[number_y]+"_"+label_centbin[icent]);
          h_mc_match_signal1 = (TH1D*)file0->Get("hist_D0_mass_mc_signal_"+label_pTbin[i+1]+"_"+label_ybin[number_y]+"_"+label_centbin[icent]);
          h_mc_match_all = (TH1D*)file0->Get("hist_D0_mass_mc_signalPlusSwap_"+label_pTbin[i]+"_"+label_ybin[number_y]+"_"+label_centbin[icent]);
          h_mc_match_all1 = (TH1D*)file0->Get("hist_D0_mass_mc_signalPlusSwap_"+label_pTbin[i+1]+"_"+label_ybin[number_y]+"_"+label_centbin[icent]);
	  if(isForward){
            TH1D* aux_h_mc_match_signal = (TH1D*)file0->Get("hist_D0_mass_mc_signal_"+label_pTbin[i]+"_y-2.0to-1.0_"+label_centbin[icent]);
            TH1D* aux_h_mc_match_all = (TH1D*)file0->Get("hist_D0_mass_mc_signalPlusSwap_"+label_pTbin[i]+"_y-2.0to-1.0_"+label_centbin[icent]);
            h_mc_match_signal->Sumw2();
            h_mc_match_all->Sumw2();
            h_mc_match_signal->Add(aux_h_mc_match_signal);
            h_mc_match_all->Add(aux_h_mc_match_all);
          }
          h_mc_match_signal->Sumw2();
          h_mc_match_signal1->Sumw2();
          h_mc_match_signal->Add(h_mc_match_signal1);
          cout << "h_mc_match_signal" << " " << h_mc_match_signal->GetEntries() << endl;
          
          h_mc_match_all->Sumw2();
          h_mc_match_all1->Sumw2();
          h_mc_match_all->Add(h_mc_match_all1);
          cout << "h_mc_match_all" << " " << h_mc_match_all->GetEntries() << endl;
         }
        
        TH1D* h_data = (TH1D*)file1->Get("hist_D0_mass_"+label_centbin[icent]+"_"+label_pTbin[i]+"_"+label_ybin[number_y]);
        TH1D* h_data1 = (TH1D*)file1->Get("hist_D0_mass_"+label_centbin[icent]+"_"+label_pTbin[i+1]+"_"+label_ybin[number_y]);
        TH1D* h_dataD0bar = (TH1D*)file1->Get("hist_D0bar_mass_"+label_centbin[icent]+"_"+label_pTbin[i]+"_"+label_ybin[number_y]);
        TH1D* h_dataD0bar1 = (TH1D*)file1->Get("hist_D0bar_mass_"+label_centbin[icent]+"_"+label_pTbin[i+1]+"_"+label_ybin[number_y]);
        if(isForward){//should also include negative side, since in script can only use one side - to be improved in the future...
          TH1D* aux_h_data = (TH1D*)file1->Get("hist_D0_mass_"+label_centbin[icent]+"_"+label_pTbin[i]+"_y-2.0to-1.0");
          TH1D* aux_h_dataD0bar = (TH1D*)file1->Get("hist_D0bar_mass_"+label_centbin[icent]+"_"+label_pTbin[i]+"_y-2.0to-1.0");
          h_data->Sumw2();
          h_dataD0bar->Sumw2();
          h_data->Add(aux_h_data);
          h_dataD0bar->Add(aux_h_dataD0bar);
        }
        
        h_data->Sumw2();
        h_data1->Sumw2();
        h_dataD0bar->Sumw2();
        h_dataD0bar1->Sumw2();

        h_data->Add(h_data1);
        h_data->Add(h_dataD0bar);
        h_data->Add(h_dataD0bar1);

        h_data->SetMinimum(0);
        h_data->SetMarkerSize(0.8);
        h_data->SetMarkerStyle(20);
        h_data->SetLineWidth(1);
        h_data->SetOption("e");
        h_data->GetXaxis()->SetRangeUser(1.74,2.0);
        h_data->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
        h_data->GetYaxis()->SetTitle("Entries / 2.5 MeV");
        h_data->GetXaxis()->CenterTitle();
        h_data->GetYaxis()->CenterTitle();
        h_data->GetXaxis()->SetTitleOffset(1.3);
        h_data->GetYaxis()->SetTitleOffset(2);
        h_data->GetXaxis()->SetLabelOffset(0.007);
        h_data->GetYaxis()->SetLabelOffset(0.007);
        h_data->GetXaxis()->SetTitleSize(0.045);
        h_data->GetYaxis()->SetTitleSize(0.045);
        h_data->GetXaxis()->SetTitleFont(42);
        h_data->GetYaxis()->SetTitleFont(42);
        h_data->GetXaxis()->SetLabelFont(42);
        h_data->GetYaxis()->SetLabelFont(42);
        h_data->GetXaxis()->SetLabelSize(0.04);
        h_data->GetYaxis()->SetLabelSize(0.04);
        h_data->SetStats(kFALSE);
        
        h_data->GetXaxis()->SetNoExponent(true);
        ((TGaxis*)h_data->GetXaxis())->SetMaxDigits(7);
        
        h_data->SetMaximum(h_data->GetMaximum()*1.5);

        c1[icent][i]->cd(1);
        //The full fitting function is constructed as follow
        //[0] is signal + swap yield;
        //[1] is common mean of double gaussian;
        //[2] is signal gaussian 1 sigma;
        //[3] is signal gaussian 2 sigma;
        //[4] is fractional signal gaussian 1 yield; 1-[4] is fractional signal gaussian 2 yield;
        //[5] is fractional double gaussian signal yield, 1-[5] is fractional swap yield;
        //[6] is a factor to let width of the gaussians to vary in data;
        //[7] is swap gaussian sigma;
        //[8] is swap gaussian mean;
        //[9-12] is 3rd order poly parameters
        
        
	TF1* f = new TF1(Form("f_%d_%d", icent, i),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x", fit_range_low, fit_range_high);
 
        f->SetLineColor(2);
        f->SetLineWidth(1);

        //first fit MC signal, swap and poly bkg set to 0
        
        f->SetParameter(0,100.);
        f->SetParameter(1,D0_mass);
        f->SetParameter(2,0.03);
        f->SetParameter(3,0.005);
        f->SetParameter(4,0.1);
        
        f->FixParameter(5,1);
        f->FixParameter(6,0); //always 0 in MC
        f->FixParameter(7,0.1); //does not really mater here as yield is fix to 0
        f->FixParameter(8,D0_mass); //does not really mater here as yield is fix to 0
        f->FixParameter(9,0);
        f->FixParameter(10,0);
        f->FixParameter(11,0);
        f->FixParameter(12,0);
       
        f->SetParLimits(2,0.01,0.1);//0.1
        f->SetParLimits(3,0.003,0.05);//0.05
        f->SetParLimits(4,0,1);//1,0
       
        f->FixParameter(1,1.8648); //for first few attempt fix mean of gaussian to get reasonable estimation of other pars; later open it up
        h_mc_match_signal->Fit(Form("f_%d_%d",icent, i),"q","",fit_range_low,fit_range_high);
        h_mc_match_signal->Fit(Form("f_%d_%d",icent, i),"q","",fit_range_low,fit_range_high);
        f->ReleaseParameter(1);
        h_mc_match_signal->Fit(Form("f_%d_%d",icent, i),"L q","",fit_range_low,fit_range_high);
        h_mc_match_signal->Fit(Form("f_%d_%d",icent, i),"L q","",fit_range_low,fit_range_high);
        h_mc_match_signal->Fit(Form("f_%d_%d",icent, i),"L q","",fit_range_low,fit_range_high);
        h_mc_match_signal->Fit(Form("f_%d_%d",icent, i),"L m","",fit_range_low,fit_range_high);
        
        //now fix signal double gaussian mean, sigma and gaus1,gaus2 yield ratio
        f->FixParameter(1,f->GetParameter(1));
        f->FixParameter(2,f->GetParameter(2));
        f->FixParameter(3,f->GetParameter(3));
        f->FixParameter(4,f->GetParameter(4));
        
        //now release swap bkg parameters to fit signal+swap MC
        f->ReleaseParameter(5);
        f->SetParLimits(5,0.,1);//0, 1
        f->ReleaseParameter(7);
        f->ReleaseParameter(8);
        
        f->SetParameter(7,0.1);
        f->SetParameter(8,D0_mass);
        
        //fit signal+swap MC
        h_mc_match_all->Fit(Form("f_%d_%d",icent, i),"L q","",fit_range_low,fit_range_high);
        h_mc_match_all->Fit(Form("f_%d_%d",icent, i),"L q","",fit_range_low,fit_range_high);
        h_mc_match_all->Fit(Form("f_%d_%d",icent, i),"L q","",fit_range_low,fit_range_high);
        h_mc_match_all->Fit(Form("f_%d_%d",icent, i),"L m","",fit_range_low,fit_range_high);
        
        //now fix swap bkg parameters to fit data
        f->FixParameter(5,f->GetParameter(5));
        f->FixParameter(7,f->GetParameter(7));
        f->FixParameter(8,f->GetParameter(8));
        
        //now release poly bkg pars
        f->ReleaseParameter(9);
        f->ReleaseParameter(10);
        f->ReleaseParameter(11);
        f->ReleaseParameter(12);
        
        f->SetParameter(9, 5747.53);
        f->SetParameter(10, -2811.87);
        f->SetParameter(11, -1516.12);
        f->SetParameter(12, 764.207);
            
        //now fit data
        h_data->Fit(Form("f_%d_%d",icent, i),"q","",fit_range_low,fit_range_high);
        h_data->Fit(Form("f_%d_%d",icent, i),"q","",fit_range_low,fit_range_high);
          
        cout << "par9" << "=" << f->GetParameter(9) << endl;
        cout << "par10" << "=" << f->GetParameter(10) << endl;
        cout << "par11" << "=" << f->GetParameter(11) << endl;
        cout << "par12" << "=" << f->GetParameter(12) << endl;

        f->ReleaseParameter(1); //allow data to have different mass peak mean than MC
        f->ReleaseParameter(6); //allow data to have different peak width than MC
       
        //f->SetParameter(1, D0_mass);//Liuyao
        f->SetParameter(6,0);
        f->SetParLimits(6,-1,1);
        h_data->Fit(Form("f_%d_%d",icent, i),"L q","",fit_range_low,fit_range_high);
        h_data->Fit(Form("f_%d_%d",icent, i),"L q","",fit_range_low,fit_range_high);
        h_data->Fit(Form("f_%d_%d",icent, i),"L q","",fit_range_low,fit_range_high);
        h_data->Fit(Form("f_%d_%d",icent, i),"L m","",fit_range_low,fit_range_high);
        
	h_data->GetFunction(Form("f_%d_%d",icent, i))->SetBit(TF1::kNotDraw); //cesar
        
        //draw D0 signal separately
        TF1* f1 = new TF1(Form("f_sig_%d_%d",icent, i),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
        f1->SetLineColor(kOrange-3);
        f1->SetLineWidth(1);
        f1->SetLineStyle(2);
        f1->SetFillColorAlpha(kOrange-3,0.3);
        f1->SetFillStyle(1001);
        f1->FixParameter(0,f->GetParameter(0));
        f1->FixParameter(1,f->GetParameter(1));
        f1->FixParameter(2,f->GetParameter(2));
        f1->FixParameter(3,f->GetParameter(3));
        f1->FixParameter(4,f->GetParameter(4));
        f1->FixParameter(5,f->GetParameter(5));
        f1->FixParameter(6,f->GetParameter(6));
        
        fmasssig[i] = (TF1*)f1->Clone();
        fmasssig[i]->SetName(Form("masssigfcn_%d_%d",icent, i));
        fmasssig[i]->Write();
        
        f1->Draw("LSAME");
         
        //draw swap separately
        ///TF1* f2 = new TF1(Form("f_swap_%d",i),"[0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))", fit_range_low, fit_range_high);
        TF1* f2 = new TF1(Form("f_swap_%d_%d",icent, i),"[0]*((1-[1])*TMath::Gaus(x,[4],[3]*(1.0 +[2]))/(sqrt(2*3.14159)*[3]*(1.0 +[2])))", fit_range_low, fit_range_high);
        f2->SetLineColor(kGreen+4);
        f2->SetLineWidth(1);
        f2->SetLineStyle(1);
        f2->SetFillColorAlpha(kGreen+4,0.3);
        f2->SetFillStyle(1001);
        f2->FixParameter(0,f->GetParameter(0));
        f2->FixParameter(1,f->GetParameter(5));
        f2->FixParameter(2,f->GetParameter(6));
        f2->FixParameter(3,f->GetParameter(7));
        f2->FixParameter(4,f->GetParameter(8));
        
        fmassswap[i] = (TF1*)f2->Clone();
        fmassswap[i]->SetName(Form("massswapfcn_%d_%d",icent, i));
        fmassswap[i]->Write();
        
        f2->Draw("LSAME");
    
        //draw poly bkg separately
        TF1* f3 = new TF1(Form("f_bkg_%d_%d",icent, i),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", fit_range_low, fit_range_high);
        f3->SetLineColor(4);
        f3->SetLineWidth(1);
        f3->SetLineStyle(2);
        f3->FixParameter(0,f->GetParameter(9));
        f3->FixParameter(1,f->GetParameter(10));
        f3->FixParameter(2,f->GetParameter(11));
        f3->FixParameter(3,f->GetParameter(12));
        
        fmassbkg[i] = (TF1*)f3->Clone();
        fmassbkg[i]->SetName(Form("massbkgfcn_%d_%d",icent,i));
        fmassbkg[i]->Write();
        
        f3->Draw("LSAME");
        
        tex->DrawLatex(0.22,0.86,range_centbin[icent]);
        tex->DrawLatex(0.22,0.80,Form("%.1f < p_{T} < %.1f GeV/c",momentumbin[i],momentumbin[i+1]));
        tex->DrawLatex(0.22,0.74,range_ybin[number_y]);
        
        texCMS->DrawLatex(.18,.97,"#font[61]{CMS} #it{Preliminary}");
        texCMS->DrawLatex(0.62,0.97, "#scale[0.8]{PbPb #sqrt{s_{NN}} = 5.02 TeV}");
        
        TLegend* leg = new TLegend(0.65,0.58,0.81,0.9,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetTextSize(0.045);
        leg->SetTextFont(42);
        leg->SetFillStyle(0);
        leg->AddEntry(h_data," Data","p");
        leg->AddEntry(f," Fit","L");
        leg->AddEntry(f1," D^{0}+#bar{D^{#lower[0.2]{0}}} Signal","f");
        leg->AddEntry(f2," K-#pi swap","f");
        leg->AddEntry(f3," Combinatorial","l");
        leg->Draw("SAME");
         
        
        //c->Print(Form("plots/massfit_pt%d.pdf",i));
        
        //fit vn
        //[13] id vn (D0 + D0bar)
        ///vn denominator - uncertainty is negligible for v1,v2, and v3
        Double_t vn_den = -99999999.9;
        Double_t v24_ref = 0; 
          
        if(isV24tracker){
          TH1D* hQ2Q2Q2Q2_TrkmTrkmTrkpTrkp = (TH1D*)file1->Get("hist_AveNoAutoCorr_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_Re_"+label_centbin[icent]);
          TH1D* h_product_4Q2 = (TH1D*) hQ2Q2Q2Q2_TrkmTrkmTrkpTrkp->Clone("h_product_4Q2");
          TH1D* hQ2Q2_TrkmTrkp = (TH1D*)file1->Get("hist_Q2Q2_TrkmTrkp_Re_"+label_centbin[icent]);

                
          Double_t aux_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_mean = h_product_4Q2->GetMean();
          Double_t aux_Q2Q2_TrkmTrkp_mean = hQ2Q2_TrkmTrkp->GetMean();
          vn_den = TMath::Power(fabs(aux_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_mean - 2.0*aux_Q2Q2_TrkmTrkp_mean*aux_Q2Q2_TrkmTrkp_mean),3./4.);
          cout << "======v24_ref======:"<<TMath::Power(fabs(aux_Q2Q2Q2Q2_TrkmTrkmTrkpTrkp_mean - 2.0*aux_Q2Q2_TrkmTrkp_mean*aux_Q2Q2_TrkmTrkp_mean),1./4.)<< endl; 
         }
        else if(isV24HF){
          TH1D* hQ2Q2Q2Q2_HFmHFmHFpHFp = (TH1D*)file1->Get("hist_AveNoAutoCorr_Q2Q2Q2Q2_HFmHFmHFpHFp_Re_"+label_centbin[icent]);
          TH1D* h_product_4Q2 = (TH1D*) hQ2Q2Q2Q2_HFmHFmHFpHFp->Clone("h_product_4Q2");
          TH1D* hQ2Q2_HFmHFp = (TH1D*)file1->Get("hist_Q2Q2_HFmHFp_Re_"+label_centbin[icent]);      
           
          Double_t aux_Q2Q2Q2Q2_HFmHFmHFpHFp_mean = h_product_4Q2->GetMean();
          Double_t aux_Q2Q2_HFmHFp_mean = hQ2Q2_HFmHFp->GetMean();
          vn_den = TMath::Power(fabs(aux_Q2Q2Q2Q2_HFmHFmHFpHFp_mean - 2.0*aux_Q2Q2_HFmHFp_mean*aux_Q2Q2_HFmHFp_mean),3./4.);
          cout << "=======v24_ref====:"<<TMath::Power(fabs(aux_Q2Q2Q2Q2_HFmHFmHFpHFp_mean - 2.0*aux_Q2Q2_HFmHFp_mean*aux_Q2Q2_HFmHFp_mean),1./4.)<<endl;
        }
        
       
        ///Delta vn numerator
        //D0
        TProfile *tp_vn_data;
        TH1D *vn_data;
        TH1D *hvn_data_Q2Q2;
        if(isV24tracker){
          TH2D* hQ2Q2Q2Q2_D0plusD0barTrkTrkTrk = (TH2D*)file1->Get("hist_AveNoAutoCorr_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_"+label_centbin[icent]+"_"+label_pTbin[i]+"_"+label_ybin[number_y]);
          TH2D* hQ2Q2Q2Q2_D0plusD0barTrkTrkTrk1 = (TH2D*)file1->Get("hist_AveNoAutoCorr_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_"+label_centbin[icent]+"_"+label_pTbin[i+1]+"_"+label_ybin[number_y]);
          hQ2Q2Q2Q2_D0plusD0barTrkTrkTrk->Add(hQ2Q2Q2Q2_D0plusD0barTrkTrkTrk1);
          
          TH2D* hQ2Q2_D0plusD0barTrk = (TH2D*)file1->Get("hist_Q2Q2_D0plusD0barTrk_Re_"+label_centbin[icent]+"_"+label_pTbin[i]+"_"+label_ybin[number_y]);
          TH2D* hQ2Q2_D0plusD0barTrk1 = (TH2D*)file1->Get("hist_Q2Q2_D0plusD0barTrk_Re_"+label_centbin[icent]+"_"+label_pTbin[i+1]+"_"+label_ybin[number_y]);
          hQ2Q2_D0plusD0barTrk->Add(hQ2Q2_D0plusD0barTrk1);

          if(isForward){
            TH2D* aux_hQ2Q2Q2Q2_D0plusD0barTrkTrkTrk = (TH2D*)file1->Get("hist_AveNoAutoCorr_Q2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_"+label_centbin[icent]+"_"+label_pTbin[i]+"_y-2.0to-1.0");
            TH2D* aux_hQ2Q2_D0plusD0barTrk = (TH2D*)file1->Get("hist_Q2Q2_D0plusD0barTrk_Re_"+label_centbin[icent]+"_"+label_pTbin[i]+"_y-2.0to-1.0");
	    hQ2Q2Q2Q2_D0plusD0barTrkTrkTrk->Sumw2();
            hQ2Q2Q2Q2_D0plusD0barTrkTrkTrk->Add(aux_hQ2Q2Q2Q2_D0plusD0barTrkTrkTrk);
             
            hQ2Q2_D0plusD0barTrk->Sumw2();
            hQ2Q2_D0plusD0barTrk->Add(aux_hQ2Q2_D0plusD0barTrk);
          }
          TH2D* hproduct_4D0plusD0bar  = (TH2D*)hQ2Q2Q2Q2_D0plusD0barTrkTrkTrk->Clone();
          TProfile* tpQ2Q2Q2Q2_D0plusD0barTrkTrkTrk = (TProfile*)hproduct_4D0plusD0bar->ProfileX(Form("tpQ2Q2Q2Q2_D0plusD0barTrkTrkTrk_Re_%d", icent)+label_centbin[icent]+"_"+label_pTbin[i]+"_"+label_ybin[number_y],1,-1,"i");
	  vn_data = tpQ2Q2Q2Q2_D0plusD0barTrkTrkTrk->ProjectionX("");
	  TProfile* tpQ2Q2_D0plusD0barTrk = (TProfile *)hQ2Q2_D0plusD0barTrk->ProfileX(Form("tpQ2Q2_D0plusD0barTrk_Re_%d", icent)+label_centbin[icent]+"_"+label_pTbin[i]+"_"+label_ybin[number_y],1,-1,"i");
	  TH1D *hQ2Q2_TrkmTrkp = (TH1D*)file1->Get("hist_Q2Q2_TrkmTrkp_Re_"+label_centbin[icent]);
	  Double_t aux_x = -2.0*(hQ2Q2_TrkmTrkp->GetMean());

	  tpQ2Q2_D0plusD0barTrk->Scale(aux_x);
	  hvn_data_Q2Q2 = tpQ2Q2_D0plusD0barTrk->ProjectionX("");
	 
          vn_data->Add(hvn_data_Q2Q2);
	  vn_data->Scale(-1.0/vn_den);
	}
          
          if(isV24HF){
          TH2D* hQ2Q2Q2Q2_D0plusD0barHFHFHF = (TH2D*)file1->Get("hist_AveNoAutoCorr_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_"+label_centbin[icent]+"_"+label_pTbin[i]+"_"+label_ybin[number_y]);
          TH2D* hQ2Q2Q2Q2_D0plusD0barHFHFHF1 = (TH2D*)file1->Get("hist_AveNoAutoCorr_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_"+label_centbin[icent]+"_"+label_pTbin[i+1]+"_"+label_ybin[number_y]);
          hQ2Q2Q2Q2_D0plusD0barHFHFHF->Add(hQ2Q2Q2Q2_D0plusD0barHFHFHF1);
          
          TH2D* hQ2Q2_D0plusD0barHF = (TH2D*)file1->Get("hist_Q2Q2_D0plusD0barHF_Re_"+label_centbin[icent]+"_"+label_pTbin[i]+"_"+label_ybin[number_y]);
          TH2D* hQ2Q2_D0plusD0barHF1 = (TH2D*)file1->Get("hist_Q2Q2_D0plusD0barHF_Re_"+label_centbin[icent]+"_"+label_pTbin[i+1]+"_"+label_ybin[number_y]);
          hQ2Q2_D0plusD0barHF->Add(hQ2Q2_D0plusD0barHF1);

          if(isForward){//should also include negative side, since in script can only use one side ...
            TH2D* aux_hQ2Q2Q2Q2_D0plusD0barHFHFHF = (TH2D*)file1->Get("hist_AveNoAutoCorr_Q2Q2Q2Q2_D0plusD0barHFHFHF_Re_"+label_centbin[icent]+"_"+label_pTbin[i]+"_y-2.0to-1.0"); 
            TH2D* aux_hQ2Q2_D0plusD0barHF = (TH2D*)file1->Get("hist_Q2Q2_D0plusD0barHF_Re_"+label_centbin[icent]+"_"+label_pTbin[i]+"_y-2.0to-1.0");
            hQ2Q2Q2Q2_D0plusD0barHFHFHF->Sumw2();
            hQ2Q2Q2Q2_D0plusD0barHFHFHF->Add(aux_hQ2Q2Q2Q2_D0plusD0barHFHFHF);
           
            hQ2Q2_D0plusD0barHF->Sumw2();
            hQ2Q2_D0plusD0barHF->Add(aux_hQ2Q2_D0plusD0barHF);
          }
          TH2D* hproduct_4D0plusD0bar = (TH2D*) hQ2Q2Q2Q2_D0plusD0barHFHFHF->Clone("hproduct_4D0plusD0bar");
          TProfile* tpQ2Q2Q2Q2_D0plusD0barHFHFHF = (TProfile*)hproduct_4D0plusD0bar->ProfileX(Form("tpQ2Q2Q2Q2_D0plusD0barHFHFHF_Re_%d", icent)+label_centbin[icent]+"_"+label_pTbin[i]+"_"+label_ybin[number_y],1,-1,"i");
          vn_data = tpQ2Q2Q2Q2_D0plusD0barHFHFHF->ProjectionX("");
          TProfile* tpQ2Q2_D0plusD0barHF = (TProfile *)hQ2Q2_D0plusD0barHF->ProfileX(Form("tpQ2Q2_D0plusD0barHF_Re_%d", icent)+label_centbin[icent]+"_"+label_pTbin[i]+"_"+label_ybin[number_y],1,-1,"i");
          TH1D *hQ2Q2_HFmHFp = (TH1D*)file1->Get("hist_Q2Q2_HFmHFp_Re_"+label_centbin[icent]);
          Double_t aux_x = -2.0*(hQ2Q2_HFmHFp->GetMean());
          tpQ2Q2_D0plusD0barHF->Scale(aux_x);
          hvn_data_Q2Q2 = tpQ2Q2_D0plusD0barHF->ProjectionX("");
 
          vn_data->Add(hvn_data_Q2Q2);
          vn_data->Scale(-1.0/vn_den);
        }
        vn_data->SetStats(kTRUE);

        c1[icent][i]->cd(2);
        hist111->Draw();
    
	TF1* fmass_combinemassvnfit = new TF1(Form("fmass_combinemassvnfit_%d_%d",icent, i),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x", fit_range_low, fit_range_high);
        

        TF1* fvn_combinemassvnfit = new TF1(Form("fvn_combinemassvnfit_%d_%d",icent, i),"[13]*( (([0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))) + ([0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))))/([0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x) ) + (1 - (([0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))) + ([0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))))/([0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x) )*( [14] + [15] * x )");//[13] is vn
 

        fmass_combinemassvnfit->SetLineColor(2);
        fmass_combinemassvnfit->SetLineWidth(1);
        
        fvn_combinemassvnfit->SetLineColor(4);
        fvn_combinemassvnfit->SetLineWidth(1);
        
        ROOT::Math::WrappedMultiTF1 wfmass_combinemassvnfit(*fmass_combinemassvnfit, fmass_combinemassvnfit->GetNdim()); //1 is dimension of function f(x)
        ROOT::Math::WrappedMultiTF1 wfvn_combinemassvnfit(*fvn_combinemassvnfit, fvn_combinemassvnfit->GetNdim());
        
        ROOT::Fit::DataOptions opt;
        ROOT::Fit::DataRange range_massfit;

        range_massfit.SetRange(fit_range_low,fit_range_high);
        ROOT::Fit::BinData datamass(opt,range_massfit);
        ROOT::Fit::FillData(datamass, h_data);
        
        ROOT::Fit::DataRange range_vnfit;
        range_vnfit.SetRange(fit_range_low,fit_range_high);
        ROOT::Fit::BinData datavn(opt,range_vnfit);
        ROOT::Fit::FillData(datavn, vn_data);
        
        ROOT::Fit::Chi2Function chi2_B(datamass, wfmass_combinemassvnfit);
        ROOT::Fit::Chi2Function chi2_SB1(datavn, wfvn_combinemassvnfit);
        
	GlobalChi2_poly3bkg_floatwidth globalChi2(chi2_B, chi2_SB1); 

        ROOT::Fit::Fitter fitter;
        
	const int Npar = 16; 
        double par0[Npar];
        for( int ipar = 0; ipar < f->GetNpar(); ipar++ ) par0[ipar] = f->GetParameter(ipar);
        par0[13] = 0.15; 
        par0[14] = 0.25; 
        par0[15] = 0.00; 
        
        
        fitter.Config().SetParamsSettings(Npar,par0);
        // fix parameter
        fitter.Config().ParSettings(2).Fix();
        fitter.Config().ParSettings(3).Fix();
        fitter.Config().ParSettings(4).Fix();
        fitter.Config().ParSettings(5).Fix();
        fitter.Config().ParSettings(7).Fix();
        fitter.Config().ParSettings(8).Fix();
        //fitter.Config().ParSettings(1).SetLimits(1.74, 2.0);
        fitter.Config().ParSettings(0).Fix();
        fitter.Config().ParSettings(1).Fix();
        fitter.Config().ParSettings(6).Fix();
        fitter.Config().ParSettings(9).Fix();
        fitter.Config().ParSettings(10).Fix();
        fitter.Config().ParSettings(11).Fix();
        fitter.Config().ParSettings(12).Fix();
        
        fitter.Config().MinimizerOptions().SetPrintLevel(0);
        fitter.Config().SetMinimizer("Minuit2","Migrad");

        fitter.FitFCN(Npar,globalChi2,0,datamass.Size()+datavn.Size(),true);	 
        ROOT::Fit::FitResult result = fitter.Result();
        result.Print(std::cout);
        
        fmass_combinemassvnfit->SetFitResult( result, iparmassfit_poly3bkg_floatwidth);
        fmass_combinemassvnfit->SetRange(range_massfit().first, range_massfit().second);
        fmass_combinemassvnfit->SetLineColor(kRed);
        h_data->GetListOfFunctions()->Add(fmass_combinemassvnfit);
        
        fvn_combinemassvnfit->SetFitResult( result, iparvnfit1_poly3bkg_floatwidth);
        fvn_combinemassvnfit->SetRange(range_vnfit().first, range_vnfit().second);
        fvn_combinemassvnfit->SetLineColor(4); 
        vn_data->GetListOfFunctions()->Add(fvn_combinemassvnfit);
        vn_data->SetTitle(Form("vn_data_cent%d_pT%d", icent, i));
        vn_data->SetName(Form("vn_data_cent%d_pT%d", icent, i));
        vn_data->SetMarkerSize(0.8);
        vn_data->SetMarkerStyle(20);
        vn_data->SetMarkerColor(kRed);
        vn_data->SetLineColor(kRed);
        vn_data->SetLineWidth(1);
        //c1->cd();
        //hist->Draw();

        vn_data->Draw("PESAME");
        vn_data->Write();
        
        fvn[i] = (TF1*)fvn_combinemassvnfit->Clone();
        fvn[i]->SetName(Form("vnfit1_cent%d_pt%d", icent,i));
        fvn[i]->Write();
        
        tex->DrawLatex(0.22,0.86,range_centbin[icent]);
        tex->DrawLatex(0.22,0.80,Form("%.1f < p_{T} < %.1f GeV/c",momentumbin[i],momentumbin[i+1]));
        tex->DrawLatex(0.22,0.74,range_ybin[number_y]);
        tex->DrawLatex(0.60,0.86,Form(plots_label1,fvn_combinemassvnfit->GetParameter(13),fvn_combinemassvnfit->GetParError(13)));
        
        texCMS->DrawLatex(.18,.97,"#font[61]{CMS} #it{Preliminary}");
        texCMS->DrawLatex(0.62,0.97, "#scale[0.8]{PbPb #sqrt{s_{NN}} = 5.02 TeV}");
        
        fmasstotal[i] = (TF1*)fmass_combinemassvnfit->Clone();
        fmasstotal[i]->SetName(Form("masstotalfcn_cent%d_pt%d",icent,i));
        fmasstotal[i]->Write();
        
        TLegend* leg1 = new TLegend(0.65,0.78,0.95,0.9,NULL,"brNDC");
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.045);
        leg1->SetTextFont(42);
        leg1->SetFillStyle(0);
        leg1->AddEntry(h_data,"data","p");

        TF1* falpha = new TF1(Form("falpha_%d_%d",icent, i),"( [0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) )/( [0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x )", fit_range_low,fit_range_high);

        for(int j=0;j<13;j++)
        {
            falpha->FixParameter(j,fmass_combinemassvnfit->GetParameter(j));
        }
        
        cout <<"=================================" << endl;
        cout << "fmass_combinemassvnfit->GetParameter(j)" << " " << fmass_combinemassvnfit->GetParameter(0)<< endl; 
        //ofparameter << fmass_combinemassvnfit->GetParameter(0)<< endl; 
        falpha->SetName(Form("sigfrac_fcn_cent%d_pt%d",icent, i));
        falpha->Write();
        
        double xmass[10000];
        double pullmass[10000];
        
        float Chi2=0;
        int ndf = 100 - 13;
        
        for(int k=0;k<h_data->GetNbinsX();k++)
        {
            xmass[k] = h_data->GetXaxis()->GetBinCenter(k+1);
            pullmass[k] = (h_data->GetBinContent(k+1) - fmass_combinemassvnfit->Eval(xmass[k]))/h_data->GetBinError(k+1);
            //cout<<pullmass[k]<<endl;
            Chi2 += pullmass[k]*pullmass[k];
        }

        c1[icent][i]->cd(1);
        tex->DrawLatex(0.22,0.68,Form("Chi2/ndf = %.0f/%d",Chi2,ndf));
	Double_t prob_massfit = TMath::Prob(Chi2,ndf);
	tex->DrawLatex(0.22,0.62,Form("Prob. = %.4f",prob_massfit));

        double xvn[200];
        double pullvn[200];
        double vny[200];
        
        float Chi2vn=0.;
        int ndfvn = 13 - 3;
        

        for(int k=0;k<vn_data->GetNbinsX();k++)
        {
            vny[k] = vn_data->GetBinContent(k+1);
            xvn[k] = vn_data->GetBinCenter(k+1);
            pullvn[k] = (vny[k] - fvn_combinemassvnfit->Eval(xvn[k]))/vn_data->GetBinError(k+1);
            cout<<pullvn[k]<<endl;
            //cout<<pullmass[k]<<endl;
            Chi2vn += pullvn[k]*pullvn[k];
        }

        c1[icent][i]->cd(2);
        tex->DrawLatex(0.22,0.68,Form("Chi2/ndf = %.0f/%d",Chi2vn,ndfvn));
	Double_t prob_vnfit = TMath::Prob(Chi2vn,ndfvn);
	tex->DrawLatex(0.22,0.62,Form("Prob. = %.4f",prob_vnfit));
     
        sumyvn_icent += fmass_combinemassvnfit->GetParameter(0)*fvn_combinemassvnfit->GetParameter(13);
        sumyvne_icent += fmass_combinemassvnfit->GetParameter(0)*fvn_combinemassvnfit->GetParError(13);
        
        sumy_icent += fmass_combinemassvnfit->GetParameter(0);
        
        c1[icent][i]->Print(Form(folder_plots+label2_ybin+"_%d_%d.pdf",icent, i));
       ofout << fvn_combinemassvnfit->GetParameter(13) << "," <<  fmass_combinemassvnfit->GetParameter(0) << endl; 
   }//end pT loop; 
    vn[icent] = sumyvn_icent/sumy_icent; 
    vne[icent] = sumyvne_icent/sumy_icent; 
    if(icent==0)   ofout  << vn[icent] << ", "<< icent << endl; 
    if(icent==1)   ofout  << vn[icent] << ", "<< icent << endl; 
    if(icent==2)   ofout  << vn[icent] << ", "<< icent << endl; 
    if(icent==3)   ofout  << vn[icent] << ", "<< icent << endl; 
    if(icent==4)   ofout  << vn[icent] << ", "<< icent << endl; 
    }//centrality 
    
    TGraphErrors* vnplot = new TGraphErrors(n_centbin,centbin,vn,centbinwidth,vne); 
     
    vnplot->Print("all");
    vnplot->SetName("vnvspt");
    vnplot->GetYaxis()->SetTitle(plots_label2);
    vnplot->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    vnplot->SetTitle("D^{0} + #bar{D^{0}}");
    
    TCanvas *c3 = new TCanvas();
    c3->cd();
    vnplot->Draw("ApE");
    vnplot->SetTitle(""); 
    vnplot->SetMarkerStyle(20);
    vnplot->SetMarkerColor(1);
    vnplot->SetLineColor(1);
    ///vnplot->GetYaxis()->SetRangeUser(-0.07,0.27);//v2
    ///vnplot->GetYaxis()->SetRangeUser(-0.14,0.17); //v3
    //vnplot->GetYaxis()->SetRangeUser(-0.14,0.27); //v2 and v3
    vnplot->GetYaxis()->SetRangeUser(-0.05,0.25); //v2 and v3
    vnplot->GetXaxis()->SetLimits(0.,70.0); 
    texCMS->DrawLatex(.10,.93,"#font[61]{CMS} #it{Preliminary}");
    texCMS->DrawLatex(0.64,0.93, "#scale[0.8]{PbPb #sqrt{s_{NN}} = 5.02 TeV}");
   
    auto leg_vnplot = new TLegend(0.12,0.7,0.6,0.9);
    //leg_vnplot->SetHeader(range_centbin[icent]+", "+range_ybin[number_y],""); // option "C" allows to center the header
    leg_vnplot->AddEntry(vnplot,"D^{0} + #bar{D^{0}}","pl");
    leg_vnplot->SetFillStyle(0);
    leg_vnplot->Draw();
    c3->Print(folder_plots+label2_ybin+"_Vn.pdf");

    vnplot->Write();
    ofile.Close();   
}
