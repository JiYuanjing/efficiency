#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TPaveStats.h"

double myGaus(Double_t *x,Double_t *par){ 
	Double_t fitval=0;
	Double_t arg = 0;
	if (par[2]!=0) arg = (x[0] - par[1])/par[2];
	Double_t bin_width = 0.0005;
	Double_t constant = bin_width*par[0]/sqrt(2*TMath::Pi())/par[2];
	fitval = constant*TMath::Exp(-0.5*arg*arg)+par[3]+par[4]*x[0];
	//fitval = bin_width*par[0]*TMath::Gaus(x[0],par[1],par[2],kTRUE)+par[3]+par[4]*x[0];
	return fitval;
} 

void mixedmassDstarCent(int option=0,Double_t lpt = 1.5, Double_t hpt =2, int centlow =6, int centhigh = 8,char* outname="fitparametercent.txt", 
								  TString rootname="../D0withrapiditybin.root")
{
gStyle->SetOptFit(111);
//TString rootname = "all_sample_hists.root";
TFile *f = TFile::Open(rootname);
TH3F* imp = (TH3F*)f->Get("mh3InvariantMassVsPtVsCentDstar");
TH3F* implk = (TH3F*)f->Get("mh3InvariantMassVsPtVsCentLikeDstar");
TH3F* impsb = (TH3F*)f->Get("mh3InvariantMassVsPtVsCentSBDstar;1");
TH3F* impmixed = (TH3F*)f->Get("mh3InvariantMassVsPtVsCentDstarMixedEvent;1");
TH1F* ref = (TH1F*)f->Get("mh1CentWg")->Clone("mh1CentWg");
ref->SetDirectory(0);

imp->Sumw2();
implk->Sumw2();
impsb->Sumw2();
impmixed->Sumw2();
TCanvas* c = new TCanvas("c","c");

//lastbin=120;

int firstbin_pt = imp->GetXaxis()->FindBin(lpt+1e-6);
int lastbin_pt = imp->GetXaxis()->FindBin(hpt-1e-6);
int firstbin_cent = imp->GetYaxis()->FindBin(centlow+1e-6);
int lastbin_cent = imp->GetYaxis()->FindBin(centhigh-1e-6);
cout<<firstbin_cent<<"\t"<<lastbin_cent<<endl;
TH1* h = imp->ProjectionZ("h_pz",firstbin_pt,lastbin_pt, firstbin_cent, lastbin_cent,"e");
TH1* hlk = implk->ProjectionZ("hlk_pz",firstbin_pt,lastbin_pt,  firstbin_cent, lastbin_cent,"e");
TH1* hsb = impsb->ProjectionZ("hsb_pz",firstbin_pt,lastbin_pt, firstbin_cent, lastbin_cent ,"e");
TH1* hmixed = impmixed->ProjectionZ("hmixed_pz",firstbin_pt,lastbin_pt,  firstbin_cent, lastbin_cent,"e");

/*
h->Sumw2();
hlk->Sumw2();
hsb->Sumw2();
hmixed->Sumw2();

if (lpt<){
h->Rebin();
hlk->Rebin();
hmixed->Rebin();
}
*/

//scale
Double_t  bwd_hsb=hsb->GetBinWidth(1);
Double_t scale_sb = h->Integral(h->FindBin(0.152+1e-6), h->FindBin(0.175-1e-6))/hsb->Integral(hsb->FindBin(0.152+1e-6), hsb->FindBin(0.175-1e-6));
cout<<"sideband scale: "<<1.0/scale_sb<<endl;
hsb->Scale(scale_sb);
Double_t  bwd_hmixed=hmixed->GetBinWidth(1);
Double_t scale = h->Integral(h->FindBin(0.152+1e-6), h->FindBin(0.175-1e-6))/hmixed->Integral(hmixed->FindBin(0.152+1e-6), hmixed->FindBin(0.175-1e-6));
cout<<"mixed event scale: "<<1.0/scale<<endl;
hmixed->Scale(scale);

if (centhigh<3) hmixed=(TH1*)hsb->Clone("sb");

//substract background
TH1* hdelta = (TH1*)h->Clone("delta");
hdelta->Sumw2();
if (option==0) hdelta->Add(h, hmixed, 1, -1);
else if (option==1) hdelta->Add(h, hlk, 1, -1);

//set the name and the 
char Ytitle[50];
sprintf(Ytitle,"Counts/%.2fMeV/c^{2}",bwd_hmixed*1000);
h->GetYaxis()->SetTitle(Ytitle);
hlk->GetYaxis()->SetTitle(Ytitle);
hsb->GetYaxis()->SetTitle("Counts");
hmixed->GetYaxis()->SetTitle(Ytitle);
h->GetXaxis()->SetTitle("m_{K#pi#pi}-m_{K#pi}(GeV/c^{2})");
hdelta->GetYaxis()->SetTitle(Ytitle);
hdelta->GetXaxis()->SetTitle("m_{K#pi#pi}-m_{K#pi}(GeV/c^{2})");
h->SetLineColor(kBlue);
h->SetMarkerColor(kBlue);
h->SetMarkerStyle(8);
hlk->SetLineColor(kCyan);
hlk->SetMarkerStyle(1);
hlk->SetMarkerColor(kCyan);
hsb->SetLineColor(kGreen);
hsb->SetMarkerStyle(1);
hsb->SetMarkerColor(kGreen);
hmixed->SetLineColor(kRed);
hmixed->SetMarkerStyle(1);
hmixed->SetMarkerColor(kRed);
hdelta->SetLineColor(kBlack);
hdelta->SetMarkerColor(kBlack);
hdelta->SetMarkerStyle(8);
hdelta->GetXaxis()->SetTitleOffset(1);
hdelta->Draw();
//gPad->Update();
//TPaveStats *ps = (TPaveStats*)hdelta->GetListOfFunctions()->FindObject("stats");
//ps->SetX1NDC(0.65); ps->SetX2NDC(0.85);
hmixed->Draw("same");
h->Draw("same");
hlk->Draw("same");
hmixed->Draw("same");

//fit
Double_t fitrange1=0.139,fitrange2=0.175 ,fitrange3=0.175;
if (hpt<=2.5) { fitrange1=0.14;fitrange2=0.148;}
else {fitrange1=0.139;fitrange2=0.155;}
if (centhigh<1&&lpt>=6) {fitrange1=0.141;fitrange2=0.150;fitrange3=0.165;}
Double_t para[5];
TF1* fit = new TF1("fit",myGaus,fitrange1,fitrange3,5);
TF1* fitgaus = new TF1("fitgaus","gausn",fitrange1,fitrange2);
TF1*fitbk = new TF1("fitbk","pol1",0.155,0.175);
hdelta->Fit(fitgaus,"R");
hdelta->Fit(fitbk,"R");
fitgaus->GetParameters(para);
fitbk->GetParameters(para+3);
Double_t bin_width = hdelta->GetBinWidth(1);
cout<<"bin width"<<bin_width<<endl;
fit->SetParNames ("D* Counts","#DeltaM(GeV/c^{2})","#sigma(GeV/c^{2})","p0","p1");
fit->SetParameter(0,hdelta->Integral());
fit->SetParameter(1,0.1454);
fit->SetParLimits(1,0.145,0.146);
fit->SetParameter(2,para[2]);
fit->SetParLimits(2,0.0005,0.00095);
fit->SetParameter(3,para[3]);
fit->SetParameter(4,para[4]);
fit->SetLineColor(kRed);
hdelta->Fit(fit,"R");
fit->SetNpx(200);
fit->Draw("same");
fit->GetParameters(para);
Double_t significance = para[0] / fit->GetParError(0);
Double_t peakerr = fit->GetParError(1);
Double_t firstb = 0.135;
cout<<"firstbin "<<firstb<<endl;
Double_t hb=0.15,lb=0.14;
Int_t hrange = (hb-firstb)/(Double_t)bin_width+1;
Int_t lrange = (lb-firstb)/(Double_t)bin_width+1;
cout<<hrange<<" "<<lrange<<endl;
cout<<hdelta->Integral(lrange,hrange)<<endl;
double y_max = h->GetMaximum();
Double_t drawmax=1;
if (hpt<=2) { 
	drawmax=0.4;
	hdelta->GetXaxis()->SetRangeUser(0.138,0.160);
}
else if (hpt<=3&&hpt>2) drawmax=1;
else drawmax=1.6;
hdelta->GetYaxis()->SetRangeUser(-0.08*y_max,drawmax*y_max);
TLegend* leg = new TLegend(0.62,0.3,0.88,0.5);
leg->AddEntry(h,"right sign(RS)","lpfe");
//leg->AddEntry(hsb,"side band(SB)","lpfe");
leg->AddEntry(hlk,"wrong sign(WS)","lpfe");
leg->AddEntry(hmixed,"Mixed Event(ME)","lpfe");
if (option==0) leg->AddEntry(hdelta,"RS-ME","lpfe");
else if (option==1) leg->AddEntry(hdelta,"RS-WS","lpfe");
leg->SetBorderSize(0);
char ptrange[50];
sprintf(ptrange,"%2.2f<p_{T}<%2.2f (GeV/c)",lpt,hpt);
char peak[50];
char s[50];
char counts[50];
sprintf(s,"S/#DeltaS = %.0f",significance);
sprintf(counts,"D* Counts=%.0f", para[0]);
sprintf(peak,"#DeltaM=%.2f #pm %.2fMeV/c^{2}",para[1]*1000,peakerr*1000);
TLatex latex;
latex.SetTextSize(0.04);
Double_t ymax = hdelta->GetMaximum();
Double_t xmax;
char sigma[50];
sprintf(sigma,"#sigma=%.2f#times10^{-4}GeV/c^{2}",para[2]*pow(10,4));
Double_t ltxp0=0.136;
if (hpt<=2) ltxp0=0.139;
int centedge[10]={80,70,60,50,40,30,20,10,5,0};
latex.DrawLatex(ltxp0,ymax*0.91,Form("%d~%d % MB",centedge[centhigh+1],centedge[centlow]));
latex.DrawLatex(ltxp0,ymax*0.84,ptrange);
latex.DrawLatex(ltxp0,ymax*0.77,peak);
//latex.DrawLatex(0.136,ymax*0.70,sigma);
//latex.DrawLatex(0.136,ymax*0.63,counts);
latex.DrawLatex(ltxp0,ymax*0.70,s);

char events[50];
sprintf(events,"%.0fM Events",ref->Integral( ref->FindBin(centlow), ref->FindBin(centhigh))*pow(10,-6));
latex.DrawLatex(0.147,ymax*0.91,events);
leg->Draw();
char gifname[50];
sprintf(gifname,"pic/%.1fto%.1fptcent_%d_%d.gif",lpt,hpt,centlow,centhigh);
c->SaveAs(gifname);

//output the parameter
ofstream outf;
outf.open(outname, ios::app);
if (lpt==1111&&hpt==10) {
outf<<para[0]/2.0<<"\t"<<fit->GetParError(0)/2.0<<"\t"<<
	    para[1]<<"\t"<<fit->GetParError(1)<<"\t"<<
		para[2]<<"\t"<<fit->GetParError(2)<<endl;}
else {
	outf<<para[0]<<"\t"<<fit->GetParError(0)<<"\t"<<
	    para[1]<<"\t"<<fit->GetParError(1)<<"\t"<<
		para[2]<<"\t"<<fit->GetParError(2)<<endl;
	}
}


