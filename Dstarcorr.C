#include "drawstyle.h"
#include "anaconst.h"
TString effversion = "effDstar_1218.root";
// TString effversion = "effDstar.root";
//part 1: N binary
 const float Nbin[9] = {12.91047, 29.32408, 62.25545, 120.65795, 215.95050, 360.58912, 579.89409, 812.73278, 1042.75372};
void makeefficiency(int lcent=4,int hcent=6){
TFile* f = TFile::Open(effversion,"update");
TH1D* hnocuts[9];
TH1D* hwithallcuts[9];
TH1D* heff_rw[9];
TH1D* heff;
TH1D* sumnocut;
TH1D* sumwithcut;
TH1D* hefftot;
double totnum=0;
f->Delete("rebinsumnocuts;*");
f->Delete("sumnocuts;*");
f->Delete("rebinsumwithcuts;*");
f->Delete("sumwithcut;*");
f->Delete("rebineff;*");
f->Delete(Form("eff_%d~%d;*",centedge[hcent+1],centedge[lcent]));
f->Delete(Form("htoteff%d~%d%s;*",centedge[hcent+1],centedge[lcent],"%"));

for (int cent=lcent;cent<hcent+1;cent++) {
hnocuts[cent]= (TH1D*)f->Get(Form("hNoCuts_cent%i",cent))->Clone(Form("hNoCuts_cent%d",cent));
hwithallcuts[cent] = (TH1D*)f->Get(Form("hafterCuts_cent%i",cent))->Clone(Form("hafterCuts_cent%i",cent));
hwithallcuts[cent]=(TH1D*)hwithallcuts[cent]->Rebin(nPtBins,Form("rebinhwithcutscent%d", cent), PtEdge);
hnocuts[cent]=(TH1D*)hnocuts[cent]->Rebin(nPtBins,Form("rebinhnocutscent%d", cent), PtEdge);
heff_rw[cent]=(TH1D*)hwithallcuts[cent]->Clone(Form("heff_rw%d",cent));
heff_rw[cent]->Divide(hnocuts[cent]);
}

sumnocut=(TH1D*)hnocuts[lcent]->Clone("sumnocuts");
sumwithcut=(TH1D*)hwithallcuts[lcent]->Clone("sumwithcut");
hefftot=(TH1D*)heff_rw[lcent]->Clone(Form("htoteff%d~%d%s",centedge[hcent+1],centedge[lcent],"%"));
sumnocut->Scale(0);
sumwithcut->Scale(0);
hefftot->Scale(0);

for (int cent=lcent;cent<hcent+1;cent++){
	sumnocut->Add(hnocuts[cent], Nbin[cent]);
	sumwithcut->Add(hwithallcuts[cent], Nbin[cent]);
	hefftot->Add(heff_rw[cent], Nbin[cent]);
	totnum+=Nbin[cent];	
}

hefftot->Scale(1.0/totnum);
hefftot->SetDirectory(0);
heff=(TH1D*)sumwithcut->Clone(Form("eff_%d~%d",centedge[hcent+1],centedge[lcent]));
heff->Divide(sumnocut);
heff->SetTitle(Form("%d~%d%s Eff.",centedge[hcent+1],centedge[lcent],"%"));
heff->GetXaxis()->SetTitle("p_{T}");
heff->GetYaxis()->SetTitle("Eff");
heff->Draw();
heff->Write();
hefftot->GetXaxis()->SetTitle("p_{T}");
hefftot->GetYaxis()->SetTitle("Eff");
hefftot->Write();
sumnocut->Write();
sumwithcut->Write();
f->Close();
//return hefftot;
}

void Dstarcorr(int lcent=0,int hcent=3){
	globleStyle();
	myStyle();
	//TH1D* heff = makeefficiency(lcent,hcent);
	makeefficiency(lcent,hcent);
	int const centedge[10]={80,70,60,50,40,30,20,10,5,0};

	TFile* f = TFile::Open(effversion);

	// TH1D* heff = (TH1D*)f->Get(Form("eff_%d~%d",centedge[hcent+1],centedge[lcent]))->Clone("eff");
	TH1D* heff = (TH1D*)f->Get(Form("htoteff%d~%d%s",centedge[hcent+1],centedge[lcent], "%"))->Clone("eff");
	heff->SetDirectory(0);
	f->Close();

	TCanvas* c = new TCanvas("spectra","spectra",600,800);
	c->Divide(1,2); 
  c->cd(1);
  gPad->SetPad(0, 0.35,1,1);
  setPad(0.15,0.05,0.1,0);
  c->cd(2);
  gPad->SetPad(0,0,1,0.35);
  setPad(0.15,0.05,0,0.2);

	TFile* ff = TFile::Open(Form("Dstar/data/Dstar_hists_%d_%d.root", centedge[hcent+1],centedge[lcent]),"update");
	//f->Delete("Correct;*");
	TH1D* hraw = (TH1D*)ff->Get("rawspectra")->Clone("rawspectra");
	hraw->SetDirectory(0);
	hraw->GetYaxis()->SetTitle("#frac{d^{2}N}{2#pi p_{T}dp_{T}dy}");
	TH1D* hcorr = (TH1D*)hraw->Clone(Form("hCorrect_%d_%d",centedge[hcent+1],centedge[lcent]));
	hcorr->SetDirectory(0);
	double branchratio=DstarBranch*D0Branch;
	hcorr->Divide(hcorr,heff,0.5,branchratio);
  TH1D* hrawyield = (TH1D*)ff->Get("Counts")->Clone("rawyield");
	hrawyield->SetDirectory(0);
	TH1D* heffyield = (TH1D*)hrawyield->Clone(Form("heffyield_%d_%d",centedge[hcent+1],centedge[lcent]));
	heffyield->SetDirectory(0);
	heffyield->Divide(heff);
	c->cd(1);
	hcorr->Draw();
	//hcorr->Write();
	gPad->SetLogy();
	TF1* fitfun = new TF1("levy","(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+2.01026*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+2.01026*2.01026)-2.01026)/([2]*[1]),-[2])*x[0]", 0, 12);
	double p[3]={0.5,0.33,15};
	fitfun->SetParameters(p);
	fitfun->SetLineColor(kGreen);
	//hcorr->Fit(fitfun,"R+");

	TFile* fD0 = TFile::Open(Form("D0/data/hD0_myspectra_%d_%d.root",centedge[hcent+1],centedge[lcent]));
	TH1* hD0 = (TH1*)fD0->Get("D0spectra")->Clone(Form("D0spectra_%d_%d",centedge[hcent+1],centedge[lcent]));
	hD0->Draw();
	hD0->SetMarkerColor(kRed);
	hcorr->Draw("same");
//TFile* fD0spectrum = TFile::Open("D0_Spectra_Run14_HFT_beforePtShift.root");
//TGraphErrors* g0_10 = (TGraphErrors*)fD0spectrum->Get(Form("gD0_sys_%d_%d",centedge[hcent+1],centedge[lcent]));
//g0_10->Draw("psame");
//g0_10->SetMarkerColor(kRed);
//TF1* flvey0_10 = (TF1*)fD0spectrum->Get(Form("flevy_%d_%d",centedge[hcent+1],centedge[lcent]));
//flevy_0_10->Draw("same");
//g0_10->Fit(flvey0_10);
//flevy_0_10->SetLineColor(kRed);

/*
TFile* fD0my = TFile::Open( Form("hD0_myspectra_%d_%d",centedge[hcent+1],centedge[lcent]));
TH1D* hD0my = fD0my->Get("D0spectra")->Clone("myspectra");
hD0my->Draw("samep");
hD0my->SetDirectory(0);
hD0my->SetMarkerColor(kViolet);
*/

	TLegend* l =new TLegend(0.65,0.7,0.9,0.9);
	l->SetHeader(Form("Au+Au %d~%d%s",centedge[hcent+1],centedge[lcent],"%"));
	l->AddEntry(hcorr, "D*","lpe" );
	l->AddEntry(hD0,"D^{0}", "lpe");
	//l->AddEntry(hD0my,"my D^{0}","lpe");
	l->Draw("same");

	c->cd(2);
	TH1* hratio = (TH1*)hcorr->Clone("ratio");
	hratio->Divide(hD0);
	hratio->GetYaxis()->SetNdivisions(50205);
	hratio->GetYaxis()->SetLabelSize(0.07);
	hratio->GetXaxis()->SetLabelSize(0.07);
	hratio->GetYaxis()->SetRangeUser(0,1.1);
  hratio->SetTitle("");
  hratio->GetYaxis()->SetTitle("D*/D^{0}");
  hratio->SetDirectory(0);
	hratio->Draw();

  fstream fratio;
  fratio.open(Form("data/ratio_%d_%d.txt",centedge[hcent+1],centedge[lcent]), ios::out | ios::trunc);

  for (int i=0;i<nPtBins;i++){
    fratio<<hratio->GetBinCenter(i+1)<<"\t"<<hratio->GetBinContent(i+1)<<"\t"<<hratio->GetBinError(i+1)<<endl;
  } 
  fratio.close();

	drawLine(PtEdge[0], 1, PtEdge[nPtBins], 1, 0.5, 2 ,kBlue);
	c->SaveAs(Form("pic/ratio_%d_%d.gif",centedge[hcent+1],centedge[lcent]));
	TFile* after = new TFile(Form("data/correct_%d_%d.root",centedge[hcent+1],centedge[lcent]),"recreate");
	heffyield->Write();
  hrawyield->Write();
  hcorr->Write();
  hratio->Write();
}
