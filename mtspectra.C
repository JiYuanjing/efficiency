#include "drawstyle.h"
#include "anaconst.h"
TString effversion = "effDstar_1218.root";
// TString effversion = "effDstar.root";
//part 1: N binary
 const float Nbin[9] = {12.91047, 29.32408, 62.25545, 120.65795, 215.95050, 360.58912, 579.89409, 812.73278, 1042.75372};

void mtspectra(int lcent=0,int hcent=3){
	globleStyle();
	myStyle();

	int const centedge[10]={80,70,60,50,40,30,20,10,5,0};

	TCanvas* c = new TCanvas("spectra","spectra",600,800);
	c->Divide(1,2); 
	c->cd(1);
	gPad->SetPad(0, 0.35,1,1);
	setPad(0.15,0.05,0.1,0);
	c->cd(2);
	gPad->SetPad(0,0,1,0.35);
	setPad(0.15,0.05,0,0.2);

	double mtEdge[nPtBins+1];
	for (int i=0;i<nPtBins+1;i++){
		mtEdge[i] = sqrt(PtEdge[i]*PtEdge[i]+mDstar*mDstar);
	}
	
	TH1D* hDstarmt = new TH1D(Form("Dstarmt_%d_%d", centedge[hcent+1],centedge[lcent]), "D* mtspectra;m_{T}(GeV/c^{2});#frac{d^{2}N}{2#pi m_{T}dm_{T}dy}",nPtBins, mtEdge);
	
	TFile* fD0 = TFile::Open(Form("D0/data/hD0_myspectra_%d_%d.root",centedge[hcent+1],centedge[lcent]));
	TH1* hD0 = (TH1*)fD0->Get("D0spectra")->Clone(Form("D0spectra_%d_%d",centedge[hcent+1],centedge[lcent]));
	hD0->SetDirectory(0);
	fD0->Close();
	
	TFile* f = TFile::Open(Form("data/correct_%d_%d.root",centedge[hcent+1],centedge[lcent]),"update");
	f->Delete(Form("Dstarmt_%d_%d;*", centedge[hcent+1],centedge[lcent]));
	f->Delete(Form("D0mt_%d_%d;*", centedge[hcent+1],centedge[lcent]));
	TH1D* hptspectra = (TH1D*)f->Get(Form("hCorrect_%d_%d",centedge[hcent+1],centedge[lcent]));
	double tmpyield, tmperr;
	for (int i=0;i<nPtBins;i++){
		tmpyield = hptspectra->GetBinContent(i+1);
		tmperr = hptspectra->GetBinError(i+1);
		hDstarmt->SetBinContent(i+1, tmpyield);
		hDstarmt->SetBinError(i+1, tmperr);
	}

	c->cd(1);
	hDstarmt->Draw();
	gPad->SetLogy();
	TF1* fitfunDstar = new TF1("levyDstar","(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+2.01026*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+2.01026*2.01026)-2.01026)/([2]*[1]),-[2])*x[0]", 0, 12);
	TF1* fitfunD0 = new TF1("levyD0","(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+1.86483*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+1.86483*1.86483)-1.86483)/([2]*[1]),-[2])*x[0]", 0, 12);
	
	double D0mtEdge[nPtBins+1];
	for (int i=0;i<nPtBins+1;i++){
		D0mtEdge[i] = sqrt(PtEdge[i]*PtEdge[i]+mD0*mD0);
	}
	TH1D* hD0mt = new TH1D(Form("D0mt_%d_%d", centedge[hcent+1],centedge[lcent]), "D0 mtspectra;m_{T}(GeV/c^{2});#frac{d^{2}N}{2#pi m_{T}dm_{T}dy}", nPtBins, D0mtEdge);

	for (int i=0;i<nPtBins;i++){
		tmpyield = hD0->GetBinContent(i+1);
		tmperr = hD0->GetBinError(i+1);
		hD0mt->SetBinContent(i+1, tmpyield);
		hD0mt->SetBinError(i+1, tmperr);
	}
	
	hD0mt->SetMarkerColor(kRed);
	hD0mt->Draw();
	double p[3]={0.5,0.33,15};
	fitfunDstar->SetParameters(p);
	fitfunDstar->SetLineColor(kGreen);
	fitfunD0->SetParameters(p);
	fitfunD0->SetLineColor(kBlue);
	hD0mt->Fit(fitfunD0, "R+");
	hDstarmt->Draw("same");
	hDstarmt->Fit(fitfunDstar,"R+");
	TLegend* l =new TLegend(0.65,0.7,0.9,0.9);
	l->SetHeader(Form("Au+Au %d~%d%s",centedge[hcent+1],centedge[lcent],"%"));
	l->AddEntry(hDstarmt, "D*","lpe" );
	l->AddEntry(hD0mt,"D^{0}", "lpe");
	//l->AddEntry(hD0my,"my D^{0}","lpe");
	l->Draw("same");

	c->cd(2);
	fitfunD0->SetNpx(100);
	fitfunDstar->SetNpx(100);
	TH1* hfitDstar = (TH1*)fitfunDstar->GetHistogram();
	TH1* hfitD0 = (TH1*)fitfunD0->GetHistogram();
	TH1* hratio = (TH1*)hfitDstar->Clone("mtratio");
	hratio->Divide(hfitD0);
	// for (int i=0;i<nPtBins;i++){
		// double x, y;
		// x=hratio->GetBinCenter(i+1);
		// y=hratio->GetBinContent(i+1);
		// hratio->SetBinContent(i+1, y/fitfunD0->Eval(x));
	// }
	
	hratio->GetYaxis()->SetNdivisions(50205);
	hratio->GetYaxis()->SetLabelSize(0.07);
	hratio->GetXaxis()->SetLabelSize(0.07);
	hratio->GetYaxis()->SetRangeUser(0,1.1);
	hratio->GetXaxis()->SetRangeUser(D0mtEdge[0],mtEdge[nPtBins]);
	hratio->SetTitle("");
	hratio->GetYaxis()->SetTitle("D*/D^{0}");
	hratio->SetDirectory(0);
	hratio->Draw();

	// fstream fratio;
	// fratio.open(Form("data/mtratio_%d_%d.txt",centedge[hcent+1],centedge[lcent]), ios::out | ios::trunc);

	// for (int i=0;i<nPtBins;i++){
		// fratio<<hratio->GetBinCenter(i+1)<<"\t"<<hratio->GetBinContent(i+1)<<"\t"<<hratio->GetBinError(i+1)<<endl;
	// } 
	// fratio.close();

	TF1* ratio = new TF1("ratio", "/((1/(2*TMath::Pi()))*[3]*([5]-1)*([5]-2)/([5]*[4]*([5]*[4]+1.86483*([5]-2)))*TMath::Power(([5]*[4]+TMath::Sqrt(x[0]*x[0]+1.86483*1.86483)-1.86483)/([5]*[4]),-[5])*x[0])");
	
	drawLine(D0mtEdge[0], 1, mtEdge[nPtBins], 1, 0.5, 2 ,kBlue);
	c->SaveAs(Form("pic/mtratio_%d_%d.gif",centedge[hcent+1],centedge[lcent]));
	hDstarmt->Write();
	hD0mt->Write();
	// hratio->Write();
}
