#include "anaconst.h"
#include "drawstyle.h"

void sethist(TH1* h, char* xtitle, char* ytitle);
void drawratio(){
  globleStyle();
  myStyle();
	TH1D* hratio[nCentBins];
  TFile* f=NULL;
	for (int i=0;i<nCentBins;i++){
		int hcent = hCentEdge[i], lcent = lCentEdge[i];
		f = TFile::Open(Form("data/correct_%d_%d.root", centedge[hcent+1], centedge[lcent]));
		hratio[i] = (TH1D*)f->Get("ratio")->Clone(Form("ratio_%d_%d",  centedge[hcent+1], centedge[lcent]));
		hratio[i]->SetDirectory(0);
    f->Close();
		sethist(hratio[i], "p_{T}/GeV","D*/D^{0}");
    hratio[i]->SetMarkerColor(2+i);
    hratio[i]->SetLineColor(2+i);
    hratio[i]->GetXaxis()->SetRangeUser(1.5,6);
		hratio[i]->Draw("same");
	}
	gPad->SetGridy(1);
  TLegend* l = new TLegend(0.7,0.2,0.9,0.5);
  for (int i=0;i<nCentBins;i++){
		int hcent = hCentEdge[i], lcent = lCentEdge[i];
    l->AddEntry(hratio[i],Form("%d%s~%d%s",centedge[hcent+1],"%",centedge[lcent],"%"),"lep");
  }

  TFile* fphy = TFile::Open("DstarD0ratio_beforescale.root");
  TH1* hpy = (TH1*)fphy->Get("DstarD0ratio");
  hpy->SetLineColor(1);
  hpy->SetLineWidth(2);
  hpy->Draw("same");

  l->AddEntry(hpy,"pythia","lep");
  l->Draw();
	gPad->SaveAs("pic/ratio.gif");
}

void sethist(TH1* h, char* xtitle, char* ytitle){
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->SetLineWidth(2);
  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetNdivisions(515);
  h->GetZaxis()->SetNdivisions(505);
  double titleSize = 0.05;
	h->GetXaxis()->SetTitleSize(titleSize);
	h->GetYaxis()->SetTitleSize(titleSize);
  h->GetXaxis()->SetTitleOffset(1);
  h->GetYaxis()->SetTitleOffset(1);
}

