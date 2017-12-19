#include "drawstyle.h"
#include "anaconst.h"
TString effversion = "effDstar.root";
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
