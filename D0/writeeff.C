void writeeff(){ 
//TFile* f = TFile::Open("D0_eff_from_xiaolong.root");
//TFile* f = TFile::Open("D0_eff_combine_newPID_newCuts.root");
// TFile* f = TFile::Open("D0.extend.root");
TFile* f = TFile::Open("D0.extend_1207.root");
//TFile* f = TFile::Open("D0_eff_combine_newPID_newCuts_1027DoubleStat.root");
//TFile* f = TFile::Open("D0_eff_combine_newPID_newCuts_extendpT_Rebin_1113_ExtractDiff.root");
//TFile* fvtx = TFile::Open("vtxCorr_default.root");
//TFile* fvtx = TFile::Open("vtxCorr_pt0_3.root");
TH3* h3nocuts = (TH3*)f->Get("h3PtCentY")->Clone("h3PtCentY");
//TH2* h2allcuts = (TH2*)f->Get("h2PtCut_pt1")->Clone("h3PtCentYCut_pt1");
TH3* h3allcuts = (TH3*)f->Get("h3PtCentYCut_pt1")->Clone("h3PtCentCutY");
h3nocuts->SetDirectory(0);
h3allcuts->SetDirectory(0);
f->Close();
TFile* effout = new TFile("D0_newPID_extend_eff.root","recreate");
const int nybin = 10;
const int ncentbin = 9;
TH1D* hnocuts[ncentbin][nybin];
TH1D* hwithallcuts[ncentbin][nybin];
TH1D* heff[ncentbin][nybin];
int const centedge[10]={80,70,60,50,40,30,20,10,5,0};
int centbin=-1;
for (int cent=0;cent<ncentbin;cent++){
	centbin = h3nocuts->GetYaxis()->FindBin(cent+0.5); 
	for (int iy=0;iy<nybin;iy++){
		int ylow = h3nocuts->GetZaxis()->FindBin(2.0/(double)nybin*iy-1+1e-6);
		int yhigh = h3nocuts->GetZaxis()->FindBin(2.0/(double)nybin*(iy+1)-1-1e-6);
		cout<<"cenbin "<<centbin<<" iy "<<2.0/(double)nybin*iy-1<<"\t"<<2.0/(double)nybin*(iy+1)-1<<endl;
		hnocuts[cent][iy] = h3nocuts->ProjectionX(Form("hnocent%dy%d",cent), centbin, centbin, ylow,yhigh ,"e");
		hwithallcuts[cent][iy] = h3allcuts->ProjectionX(Form("hallcent%dy%d",cent), centbin, centbin, ylow, yhigh,"e");
		hnocuts[cent][iy]->SetDirectory(0);
		hwithallcuts[cent][iy]->SetDirectory(0);
		heff[cent][iy]=(TH1D*)hwithallcuts[cent][iy]->Clone(Form("heffcent%dy%d",cent,iy));
		heff[cent][iy]->Divide(hnocuts[cent][iy]);
		heff[cent][iy]->Write();
	}
}
effout->Close();
}
