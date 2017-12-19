#include "../drawstyle.h"
#include "../anaconst.h"
void D0spectra(int lcent=7, int hcent=8){
	
	globleStyle();
	myStyle();

	float ptedge[nPtBins+1];
	for (int i=0;i<nPtBins+1;i++){
		ptedge[i]=PtEdge[i];
  }
	
	TH1D* heff = makeefficiencyxiaolong(hcent,lcent, PtEdge, nPtBins);
// TH1D* heff = makeefficiencycent(hcent,lcent, PtEdge, nPtBins);
	heff->SetDirectory(0);
	heff->Draw();
	gPad->Update();
// return;

	double D0counts[nPtBins] = {0};
	double D0error[nPtBins] = {1.e-6};
	ifstream fD0counts;
	fD0counts.open(Form("data/D0tot_%d_%d.txt", centedge[hcent+1], centedge[lcent]));
	// cout<<Form("data/D0tot_%d_%d.txt", centedge[hcent+1], centedge[lcent])<<endl;
  //fD0counts.open(Form("rawY_%d_%d.txt", centedge[hcent+1], centedge[lcent]));
	double events = -1;
	fD0counts>>events;
	cout<<"events: "<<events<<endl;
	for(int i=0;i<nPtBins;i++){
	fD0counts>>D0counts[i]>>D0error[i];
		//fD0counts>>D0counts[i];
	}
	fD0counts.close();
 
	TH1D* hrawyield = new TH1D("D0rawyield", "D0rawyield;p_{T}(GeV/c);Counts",nPtBins,PtEdge);
	TH1D* heffyield = new TH1D("D0effyield", "D0effyield;p_{T}(GeV/c);Counts",nPtBins,PtEdge);
	TH1D* h = new TH1D("D0spectra", "D0Spectra;p_{T}(GeV/c);#frac{dN^{2}}{2#pi p_{T}dp_{T}dy}",nPtBins,PtEdge);
	h->SetDirectory(0);
	//TFile* fD0spectrum = TFile::Open("D0_Spectra_Run14_HFT_beforePtShift.root");
	TFile* fD0spectrum = TFile::Open("D0_Spectra_Run14_HFT_beforePtShift.reyield.root");
	//TFile* fD0spectrum = TFile::Open("D0_Spectra_Run14HFT_4Yuanjing_1027.root");
	TGraphErrors* g = (TGraphErrors*)fD0spectrum->Get(Form("gD0_sys_%d_%d",centedge[hcent+1],centedge[lcent]));
	TGraphErrors* gerr = (TGraphErrors*)fD0spectrum->Get(Form("gD0_err_%d_%d",centedge[hcent+1],centedge[lcent]));
	g->SetMarkerColor(kRed);
	// TF1* flvey0_10 = (TF1*)fD0spectrum->Get(Form("flevy_%d_%d",centedge[hcent+1],centedge[lcent]));
	// flevy_0_10->Draw("same");
	// flevy_0_10->SetLineColor(kRed);
	fD0spectrum->Close();
  TGraphErrors* gratio;
  double ratiopt[nPtBins];
  double ratioerr[nPtBins];
  double d0ratio[nPtBins];;

	TCanvas* c = new TCanvas("spectra","spectra",600,800);
	c->Divide(1,2); 
  c->cd(1);
  gPad->SetPad(0, 0.35,1,1);
  setPad(0.15,0.05,0.1,0);
  c->cd(2);
  gPad->SetPad(0,0,1,0.35);
  setPad(0.15,0.05,0,0.2);
	
  fstream yieldtxt;
  yieldtxt.open(Form("data/yield_%d_%d.txt", centedge[hcent+1], centedge[lcent]), ios::out | ios::trunc);
  yieldtxt<<events<<endl;
	double pi = TMath::Pi();
  int j = -1, k=0;
  double x=-999, y=-999;
	for (int i =0;i<nPtBins;i++){
		double pt =h->GetBinCenter(i+1);
		double eff= heff->GetBinContent(heff->FindBin(pt));
		//double eff= 1.0;
		double width = h->GetBinWidth(i+1);
	//	cout<<pt<<"\t"<<width<<"\t"<<eff<<endl;	
		double myd0scale = 1./(eff*pt*2.0*2.0*pi*width*(double)events*D0Branch*2.0);
		double myd0 = D0counts[i]*myd0scale;
		//double myd0 = D0counts[i]/eff/pt/2.0/2.0/pi/width/(double)events/0.0388/2.0;
    double myerr = D0error[i]/eff/pt/2.0/2.0/pi/width/(double)events/D0Branch/2.0;
    yieldtxt<<D0counts[i]/eff<<endl;
		h->SetBinContent(i+1,myd0);
		h->SetBinError(i+1,myerr);
    hrawyield->SetBinContent(i+1,D0counts[i]);
    hrawyield->SetBinError(i+1,D0error[i]);
    heffyield->SetBinContent(i+1,D0counts[i]/eff);
    heffyield->SetBinError(i+1,D0error[i]/eff);
    while (x<=pt && j<g->GetN()){
      if (x==pt){ 
      	//	cout<<x<<"\t"<<y<<endl;
		  double d0err =gerr->GetErrorY(j);  //Tgraph bin start from 0
      ratiopt[k]=pt;
      d0ratio[k]=myd0/y;
      ratioerr[k]=sqrt(myerr*myerr/y/y+myd0*myd0/y/y/y/y*d0err*d0err);
      k++;
      }
	    j++;
      g->GetPoint(j,x,y);
    }
	}
	yieldtxt.close();
	c->cd(1);
	gPad->SetLogy();
	h->Draw("p");
	 g->Draw("psame");
	TLegend* l =new TLegend(0.65,0.7,0.9,0.85);
	l->SetHeader(Form("Au+Au %d~%d%s",centedge[hcent+1],centedge[lcent],"%"));
	l->AddEntry(h, "my D^{0}","lpe" );
	l->AddEntry(g,"Guannan's D^{0}", "lpe");
	l->Draw();

//draw D0 ratio
	c->cd(2);
  gratio = new TGraphErrors(k,ratiopt,d0ratio,0,ratioerr); 
	gratio->GetYaxis()->SetNdivisions(50205);
	gratio->GetYaxis()->SetLabelSize(0.07);
	gratio->GetXaxis()->SetLabelSize(0.07);
	gratio->GetYaxis()->SetRangeUser(0.87,1.23);
	gratio->GetXaxis()->SetRangeUser(PtEdge[0],PtEdge[nPtBins]);
	gratio->Draw("AP");
	drawLine(PtEdge[0], 1, PtEdge[nPtBins], 1, 0.5, 2 ,kBlue);
	TFile* fout = new TFile(Form("data/hD0_myspectra_%d_%d.root",centedge[hcent+1],centedge[lcent]),"recreate");
	h->Write();
	hrawyield->Write();
	heffyield->Write();
//	gPad->SaveAs(Form("D0_myspectra_%d_%d.gif", centedge[hcent+1],centedge[lcent]));

	c->SaveAs(Form("pic/D0_myspectra_%d_%d.gif",centedge[hcent+1],centedge[lcent]));
  fout->Close();
}

TH1D* makeefficiencycent(int hcent=8,int lcent=7,const double* PtEdge,const int nPtBins){
	 const float Nbin[9] = {12.91047, 29.32408, 62.25545, 120.65795, 215.95050, 360.58912, 579.89409, 812.73278, 1042.75372};
  TFile* f = TFile::Open("eff1026.root");
  TFile* fvtx = TFile::Open("vtxCorr_default.root");
  TH1D* hnocuts[9][10];
  TH1D* hwithallcuts[9][10];
  TH1D* hvtx[9];
  TH1D* heff;
  TH1D* hnocutsadd;
  TH1D* hallcutsadd;
  TH1D* hvtxadd;
  int const centedge[10]={80,70,60,50,40,30,20,10,5,0};
  for (int cent=0;cent<9;cent++){
	for(int y=0;y<10;y++){
		hnocuts[cent][y] = (TH1D*)f->Get(Form("minPt300MeV/hNoCuts_minPt300_cent%i_y%i",cent,y))->Clone(Form("minPt300MeV/hNoCuts_minPt300_cent%i_y%i",cent,y));
		hwithallcuts[cent][y] = (TH1D*)f->Get(Form("minPt300MeV/hTpcHftTopo_minPt300_cent%i_y%i",cent,y))->Clone(Form("minPt300MeV/hTpcHftTopo_minPt300_cent%i_y%i",cent,y));
		hnocuts[cent][y]->SetDirectory(0);
		hwithallcuts[cent][y]->SetDirectory(0);
	}
	hvtx[cent] = (TH1D*)fvtx->Get(Form("mean_%d_%d",centedge[cent+1],centedge[cent]))->Clone(Form("mean_%d_%d",centedge[cent+1],centedge[cent]));
	hvtx[cent]->SetDirectory(0);
  }
  hnocutsadd=(TH1D*)hnocuts[0][0]->Clone("nocutsadd");
  hallcutsadd=(TH1D*)hwithallcuts[0][0]->Clone("allcutsadd");
  hvtxadd=(TH1D*)hvtx[0]->Clone("vtxadd");
  hnocutsadd->Scale(0);
  hallcutsadd->Scale(0);
  hvtxadd->Scale(0);
  double tot=0;
  for (int cent=lcent;cent<hcent+1;cent++)
  {
	for(int y=0;y<10;y++){
		hnocutsadd->Add(hnocuts[cent][y], Nbin[cent]);
		hallcutsadd->Add(hwithallcuts[cent][y], Nbin[cent]);
	//	cout<<"add "<<cent<<" "<<y<<endl;
	}
	hvtxadd->Add( hvtx[cent], Nbin[cent]);
	tot+=Nbin[cent];
  }
  hvtxadd->Scale(1./tot);
  //rebin the hist
  hnocutsadd = (TH1D*)hnocutsadd->Rebin(nPtBins,"nocutsadd",PtEdge);
  hallcutsadd = (TH1D*)hallcutsadd->Rebin(nPtBins,"withcutsadd",PtEdge);
  heff=(TH1D*)hallcutsadd->Clone(Form("heffcent_%d_%d",centedge[hcent+1],centedge[lcent]));
  heff->SetTitle(Form("heffcent_%d_%d",centedge[hcent+1],centedge[lcent]));
  heff->Divide(hnocutsadd);
  int Nptbins=heff->GetNbinsX();
  double tmp;
  for (int i=1;i<Nptbins+1;i++){
	tmp = heff->GetBinContent(i)*hvtxadd->GetBinContent(hvtxadd->FindBin(heff->GetBinCenter(i)));
	heff->SetBinContent(i, tmp);
  }
  heff->SetDirectory(0);
  f->Close();
  return heff;
  //return hvtxadd;
}

TH1D* makeefficiencyxiaolong(int hcent=8,int lcent=7,const double* PtEdge,const int nPtBins){ 
  //TFile* f = TFile::Open("D0_eff_combine_newPID_newCuts.root");
  //TFile* f = TFile::Open("D0_eff_combine_newPID.root");
  //TFile* f = TFile::Open("D0_eff_combine_newPID_woTheta.root");
  TFile* f = TFile::Open("D0.extend.root");
  //TFile* f = TFile::Open("D0_eff_combine_newPID_newCuts_1027DoubleStat.root");
  //TFile* fvtx = TFile::Open("vtxCorr_default.root");
  TFile* fvtx = TFile::Open("vtxCorr_pt0_3.root");
  TH2* h2nocuts = (TH2*)f->Get("h2Pt")->Clone("h2PtCent");
  TH2* h2allcuts = (TH2*)f->Get("h2PtCut_pt1")->Clone("h3PtCentYCut_pt1");
  //TH2* h2allcuts = (TH2*)f->Get("h2PtCut")->Clone("h2PtCentCut");
  h2nocuts->SetDirectory(0);
  h2allcuts->SetDirectory(0);
  TH1D* hnocuts[9];
  TH1D* hwithallcuts[9];
  TH1D* hvtx[9];
  TH1D* heff;
  TH1D* hnocutsadd;
  TH1D* hallcutsadd;
  TH1D* heffrebinadd;
  TH1D* hvtxadd;
  TH1D* heffrebin[9];
  TH1D* hnocutsrebin[9];
  TH1D* hwithcutsrebin[9];
  int const centedge[10]={80,70,60,50,40,30,20,10,5,0};
  //int etalow = h3nocuts->GetZaxis()->FindBin(-1+1e-6);
  //int etahigh = h3nocuts->GetZaxis()->FindBin(1-1e-6);
  //cout<<"eta"<<etahigh<<"\t"<<etalow<<endl;
  for (int cent=0;cent<9;cent++){ 
    int centbin = h2nocuts->GetYaxis()->FindBin(cent+0.5); 
    // cout<<"cetbin"<<centbin<<endl;
    hnocuts[cent] = h2nocuts->ProjectionX(Form("hno%d",cent), centbin, centbin ,"e");
    hwithallcuts[cent] = h2allcuts->ProjectionX(Form("hwithcut%d",cent), centbin, centbin ,"e");
    //hnocuts[cent] = h3nocuts->ProjectionX(Form("hno%d",cent), centbin, centbin, etalow,etahigh ,"e");
    //hwithallcuts[cent] = h3allcuts->ProjectionX(Form("hall%d",cent), centbin, centbin, etalow, etahigh,"e");
    hnocuts[cent]->SetDirectory(0);
    hwithallcuts[cent]->SetDirectory(0);
    hnocutsrebin[cent] = (TH1D*)hnocuts[cent]->Rebin(nPtBins,Form("hnocutsrebin%d",cent),PtEdge);
    hwithallcuts[cent] = (TH1D*)hwithallcuts[cent]->Rebin(nPtBins,Form("hwthcutsrebin%d",cent),PtEdge);
    heffrebin[cent]=(TH1D*)hwithallcuts[cent]->Clone(Form("heffrebin%d",cent));
    heffrebin[cent]->Divide(hnocutsrebin[cent]);
    
    hvtx[cent] = (TH1D*)fvtx->Get(Form("mean_%d_%d",centedge[cent+1],centedge[cent]))->Clone(Form("mean_%d_%d",centedge[cent+1],centedge[cent]));
    hvtx[cent]->SetDirectory(0);
  }
  hnocutsadd=(TH1D*)hnocuts[0]->Clone("nocutsadd");
  hallcutsadd=(TH1D*)hwithallcuts[0]->Clone("allcutsadd");
  heffrebinadd=(TH1D*)heffrebin[0]->Clone("effrebinadd");
  hvtxadd=(TH1D*)hvtx[0]->Clone("vtxadd");
  hnocutsadd->Scale(0);
  hallcutsadd->Scale(0);
  hvtxadd->Scale(0);
  heffrebinadd->Scale(0);
  hallcutsadd->Draw();
  gPad->Update();
  double tot=0;
  for (int cent=lcent;cent<hcent+1;cent++)
  {
    hnocutsadd->Add(hnocuts[cent], Nbin[cent]);
    hallcutsadd->Add(hwithallcuts[cent], Nbin[cent]);
    hvtxadd->Add( hvtx[cent], Nbin[cent]);
    heffrebinadd->Add(heffrebin[cent],Nbin[cent]);
    tot+=Nbin[cent];
  }
  hvtxadd->Scale(1./tot);
  heffrebinadd->Scale(1.0/tot);
  //rebin the hist
  hnocutsadd = (TH1D*)hnocutsadd->Rebin(nPtBins,"nocutsadd",PtEdge);
  hallcutsadd = (TH1D*)hallcutsadd->Rebin(nPtBins,"withcutsadd",PtEdge);
  heff=(TH1D*)hallcutsadd->Clone(Form("heffcent_%d_%d",centedge[hcent+1],centedge[lcent]));
  heff->SetTitle(Form("heffcent_%d_%d",centedge[hcent+1],centedge[lcent]));
  heff->Divide(hnocutsadd);
  int Nptbins=heff->GetNbinsX();
  double tmp, vtxeff=0;
  TH1D* neweff = new TH1D("neweff", "neweff", nPtBins, PtEdge);
  for (int i=1;i<Nptbins+1;i++){
    if (heff->GetBinCenter(i)>7) vtxeff = hvtxadd->GetBinContent(hvtxadd->FindBin(7-1e-6));
    else vtxeff = hvtxadd->GetBinContent(hvtxadd->FindBin(heff->GetBinCenter(i))); 
    tmp = heff->GetBinContent(i)*vtxeff;
    heff->SetBinContent(i, tmp);
    neweff->SetBinContent(i, tmp);
    tmp = heffrebinadd->GetBinContent(i)*vtxeff;
    //heffrebinadd->SetBinContent(i,tmp);
  }
  heff->SetDirectory(0);
  f->Close();
  return heff;
  //return heffrebinadd;
  //return hvtxadd;
}
