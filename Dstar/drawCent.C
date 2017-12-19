#include "mixedDstar.h"
#include "../anaconst.h"
void drawCent(int lowcent=7, int highcent=8){
//	gSystem->Load("mixedmassDstarCent.C");
  TString rtname="../D0withrapiditybin.root";
  char outname[50]="fitparameterCent.txt";
//clear the old result
  fstream in;
  in.open(outname,ios::out | ios::trunc);
//  in.open(outname,ios::out | ios::app);
  in.close();

  for(int i=0;i<nPtBins;i++){
    mixedmassDstarCent(0,PtEdge[i],PtEdge[i+1],lowcent,highcent,outname,rtname);
    cout<<"finish: "<< PtEdge[i]<<"\t"<<PtEdge[i+1] <<endl;
  }

  TFile* file = TFile::Open(rtname);
  TH1F* ref = (TH1F*)file->Get("mh1CentWg")->Clone("mh1CentWg");
  ref->SetDirectory(0);
  file->Close();
  in.open(outname,ios::in);
  if (!in){
  cout<<"cannot open the file!"<<endl;
  return;
  }
  
  Double_t par[3][nPtBins]={0};
  Double_t parerr[3][nPtBins]={0};
  for(int N=0;N<nPtBins;N++){
    if (!in.good()) break;
    for (int i=0;i<3;i++){
	    in>>par[i][N]>>parerr[i][N];
	    cout<<par[i][N]<<"\t"<<parerr[i][N]<<"\t";
	  }
	  cout<<endl;
  }

  TCanvas* c[3];
  TH1D* h[3];
  char histname[3][20]={"Counts","mean(MeV/c^{2})","#sigma (MeV/c^{2})"};
  TString gifname[3]={Form("pic/countscent_%d_%d.gif", lowcent, highcent), Form("pic/meancent_%d_%d.gif", lowcent, highcent), Form("pic/sigmacent_%d_%d.gif", lowcent, highcent)};
  //hist bin is caculate from 1
  TFile* f=new TFile(Form("data/Dstar_hists_%d_%d.root", centedge[highcent+1],centedge[lowcent]),"recreate");

  for (int i=0;i<3;i++){
	  c[i] = new TCanvas(histname[i],histname[i]);
	  h[i] = new TH1D(histname[i],histname[i],nPtBins,PtEdge);
	  h[i]->GetXaxis()->SetTitle("p_{T}(GeV/c)");
	  h[i]->GetYaxis()->SetTitle(histname[i]);

  	//fill the hist
    for (int N=1;N<=nPtBins;N++){
	  if (i==1 || i==2) {
  	  par[i][N-1]=par[i][N-1]*1000;
	    parerr[i][N-1]=parerr[i][N-1]*1000;
	  }
	  h[i]->SetBinContent(N,par[i][N-1]);
	  h[i]->SetBinError(N,parerr[i][N-1]);
	  }
	
	  c[i]->cd();
	  h[i]->Draw();
	  if (i==1) {
		  TLine* l = new TLine(1.5,145.4257,12,145.4257);
		  l->SetLineColor(kGreen);
		  TBox* b = new TBox(1.5,145.4257+0.0017,10,145.4257-0.0017);
	  	//l->SetLineStyle(2);
		  //l->SetLineWidth();
		  b->Draw("same");
		  l->Draw("same");
	  }
  	c[i]->SaveAs(gifname[i]);
	  h[i]->Write();
	}
	long double events =ref->Integral(ref->FindBin(lowcent),ref->FindBin(highcent));
	cout<<events<<endl;
	const double pi =3.1415926535;
	TH1D*	hrawspectra= new TH1D("rawspectra","D* raw spectra;p_{T}(GeV/c);#frac{d^{2}N}{2#pi p_{T}dp_{T}dy}",nPtBins,PtEdge);
  fstream frawyield;
  frawyield.open(Form("data/Dstar_rawyield_%d_%d.txt", centedge[highcent+1], centedge[lowcent]), ios::out | ios::trunc);
  frawyield<<events<<endl;
	for (int N=1;N<=nPtBins;N++){
		double ptcenter =hrawspectra->GetBinCenter(N);
		double width = hrawspectra->GetBinWidth(N);
	  hrawspectra->SetBinContent(N,par[0][N-1]/2.0/width/ptcenter/2.0/pi/events);
	  hrawspectra->SetBinError(N,parerr[0][N-1]/width/2.0/ptcenter/2.0/pi/events);
	  cout<<hrawspectra->GetBinWidth(N)<<endl;
    frawyield<<par[0][N-1]<<endl;
  }
  frawyield.close();

	TCanvas* cc = new TCanvas("cc","cc");
	hrawspectra->Draw();
	cc->SetLogy();
	hrawspectra->Write();
	ref->Write();
	//f->Close();
	//gSystem->Exec(Form("cp histCent_%d_%d.root ../Dstareffi/", centedge[highcent+1],centedge[lowcent]));
}
