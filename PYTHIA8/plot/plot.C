#include <stdlib.h>
#include <stdio.h>

#include "TStyle.h"
#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"
#include "TSystem.h"
#include "TString.h"
#include "TPaveStats.h"

int plot(){
  gStyle->SetOptStat(1);
  TString files[5] = {
    "CUETP8M1",
    "CUETP8M2T4",
    "CP5",
    "Vanilla",
    "PYTHIA6"
  };
  gSystem->Exec("mkdir comp");
  for(int i=0; i<5; i++){
    gSystem->Exec(Form("mkdir %s ",files[i].Data()));
    gSystem->Exec(Form("mkdir %s ",Form("%s/hist",files[i].Data())));
  }

  // plot everything
  for(int i=0; i<5; i++){
    TString filename = "../output/" + files[i] + ".root";
    TFile* f = TFile::Open(filename);
    TIter next(f->GetListOfKeys());
    TKey* key;
    while((key = (TKey*)next())){
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1D")) continue;
      TH1D *h = (TH1D*)key->ReadObj();
      h->SetTitle(files[i].Data());
      TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
      c1->cd();
      h->Draw();
      gPad->Update();
      TPaveStats* sb=(TPaveStats*) h->FindObject("stats");
      sb->SetX1NDC(0.5);
      sb->SetX2NDC(0.7);
      sb->SetY1NDC(.75);
      sb->SetY2NDC(.9);
      gPad->Modified();
      TString name = h->GetName();
      c1->SaveAs("./" + files[i]+ "/hist/" + name + ".pdf");
      delete c1;
    }
  }


  TH1D *h_jet_spectrum[5][4] = {0};
  TH1D* jetpt_weighted_comp[5] = {0};

  //plot jet spectrum
  for(int i=0; i<5; i++){
    TString filename = "../output/" + files[i] + ".root";
    TFile* f = TFile::Open(filename);

      if(i!=4){
      f->GetObject("jetpt_flatPthat_Unweighted_h",h_jet_spectrum[i][0]);
      f->GetObject("jetpt_flatPthat_Weighted_h",h_jet_spectrum[i][1]);
    }
    f->GetObject("jetpt_stagPthat_Unweighted_h",h_jet_spectrum[i][2]);
    f->GetObject("jetpt_stagPthat_Weighted_h",h_jet_spectrum[i][3]);
    for(int j=0; j<4; j++) {
      if(i==4&&(j==0||j==1)) continue;
      h_jet_spectrum[i][j]->Rebin(5);
      h_jet_spectrum[i][j]->SetTitle(files[i].Data());
      h_jet_spectrum[i][j]->Sumw2();
      double binsize = h_jet_spectrum[i][j]->GetXaxis()->GetBinWidth(1);
      double integral = h_jet_spectrum[i][j]->Integral();
      h_jet_spectrum[i][j]->Scale(1/binsize/integral);
      TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
      c1->cd();
      gPad->SetLogy();
      h_jet_spectrum[i][j]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{dp_{T}}");
      h_jet_spectrum[i][j]->GetYaxis()->SetTitleOffset(1.3);
      h_jet_spectrum[i][j]->Draw();
      TString name = h_jet_spectrum[i][j]->GetName();
      c1->SaveAs("./" + files[i]+ "/" + name + "spectrum_.pdf");
      delete c1;
    }
    if(i==4) continue;
    TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
    c2->cd();
    gStyle->SetOptStat(0);
    jetpt_weighted_comp[i] = (TH1D*) h_jet_spectrum[i][1]->Clone();
    jetpt_weighted_comp[i]->Sumw2();
    jetpt_weighted_comp[i]->Divide(h_jet_spectrum[i][3]);
    jetpt_weighted_comp[i]->GetYaxis()->SetRangeUser(0.7,2);
    jetpt_weighted_comp[i]->GetYaxis()->SetTitle("Flat/Stag");
    jetpt_weighted_comp[i]->Draw();
    c2->SaveAs("./" + files[i]+ "/" + "jetpt_spectrum_flatoverstag_.pdf");
    delete c2;
  }
  TH1D* h_dijet[5][4] = {0};
  TH1D* dijet_flatoverstag[4] = {0};
  TH1D* dphi_flatoverstag[4] = {0};

  //plot dphi and xj distribution
  gStyle->SetOptStat(0);
  for(int i=0; i<5; i++){
    TString filename = "../output/" + files[i] + ".root";
    TFile* f = TFile::Open(filename);
    f->GetObject("dphi_flatPthat_Weighted_h",h_dijet[i][0]);
    f->GetObject("dphi_stagPthat_Weighted_h",h_dijet[i][1]);
    f->GetObject("dijet_flatPthat_Weighted_h",h_dijet[i][2]);
    f->GetObject("dijet_stagPthat_Weighted_h",h_dijet[i][3]);
    for(int j=0; j<4; j++) {
      if(i==4&&(j==0||j==2)) continue;
      h_dijet[i][j]->SetTitle(files[i].Data());
      h_dijet[i][j]->Sumw2();
      double binsize = h_dijet[i][j]->GetXaxis()->GetBinWidth(1);
      double integral = h_dijet[i][j]->Integral();
      h_dijet[i][j]->Scale(1/binsize/integral);
      TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
      c1->cd();
      if(j==0 || j==1) {
      gPad->SetLogy();
      h_dijet[i][j]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{d#phi}");
      h_dijet[i][j]->GetYaxis()->SetRangeUser(5E-3,5);
      }
      if(j==2 || j==3) {
      h_dijet[i][j]->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{dx_{J}}");
      h_dijet[i][j]->GetYaxis()->SetRangeUser(-0.2,3.2);
      }
      h_dijet[i][j]->GetYaxis()->SetTitleOffset(1.3);
      h_dijet[i][j]->Draw();
      TString name = h_dijet[i][j]->GetName();
      c1->SaveAs("./" + files[i]+ "/" + name + "spectrum_.pdf");
      delete c1;
    }

    // plot flat over stag
    if(i==4) continue;
    TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
    c2->cd();
    dijet_flatoverstag[i] = (TH1D*) h_dijet[i][2]->Clone();
    dijet_flatoverstag[i]->Sumw2();
    dijet_flatoverstag[i]->Divide(h_dijet[i][3]);
    dijet_flatoverstag[i]->GetYaxis()->SetRangeUser(0.6,1.4);
    dijet_flatoverstag[i]->Draw();
    c2->SaveAs("./" + files[i]+ "/" + "dijet_flatoverstag_.pdf");
    delete c2;

    TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
    c3->cd();
    dphi_flatoverstag[i] = (TH1D*) h_dijet[i][0]->Clone();
    dphi_flatoverstag[i]->Sumw2();
    dphi_flatoverstag[i]->Divide(h_dijet[i][1]);
    dphi_flatoverstag[i]->GetYaxis()->SetRangeUser(0.6,1.4);
    dphi_flatoverstag[i]->Draw();
    c3->SaveAs("./" + files[i]+ "/" + "dphi_flatoverstag_.pdf");
    delete c3;
  }

    // plot in one frame
    for(int j=0; j<4; j++){
      TCanvas* c_dijet = new TCanvas("c_dijet", "c_dijet", 800, 600);
      TLegend* legend = new TLegend(0.15,0.75,0.35,0.9);
      // gPad->SetLogy();
      for(int i=0; i<5; i++){
        if(i==4&&(j==0||j==2)) continue;
        if(j==0||j==1) gPad->SetLogy();
        c_dijet->cd();
        h_dijet[i][j]->SetTitle("");
        h_dijet[i][j]->Draw("same");
        h_dijet[i][j]->SetMarkerStyle(20+i);
        h_dijet[i][j]->SetMarkerColor(i==4?i+2:1+i);
        h_dijet[i][j]->SetLineColor(i==4?i+2:1+i);
        legend->AddEntry(h_dijet[i][j],Form("%s",files[i].Data()));
      }
      legend->Draw("same");
      c_dijet->SaveAs(Form("./comp/%d.pdf",j));

      delete c_dijet;
    }


    // jet specrrum compare

    // i
    //CUETP8M1 0 -- base
    //CUETP8M2T4 1
    // CP5 2
    //Vanilla 3
    // pythia6 -- no flat

    // j
    //flat weighted 1
    //stag weighted 3

    gStyle->SetOptStat(0);
    TCanvas* c_spectrum_flat = new TCanvas("c_spectrum_flat", "c_spectrum_flat", 800, 600);
    TCanvas* c_spectrum_stag = new TCanvas("c_spectrum_stag", "c_spectrum_stag", 800, 600);
    TCanvas* c_dphi_flat = new TCanvas("c_dphi_flat", "c_dphi_flat", 800, 600);
    TCanvas* c_dphi_stag = new TCanvas("c_dphi_stag", "c_dphi_stag", 800, 600);
    TCanvas* c_dijet_flat = new TCanvas("c_dijet_flat", "c_dijet_flat", 800, 600);
    TCanvas* c_dijet_stag = new TCanvas("c_dijet_stag", "c_dijet_stag", 800, 600);
    TCanvas* c_dijet_stag_pythia6 = new TCanvas("c_dijet_stag_pythia6", "c_dijet_stag_pythia6", 800, 600);
    TLegend* l_spectrum = new TLegend(0.6,0.75,0.9,0.9);
    TLegend* l_spectrum2 = new TLegend(0.6,0.75,0.9,0.9);
    TLegend* l_spectrum3 = new TLegend(0.6,0.75,0.9,0.9);
    TLine *line1 = new TLine(0,1,1,1);
    line1->SetLineStyle(2);
    TLine *line2 = new TLine(0,1,3.141592,1);
    line2->SetLineStyle(2);

    // jet spectrum over CUETP8M1 flat in one frame -- no pythia6
    for(int i=1; i<5; i++){
      if(i!=4){
      TCanvas* c4 = new TCanvas("c4", "c4", 800, 600);
      c4->cd();
      TH1D* jetpt_flat = (TH1D*) h_jet_spectrum[i][1]->Clone();
      jetpt_flat->Sumw2();
      jetpt_flat->Divide(h_jet_spectrum[0][1]);
      jetpt_flat->GetYaxis()->SetRangeUser(0.6,1.4);
      jetpt_flat->GetXaxis()->SetTitle("inclusive jet p_{T} (Flat Weighted)");
      jetpt_flat->GetYaxis()->SetTitle(Form("%s/CUETP8M1",files[i].Data()));
      jetpt_flat->Draw();
      c4->SaveAs("./comp/jetpt_flat_" + files[i]+ "_CUETP8M1.pdf");
      c_spectrum_flat->cd();
      jetpt_flat->GetYaxis()->SetTitle("Tune/CUETP8M1");
      jetpt_flat->SetTitle("");
      jetpt_flat->SetMarkerStyle(20+i);
      jetpt_flat->SetMarkerColor(i==4?i+2:1+i);
      jetpt_flat->SetLineColor(i==4?i+2:1+i);
      jetpt_flat->Draw("same");
      l_spectrum->AddEntry(jetpt_flat, Form("%s/CUETP8M1",files[i].Data()));
      l_spectrum2->AddEntry(jetpt_flat, Form("%s/CUETP8M1",files[i].Data()));
      l_spectrum->Draw("same");
    }

    // xj distribution over CUETP8M1 flat in one frame -- no pythia6
    if(i!=4){
      TCanvas* c6 = new TCanvas("c6", "c6", 800, 600);
      c6->cd();
      TH1D* dijet_flat = (TH1D*) h_dijet[i][2]->Clone();
      dijet_flat->Sumw2();
      dijet_flat->Divide(h_dijet[0][2]);
      dijet_flat->GetYaxis()->SetRangeUser(0.6,1.4);
      dijet_flat->GetXaxis()->SetTitle("x_{J} (Flat Weighted)");
      dijet_flat->GetYaxis()->SetTitle(Form("%s/CUETP8M1",files[i].Data()));
      dijet_flat->Draw();
      TLine *line1 = new TLine(0,1,1,1);
      line1->SetLineStyle(2);
      line1->Draw("same");
      c6->SaveAs("./comp/dijet_flat_" + files[i]+ "_CUETP8M1.pdf");
      c_dijet_flat->cd();
      dijet_flat->GetYaxis()->SetTitle("Tune/CUETP8M1");
      dijet_flat->SetTitle("");
      dijet_flat->SetMarkerStyle(20+i);
      dijet_flat->SetMarkerColor(i==4?i+2:1+i);
      dijet_flat->SetLineColor(i==4?i+2:1+i);
      dijet_flat->Draw("same");
      l_spectrum->Draw("same");
      line1->Draw("same");
    }



    // dphi distribution over CUETP8M1 flat in one frame -- no pythia6
    if(i!=4){
      TCanvas* c9 = new TCanvas("c9", "c9", 800, 600);
      c9->cd();
      gPad->SetLogy();
      TH1D* dphi_flat = (TH1D*) h_dijet[i][0]->Clone();
      dphi_flat->Sumw2();
      dphi_flat->Divide(h_dijet[0][0]);
      dphi_flat->GetYaxis()->SetRangeUser(0.6,1.4);
      dphi_flat->GetXaxis()->SetTitle("|#Delta#phi| (Flat Weighted)");
      dphi_flat->GetYaxis()->SetTitle(Form("%s/CUETP8M1",files[i].Data()));
      dphi_flat->Draw();
      line2->Draw("same");
      c9->SaveAs("./comp/dphi_flat_" + files[i]+ "_CUETP8M1.pdf");
      c_dphi_flat->cd();
      dphi_flat->GetYaxis()->SetTitle("Tune/CUETP8M1");
      dphi_flat->SetTitle("");
      dphi_flat->SetMarkerStyle(20+i);
      dphi_flat->SetMarkerColor(i==4?i+2:1+i);
      dphi_flat->SetLineColor(i==4?i+2:1+i);
      dphi_flat->Draw("same");
      l_spectrum->Draw("same");
    }

    // jet spectrum over CUETP8M1 stag in one frame
    TCanvas* c5 = new TCanvas("c5", "c5", 800, 600);
    c5->cd();
    TH1D* jetpt_stag = (TH1D*) h_jet_spectrum[i][3]->Clone();
    jetpt_stag->Sumw2();
    jetpt_stag->Divide(h_jet_spectrum[0][3]);
    jetpt_stag->GetYaxis()->SetRangeUser(0.6,1.4);
    jetpt_stag->GetXaxis()->SetTitle("inclusive jet p_{T} (Stag Weighted)");
    jetpt_stag->GetYaxis()->SetTitle(Form("%s/CUETP8M1",files[i].Data()));
    jetpt_stag->Draw();
    c5->SaveAs("./comp/jetpt_stag_" + files[i]+ "_CUETP8M1.pdf");
    c_spectrum_stag->cd();
    jetpt_stag->GetYaxis()->SetTitle("Tune/CUETP8M1");
    jetpt_stag->SetTitle("");
    jetpt_stag->SetMarkerStyle(20+i);
    jetpt_stag->SetMarkerColor(i==4?i+2:1+i);
    jetpt_stag->SetLineColor(i==4?i+2:1+i);
    jetpt_stag->Draw("same");
    if(i==4) l_spectrum2->AddEntry(jetpt_stag, Form("%s/CUETP8M1",files[4].Data()));
    l_spectrum2->Draw("same");



    // dphi distribution over CUETP8M1 stag in one frame
    TCanvas* c8 = new TCanvas("c8", "c8", 800, 600);
    c8->cd();
    gPad->SetLogy();
    TH1D* dphi_stag = (TH1D*) h_dijet[i][1]->Clone();
    dphi_stag->Sumw2();
    dphi_stag->Divide(h_dijet[0][1]);
    dphi_stag->GetYaxis()->SetRangeUser(0.6,1.4);
    dphi_stag->GetXaxis()->SetTitle("|#Delta#phi| (Stag Weighted)");
    dphi_stag->GetYaxis()->SetTitle(Form("%s/CUETP8M1",files[i].Data()));
    dphi_stag->Draw();
    line2->Draw("same");
    c8->SaveAs("./comp/dphi_stag_" + files[i]+ "_CUETP8M1.pdf");
    c_dphi_stag->cd();
    dphi_stag->GetYaxis()->SetTitle("Tune/CUETP8M1");
    dphi_stag->SetTitle("");
    dphi_stag->SetMarkerStyle(20+i);
    dphi_stag->SetMarkerColor(i==4?i+2:1+i);
    dphi_stag->SetLineColor(i==4?i+2:1+i);
    dphi_stag->Draw("same");
    l_spectrum2->Draw("same");



    // xj distribution over CUETP8M1 stag in one frame
    TCanvas* c7 = new TCanvas("c7", "c7", 800, 600);
    c7->cd();
    TH1D* dijet_stag = (TH1D*) h_dijet[i][3]->Clone();
    dijet_stag->Sumw2();
    dijet_stag->Divide(h_dijet[0][3]);
    dijet_stag->GetYaxis()->SetRangeUser(0.6,1.4);
    dijet_stag->GetXaxis()->SetTitle("x_{J} (Stag Weighted)");
    dijet_stag->GetYaxis()->SetTitle(Form("%s/CUETP8M1",files[i].Data()));
    dijet_stag->Draw();
    line1->Draw("same");
    c7->SaveAs("./comp/dijet_stag_" + files[i]+ "_CUETP8M1.pdf");
    c_dijet_stag->cd();
    dijet_stag->GetYaxis()->SetTitle("Tune/CUETP8M1");
    dijet_stag->SetTitle("");
    dijet_stag->SetMarkerStyle(20+i);
    dijet_stag->SetMarkerColor(i==4?i+2:1+i);
    dijet_stag->SetLineColor(i==4?i+2:1+i);
    dijet_stag->Draw("same");
    l_spectrum2->Draw("same");
    line1->Draw("same");
    }

    // xj distribution over PYTHIA6 stag in one frame
    for(int i=0; i<4; i++){
      TCanvas* c7 = new TCanvas("c7", "c7", 800, 600);
      c7->cd();
      TH1D* dijet_stag = (TH1D*) h_dijet[i][3]->Clone();
      dijet_stag->Sumw2();
      dijet_stag->Divide(h_dijet[4][3]);
      dijet_stag->GetYaxis()->SetRangeUser(0.6,1.4);
      dijet_stag->GetXaxis()->SetTitle("x_{J} (Stag Weighted)");
      dijet_stag->GetYaxis()->SetTitle(Form("%s/PYTHIA6",files[i].Data()));
      dijet_stag->Draw();
      line1->Draw("same");
      c7->SaveAs("./comp/dijet_stag_" + files[i]+ "_pythia6.pdf");
      c_dijet_stag_pythia6->cd();
      dijet_stag->GetYaxis()->SetTitle("Tune/CUETP8M1");
      dijet_stag->SetTitle("");
      dijet_stag->SetMarkerStyle(20+i);
      dijet_stag->SetMarkerColor(i==4?i+2:1+i);
      dijet_stag->SetLineColor(i==4?i+2:1+i);
      dijet_stag->Draw("same");
      l_spectrum3->AddEntry(dijet_stag, Form("%s/PYTHIA6",files[i].Data()));
      l_spectrum3->Draw("same");
      line1->Draw("same");
    }

    c_spectrum_flat->SaveAs("./comp/spectrum_flat.pdf");
    c_spectrum_stag->SaveAs("./comp/spectrum_stag.pdf");
    c_dijet_flat->SaveAs("./comp/dijet_flat.pdf");
    c_dijet_stag->SaveAs("./comp/dijet_stag.pdf");
    c_dphi_flat->SaveAs("./comp/dphi_flat.pdf");
    c_dphi_stag->SaveAs("./comp/dphi_stag.pdf");
    c_dijet_stag_pythia6->SaveAs("./comp/dijet_stag_pythia6.pdf");



  return 0;
}
