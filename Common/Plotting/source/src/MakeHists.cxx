
/// this is the heavily lifted from Oleksandr Borysov



#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm> 
#include <iterator>
#include "stdlib.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TClass.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TPad.h"

#include "TCollection.h"
#include "TList.h"
#include "TF1.h"

#include "TROOT.h"
#include "TMath.h"

#include "MakeHists.h"




MakeHists::~MakeHists()
{
  for (std::map< std::string, std::vector <TH1*> >::iterator hmitr = fhist_map.begin(); hmitr != fhist_map.end(); ++hmitr) {
    for (std::vector <TH1*>::iterator hvitr = hmitr->second.begin(); hvitr != hmitr->second.end(); ++hvitr) {
      if (*hvitr) delete *hvitr;
    }
    hmitr->second.clear();
  }
  fhist_map.clear();

  for (std::map< std::string, std::vector <TLegend*> >::iterator hmitr = flgnd_map.begin(); hmitr != flgnd_map.end(); ++hmitr) {
    for (std::vector <TLegend*>::iterator hvitr = hmitr->second.begin(); hvitr != hmitr->second.end(); ++hvitr) {
      if (*hvitr) delete *hvitr;
    }
    hmitr->second.clear();
  }
  flgnd_map.clear();
}



Int_t MakeHists::SaveHists (const std::string &fname)
{
  if (fname.empty()) {
    std::cout << "Warning: MakeHists::SaveHists: File name is empty. Return.\n";
    return 0;
  }
  TFile *pfile = TFile::Open(fname.c_str(), "RECREATE");
  if (!pfile) {
    std::cout << "Warning: MakeHists::SaveHists: Can not open file " << fname << "! Return!\n";
    return 0;
  }
  TH1D  *hinfo = new TH1D("fhinfo", "fhinfo", fhist_map.size(), 0, fhist_map.size());
  Int_t nh = 0;
  for (std::map< std::string, std::vector <TH1*> >::iterator hmitr = fhist_map.begin(); hmitr != fhist_map.end(); ++hmitr) {
    hinfo->Fill(hmitr->first.c_str(), hmitr->second.size());
    for (std::vector <TH1*>::iterator hvitr = hmitr->second.begin(); hvitr != hmitr->second.end(); ++hvitr) {
      if (*hvitr) (*hvitr)->Write();
    }
    nh+=hmitr->second.size();
  }
  hinfo->Write();
  pfile->Close();
  delete pfile;
  return nh;
}





Int_t MakeHists::AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1)
{ // 1D
  if (nhists < 0) return 0;
  int hfail = 0;
  for (int ii = 0; ii < nhists; ++ii) {
    std::stringstream  hname("");
    hname << hmapid << "_" << ii;
    TH1 *hh = new TH1D(hname.str().c_str(), hname.str().c_str(), n1, a1, b1);
    if (hh) {
      hh->SetDirectory(0);
      hh->Sumw2();
      fhist_map[hmapid].push_back(hh);
    } else  { ++hfail; std::cerr << "Error:  MakeHists::AddHists: Problem to create histogram: " << hname.str() << ". Continue.\n"; }
    if (debugl > 1) std::cout << "Histogram " << hname.str() << "  " << " has been created.\n";
  }
  return nhists - hfail;
}



Int_t MakeHists::AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1,  
                                                                     Int_t n2, Double_t a2, Double_t b2)
{ // 2D
  if (nhists < 0) return 0;
  int hfail = 0;
  for (int ii = 0; ii < nhists; ++ii) {
    std::stringstream  hname("");
    hname << hmapid << "_" << ii;
    TH1 *hh = new TH2D(hname.str().c_str(), hname.str().c_str(), n1, a1, b1, n2, a2, b2);
    if (hh) {
      hh->SetDirectory(0);
      hh->Sumw2();
      fhist_map[hmapid].push_back(hh);
    } else  { ++hfail; std::cerr << "Error:  MakeHists::AddHists: Problem to create histogram: " << hname.str() << ". Continue.\n"; }
  }
  return nhists - hfail;
}



Int_t MakeHists::AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1,  
                                                                     Int_t n2, Double_t a2, Double_t b2,
                                                                     Int_t n3, Double_t a3, Double_t b3)
{ // 3D
  if (nhists < 0) return 0;
  int hfail = 0;
  for (int ii = 0; ii < nhists; ++ii) {
    std::stringstream  hname("");
    hname << hmapid << "_" << ii;
    TH1 *hh = new TH3D(hname.str().c_str(), hname.str().c_str(), n1, a1, b1, n2, a2, b2, n3, a3, b3);
    if (hh) {
      hh->SetDirectory(0);
      hh->Sumw2();
      fhist_map[hmapid].push_back(hh);
    } else  { ++hfail; std::cerr << "Error:  MakeHists::AddHists: Problem to create histogram: " << hname.str() << ". Continue.\n"; }
  }
  return nhists - hfail;
}



Int_t MakeHists::AddHistsLogBin(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t *a1)
{ // 1D log axis
  if (nhists < 0) return 0;
  int hfail = 0;
  for (int ii = 0; ii < nhists; ++ii) {
    std::stringstream  hname("");
    hname << hmapid << "_" << ii;
    TH1 *hh = new TH1D(hname.str().c_str(), hname.str().c_str(), n1, a1);
    if (hh) {
      hh->SetDirectory(0);
      hh->Sumw2();
      fhist_map[hmapid].push_back(hh);
    } else  { ++hfail; std::cerr << "Error:  MakeHists::AddHists: Problem to create histogram: " << hname.str() << ". Continue.\n"; }
  }
  return nhists - hfail;
}



Int_t MakeHists::AddHistsLogBin(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1, Int_t n2, Double_t *a2)
{ // 2D log axis
  if (nhists < 0) return 0;
  int hfail = 0;
  for (int ii = 0; ii < nhists; ++ii) {
    std::stringstream  hname("");
    hname << hmapid << "_" << ii;
    TH1 *hh = new TH2D(hname.str().c_str(), hname.str().c_str(), n1, a1, b1, n2, a2);
    if (hh) {
      hh->SetDirectory(0);
      hh->Sumw2();
      fhist_map[hmapid].push_back(hh);
    } else  { ++hfail; std::cerr << "Error:  MakeHists::AddHists: Problem to create histogram: " << hname.str() << ". Continue.\n"; }
  }
  return nhists - hfail;
}


Int_t MakeHists::AddHists(const std::string hmapid, const Int_t nhists)
{
  if (nhists < 0) return 0;
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    std::cout << "Warning:  MakeHists::AddHists: " << hmapid << " already exists. Nothing was done.\n";
    return 0;
  }
  fhist_map[hmapid].reserve(nhists);
  fhist_map[hmapid].assign(nhists, (TH1D*)0);
  return nhists;
}



Int_t MakeHists::AddHists(const std::string hmapid, TH1 *hist, const Int_t posid)
{
  if (posid < 0 || !hist) return 0;
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    if ( static_cast<unsigned int>(posid) < fhist_map[hmapid].size() ) {
      fhist_map[hmapid][posid] = hist;
      return 1;
    }
  } else {
    std::cout << "Warning:  MakeHists::AddHists: Cannnot assign histogram, collection " << hmapid << " does not exist. Nothing was done.\n";
  }
  return 0;
}




TH1* MakeHists::GetHist(const std::string hmapid, Int_t histid)
{
  std::map< std::string, std::vector <TH1*> >::iterator  hitr = fhist_map.find(hmapid);
  if ( hitr != fhist_map.end() ) {
    if (histid >= 0 && hitr->second.size() > static_cast<unsigned int>(histid))
      return hitr->second[histid];
  } 
  return 0;  
}


void MakeHists::FillHist(const std::string hmapid, Int_t histid, Double_t val) 
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH1") != std::string::npos) hh->Fill(val);
    else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 1D. Continue.\n";
  } else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}

void MakeHists::SetBinHist(const std::string hmapid, Int_t histid, Int_t binVal, Double_t val) 
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH1") != std::string::npos) hh->SetBinContent(binVal, val);
    else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 1D. Continue.\n";
  } else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}

void MakeHists::FillHist(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2)
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH2") != std::string::npos) hh->Fill(val1, val2);
    else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 2D. Continue.\n";
  } else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}



void MakeHists::FillHistW(const std::string hmapid, Int_t histid, Double_t val, Double_t weight) 
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH1") != std::string::npos) hh->Fill(val, weight);
    else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 1D. Continue.\n";
  } else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}


void MakeHists::FillHistW(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2, Double_t weight)
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH2") != std::string::npos) dynamic_cast<TH2*>(hh)->Fill(val1, val2, weight);
    else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 2D. Continue.\n";
  } else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}


void MakeHists::FillHistW(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2, Double_t val3, Double_t weight)
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH3") != std::string::npos) dynamic_cast<TH3*>(hh)->Fill(val1, val2, val3, weight);
    else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 3D. Continue.\n";
  } else  std::cerr << "Error:  MakeHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}



void MakeHists::SetRange(const std::string hmapid, const Double_t val1, const Double_t val2, const Int_t axis)
{
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1D *hh = dynamic_cast<TH1D*>(*itr);
      if (!hh) continue;
      switch (axis) {
        case 0: hh->GetXaxis()->SetRangeUser(val1, val2); break;
        case 1: hh->GetYaxis()->SetRangeUser(val1, val2); break;
        case 2: hh->GetZaxis()->SetRangeUser(val1, val2); break;
      }
    }
  }
}


void MakeHists::SetLineColor(const std::string hmapid, const Int_t clr)
{
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1D *hh = dynamic_cast<TH1D*>(*itr);
      if (!hh) continue;
      hh->SetLineColor(clr);
      hh->SetMarkerColor(clr);
    }
  }
}


void MakeHists::SetLineWidth(const std::string hmapid, const Int_t wdth)
{
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1D *hh = dynamic_cast<TH1D*>(*itr);
      if (!hh) continue;
      hh->SetLineWidth(wdth);
    }
  }
}


void MakeHists::SetMarkerSize(const std::string hmapid, const Double_t msize)
{
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1D *hh = dynamic_cast<TH1D*>(*itr);
      if (!hh) continue;
      hh->SetMarkerSize(msize);
    }
  }
}


void MakeHists::Scale(const std::string hmapid, const Double_t x)
{
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1 *hh = dynamic_cast<TH1*>(*itr);
      if (!hh) continue;
      hh->Scale(x);
    }
  }
}


double MakeHists::Integral(const std::string hmapid)
{ 
  double integral = 0;
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1 *hh = dynamic_cast<TH1*>(*itr);
      if (!hh) continue;
      integral = hh->Integral();
    }
  }
  return integral;
}
