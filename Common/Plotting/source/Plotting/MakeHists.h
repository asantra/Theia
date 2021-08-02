// This class allows to create and manage histograms of similar type.
// E.g. create 128 identical histograms to study signal of each channel of one APV chip.

#ifndef MAKEHISTS_h
#define MAKEHISTS_h

// #if !defined(__CINT__) || defined(__MAKECINT__)
//#if !defined(__CLING__) || defined(__ROOTCLING__)  // for root 6

#include <string>
#include <vector>
#include <map>

#include <Rtypes.h>

//#endif
 


class TH1;
class TLegend;
class TCanvas;


class MakeHists{
public:  
  MakeHists();
  virtual ~MakeHists();
  Int_t AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1); // 1D
  Int_t AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1,  // 2D
                                                               Int_t n2, Double_t a2, Double_t b2);
  Int_t AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1,
                                                               Int_t n2, Double_t a2, Double_t b2,
                                                               Int_t n3, Double_t a3, Double_t b3); // 3D
  
  Int_t AddHists(const std::string hmapid, const Int_t nhists);
  Int_t AddHists(const std::string hmapid, TH1 *hist, const Int_t posid);
  Int_t AddHistsLogBin(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t *a1); // logBin 1D
  Int_t AddHistsLogBin(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1, Int_t n2, Double_t *a2); // logBin 2D

  TH1* GetHist(const std::string hmapid, const Int_t histid);
  std::vector<TH1*>& GetHist(const std::string hmapid) { return fhist_map[hmapid]; }
  
  void FillHist(const std::string hmapid, Int_t histid, Double_t val);
  void SetBinHist(const std::string hmapid, Int_t histid, Int_t binVal, Double_t val);
  void FillHist(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2);

  void FillHistW(const std::string hmapid, Int_t histid, Double_t val, Double_t weight);
  void FillHistW(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2, Double_t weight);
  void FillHistW(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2, Double_t val3, Double_t weight);

  Int_t SaveHists (const std::string &fname);
  Int_t SaveHists (const std::vector<std::string> hmapids, const std::string &fname);
  Int_t LoadFromFile(const std::string &fname);
  
  void SetRange(const std::string hmapid, const Double_t val1, const Double_t val2, const Int_t axis = 0);
  void SetLineColor(const std::string hmapid, const Int_t clr = 0);
  void SetLineWidth(const std::string hmapid, const Int_t wdth = 1);
  void SetMarkerSize(const std::string hmapid, const Double_t msize = 1.0);
  void Scale(const std::string hmapid, const Double_t x);
  double Integral(const std::string hmapid);
  static int SetLog (TCanvas *cc, const int axis = 1);
  
  
protected:

  std::map< std::string, std::vector <TH1*> > fhist_map;
  
  std::map< std::string, std::vector<TLegend*> > flgnd_map;
  std::vector<Int_t> clrarr;
  
  Int_t debugl;
  
};


#endif
