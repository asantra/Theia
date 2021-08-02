#ifndef DETECTOR_H
#define DETECTOR_H

#include <TArrayD.h>
#include <TList.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include "Utils.h"
#include "Probe.h"
#include "Cluster.h"
#include "Layer.h"
#include "Material.h"
#include <Riostream.h>
#include <TMaterial.h>
#include <TSystem.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <stdio.h>


class NaCardsInput;


class Detector : public TNamed {
 public:
  static const Double_t kMassP;
  static const Double_t kMassK;
  static const Double_t kMassPi;
  static const Double_t kMassMu;
  static const Double_t kMassE;
 public:
  enum {kMagAlice=-1};
  enum {kNoDecay,kDoRealDecay,kApplyDecay};
  Detector(const char *name="test_detector", const char *title="detector");
  Detector();
  ~Detector();
  //
  void         ReadMaterials(const char* fname);
  void         ReadSetup(const char* setup, const char* materials);
  TObjArray*   GetMaterials() const {return (TObjArray*)&fMaterials;}
  Material*  GetMaterial(const char* name) const {return (Material*)fMaterials.FindObject(name);}
  //
  Layer* AddLayer(const char *type, const char *name, Float_t zPos, Float_t radL, Float_t density, Float_t thickness, Float_t xRes=999999, Float_t yRes=999999, Float_t eff=1,Material* mat=0);
  void         AddBeamPipe(Float_t r, Float_t dr, Float_t radL, Float_t density, Material* mat=0);
  BeamPipe*    GetBeamPipe() const {return fBeamPipe;}

  void         ClassifyLayers();
  void         ResetMCTracks(Int_t maxLr);
  //
  void     SetApplyBransonPCorrection(float v = 0.1) {fApplyBransonPCorrection = v;} // set to negative to not apply BP correction
  void     SetIncludeVertex(Bool_t v=kTRUE)       {fIncludeVertex = v;}
  Bool_t   GetIncludeVertex()               const {return fIncludeVertex;}
  //-------------------
  void         CreateProbe(Probe* adr, double pt, double yrap, double phi, double mass, int charge, double x,double y, double z);
  Probe* CreateProbe(double pt, double yrap, double phi, double mass, int charge, double x,double y, double z);
  Probe* PrepareProbe(double pt, double yrap, double phi, double mass, int charge, double x,double y, double z);
  Bool_t SolveSingleTrackViaKalman(double pt, double yrap, double phi,double mass, int charge, double x,double y, double z);
  Bool_t SolveSingleTrack(double pt, double yrap, double phi, double mass, int charge, double x=0,double y=0, double z=0, TObjArray* sumArr=0,int nMC=1, int offset=-1);
  Int_t  PropagateToLayer(Probe* trc, Layer* lrFrom, Layer* lrTo, int dir, Bool_t modeMC=kFALSE);
  Bool_t PropagateToZBxByBz(Probe* trc,double z, double maxDZ=1.0, Double_t xOverX0=0., Double_t xTimesRho=0., Bool_t modeMC=kFALSE);
  Bool_t SolveSingleTrackViaKalmanMC(int offset);
  Bool_t SolveSingleTrackViaKalmanMC_Noam(double pt, double yrap, double phi, double mass, int charge, double x=0,double y=0, double z=0, int offset=-1);
  Bool_t SolveSingleTrackViaKalmanMC_Noam_multiseed(std::vector<TLorentzVector>& pseed, double mass, int charge, int offset=-1, bool doPrint=false);
  Bool_t TransportKalmanTrackWithMS(Probe *probTr, int maxLr, Bool_t bg=kFALSE);
  Int_t GetFieldReg(double z);
  //-------------------
  TList*       GetLayers()                 const {return (TList*)&fLayers;}
  Layer* GetLayer(Int_t i)           const {return (Layer*)fLayers.At(i);}
  Layer* GetLayer(const char *name)  const {return (Layer*)fLayers.FindObject(name);}
  Layer* GetLayerITS(int i)          const {return (Layer*)fLayersITS[i];}
  Layer* GetLayerMS(int i)           const {return (Layer*)fLayersMS[i];}
  Layer* GetLayerTR(int i)           const {return (Layer*)fLayersTR[i];}
  //
  void         SetDefStepAir(Double_t v=1) {fDefStepAir = v;}
  void         SetDefStepMat(Double_t v=1) {fDefStepMat = v;}
  //
  Probe* GetProbe()                  const {return (Probe*)&fProbe;}
  Double_t     GetZDecay()                 const {return fZDecay;}

  Probe* GetMuBransonCorrVtx()    const {
    return fMuTrackBCVertex.GetUniqueID()==0 ? (Probe*)&fMuTrackBCVertex : 0;
  }
  Probe* GetMuBransonCorrLastITS()    const {
    return fMuTrackBCLastITS.GetUniqueID()==0 ? (Probe*)&fMuTrackBCLastITS : 0;
  }
  Probe* GetMuLastITS()    const {
    return fMuTrackLastITS.GetUniqueID()==0 ? (Probe*)&fMuTrackLastITS : 0;
  }
  Probe* GetMuVtx()    const {
    return fMuTrackVertex.GetUniqueID()==0 ? (Probe*)&fMuTrackVertex : 0;
  }
  
  //
  Bool_t       UpdateTrack(Probe* trc, const Layer* lr, const Cluster* cl) const;
  Bool_t       NeedToKill(Probe* probe) const;
  Bool_t       GetUseBackground()          const {return fUseBackground;}
  void         SetUseBackground(Bool_t v=kTRUE)  {fUseBackground = v;}
  void         CheckClusters(int i1, int i2, int i3, int i4);
  void         CheckTrackProlongations(Probe *probe, Layer* lrP, Layer* lr, bool doPrint=false);
  Bool_t       IsCorrect(Probe *probe);
  void         RequirePattern(UInt_t patt);
  void         SetMinITSHits(int n);
  void         SetMinMSHits(int n);
  void         SetMinTRHits(int n);
  void         SetMaxSeedToPropagate(int n) {fMaxSeedToPropagate = n;}
  void         SetErrorScale(double v=500)  {fErrScale = v;}
  void         SetMaxChi2Cl(double v=10)  {fMaxChi2Cl = v;}
  void         SetMaxChi2NDF(double v=10) {fMaxNormChi2NDF = v;}
  void         SetMaxChi2Vtx(double v=20) {fMaxChi2Vtx = v;}
  void         SetMinP2Propagate(double val=1) {fMinP2Propagate = val;}
  void         BookControlHistos();
  //
  void         SetDecayZProfile(TH1* histo)       {fDecayZProf = histo;}
  TH1*         GetDecayZProfile()           const {return fDecayZProf;}
  void         InitDecayZHisto(double absorberLambda=50.0);
  //
  void         GenBgEvent(double x, double y, double z, Int_t offset=0);
  void         InitBgGeneration(int dndeta=300, double y0=2., double sigy=1.2, double ymin=1.8,double ymax=4.,
				double T=0.17, double ptmin=0.01, double ptmax=5);
  void         InitBgGenerationPart(double npi=300, double nKplus=300, double nKminus=300,double nP=300,double ratio=1, 
				    double y0=2., double y0BGPi= 0.8,double y0BGKplus =0.8, double y0BGKminus =0.8,double y0BGP=0.8, 
				    double sigyBGPi=1.2, double sigmayBGKplus=0.8, double sigmayBGKminus=0.8, double sigmaBGP =0.8, double ymin=1.8,double ymax=4.,
				    double Tpi=0.17, double TK=0.23, double TP=0.25, double ptmin=0.01, double ptmax=5);
  void         InitBkg(double beamenergy);

  void         SetExternalInput(Bool_t v = kFALSE)    {fExternalInput = v;}
  Bool_t       GetExternalInput()           const  {return fExternalInput;}

  void         ImposeVertex(float x=0.,float y=0., float z=0.) {
    fImposeVertexPosition = kTRUE;
    double tmpLab[3];
    tmpLab[0] = x;
    tmpLab[1] = y;
    tmpLab[2] = z;
    // assign in tracking frame
    Probe::Lab2Trk(tmpLab, fRefVtx);
  }
  void         UseTrackOriginAsVertex() {
    fImposeVertexPosition = kFALSE;
    fRefVtx[0] = fRefVtx[1] = fRefVtx[2] = 0.;    
  }
  //
  void Print(const Option_t* opt = "") const; 
  //
  Int_t GetNumberOfActiveLayers()    const {return fNActiveLayers;}  
  Int_t GetNumberOfActiveLayersITS() const {return fNActiveLayersITS;}
  Int_t GetNumberOfActiveLayersMS()  const {return fNActiveLayersMS;}
  Int_t GetNumberOfActiveLayersTR()  const {return fNActiveLayersTR;}
  void   SetLastActiveLayerTracked(int lr)   {fLastActiveLayerTracked = lr;}
  Int_t  GetLastActiveLayerTracked() const {return fLastActiveLayerTracked;}
  Int_t  GetLastActiveLayer()        const {return fLastActiveLayer;}
  Int_t  GetLastActiveLayerITS()     const {return fLastActiveLayerITS;}
  void  ForceLastActiveLayer(int lr);

  Double_t GetNCh() const {return fNCh;}
  Double_t GetNChPi() const {return fNChPi;}
  Double_t GetNChK() const {return fNChK;}
  Double_t GetNChP() const {return fNChP;}
  TF1*     GetdNdYPi() const {return fdNdYPi;}
  TF1*     GetdNdYK() const {return fdNdYK;}
  TF1*     GetdNdYP() const {return fdNdYP;}
  TF1*     GetdNdPtPi() const {return fdNdPtPi;}
  TF1*     GetdNdPtK() const {return fdNdPtK;}
  TF1*     GetdNdPtP() const {return fdNdPtP;}

  //
  // Helper functions
  Double_t ThetaMCS                 ( Double_t mass, Double_t RadLength, Double_t momentum ) const;

  Double_t HitDensity(double xLab,double ylab,double zlab) const;
  void     PerformDecay(Probe* trc);
  //
  float GetChi2MuAtVtx() const {return fChi2MuVtx;}
  
  TH1*     GetHChi2Branson()    const {return fHChi2Branson;}
  TH2*     GetHChi2LrCorr()    const {return fHChi2LrCorr;}
  TH2*     GetHChi2NDFCorr()   const {return fHChi2NDFCorr;}
  TH2*     GetHChi2NDFFake()   const {return fHChi2NDFFake;}
  TH2*     GetHChi2VCorr()     const {return fHChi2VtxCorr;}
  TH2*     GetHChi2VFake()     const {return fHChi2VtxFake;}
  TH2*     GetHNCand()         const {return fHNCand;}
  TH2*     GetHCandCorID()     const {return fHCandCorID;}
  TH2*     GetHChi2MS()        const {return fHChi2MS;}

  // Howard W. hit distribution and convolution integral
  //  Double_t Dist              ( Double_t Z, Double_t radius ) ;  

  //  Double_t UpcHitDensity     ( Double_t zPos )   ;
  //  Double_t IntegratedHitDensity  ( Double_t multiplicity, Double_t zPos )   ;
  //  Double_t OneEventHitDensity    ( Double_t multiplicity, Double_t zPos ) const   ;

  //
  // ------------------------------
  //
 protected:
  Int_t  fNLayers;               // total number of layers in the model
  Int_t  fNActiveLayers;         // number of active layers in the model
  Int_t  fNActiveLayersITS;      // number of active ITS layers in the model
  Int_t  fNActiveLayersMS;      // number of active MS layers in the model
  Int_t  fNActiveLayersTR;      // number of active Trigger layers in the model
  Int_t  fLastActiveLayerITS;    // id of last active ITS layer
  Int_t  fLastActiveLayer;       // id of last active layer
  Int_t  fLastActiveLayerTracked;    // id of last active layer really used for tracking of given pt
  //-------------------------
  TList  fLayers;                // List of layer pointers
  TObjArray fLayersITS;          // vertex tracker layers
  TObjArray fLayersMS;           // MS layers
  TObjArray fLayersTR;           // Trigger layers
  BeamPipe* fBeamPipe;           // special object - beam pipe
  Layer* fVtx;             // special layer: vertex
  TObjArray fMaterials;                                 // defined materials
  Int_t     fMagFieldID;                                // type of mag field
  //
  //-------------RS---------------------------
  Probe fProbe;
  Probe fMuTrackVertex;
  Probe fMuTrackLastITS;  
  Probe fMuTrackBCVertex;
  Probe fMuTrackBCLastITS;  
  
  Bool_t   fExternalInput; // MC particles are set externally
  Bool_t   fIncludeVertex;
  float    fApplyBransonPCorrection; // if >=0, apply BP correction with additional error on the vertex
  Bool_t   fUseBackground; // do we want to simulate background?
  // reconstruction settings
  Double_t fErrScale;   // parameter defining the initial cov.matrix error wrt sensor resolution 
  Double_t fMaxChi2Cl;   // max cluster-track chi2 
  Double_t fMaxNormChi2NDF;// max chi2/NDF to accept
  Double_t fMaxChi2Vtx; // cut on chi2 to vtx
  Int_t    fMinITSHits;  // min ITS hits in track to accept
  Int_t    fMinMSHits;  // min MS hits in track to accept
  Int_t    fMinTRHits;  // min Trigger hits in track to accept
  Double_t fMinP2Propagate; // min p to continue propagation
  //
  Double_t fMaxChi2ClSQ; // precalculated sqrt(chi2);
  //
  Int_t    fMaxSeedToPropagate; // take 1st fMaxSeedToPropagate seeds from each layer
  //
  TH1*     fDecayZProf; // optional decay Z profile
  Double_t fZDecay;     // impose decay here
  Int_t    fDecMode;    // decay mode
  float fChi2MuVtx; // chi2 of muon at vertex, if Branson correction is on
  //
  // field stepping optimization
  Int_t     fFldNReg;   // number of field regions
  const Double_t* fFldZMins;  // z of field regions beginning
  const Double_t* fFldZMaxs;  // z of field regions end
  Double_t* fFieldBndZ;  // Z boundaries of the field
  Int_t     fFieldBndZN; // number of boundaries
  Double_t  fDefStepAir;    // default step size in air
  Double_t  fDefStepMat;    // default step size in material
  Double_t  fGenPnt[3];     // particle generation point
  Double_t  fRefVtx[3];     // reference vertex (as all constraints, stored in tracking frame)
  Bool_t fImposeVertexPosition; // impose vertex position externally
  //
  TArrayI  fPattITS;                     // bit pattern of each group of active layers where hit is required
  Double_t fNCh;        // rapidity density
  TF1*     fdNdY;       // Y profile of bg
  TF1*     fdNdPt;      // Pt profile of bg
  Double_t fNChPi;        // pion multiplicity
  Double_t fNChK;        // kaon multiplicity 
  Double_t fNChP;        // proton multiplicity 
  TF1*     fdNdYPi;      // pion y
  TF1*     fdNdYK;       // kaon y
  TF1*     fdNdYP;       // proton p
  TF1*     fdNdPtPi; // pion pT
  TF1*     fdNdPtK;  // kaon pT
  TF1*     fdNdPtP;  // proton pT

  //
  // control histos
  TH1F*    fHChi2Branson;    // chi2 of muon at the vertex
  TH2F*    fHChi2LrCorr;    // chi2 of correct cluster
  TH2F*    fHChi2NDFCorr;   // total chi2 of correct tracks
  TH2F*    fHChi2NDFFake;   // total chi2 of fake tracks
  TH2F*    fHChi2VtxCorr;   // chi2 to vtx, correct tracks
  TH2F*    fHChi2VtxFake;   // chi2 to vtx, fake tracks
  TH2F* fHNCand;         // candidates per layer
  TH2F* fHChi2MS;
  TH2F* fHCandCorID;     // ID of correct in the candidates
  //
  static const Double_t fgkFldEps; // offset for field region boundaries crossing
  Double_t *fZSteps;
  //
 protected:

  ClassDef(Detector,1);
};

//====================================================

const double kVeryLarge = 1e16;
//==========================================================================

class NaCardsInput {
 public:
  NaCardsInput();
  NaCardsInput(const char* fname);
  virtual ~NaCardsInput();
  FILE* OpenFile(const char* fname);
  //
  virtual int FindEntry(const char* key, const char* mod="", 
			const char* form="", int rew=0, int errmsg=1);
  virtual int NextEntry(const char* key, const char* mod="",const char* form="");
  virtual int NextEntry();
  virtual int GetNArgs() {return fNArgs;}
  void  Rewind() {if (fCardsFile) rewind(fCardsFile);}
  void  StepBack() {if (fCardsFile) fseek(fCardsFile, fLastPos, SEEK_SET);}
  char* GetKey() {return fKey;}
  char* GetModifier() {return fModifier;}
  char* GetArg(const int iarg, const char* option="", int *err=0);
  char* GetArg(char* &dest, const int iarg, const char* option="", int *err=0);
  float GetArgF(const int iarg, int *err=0);
  int   GetArgD(const int iarg, int *err=0);
  unsigned int GetArgB(const int iarg, int *err=0);
  char  **GetArgs() {return fArgs;} 
  int   CompareKey(const char *key);
  int   CompareModifier(const char *mod);
  int   CompareArgList(const char *form);
  char* GetLastComment() const {return (char*)fLastComment.Data();}
  char* GetLastBuffer()  const {return (char*)fbuffer;}
  virtual void Print();
  //
 protected:
  virtual void ClearArgs(); 
 protected:
  static const char fgkComment='#';   // comment identifier
  static const char fgkDelimiter=':'; // delimiter between keyword and modifier
  static const char fgkContinuation='\\'; // delimiter between keyword and modifier
  static const int  fgkMaxLen = 2048;  // max. length of the entry
  //
  FILE* fCardsFile;        // pointer on the opened file
  long  fLastPos;            // position in the stream of last record read
  char  fbuffer[fgkMaxLen];// string buffer for current line
  char* fPoint;            // pointer on the beginning of data in the line
  int   fNArgs;            // number of the arguments in the entry
  char** fArgs;            // list of the arguments
  char* fKey;              // current Key
  char* fModifier;         // current Modifier
  TString fLastComment;     // last comments block read
  //
};


#endif
