#include <TStopwatch.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include "TrkDetector.h"
#include "LUXELog.h"

ClassImp(TrkDetector)

const Double_t TrkDetector::kMassP  = 0.938;
const Double_t TrkDetector::kMassK  = 0.4937;
const Double_t TrkDetector::kMassPi = 0.1396;
const Double_t TrkDetector::kMassMu = 0.1057;
const Double_t TrkDetector::kMassE  = 0.0005;

const Double_t TrkDetector::fgkFldEps = 1.e-4;

TrkDetector::TrkDetector(const char *name, const char *title) 
:  TNamed(name,title)
  ,fNTrkLayers(0)
  ,fNActiveTrkLayers(0)
  ,fNActiveTrkLayersITS(0)
  ,fNActiveTrkLayersMS(0)
  ,fNActiveTrkLayersTR(0)
  ,fLastActiveTrkLayerITS(0)
  ,fLastActiveTrkLayer(0)
  ,fLastActiveTrkLayerTracked(0)
  ,fTrkLayers()
  ,fVtx(0)
  ,fMaterials()
  ,fMagFieldID(kMagAlice)
  ,fTrkProbe()
  ,fExternalInput(kFALSE)
  ,fIncludeVertex(kTRUE)
  ,fUseBackground(kTRUE)
  ,fMaxChi2Cl(10.)
  ,fMaxNormChi2NDF(10)
  ,fMaxChi2Vtx(20)
  ,fMinITSHits(0)
  ,fMinMSHits(0)
  ,fMinTRHits(0)
  ,fMinP2Propagate(1.)
  ,fMaxChi2ClSQ(0)
  ,fMaxSeedToPropagate(0)
   //
,fDecayZProf(0)
,fZDecay(0)
,fDecMode(kNoDecay)
,fFldNReg(0)
,fFldZMins(0)
,fFldZMaxs(0)
,fDefStepAir(1.0)
,fDefStepMat(1.0)
,fImposeVertexPosition(kFALSE)
,fPattITS(0)
,fNCh(-1)
,fdNdY(0)
,fdNdPt(0)
,fHChi2LrCorr(0)
,fHChi2NDFCorr(0)
,fHChi2NDFFake(0)
,fHChi2VtxCorr(0)
,fHChi2VtxFake(0)
,fHNCand(0)
,fHChi2MS(0)
,fHCandCorID(0)
,fZSteps(0)
{
  //
  // default constructor
  //
  //  fTrkLayers = new TObjArray();
  fRefVtx[0] = fRefVtx[1] = fRefVtx[2] = 0;
}

TrkDetector::TrkDetector() {}
TrkDetector::~TrkDetector() { // 
  // virtual destructor
  //
  //  delete fTrkLayers;
  delete fdNdY;
  delete fdNdPt;
}


void TrkDetector::Print(const Option_t *opt) const
{
  // Prints the detector layout
  printf("TrkDetector %s: \"%s\"\n",GetName(),GetTitle());
  for (Int_t i = 0; i<fTrkLayers.GetEntries(); i++) GetTrkLayer(i)->Print(opt);
}

Double_t TrkDetector::ThetaMCS ( Double_t mass, Double_t radLength, Double_t momentum ) const
{
  //
  // returns the Multiple Couloumb scattering angle (compare PDG boolet, 2010, equ. 27.14)
  //
  Double_t beta  =  momentum / TMath::Sqrt(momentum*momentum+mass*mass)  ;
  Double_t theta =  0.0 ;    // Momentum and mass in GeV
  // if ( RadLength > 0 ) theta  =  0.0136 * TMath::Sqrt(RadLength) / ( beta * momentum );
  if ( radLength > 0 ) theta  =  0.0136 * TMath::Sqrt(radLength) / ( beta * momentum ) * (1+0.038*TMath::Log(radLength)) ;
  return (theta) ;
}

//__________________________________________________________________________
TrkLayer* TrkDetector::AddTrkLayer(const char* type, const char *name, Float_t zPos, Float_t radL, Float_t density, 
				      Float_t thickness, Float_t xRes, Float_t yRes, Float_t eff, Material* mat) 
{
  //
  // Add additional layer to the list of layers (ordered by z position)
  // 
  TrkLayer *newTrkLayer = GetTrkLayer(name);
  //
  if (!newTrkLayer) {
    TString types = type;
    types.ToLower();
    newTrkLayer = new TrkLayer(name);
    newTrkLayer->SetZ(zPos);
    newTrkLayer->SetThickness(thickness);
    newTrkLayer->SetX2X0( radL>0 ? thickness*density/radL : 0);
    newTrkLayer->SetXTimesRho(thickness*density);
    newTrkLayer->SetXRes(xRes);
    newTrkLayer->SetYRes(yRes);
    newTrkLayer->SetTrkLayerEff(eff);
    if      (types=="vt")   newTrkLayer->SetType(TrkLayer::kITS);
    else if (types=="ms")   newTrkLayer->SetType(TrkLayer::kMS);
    else if (types=="tr")   newTrkLayer->SetType(TrkLayer::kTRIG);
    else if (types=="vtx")  {newTrkLayer->SetType(TrkLayer::kVTX); }
    else if (types=="abs")  {newTrkLayer->SetType(TrkLayer::kABS); newTrkLayer->SetDead(kTRUE); }
    else if (types=="dummy")  {newTrkLayer->SetType(TrkLayer::kDUMMY); newTrkLayer->SetDead(kTRUE); }
    //
    if (!newTrkLayer->IsDead()) newTrkLayer->SetDead( xRes==kVeryLarge && yRes==kVeryLarge);
    //
    if (fTrkLayers.GetEntries()==0) fTrkLayers.Add(newTrkLayer);
    else {      
      for (Int_t i = 0; i<fTrkLayers.GetEntries(); i++) {
	TrkLayer *l = GetTrkLayer(i);
	if (zPos<l->GetZ()) { fTrkLayers.AddBefore(l,newTrkLayer); break; }
	if (zPos>l->GetZ() && (i+1)==fTrkLayers.GetEntries() ) { fTrkLayers.Add(newTrkLayer); } // even bigger then last one
      }      
    }
    //
  } else printf("TrkLayer with the name %s does already exist\n",name);
  newTrkLayer->SetMaterial(mat);
  //
  return newTrkLayer;
}



//________________________________________________________________________________
void TrkDetector::ClassifyTrkLayers()
{
  // assign active Id's, etc
  fLastActiveTrkLayer = -1;
  fLastActiveTrkLayerITS = -1;
  fNActiveTrkLayers = 0;
  fNActiveTrkLayersITS = 0;
  fNActiveTrkLayersMS = 0;
  fNActiveTrkLayersTR = 0;  
  //
  fNTrkLayers = fTrkLayers.GetEntries();
  for (int il=0;il<fNTrkLayers;il++) {
    TrkLayer* lr = GetTrkLayer(il);
    lr->SetID(il);
    if (!lr->IsDead()) {
      fLastActiveTrkLayer = il; 
      lr->SetActiveID(fNActiveTrkLayers++);
      if (lr->IsITS()) {
	fLastActiveTrkLayerITS = il;
	fNActiveTrkLayersITS++;
	fTrkLayersITS.AddLast(lr);
      }
      if (lr->IsMS())   {
	fNActiveTrkLayersMS++;
	fTrkLayersMS.AddLast(lr);
      }
      if (lr->IsTrig()) {
	fNActiveTrkLayersTR++;
	fTrkLayersTR.AddLast(lr);
      }
    }
  }
  //
  TVirtualMagField* fld = TGeoGlobalMagField::Instance()->GetField();
  if (fld->IsA() == MagField::Class()) {
    MagField* fldm = (MagField*) fld;
    fFldNReg = fldm->GetNReg();
    fFldZMins = fldm->GetZMin();
    fFldZMaxs = fldm->GetZMax();
    fDefStepAir = 100;
    fDefStepMat = 5;
  }
  else {
    fDefStepAir = 1;
    fDefStepMat = 1;
  }
  fZSteps = new Double_t[fFldNReg+2];
  //
  printf("DefStep: Air %f Mat: %f | N.MagField regions: %d\n",fDefStepAir,fDefStepMat,fFldNReg);
  //
  TrkProbe::SetNITSTrkLayers(fNActiveTrkLayersITS + ((fVtx && !fVtx->IsDead()) ? 1:0));
  printf("TrkProbe is initialized with %d slots\n",TrkProbe::GetNITSTrkLayers());
}

//_________________________________________________________________
TrkProbe* TrkDetector::CreateTrkProbe(double pt, double yrap, double phi, double mass, int charge, double x,double y, double z)
{
  // create track of given kinematics
  double xyz[3] = {x,y,z};
  double pxyz[3] = {pt*TMath::Cos(phi),pt*TMath::Sin(phi),TMath::Sqrt(pt*pt+mass*mass)*TMath::SinH(yrap)};
  TrkProbe* probe = new TrkProbe(xyz,pxyz,charge);
  probe->SetMass(mass);
  probe->SetTrID(-1);
  return probe;
}

//_________________________________________________________________
void TrkDetector::CreateTrkProbe(TrkProbe* adr, double pt, double yrap, double phi, double mass, int charge, double x,double y, double z)
{
  // create track of given kinematics
  double xyz[3] = {x,y,z};
  double pxyz[3] = {pt*TMath::Cos(phi),pt*TMath::Sin(phi),TMath::Sqrt(pt*pt+mass*mass)*TMath::SinH(yrap)};
  TrkProbe* probe = new(adr) TrkProbe(xyz,pxyz,charge);
  probe->SetTrID(-1);
  probe->SetMass(mass);
}

//_________________________________________________________________
TrkProbe* TrkDetector::PrepareTrkProbe(double pt, double yrap, double phi, double mass, int charge, double x,double y, double z)
{
  // Prepare trackable Kalman track at the farthest position
  //
  fGenPnt[0]=x;
  fGenPnt[1]=y;  
  fGenPnt[2]=z;  
  //
  if (fDecayZProf) { // decay is requested
    fZDecay  = fDecayZProf->GetRandom();
    fDecMode = kDoRealDecay;
    LUXEDebug(2,Form("Selected %.2f as decay Z",fZDecay));
  }
  // track parameters
  // Assume track started at (0,0,0) and shoots out on the X axis, and B field is on the Z axis
  fTrkProbe.Reset();
  TrkProbe* probe = CreateTrkProbe(pt,yrap,phi,mass,charge,x,y,z);
  fTrkProbe = *probe;     // store original track
  //
  // propagate to last layer
  fLastActiveTrkLayerTracked = 0;
  int resp=0;
  TrkLayer* lr=0,*lrP=0;
  for (Int_t j=0; j<=fLastActiveTrkLayer; j++) {
    lrP = lr;
    lr = GetTrkLayer(j);
    lr->Reset();
    if (!lrP) continue;
    //
    if (!(resp=PropagateToTrkLayer(probe,lrP,lr,1))) return 0;
    Cluster* cl = lr->GetCorCluster();
    double r = probe->GetR();
    //    printf("L%2d %f %f %f\n",j,r, lr->GetRMin(),lr->GetRMax());
    if (r<lr->GetRMax() && r>lr->GetRMin()) {
      if (resp>0) cl->Set(probe->GetXLoc(),probe->GetYLoc(), probe->GetZLoc(),probe->GetTrID());
      else cl->Kill();
    }
    else cl->Kill();
    //
    if (!lr->IsDead()) fLastActiveTrkLayerTracked = j;
  }
  probe->ResetCovariance();// reset cov.matrix
  //  printf("Last active layer trracked: %d (out of %d)\n",fLastActiveTrkLayerTracked,fLastActiveTrkLayer);
  //
  return probe;
}

//____________________________________________________________________________
Int_t TrkDetector::GetFieldReg(double z) 
{
  // return field region * 2 + 1
  int ir = 0;
  for (int i=0;i<fFldNReg;i++) {
    if (z<=fFldZMins[i]) return ir;
    ir++;
    if (z<=fFldZMaxs[i]) return ir;
    ir++;
  }
  return ir;
}

//____________________________________________________________________________
Bool_t TrkDetector::PropagateToZBxByBz(TrkProbe* trc,double z,double maxDZ,Double_t xOverX0,Double_t xTimesRho,Bool_t modeMC) 
{
  // propagate to given Z, checking for the field boundaries
  //
  double curZ = trc->GetZ();
  double dza = curZ-z;
  if (TMath::Abs(dza)<fgkFldEps) return kTRUE;
  // even id's correspond to Z between the field regions, odd ones - inside
  int ib0 = GetFieldReg(curZ); // field region id of start point
  int ib1 = GetFieldReg(z);    // field region id of last point
  int nzst = 0;
  //  LUXEDebug(2,Form("FldRegID: %d %d (%f : %f)",ib0,ib1, curZ,z));
  if (ib1>ib0) { // fwd propagation with field boundaries crossing
    for (int ib=ib0;ib<ib1;ib++) {
      if ( ib&0x1 ) { // we are in the odd (field ON) region, go till the end of field reg.
	//	printf("Here00 | %d %f\n",ib>>1,fFldZMaxs[ib>>1]);
	fZSteps[nzst++] = fFldZMaxs[ib>>1] + fgkFldEps;
      }
      else { // we are in even (field free) region, go till the beginning of next field reg.
	//	printf("Here01 | %d %f\n",ib>>1,fFldZMins[ib>>1]);
	fZSteps[nzst++] = fFldZMins[ib>>1] + fgkFldEps;
      }
    }
  }
  else if (ib1<ib0) { // bwd propagation
    for (int ib=ib0;ib>ib1;ib--) {
      if ( ib&0x1 ) { // we are in the odd (field ON) region, go till the beginning of field reg.
	//	printf("Here10 | %d %f\n",(ib-1)>>1,fFldZMins[(ib-1)>>1]);
	fZSteps[nzst++] = fFldZMins[(ib-1)>>1] - fgkFldEps;
      }
      else { // we are in even (field free) region, go till the beginning of next field reg.
	//	printf("Here11 | %d %f\n",(ib-1)>>1,fFldZMaxs[(ib-1)>>1]);
	fZSteps[nzst++] = fFldZMaxs[(ib-1)>>1] - fgkFldEps;
      }
    }
  }
  fZSteps[nzst++] = z; // same field region, do nothing
  //
  //  printf("ZSteps: "); for (int ist=0;ist<nzst;ist++) printf("%+.5f ",fZSteps[ist]); printf("\n");
  for (int ist=0;ist<nzst;ist++) {
    double frc = (trc->GetZ()-fZSteps[ist])/dza;
    if (!trc->PropagateToZBxByBz(fZSteps[ist], maxDZ, frc*xOverX0, frc*xTimesRho, modeMC)) return kFALSE;
  }
  return kTRUE;
  //
}

//____________________________________________________________________________
Int_t TrkDetector::PropagateToTrkLayer(TrkProbe* trc, TrkLayer* lrFrom, TrkLayer* lrTo, int dir, Bool_t modeMC)
{
  // bring the track to lrTo, moving in direction dir (1: forward, -1: backward)
  // if relevant, account for the materials on lr
  //
  LUXEDebug(2,Form("From %d to %d, dir: %d",lrFrom? lrFrom->GetUniqueID():-1, lrTo->GetUniqueID(), dir));
  if (trc->GetTrack()->GetAlpha()<0) return -1;
  if ( lrFrom && dir*(lrTo->GetZ()-lrFrom->GetZ())<0 ) LUXEFatal(Form("Dir:%d Zstart: %f Zend: %f\n",dir,lrFrom->GetZ(),lrTo->GetZ()));
  //
  if      (dir>0 && (lrTo->GetZ()-0.5*lrTo->GetThickness()) < trc->GetZ() ) return -1;
  else if (dir<0 && (lrTo->GetZ()+0.5*lrTo->GetThickness()) > trc->GetZ() ) return -1;
  double dstZ;
  //
  if (lrFrom) {
    //
    if (!lrFrom->IsDead()) { // active layers are thin, no need for step by step tracking. The track is always in the middle
      LUXEDebug(2,Form("Correcting for mat.in active layer: X/X0: %f X*rho:%f ", lrFrom->GetX2X0(), dir*lrFrom->GetXTimesRho()));
      // note: for thin layer we ignore difference between the real BB and ETP eloss params
      if (!trc->CorrectForMeanMaterial(lrFrom->GetX2X0(), -dir*lrFrom->GetXTimesRho(), modeMC)) return 0; 
    }
    else {
      //
      dstZ = lrFrom->GetZ()+0.5*dir*lrFrom->GetThickness(); // go till the end of starting layer applying corrections
      if (dir==1 && trc->GetZ()<=fZDecay && dstZ>fZDecay) { // need to perform or to apply decay
	double frac = (fZDecay-trc->GetZ())/lrFrom->GetThickness();
	// account for the difference between real BB and ETP param eloss
	double corrELoss = lrFrom->GetELoss2ETP(trc->GetP(), trc->GetMass() );
	if (!PropagateToZBxByBz(trc,fZDecay, fDefStepMat, frac*lrFrom->GetX2X0(), -frac*lrFrom->GetXTimesRho()*corrELoss, modeMC)) return 0;
	PerformDecay(trc);
	frac = 1.-frac;
	corrELoss = lrFrom->GetELoss2ETP(trc->GetP(), trc->GetMass() );
	if (!PropagateToZBxByBz(trc,dstZ, fDefStepMat, frac*lrFrom->GetX2X0(), -frac*lrFrom->GetXTimesRho()*corrELoss, modeMC)) return 0;
      }
      else {
	// account for the difference between real BB and ETP param eloss
	double corrELoss = lrFrom->GetELoss2ETP(trc->GetP(), trc->GetMass() );
	if (!PropagateToZBxByBz(trc,dstZ, fDefStepMat, lrFrom->GetX2X0(), -dir*lrFrom->GetXTimesRho()*corrELoss, modeMC)) return 0;
      }
    }
  }
  //
  dstZ = lrTo->GetZ();
  if (lrTo->IsDead()) dstZ += -dir*lrTo->GetThickness()/2; // for thick dead layers go till entrance
  //
  if (dir==1 && trc->GetZ()<=fZDecay && dstZ>fZDecay) { // need to perform or to apply decay
    if (!PropagateToZBxByBz(trc,fZDecay, fDefStepAir)) return 0;
    PerformDecay(trc);
  }
  if (!PropagateToZBxByBz(trc,dstZ, fDefStepAir)) return 0;
  //
  // if (LUXELog::GetGlobalDebugLevel()>=2) trc->GetTrack()->Print();

  return 1;
}

//________________________________________________________________________________
Bool_t TrkDetector::SolveSingleTrackViaKalman(double pt, double yrap, double phi, 
						 double mass, int charge, double x,double y, double z)
{
  // analytical estimate of tracking resolutions
  //  fTrkProbe.SetUseLogTermMS(kTRUE);
  //
  for (int i=3;i--;) fGenPnt[i] = 0;
  if (fMinITSHits>fNActiveTrkLayersITS) {
    fMinITSHits = fLastActiveTrkLayerITS; 
    printf("Redefined request of min N ITS hits to %d\n",fMinITSHits);
  }
  //
  TrkProbe* probe = PrepareTrkProbe(pt,yrap,phi,mass,charge,x,y,z);
  if (!probe) return kFALSE;
  //
  TrkLayer *lr = 0;
  //
  // Start the track fitting --------------------------------------------------------
  //
  // Back-propagate the covariance matrix along the track. 
  // Kalman loop over the layers
  //
  TrkProbe* currTr = 0;
  lr = GetTrkLayer(fLastActiveTrkLayerTracked);
  lr->SetAnTrkProbe(*probe);
  delete probe; // rethink...
  //
  int nupd = 0;
  for (Int_t j=fLastActiveTrkLayerTracked; j--; ) {  // TrkLayer loop
    //
    TrkLayer *lrP = lr;
    lr = GetTrkLayer(j);
    //
    lr->SetAnTrkProbe( *lrP->GetAnTrkProbe() );
    currTr = lr->GetAnTrkProbe();
    currTr->ResetHit(lrP->GetActiveID());
    //
    // if there was a measurement on prev layer, update the track
    if (!lrP->IsDead()) { // include "ideal" measurement
      //      printf("Before update on %d : ",j); currTr->Print("etp");
      Cluster* cl = lrP->GetCorCluster();
      if (!cl->IsKilled()) {
	if (!UpdateTrack(currTr,lrP,cl))  return kFALSE;
	nupd++;
      }
      //      printf("After update on %d (%+e %+e) : ",j, lrP->GetXRes(),lrP->GetYRes()); currTr->Print("etp");

    }

    if (!PropagateToTrkLayer(currTr,lrP,lr,-1)) return kFALSE;      // propagate to current layer
    //
  } // end loop over layers
  // is MS reco ok?
  // check trigger
  int nhMS=0,nhTR=0,nhITS=0;
  for (int ilr=fNTrkLayers;ilr--;) {
    TrkLayer *lrt = GetTrkLayer(ilr);
    if (lrt->IsTrig()) if (!lrt->GetCorCluster()->IsKilled()) nhTR++;
    if (lrt->IsMS())   if (!lrt->GetCorCluster()->IsKilled()) nhMS++;
    if (lrt->IsITS())  if (!lrt->GetCorCluster()->IsKilled()) nhITS++;
  }
  //  printf("ITS: %d MS: %d TR: %d\n",nhITS,nhMS,nhTR);
  if (nhTR<fMinTRHits) return kFALSE;
  if (nhMS<fMinMSHits) return kFALSE;
  //
  return kTRUE;
}

//____________________________________________________________________________
Bool_t TrkDetector::UpdateTrack(TrkProbe* trc, const TrkLayer* lr, const Cluster* cl) const
{
  // update track with measured cluster
  // propagate to cluster
  if (cl->IsKilled()) return kTRUE;
  double meas[2] = {cl->GetY(),cl->GetZ()}; // ideal cluster coordinate, tracking (AliExtTrParam frame)
  double rcl = TMath::Sqrt(cl->GetY()*cl->GetY()+cl->GetZ()*cl->GetZ());
  double sgY = lr->GetYRes(rcl), sgX = lr->GetXRes(rcl);
  double measErr2[3] = {sgY*sgY,0,sgX*sgX}; // !!! Lab 
  //
  double chi2 = trc->GetTrack()->GetPredictedChi2(meas,measErr2);
  //    printf("Update for lr:%s -> chi2=%f\n",lr->GetName(), chi2);
  //    printf("cluster was :"); cl->Print("lc");
  //    printf("track   was :"); trc->Print("etp");  
    if (chi2>fMaxChi2Cl) return kTRUE; // chi2 is too large
    
  if (!trc->Update(meas,measErr2)) {
    // LUXEDebug(2,Form("layer %s: Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}", lr->GetName(),meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]));
    // if (LUXELog::GetGlobalDebugLevel()>1) trc->Print("l");
    return kFALSE;
  }
  trc->AddHit(lr, chi2, cl->GetTrID());
  //
  return kTRUE;
}

//________________________________________________________________________________
Bool_t TrkDetector::SolveSingleTrack(double pt, double yrap, double phi, 
					double mass, int charge, double x,double y, double z, 
					TObjArray* sumArr,int nMC, int offset)
{
  // analityc and fullMC (nMC trials) evaluaion of tracks with given kinematics.
  // the results are filled in KMCTrackSummary objects provided via summArr array
  //
  //
  if (!SolveSingleTrackViaKalman(pt,yrap,phi,mass,charge,x,y,z)) return kFALSE;
  //
  /*
  int nsm = sumArr ? sumArr->GetEntriesFast() : 0;
  TrkLayer* vtx = GetTrkLayer(0);
  */
  //
  /*RS
  for (int i=0;i<nsm;i++) {
    KMCTrackSummary* tsm = (KMCTrackSummary*)sumArr->At(i);
    if (!tsm) continue;
    tsm->SetRefTrkProbe( GetTrkProbeTrack() ); // attach reference track (generated)
    tsm->SetAnTrkProbe( vtx->GetAnTrkProbe() ); // attach analitycal solution
  }
  */
  //
  if (offset<0) offset = fNTrkLayers;

  TStopwatch sw;
  sw.Start();
  //
  sw.Stop();
  //  printf("Total time: "); sw.Print();
  return kTRUE;
}

//____________________________________________________________________________
void TrkDetector::ResetMCTracks(Int_t maxLr)
{
  if (maxLr<0 || maxLr>=fNTrkLayers) maxLr = fNTrkLayers-1;
  for (int i=maxLr+1;i--;) GetTrkLayer(i)->ResetMCTracks();
}




//________________________________________________________________________________
// TODO: Noam code start
Bool_t TrkDetector::SolveSingleTrackViaKalmanMC_Noam(double pt, double yrap, double phi, 
						 double mass, int charge, double x,double y, double z, int offset)
{
  TrkProbe* probe = PrepareTrkProbe(pt,yrap,phi,mass,charge,x,y,z); /// delete later
  if(!probe) return kFALSE;

  fMuTrackVertex.SetUniqueID(999); // invalidate
  fMuTrackBCVertex.SetUniqueID(999); // invalidate
  fMuTrackBCLastITS.SetUniqueID(999); // invalidate
  fMuTrackLastITS.SetUniqueID(999); // invalidate

  TrkProbe *currTrP=0,*currTr=0;
  static TrkProbe trcConstr;
  int maxLr = fLastActiveTrkLayerITS;
  if (offset>0) maxLr += offset;
  if (maxLr>fLastActiveTrkLayer) maxLr = fLastActiveTrkLayer;
  if (fExternalInput) maxLr = fLastActiveTrkLayerTracked;
  if (maxLr<0) return kFALSE;

  if(fVtx && !fImposeVertexPosition)
  {
    double tmpLab[3] = {fTrkProbe.GetX(),fTrkProbe.GetY(),fTrkProbe.GetZ()};
    TrkProbe::Lab2Trk(tmpLab, fRefVtx); // assign in tracking frame
    fVtx->GetMCCluster()->Set(fRefVtx[0],fRefVtx[1],fRefVtx[2]);
  }

  TrkLayer* lr = GetTrkLayer(maxLr);
  currTr = lr->AddMCTrack(probe); // start with seed track at vertex
  // probe-sPrint("etp");

  // randomize the starting point
  // const float kErrScale = 500.; // this is the parameter defining the initial cov.matrix error wrt sensor resolution
  double r = currTr->GetR();
  currTr->ResetCovariance( fErrScale*TMath::Sqrt(lr->GetXRes(r)*lr->GetYRes(r)) ); // this is the coeff to play with

  int fst = 0;
  const int fstLim = -1;

  for(Int_t j=maxLr; j--; )  // TrkLayer loop
  {
    int ncnd=0;
	 int cndCorr=-1;
    TrkLayer *lrP = lr;
    lr = GetTrkLayer(j);
	 // printf("LR   |" );  lr->Print("cl");
	 // printf("LRP |" );   lrP->Print("cl");
	 
    int ntPrev = lrP->GetNMCTracks();
	 // LUXEInfo(Form("Got %d ntPrev",ntPrev));
    if(lrP->IsDead()) // for passive layer just propagate the copy of all tracks of prev layer >>>
	 {
      for(int itrP=ntPrev;itrP--;) // loop over all tracks from previous layer
		{
        currTrP = lrP->GetMCTrack(itrP); 
        if(currTrP->IsKilled())
        {
          continue;
        }
        currTr = lr->AddMCTrack( currTrP );
        if(fst<fstLim)
        {
          fst++;
          // currTr->Print("etp");
        }
	     if(!PropagateToTrkLayer(currTr,lrP,lr,-1)) // propagate to current layer:
        {
           currTr->Kill();
			  lr->GetMCTracks()->RemoveLast();
			  continue;
		  }
      }
      continue;
    } // end of treatment of dead layer <<<
	 
    // LUXEInfo(Form("From Lr: %d | %d seeds, %d bg clusters",j+1,ntPrev,lrP->GetNBgClusters()));
    for(int itrP=0;itrP<ntPrev;itrP++) // loop over all tracks from previous layer
	 {	 
      currTrP = lrP->GetMCTrack(itrP);
		if(currTrP->IsKilled())
      {
        continue;
      }
      currTr = lr->AddMCTrack( currTrP );
      if(fst<fstLim)
		{
         fst++;
         // currTr->Print("etp");
      }
      // LUXEInfo(Form("LastChecked before:%d",currTr->GetInnerTrkLayerChecked()));
		// printf("tr%d | ", itrP); currTr->Print("etp");
      CheckTrackProlongations(currTr, lrP,lr);

      // LUXEInfo(Form("LastChecked after:%d",currTr->GetInnerTrkLayerChecked()));
      ncnd++;
      if(currTr->GetNFakeITSHits()==0 && cndCorr<ncnd) cndCorr=ncnd;
      if(NeedToKill(currTr)) {currTr->Kill(); continue;}
    }
    if (fHNCand)     fHNCand->Fill(lrP->GetActiveID(), ncnd);
    if (fHCandCorID) fHCandCorID->Fill(lrP->GetActiveID(), cndCorr);
    lr->GetMCTracks()->Sort();
    int ntTot = lr->GetNMCTracks(); // propagate max amount of allowed tracks to current layer
    if(ntTot>fMaxSeedToPropagate && fMaxSeedToPropagate>0)
	 {
      for (int itr=ntTot;itr>=fMaxSeedToPropagate;itr--)  lr->GetMCTracks()->RemoveAt(itr);
      ntTot = fMaxSeedToPropagate;
    }
    for(int itr=ntTot;itr--;)
	 {
      currTr = lr->GetMCTrack(itr);
      if(currTr->IsKilled())
      {
        continue;
      }
      if(!PropagateToTrkLayer(currTr,lrP,lr,-1)) // propagate to current layer
      {
        currTr->Kill();
        continue;
      }
    }
    // LUXEInfo(Form("Got %d tracks on layer %s",ntTot,lr->GetName()));
  } // end loop over layers

  // do we use vertex constraint?
  if(fVtx && !fVtx->IsDead() && fIncludeVertex)
  {
    int ntr = fVtx->GetNMCTracks();
    for(int itr=0;itr<ntr;itr++)
	 {
      currTr = fVtx->GetMCTrack(itr);
      if(currTr->IsKilled())
      {
        continue;
      }
      double meas[2] = {0.,0.};
      if(fImposeVertexPosition)
		{
        meas[0] = fRefVtx[1]; // bending
        meas[1] = fRefVtx[2]; // non-bending
      }
      else
		{
         Cluster* clv = fVtx->GetMCCluster();
         meas[0] = clv->GetY(); 
         meas[1] = clv->GetZ();
      }
      double measErr2[3] = {fVtx->GetYRes()*fVtx->GetYRes(),0,fVtx->GetXRes()*fVtx->GetXRes()}; //  Lab
      if(!currTr->Update(meas,measErr2))
      {
         continue;
		}
      currTr->SetInnerLrChecked(fVtx->GetActiveID());
    }
  }
  int ntTot = lr->GetNMCTracks();
  // LUXEInfo(Form("Got %d tracks",ntTot));
  
  ntTot = TMath::Min(1,ntTot);
  for (int itr=ntTot;itr--;)
  {
    currTr = lr->GetMCTrack(itr);
    if (currTr->IsKilled()) continue;
    if (fHChi2NDFCorr&&fHChi2NDFFake)
	 {
      if (IsCorrect(currTr)) fHChi2NDFCorr->Fill(currTr->GetNITSHits(),currTr->GetNormChi2(kTRUE));
      else                   fHChi2NDFFake->Fill(currTr->GetNITSHits(),currTr->GetNormChi2(kTRUE));
    }
  }
  delete probe; // can delete it here
  return kTRUE;
} 
// TODO: Noam code end
//________________________________________________________________________________



//________________________________________________________________________________
// TODO: Noam code start
Bool_t TrkDetector::SolveSingleTrackViaKalmanMC_Noam_multiseed(std::vector<TLorentzVector>& pseeds, double mass, int charge, int offset, bool doPrint)
{
  fMuTrackVertex.SetUniqueID(999); // invalidate
  fMuTrackBCVertex.SetUniqueID(999); // invalidate
  fMuTrackBCLastITS.SetUniqueID(999); // invalidate
  fMuTrackLastITS.SetUniqueID(999); // invalidate
	
  TrkProbe *currTrP=0,*currTr=0;
  static TrkProbe trcConstr;
  int maxLr = fLastActiveTrkLayerITS;
  if (offset>0) maxLr += offset;
  if (maxLr>fLastActiveTrkLayer) maxLr = fLastActiveTrkLayer;
  if (fExternalInput) maxLr = fLastActiveTrkLayerTracked;
  if(maxLr<0) return kFALSE;
  
  if(fVtx && !fImposeVertexPosition)
  {
    double tmpLab[3] = {fTrkProbe.GetX(),fTrkProbe.GetY(),fTrkProbe.GetZ()};
    TrkProbe::Lab2Trk(tmpLab, fRefVtx); // assign in tracking frame
    fVtx->GetMCCluster()->Set(fRefVtx[0],fRefVtx[1],fRefVtx[2]);
  }

  std::vector<TrkProbe*> probes;
  for(unsigned int s=0 ; s<pseeds.size() ; ++s) 
  { 
     double x=0,y=0,z=0; 
     TrkProbe* probe = PrepareTrkProbe(pseeds[s].Pt(),pseeds[s].Rapidity(),pseeds[s].Phi(),mass,charge,x,y,z);
     if (!probe) continue;
     probes.push_back( probe );
  }
  
  /// add the tracks after caching all the probes since PrepareTrkProbe resets all layers
  for(int l=0 ; l<GetTrkLayers()->GetEntries() ; l++) GetTrkLayer(l)->Reset();
  TrkLayer* lr = GetTrkLayer(maxLr);
  for(unsigned int p=0 ; p<probes.size() ; ++p)
  {
    currTr = lr->AddMCTrack( probes[p] );
    // randomize the starting point
	 // const float kErrScale = 200.; // this is the parameter defining the initial cov.matrix error wrt sensor resolution
    double r = currTr->GetR();
    currTr->ResetCovariance( fErrScale*TMath::Sqrt(lr->GetXRes(r)*lr->GetYRes(r)) ); // this is the coeff to play with
    if(doPrint) {LUXEInfo(Form("Added track %d with %f fErrScale", p,fErrScale)); currTr->Print("etp clid");}
  }
  
  int fst = 0;
  const int fstLim = -1;
  
  for(Int_t j=maxLr; j--; )  // TrkLayer loop
  {
    int ncnd=0;
  	int cndCorr=-1;
    TrkLayer *lrP = lr;
    lr = GetTrkLayer(j);
  	if(doPrint) {printf("LR  |" ); lr->Print("cl clid");}
  	if(doPrint) {printf("LRP |" ); lrP->Print("cl clid");}
  	 
    int ntPrev = lrP->GetNMCTracks();
  	if(doPrint) LUXEInfo(Form("Got %d ntPrev",ntPrev));
    if(lrP->IsDead()) // for passive layer just propagate the copy of all tracks of prev layer >>>
  	{
		if(doPrint) LUXEInfo(Form("In lrP->IsDead()"));
      for(int itrP=ntPrev;itrP--;) // loop over all tracks from previous layer
  		{
        currTrP = lrP->GetMCTrack(itrP);
		  if(NeedToKill(currTr))
		  {
			  if(doPrint) LUXEInfo(Form("In NeedToKill(currTr)) for lrP->IsDead()==true"));
			  currTr->Kill();
			  continue;
		  }
        if(currTrP->IsKilled())
        {
			 if(doPrint) LUXEInfo(Form("In currTrP->IsKilled() for lrP->IsDead()==true"));
          continue;
        }
        currTr = lr->AddMCTrack( currTrP );
        if(fst<fstLim)
        {
          fst++;
          if(doPrint) {LUXEInfo(Form("In fst<fstLim for lrP->IsDead()==true")); currTr->Print("etp clid");}
        }
  	     if(!PropagateToTrkLayer(currTr,lrP,lr,-1)) // propagate to current layer:
        {
			  if(doPrint) LUXEInfo(Form("In !PropagateToTrkLayer(currTr,lrP,lr,-1) for lrP->IsDead()==true"));
           currTr->Kill();
  			  lr->GetMCTracks()->RemoveLast();
  			  continue;
  		  }
      }
      continue;
    } // end of treatment of dead layer <<<
    int ntPrevNew = lrP->GetNMCTracks();
  	 if(doPrint) LUXEInfo(Form("Got %d ntPrevNew",ntPrevNew));
	 
  	 
    if(doPrint) {LUXEInfo(Form("From Lr: %d | %d seeds, %d bg clusters",j+1,ntPrev,lrP->GetNBgClusters()));}
    for(int itrP=0;itrP<ntPrev;itrP++) // loop over all tracks from previous layer
  	 {	 
      currTrP = lrP->GetMCTrack(itrP);
		if(doPrint) {LUXEInfo(Form("Starting %d of %d (%d) at lr %d",itrP, ntPrev, currTrP->IsKilled(),j));}
  		if(currTrP->IsKilled())
      {
		  if(doPrint) LUXEInfo(Form("Track already killed for layer %s, in loop over ntPrev, itrP=%d",lr->GetName(),itrP));
        continue;
      }
      currTr = lr->AddMCTrack( currTrP );
      if(fst<fstLim)
  		{
         fst++;
         if(doPrint) {currTr->Print("etp clid");}
      }
      if(doPrint) {LUXEInfo(Form("LastChecked before:%d",currTr->GetInnerTrkLayerChecked()));}
  		if(doPrint) {printf("tr%d | ", itrP); currTr->Print("etp clid");}
      CheckTrackProlongations(currTr, lrP,lr,doPrint);
      if(doPrint) {LUXEInfo(Form("LastChecked after:%d",currTr->GetInnerTrkLayerChecked()));}
      ncnd++;
      if(currTr->GetNFakeITSHits()==0 && cndCorr<ncnd) cndCorr=ncnd;
      if(NeedToKill(currTr))
		{
			if(doPrint) LUXEInfo(Form("Killing track for layer %s, in loop over ntPrev (NeedToKill), itrP=%d",lr->GetName(),itrP));
			currTr->Kill();
			continue;
		}
		else { if(doPrint) LUXEInfo(Form("Track not killed after CheckTrackProlongations")); }
    }
    if (fHNCand)     fHNCand->Fill(lrP->GetActiveID(), ncnd);
    if (fHCandCorID) fHCandCorID->Fill(lrP->GetActiveID(), cndCorr);
    lr->GetMCTracks()->Sort();
    int ntTot = lr->GetNMCTracks(); // propagate max amount of allowed tracks to current layer
	 int nTotNotKilled = 0;
	 if(doPrint) {LUXEInfo(Form("ntTot: %d, fMaxSeedToPropagate: %d",ntTot,fMaxSeedToPropagate));}
    if(ntTot>fMaxSeedToPropagate && fMaxSeedToPropagate>0)
  	 {
      for (int itr=ntTot;itr>=fMaxSeedToPropagate;itr--) lr->GetMCTracks()->RemoveAt(itr);
      ntTot = fMaxSeedToPropagate;
    }
    for(int itr=ntTot;itr--;)
  	 {
      currTr = lr->GetMCTrack(itr);
      if(currTr->IsKilled())
      {
		  if(doPrint) {LUXEInfo(Form("Track is already killed layer %s (!currTr->IsKilled())",lr->GetName()));}
        continue;
      }
      if(!PropagateToTrkLayer(currTr,lrP,lr,-1)) // propagate to current layer
      {
		  if(doPrint) {LUXEInfo(Form("Killing track on layer %s (!PropagateToTrkLayer)",lr->GetName()));}
        currTr->Kill();
        continue;
      }
		nTotNotKilled++;
    }
    if(doPrint) {LUXEInfo(Form("Got %d initial tracks on layer %s with nTotNotKilled=%d",ntTot,lr->GetName(),nTotNotKilled));}
  } // end loop over layers
  
  int nTotNotKilledVtx = 0;
  int nTotKilledVtx = 0;
  // do we use vertex constraint?
  if(fVtx && !fVtx->IsDead() && fIncludeVertex)
  {
    int ntr = fVtx->GetNMCTracks();
    for(int itr=0;itr<ntr;itr++)
  	 {
      currTr = fVtx->GetMCTrack(itr);
      if(currTr->IsKilled())
      {
		  if(doPrint) LUXEInfo(Form("Already killed (VTX)"));
        continue;
      }
		else {if(doPrint) LUXEInfo(Form("Track not killed after fVtx->GetMCTrack(...)"));}
      double meas[2] = {0.,0.};
      if(fImposeVertexPosition)
  		{
        meas[0] = fRefVtx[1]; // bending
        meas[1] = fRefVtx[2]; // non-bending
      }
      else
  		{
         Cluster* clv = fVtx->GetMCCluster();
         meas[0] = clv->GetY(); 
         meas[1] = clv->GetZ();
      }
      double measErr2[3] = {fVtx->GetYRes()*fVtx->GetYRes(),0,fVtx->GetXRes()*fVtx->GetXRes()}; //  Lab
      if(!currTr->Update(meas,measErr2))
      {
			if(doPrint) LUXEInfo(Form("Failed in !currTr->Update(meas,measErr2) for itr=%d",itr));
			nTotKilledVtx++;
         continue;
  		}
		nTotNotKilledVtx++;
		if(doPrint) LUXEInfo(Form("after !currTr->Update(meas,measErr2)"));
      currTr->SetInnerLrChecked(fVtx->GetActiveID());
    }
  }
  int ntTot = lr->GetNMCTracks();
  if(doPrint) {LUXEInfo(Form("Got %d tracks after vertex constraint with nTotNotKilledVtx=%d and nTotKilledVtx=%d",ntTot,nTotNotKilledVtx,nTotKilledVtx));}
  
  int nonkiled = 0;
  ntTot = TMath::Min(1,ntTot);
  for (int itr=ntTot;itr--;)
  {
    currTr = lr->GetMCTrack(itr);
    if (currTr->IsKilled()) continue;
	 nonkiled++;
    if (fHChi2NDFCorr&&fHChi2NDFFake)
  	 {
      if (IsCorrect(currTr)) fHChi2NDFCorr->Fill(currTr->GetNITSHits(),currTr->GetNormChi2(kTRUE));
      else                   fHChi2NDFFake->Fill(currTr->GetNITSHits(),currTr->GetNormChi2(kTRUE));
    }
  }
  if(doPrint)
  {
	  LUXEInfo(Form("Got %d final tracks non killed --> ",nonkiled));
	  if(currTr and !currTr->IsKilled()) currTr->Print("etp clid");
  }
  
  // delete the probes and clear their vector here
  for(unsigned int s=0 ; s<pseeds.size() ; ++s) {if(probes[s]) delete probes[s];}
  probes.clear();
  
  return kTRUE;
}
// TODO: Noam code end
///________________________________________________________________________________




void TrkDetector::CheckClusters(int i1, int i2, int i3, int i4)
{
	TrkLayer* lr = 0;
	Cluster *cl = 0;
	if(i1>=0) {lr = GetTrkLayer(1); cl = lr->GetBgCluster(i1); if(cl->IsKilled()) LUXEInfo(Form("Cluster with id=%d and index=%d for layer %s is killed",cl->GetTrID(),i1,lr->GetName()));}
	if(i2>=0) {lr = GetTrkLayer(3); cl = lr->GetBgCluster(i2); if(cl->IsKilled()) LUXEInfo(Form("Cluster with id=%d and index=%d for layer %s is killed",cl->GetTrID(),i2,lr->GetName()));}
	if(i3>=0) {lr = GetTrkLayer(5); cl = lr->GetBgCluster(i3); if(cl->IsKilled()) LUXEInfo(Form("Cluster with id=%d and index=%d for layer %s is killed",cl->GetTrID(),i3,lr->GetName()));}
	if(i4>=0) {lr = GetTrkLayer(7); cl = lr->GetBgCluster(i4); if(cl->IsKilled()) LUXEInfo(Form("Cluster with id=%d and index=%d for layer %s is killed",cl->GetTrID(),i4,lr->GetName()));}
}

//____________________________________________________________________________
void TrkDetector::CheckTrackProlongations(TrkProbe *probe, TrkLayer* lrP, TrkLayer* lr, bool doPrint)
{
  // explore prolongation of probe from lrP to lr with all possible clusters of lrP
  // the probe is already brought to clusters frame
  // for the last ITS plane apply Branson correction
  //if (lrP->GetUniqueID()==fLastActiveTrkLayerITS) probe->Print("etp");
  /*
  if (lrP->GetUniqueID()==fLastActiveTrkLayerITS && fVtx) {
    printf("Before Branson: "); probe->Print("etp");
    double zP = probe->GetZ();
    if (!PropagateToZBxByBz(probe,fVtx->GetZ(),fDefStepAir)) return;
    double measVErr2[3] = {fVtx->GetYRes()*fVtx->GetYRes(),0,fVtx->GetXRes()*fVtx->GetXRes()}; // we work in tracking frame here!
    double measV[2] = {0,0};
    probe->Update(measV,measVErr2);
    if (!PropagateToZBxByBz(probe,zP,fDefStepAir)) return;
    printf("After Branson: "); probe->Print("etp");
  }
  */
  static TrkProbe propVtx;
  //
  int nCl = lrP->GetNBgClusters();
  double rad = probe->GetR();
  double sgy = lrP->GetYRes(rad), sgx = lrP->GetXRes(rad); 
  double measErr2[3] = { sgx*sgx, 0, sgy*sgy}; // we work in tracking frame here!
  double meas[2] = {0,0};
  double tolerY = probe->GetSigmaX2() + measErr2[0];
  double tolerZ = probe->GetSigmaY2() + measErr2[2];
  tolerY = TMath::Sqrt(fMaxChi2Cl*tolerY);
  tolerZ = TMath::Sqrt(fMaxChi2Cl*tolerZ);
  double yMin = probe->GetTrack()->GetY() - tolerY;
  double yMax = yMin + tolerY+tolerY;    
  double zMin = probe->GetTrack()->GetZ() - tolerZ;
  double zMax = zMin + tolerZ + tolerZ;
  if(doPrint) probe->GetTrack()->Print("clid");
  probe->SetInnerLrChecked(lrP->GetActiveID());
  if(doPrint) {LUXEInfo(Form("From Lr(%d) %s to Lr(%d) %s | LastChecked %d", lrP->GetActiveID(),lrP->GetName(),lr->GetActiveID(),lr->GetName(),probe->GetInnerTrkLayerChecked()));}
  for (int icl=-1;icl<nCl;icl++) {
    //
    if(gRandom->Rndm() > lrP->GetTrkLayerEff()) continue; // generate layer eff
    //
    Cluster *cl = icl<0 ? lrP->GetMCCluster() : lrP->GetBgCluster(icl);  // -1 is for true MC cluster
    if (cl->IsKilled()) {
		 if(doPrint) {LUXEInfo(Form("Skip cluster %d ",icl)); cl->Print("clid");}
      continue;
    }
    double y = cl->GetY(); // ! tracking frame coordinates
    double z = cl->GetZ(); //                             
    //
    // if(doPrint) {LUXEInfo(Form("Check against cl#%d(%d) out of %d at layer %s | y: Tr:%+8.4f Cl:%+8.4f (%+8.4f:%+8.4f) z: Tr:%+8.4f Cl: %+8.4f (%+8.4f:%+8.4f)", icl,cl->GetTrID(),nCl,lrP->GetName(), probe->GetTrack()->GetY(),y,yMin,yMax,probe->GetTrack()->GetZ(),z,zMin,zMax));}
    //
    // if (z>zMax) {if (icl==-1) continue; else break;} // all other z will be even smaller, no chance to match
    if (z>zMax) continue; // all other z will be even smaller, no chance to match
    if (z<zMin) continue;
    if (y<yMin || y>yMax) continue;
    //
    meas[0] = y; meas[1] = z;
    double chi2 = probe->GetPredictedChi2(meas,measErr2);
    //
	 if(doPrint) {LUXEInfo(Form("Seed-to-cluster chi2 = Chi2=%.2f for cl:",chi2));}
	 if(doPrint) {cl->Print("lc");}
    if (icl<0 && fHChi2LrCorr) fHChi2LrCorr->Fill(lrP->GetActiveID(), chi2);
    if (chi2>fMaxChi2Cl) continue;
	 if(doPrint) {printf("Lr%d | cl%d, chi:%.3f X:%+.4f Y:%+.4f | x:%+.4f y:%+.4f |Sg: %.4f %.4f\n", lrP->GetActiveID(),icl,chi2, (zMin+zMax)/2,(yMin+yMax)/2, z,y, tolerZ/TMath::Sqrt(fMaxChi2Cl),tolerY/TMath::Sqrt(fMaxChi2Cl));}
    // update track copy
    TrkProbe* newTr = lr->AddMCTrack( probe );
    if (!newTr->Update(meas,measErr2)) {
		if(doPrint) {LUXEInfo(Form("TrkLayer %s: Failed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}", lrP->GetName(),meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]));}
		if(doPrint) {newTr->Print("l clid");}
      newTr->Kill();
      lr->GetMCTracks()->RemoveLast();
      continue;
    }
	 else
	 {
	 	if(doPrint) { LUXEInfo(Form("TrkLayer %s: Succeed to update the track by measurement {%.3f,%3f} err {%.3e %.3e %.3e}", lrP->GetName(),meas[0],meas[1], measErr2[0],measErr2[1],measErr2[2]));}
		if(doPrint) {newTr->Print("l clid");}
	 }
    if (fMinP2Propagate>0) {
      double p = newTr->GetTrack()->GetP();
      if (p<fMinP2Propagate) {
	      newTr->Kill();
	      lr->GetMCTracks()->RemoveLast();
	      continue;
      }
    }
	 if(doPrint) { LUXEInfo(Form("Adding hit (cluster id %d) on layer %s",cl->GetTrID(),lrP->GetName()));}
    newTr->AddHit(lrP, chi2, cl->GetTrID());

    //////////////////// check chi2 to vertex
    if (fVtx && !fVtx->IsDead() && fMaxChi2Vtx>0) {
      double measVErr2[3] = {fVtx->GetXRes()*fVtx->GetXRes(),0,fVtx->GetYRes()*fVtx->GetYRes()}; // in tracking frame
      propVtx = *newTr;
      if (!PropagateToZBxByBz(&propVtx,fRefVtx[2],fDefStepAir)) { newTr->Kill(); lr->GetMCTracks()->RemoveLast();}
      double chi2V = propVtx.GetTrack()->GetPredictedChi2(fRefVtx,measVErr2);
      if (fHChi2VtxCorr && fHChi2VtxFake) {
	     if (IsCorrect(newTr)) fHChi2VtxCorr->Fill(newTr->GetNITSHits(),chi2V);
	     else                  fHChi2VtxFake->Fill(newTr->GetNITSHits(),chi2V);
      }
      if(doPrint) {LUXEInfo(Form("Chi2 to vertex: %f | y: Tr:%+8.4f Cl:%+8.4f  z: Tr:%+8.4f Cl: %+8.4f",chi2V, propVtx.GetTrack()->GetY(),fRefVtx[0], propVtx.GetTrack()->GetZ(),fRefVtx[1]));}
      if(doPrint) {LUXEInfo(Form("Chi2 to vertex: %f | y: Tr:%+8.4f Cl:%+8.4f  z: Tr:%+8.4f Cl: %+8.4f",chi2V, propVtx.GetTrack()->GetY(),fRefVtx[0], propVtx.GetTrack()->GetZ(),fRefVtx[1]));}

      if (chi2V>fMaxChi2Vtx) {
			if(doPrint) {LUXEInfo(Form("Kill due to chi2V=%8.4f>fMaxChi2Vtx=%8.4f",chi2V,fMaxChi2Vtx));}
	      newTr->Kill();
	      lr->GetMCTracks()->RemoveLast();
	      continue;
      }

      /*
      double pz = 1./TMath::Abs(propVtx.GetTrack()->Get1P());
      if (pz<1) {
	printf("LowMom: %f | chi2V=%f (%f | %f %f %f)  ",pz,chi2V, fMaxChi2Vtx, fRefVtx[0],fRefVtx[1],fRefVtx[2]);  propVtx.Print("etp");       
      }
      */
    }

    ////////////////////////////////////////

    //    if (!PropagateToTrkLayer(newTr,lrP,lr,-1)) {newTr->Kill(); continue;} // propagate to next layer
    // if (LUXELog::GetGlobalDebugLevel()>1) {
    if(doPrint) {
      LUXEInfo(Form("Cloned updated track is:"));
      newTr->Print("clid");
    }
  }
  //
}

//_________________________________________________________
//Double_t TrkDetector::HitDensity(double xLab,double ylab,double zlab ) const
Double_t TrkDetector::HitDensity(double ,double ,double  ) const
{
  // RS to do
  return 1;
}

//_________________________________________________________
Bool_t TrkDetector::NeedToKill(TrkProbe* probe) const
{
  // check if the seed should be killed
  const Bool_t kModeKillMiss = kFALSE;
  Bool_t kill = kFALSE;
  while (1) {
    int il = probe->GetInnerTrkLayerChecked();
    int nITS = probe->GetNITSHits();
    int nITSMax = nITS + il; // maximum it can have
    if (nITSMax<fMinITSHits) {
      kill = kTRUE;
		// LUXEInfo(Form("Kill due to no chance to collect enough ITS hits: nITS=%d, nITSMax=%d, fMinITSHits=%d",nITS,nITSMax,fMinITSHits));
      break;
    } // has no chance to collect enough ITS hits
    //
    int ngr = fPattITS.GetSize();
    if (ngr>0) { // check pattern
      UInt_t patt = probe->GetHitsPatt();
      // complete the layers not checked yet
      for (int i=il;i--;) patt |= (0x1<<i);
      for (int ig=ngr;ig--;) 
	   if (!(((UInt_t)fPattITS[ig]) & patt)) {
	     kill = kTRUE; 
		  // LUXEInfo(Form("Kill due to no hit pattern"));
	     break;
	   }
      //
    }
    //
    if (nITS>2) {  // check if smallest possible norm chi2/ndf is acceptable
      double chi2min = probe->GetChi2ITS();
      if (kModeKillMiss) {
	      int nMiss = fNActiveTrkLayersITS - probe->GetInnerTrkLayerChecked() - nITS; // layers already missed
	      chi2min = nMiss*probe->GetMissingHitPenalty();
      }
      chi2min /= ((nITSMax<<1)-TrkProbe::kNDOF);
      if (chi2min>fMaxNormChi2NDF) {
	      kill = kTRUE; 
			// LUXEInfo(Form("Kill due to chi2"));
	      break;
      }
    }
    //
    /*
    // loose vertex constraint
    double dst;
    if (nITS>=2) {
      probe->GetZAt(0,fBFieldG,dst);
      //printf("Zd (F%d): %f\n",probe->GetNFakeITSHits(),dst);
      if (TMath::Abs(dst)>10.) {
	kill = kTRUE; 
	break;
      }
    }
    if (nITS>=3) {
      probe->GetYAt(0,fBFieldG,dst);
      //printf("Dd (F%d): %f\n",probe->GetNFakeITSHits(),dst);
      if (TMath::Abs(dst)>10.) {
	kill = kTRUE; 
	break;
      }
    }
    */
    //
    break;
  }
  if (kill && LUXELog::GetGlobalDebugLevel()>1 && probe->GetNFakeITSHits()==0) {
    printf("Killing good seed, last upd layer was %d\n",probe->GetInnerTrkLayerChecked());
    probe->Print("l");
  }
  return kill;
  //
}

//_________________________________________________________
void TrkDetector::PerformDecay(TrkProbe* trc)
{
  // Decay track
  if (fDecMode==kNoDecay) return;
  //  printf("DecMode: %d\n",fDecMode);
  static TGenPhaseSpace decay;
  static TLorentzVector pDecParent,pDecMu,pParCM;
  //
  const double kTol = 5e-3; // 5 MeV tolerance
  double mass = trc->GetMass();  
  double pxyz[3],xyz[3]={0};
  trc->GetPXYZ(pxyz);
  LUXEDebug(2,Form(" >>Mode:%d at Z=%f, PXYZ:%+6.3f %+6.3f %+6.3f, Mass:%.3f",fDecMode,fZDecay,pxyz[0],pxyz[1],pxyz[2],mass));
  static double ctau = 1e10;
  if (fDecMode==kDoRealDecay) {
    //
    if (TMath::Abs(mass-kMassMu)<kTol) {LUXEDebug(2,Form("Decay requested but provided mass %.4f hints to muon",mass)); exit(1);}
    if (TMath::Abs(mass-kMassPi)<kTol) {
      mass = kMassPi;
      Double_t masses[2] = {kMassMu, kMassE};
      pParCM.SetXYZM(0,0,0,mass);
      decay.SetDecay(pParCM, 2, masses);
      ctau = 780.45;
    }
    else if (TMath::Abs(mass-kMassK)<kTol) {
      mass = kMassK;
      Double_t masses[3] = {kMassMu, 0};
      pParCM.SetXYZM(0,0,0,mass);
      decay.SetDecay(pParCM, 2, masses);
      ctau = 371.2;
    }
    else {LUXEDebug(2,Form("Decay requested but provided mass %.4f is not recognized as pi or K",mass)); exit(1);}
    //
    decay.Generate();
    pDecMu = *decay.GetDecay(0); // muon kinematics in parent frame
    fDecMode = kApplyDecay;
  }
  //
  pDecParent.SetXYZM(pxyz[0],pxyz[1],pxyz[2],pParCM.M());
  pDecMu = *decay.GetDecay(0);
  pDecMu.Boost(pDecParent.BoostVector());
  //
  pxyz[0] = pDecMu.Px();
  pxyz[1] = pDecMu.Py();
  pxyz[2] = pDecMu.Pz();
  //
  static TrkProbe tmpPr;
  tmpPr.Init(xyz,pxyz,trc->GetTrack()->Charge());
  double* parTr  = (double*)trc->GetTrack()->GetParameter();
  double* parNew = (double*)tmpPr.GetTrack()->GetParameter();  
  for (int i=0;i<5;i++) parTr[i] = parNew[i];
  trc->SetMass(kMassMu);
  // 
  // set decay weight
  trc->GetXYZ(xyz);
  for (int i=3;i--;) xyz[i]-=fGenPnt[i];
  double dst = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
  double ctgamma = ctau*pDecParent.Gamma();
  double exparg = dst/ctgamma;
  //  double wgh =  exparg<100 ? 1.-TMath::Exp(-exparg) : 1.0;
  double wgh =  exparg<100 ? TMath::Exp(-exparg) : 0;
  // account for the losses due to the hadron interactions >>>
  double wabs = 0;
  int nb = fDecayZProf->GetNbinsX();
  for (int ib=1;ib<nb;ib++) {
    double x = fDecayZProf->GetBinCenter(ib);
    exparg = x/ctgamma;
    if (exparg>100) break;
    double wb = fDecayZProf->GetBinContent(ib);
    wabs += TMath::Exp(-exparg)*wb*fDecayZProf->GetBinWidth(ib);
    if (wb<1e-9) break;
  }
  wgh *= wabs/ctgamma;
  // account for the losses due to the hadron interactions <<<
  trc->SetWeight(wgh);
  //  printf("Decay %.3f Z:%+7.3f ctau: %f gamma %f wgh:%e\n",mass,fZDecay,ctau, pDecMu.Gamma(),wgh);
  //
  LUXEDebug(2,Form(" <<Mode:%d at Z=%f, PXYZ:%+6.3f %+6.3f %+6.3f, Mass:%.3f, Wgh:%.3e",fDecMode,fZDecay,pxyz[0],pxyz[1],pxyz[2],kMassMu,wgh));
  //LUXEDebug(2,Form(" <<Mode:%d at Z=%f, PXYZ:%+6.3f %+6.3f %+6.3f, Mass:%.3f, Wgh:%.3e",fDecMode,fZDecay,pxyz[0],pxyz[1],pxyz[2],kMassMu,wgh));
  //
}

//_________________________________________________________
void TrkDetector::InitDecayZHisto(double absorberLambda)
{
  // prepare a profile of Zdecay: uniform till start of the absorber, then exponentially dumped
  if (absorberLambda<1) absorberLambda = 1;
  TrkLayer* labs = 0;
  for (int i=0;i<fNTrkLayers;i++) {
    if ( (labs=GetTrkLayer(i))->IsAbs()) break;
    labs = 0;
  }
  if (!labs) {
    LUXEError("Could not find beginning of the absorber, is setup loaded?");
    exit(1);
  }
  double zdmp = labs->GetZ()-labs->GetThickness()/2;
  double zmax = GetTrkLayer(fLastActiveTrkLayer)->GetZ();
  LUXEDebug(2,Form("Decay will be done uniformly till Z=%.1f, then dumped with Lambda=%.1f cm",zdmp,absorberLambda));
  //
  int nbn = int(zmax-zdmp+1);
  TH1F* hd = new TH1F("DecayZProf","Z decay profile",nbn,0,zmax);
  for (int i=1;i<=nbn;i++) {
    double z = hd->GetBinCenter(i);
    if (z<zdmp) hd->SetBinContent(i,1.);
    else {
      double arg = (z-zdmp)/absorberLambda;
      hd->SetBinContent(i, arg>100 ? 0 : TMath::Exp(-arg));
    }
  }
  SetDecayZProfile(hd);
  //
}

//_________________________________________________________
void TrkDetector::GenBgEvent(double x, double y, double z, int offset)
{
  if (fNChPi<0 && fNChK<0 && fNChP<0) return;
  // generate bg. events from simple thermalPt-gaussian Y parameterization
  if (!fdNdYPi || !fdNdYK || !fdNdYP || !fdNdPtPi || !fdNdPtK || !fdNdPtP) LUXEFatal("Background generation was not initialized");
  //
  int maxLr = fLastActiveTrkLayerITS + offset;
  if (maxLr > fLastActiveTrkLayer) maxLr = fLastActiveTrkLayer;
  //
  for (int ilr=fLastActiveTrkLayer;ilr--;) {
    TrkLayer* lr = GetTrkLayer(ilr);
    if (lr->IsDead()) continue;
    lr->ResetBgClusters();
  }
  int decMode = fDecMode;
  fDecMode = kNoDecay;
  TrkProbe bgtr;
  //
  //  double ntr = gRandom->Poisson( fNCh );
//   for (int itr=0;itr<ntr;itr++) {
//     double yrap  = fdNdY->GetRandom();
//     double pt = fdNdPt->GetRandom();
//     double phi = gRandom->Rndm()*TMath::Pi()*2;
//     int charge = gRandom->Rndm()>0.5 ? 1:-1;
//     CreateTrkProbe(&bgtr, pt, yrap, phi, kMassPi, charge, x,y,z);
//     bgtr.SetTrID(itr);
//     TransportKalmanTrackWithMS(&bgtr, maxLr,kTRUE);
//   }
  int ntrTot = 0;
  
  //
  for (int ilr=maxLr;ilr--;) {
    TrkLayer* lr = GetTrkLayer(ilr);
    if (lr->IsDead()) continue;
    lr->SortBGClusters();
  }
  fDecMode = decMode;
  //  
}

//_________________________________________________________
void TrkDetector::InitBgGeneration(int dndeta, 
				      double y0, double sigy, double ymin,double ymax,
				      double T, double ptmin, double ptmax)
{
  // initialize bg generation routines
  fNCh = dndeta*0.5*(TMath::Erf((ymax-y0)/sqrt(2.)/sigy)-TMath::Erf((ymin-y0)/sqrt(2.)/sigy));
  fdNdY = new TF1("dndy","exp( -0.5*pow( (x-[0])/[1],2) )",ymin,ymax);
  fdNdY->SetParameters(y0,sigy);
  //
  fdNdPt = new TF1("dndpt","x*exp(-sqrt(x*x+1.949e-02)/[0])",ptmin,ptmax); // assume pion
  fdNdPt->SetParameter(0,T);
}

//_________________________________________________________
void TrkDetector::InitBgGenerationPart(double NPi,double NKplus,double NKminus ,double NP,double Piratio,
					  double y0, double y0Pi,double y0Kplus,double y0Kminus,double y0P,
					  double sigyPi, double sigyKplus,double sigyKminus, double sigyP,
					  double ymin,double ymax,
					  double Tpi, double TK, double TP, double ptmin, double ptmax)
{
  // initialize bg generation routines
  fNChPi = NPi*(1+Piratio)*sigyPi*TMath::Sqrt(TMath::Pi()/2.)*(TMath::Erf((ymax-y0-y0Pi)/sqrt(2.)/sigyPi)+ TMath::Erf((ymax-y0+y0Pi)/sqrt(2.)/sigyPi)-TMath::Erf((ymin-y0-y0Pi)/sqrt(2.)/sigyPi)-TMath::Erf((ymin-y0+y0Pi)/sqrt(2.)/sigyPi));
  fNChK = NKplus*sigyKplus*TMath::Sqrt(TMath::Pi()/2.)*(TMath::Erf((ymax-y0-y0Kplus)/sqrt(2.)/sigyKplus)+ TMath::Erf((ymax-y0+y0Kplus)/sqrt(2.)/sigyKplus)-TMath::Erf((ymin-y0-y0Kplus)/sqrt(2.)/sigyKplus)-TMath::Erf((ymin-y0+y0Kplus)/sqrt(2.)/sigyKplus))+
    NKminus*sigyKminus*TMath::Sqrt(TMath::Pi()/2.)*(TMath::Erf((ymax-y0-y0Kminus)/sqrt(2.)/sigyKminus)+ TMath::Erf((ymax-y0+y0Kminus)/sqrt(2.)/sigyKminus)-TMath::Erf((ymin-y0-y0Kminus)/sqrt(2.)/sigyKminus)-TMath::Erf((ymin-y0+y0Kminus)/sqrt(2.)/sigyKminus));
  fNChP = NP*sigyP*TMath::Sqrt(TMath::Pi()/2.)*(TMath::Erf((ymax-y0-y0P)/sqrt(2.)/sigyP)+ TMath::Erf((ymax-y0+y0P)/sqrt(2.)/sigyP)-TMath::Erf((ymin-y0-y0P)/sqrt(2.)/sigyP)-TMath::Erf((ymin-y0+y0P)/sqrt(2.)/sigyP));
  //
  fdNdYPi = new TF1("dndy","exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2) )",ymin,ymax);
  fdNdYPi->SetParameters(y0,y0Pi,sigyPi);
  fdNdYK  = new TF1("dndy","exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2) )",ymin,ymax);
  fdNdYK ->SetParameters(y0,y0Kplus ,sigyKplus);
  fdNdYP  = new TF1("dndy","exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2) )",ymin,ymax);
  fdNdYP ->SetParameters(y0,y0P,sigyP);
  //
  fdNdPtPi = new TF1("dndptPi","x*exp(-sqrt(x*x+1.949e-02)/[0])",ptmin,ptmax); // pion
  fdNdPtK  = new TF1("dndptK" ,"x*exp(-sqrt(x*x+0.493*0.493)/[0])",ptmin,ptmax); // kaon
  fdNdPtP  = new TF1("dndptP" ,"x*exp(-sqrt(x*x+0.938*0.938)/[0])",ptmin,ptmax); // proton
  fdNdPtPi->SetParameter(0,Tpi);
  fdNdPtK->SetParameter(0,TK);
  fdNdPtP->SetParameter(0,TP);
}

//_____________________________________________________________________
void TrkDetector::RequirePattern(UInt_t patt)
{
  // optional pattern to satyisfy
  if (!patt) return;
  int ngr = fPattITS.GetSize();
  fPattITS.Set(ngr+1);
  fPattITS[ngr] = patt;
}

//_____________________________________________________________________
void TrkDetector::BookControlHistos()
{
  fHChi2LrCorr = new TH2F("chi2Cl","chi2 corr cluster",        fNActiveTrkLayersITS+1,0,fNActiveTrkLayersITS+1,100,0,fMaxChi2Cl);
  fHChi2NDFCorr = new TH2F("chi2NDFCorr","chi2/ndf corr tr.",  fNActiveTrkLayersITS+1,0,fNActiveTrkLayersITS+1,100,0,fMaxNormChi2NDF);
  fHChi2NDFFake = new TH2F("chi2NDFFake","chi2/ndf fake tr.",  fNActiveTrkLayersITS+1,0,fNActiveTrkLayersITS+1,100,0,fMaxNormChi2NDF);
  fHChi2VtxCorr = new TH2F("chi2VCorr","chi2 to VTX corr tr." ,fNActiveTrkLayersITS+1,0,fNActiveTrkLayersITS+1,100,0,100);
  fHChi2VtxFake = new TH2F("chi2VFake","chi2 to VTX fake tr." ,fNActiveTrkLayersITS+1,0,fNActiveTrkLayersITS+1,100,0,100);
  fHNCand     = new TH2F("hNCand","Ncand per layer",           fNActiveTrkLayersITS+1,0,fNActiveTrkLayersITS+1,200,0,-1);
  fHCandCorID = new TH2F("CandCorID","Corr.cand ID per layer", fNActiveTrkLayersITS+1,0,fNActiveTrkLayersITS+1,200,0,-1);
  //
  fHChi2MS = new TH2F("chi2ms","chi2ms",100,0,30,10,0,10);
  //
}
//_________________________________________________________
void TrkDetector::InitBkg(double beamenergy){

  int E=TMath::Nint(beamenergy);

  // default values (from 40 GeV)
  double y0BG = 2.22; // gaussian y mean - 40 GeV
  double y0BGPi = 0.666;
  double y0BGKplus = 0.694;
  double y0BGKminus = 0.569;
  double y0BGP = 0.907;
  double sigyBG = 1.2; // .. sigma
  double sigyBGPi = 0.872;
  double sigyBGKplus = 0.725;
  double sigyBGKminus = 0.635;
  double sigyBGP = 0.798;
  double yminBG = 1.5; // min y to generate
  double ymaxBG = 4.5; //
  double TBG = 0.17;   // inv.slope of thermal pt distribution
  double TBGpi = 0.17;
  double TBGK = 0.23;
  double TBGP = 0.26;
  double ptminBG = 0.01;
  double ptmaxBG = 5;
  double dndyBGPi = 615.;
  double dndyBGK = 78.;
  double dndyBGP = 150.;
  double NBGPi = 74.;
  double NBGKplus = 16.2;
  double NBGKminus = 6.03;
  double NBGP = 37.5;
  double Piratio = 0.91;

  if (E == 20){
    printf("--- Background parameters for E=20 GeV/nucleon ---\n");
    y0BG = 1.9;   // gaussian y mean - 40 GeV
    sigyBG = 1.2; // .. sigma
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    TBG = 0.17;   // inv.slope of thermal pt distribution
    ptminBG = 0.01;
    ptmaxBG = 3;
    dndyBGPi = 410.;
    dndyBGK = 51.;
    dndyBGP = 148.;
    TBGpi = 0.17;
    TBGK = 0.22;
    TBGP = 0.26;
  }else if (E == 40){ 
    // pions and Kaons from  NA49 nucl-ex/0205002 
    printf("--- Background parameters for E=40 GeV/nucleon ---\n");
    y0BG = 2.22; // gaussian y mean - 40 GeV
    y0BGPi = 0.666;
    y0BGKplus = 0.694;
    y0BGKminus = 0.569;
    y0BGP = 0.907;
    sigyBG = 1.2; // .. sigma
    sigyBGPi = 0.872;
    sigyBGKplus = 0.725;
    sigyBGKminus = 0.635;
    sigyBGP = 0.798;
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    TBG = 0.17;   // inv.slope of thermal pt distribution
    TBGpi = 0.17;
    TBGK = 0.23;
    TBGP = 0.26;
    ptminBG = 0.01;
    ptmaxBG = 5;
    dndyBGPi = 615.;
    dndyBGK = 78.;
    dndyBGP = 150.;
    NBGPi = 74.;
    NBGKplus = 16.2;
    NBGKminus = 6.03;
    NBGP = 37.5;
    Piratio = 0.91;
  }else if (E == 60){ 
    // average of values at 40 and 80
    printf("--- Background parameters for E=60 GeV/nucleon ---\n");
    y0BG = 2.42;  // gaussian y mean - 60 GeV
    y0BGPi = 0.5*(0.666+0.756);
    y0BGKplus = 0.5*(0.694+0.742);
    y0BGKminus = 0.5*(0.569+0.668);
    y0BGP = 0.5*(0.907+0.907);
    sigyBG = 1.2; // .. sigma
    sigyBGPi = 0.5*(0.872+0.974);
    sigyBGKplus = 0.5*(0.725+0.792);
    sigyBGKminus = 0.5*(0.635+0.705);
    sigyBGP = 0.5*(0.798+0.798);
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    TBG = 0.17;   // inv.slope of thermal pt distribution
    TBGpi = 0.17;
    TBGK = 0.23;
    TBGP = 0.26;
    ptminBG = 0.01;
    ptmaxBG = 5;
    dndyBGPi =  0.5*(615.+920.);
    dndyBGK =  0.5*(78.+109.);
    dndyBGP =  0.5*(150.+(30.1/41.3)*150.);
    NBGPi = 0.5*(74.+97.);
    NBGKplus = 0.5*(16.2+19.3);
    NBGKminus = 0.5*(6.03+9.16);
    NBGP = 0.5*(37.5+(30.1/41.3)*37.5);
    Piratio = 0.93;
  }else if (E == 80){
    // pions and Kaons from  NA49 nucl-ex/0205002 
    printf("--- Background parameters for E=80 GeV/nucleon ---\n");
    y0BG = 2.57;  // gaussian y mean - 80 GeV
    y0BGPi = 0.756;
    y0BGKplus = 0.742;
    y0BGKminus = 0.668;
    y0BGP = 0.907;
    sigyBGPi = 0.974;
    sigyBGKplus = 0.792;
    sigyBGKminus = 0.705;
    sigyBGP = 0.798;
    sigyBG = 1.2; // .. sigma
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    TBG = 0.18;   // inv.slope of thermal pt distribution
    TBGpi = 0.18;
    TBGK = 0.23;
    TBGP = 0.26;
    ptminBG = 0.01;
    ptmaxBG = 5;
    dndyBGPi = 920.;
    dndyBGK = 109.;
    dndyBGP = (30.1/41.3)*150.;
    NBGPi = 97.;
    NBGKplus = 19.3;
    NBGKminus = 9.16;
    NBGP = (30.1/41.3)*37.5; //ratio 80/40 from PRC73, 044910 (2006)
    Piratio = 0.94;
  }else if (E == 160){
    // pions and Kaons from  NA49 nucl-ex/0205002 
    printf("--- Background parameters for E=160 GeV/nucleon ---\n");
    y0BG = 2.9; // gaussian y mean - 160 GeV
    y0BGPi = 0.72;
    y0BGKplus = 0.839;
    y0BGKminus = 0.727;
    y0BGP = 39.8;
    sigyBGPi = 1.18;
    sigyBGKplus = 0.88;
    sigyBGKminus = 0.81;
    sigyBGP = 8.07;
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    ptminBG = 0.01;
    ptmaxBG = 5;
    dndyBGPi = 1258.;
    dndyBGK = 155.;
    dndyBGP = 292.;
    NBGPi = 107.6;
    NBGKplus = 23.4;
    NBGKminus = 12.8;
    NBGP = 2.55e+06;
    Piratio = 0.97;
    TBGpi = 0.18; // inv.slope of thermal pt distribution
    TBGK = 0.23;  // inv.slope of thermal pt distribution
    TBGP = 0.31;  // inv.slope of thermal pt distribution
  }else if (E == 400){
    y0BG = 3.37;  // gaussian y mean - 80 GeV
    sigyBG = 1.2; // .. sigma
    yminBG = 1.5; // min y to generate
    ymaxBG = 4.5; //
    TBG = 0.17;   // inv.slope of thermal pt distribution
    ptminBG = 0.01;
    ptmaxBG = 5;
    dndyBGPi = 615.;
    dndyBGK = 78.;
    dndyBGP = 150.;
    TBGpi = 0.17;
    TBGK = 0.23;
    TBGP = 0.25;
  }else{
    printf("--- Parameters not available at this energy, use those for E=40 GeV/nucleon ---\n");
  }
  
  printf("Simulation of background at %d GeV/nucleon\n", E);
  printf("pions:   total multiplicity = %f ; T = %f\n", dndyBGPi, TBGpi);
  printf("kaons:   total multiplicity = %f ; T = %f\n", dndyBGK, TBGK);
  printf("protons: total multiplicity = %f ; T = %f\n", dndyBGP, TBGP);

  InitBgGenerationPart(NBGPi, NBGKplus, NBGKminus, NBGP, Piratio, y0BG, y0BGPi, y0BGKplus, y0BGKminus, y0BGP, sigyBGPi, sigyBGKplus, sigyBGKminus, sigyBGP, yminBG, ymaxBG, TBGpi, TBGK, TBGP, ptminBG, ptmaxBG);
  return;
}

//_____________________________________________________________________
Bool_t TrkDetector::IsCorrect(TrkProbe *probTr)
{
  if (probTr->GetNFakeITSHits()) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________
void  TrkDetector::SetMinITSHits(int n)
{
  fMinITSHits = TMath::Min(n,fNActiveTrkLayersITS);
}

//_____________________________________________________________________
void  TrkDetector::SetMinMSHits(int n)
{
  fMinMSHits = TMath::Min(n,fNActiveTrkLayersMS);
}
//_____________________________________________________________________
void  TrkDetector::SetMinTRHits(int n)
{
  fMinTRHits = TMath::Min(n,fNActiveTrkLayersTR);
}

//_____________________________________________________________________
void  TrkDetector::ForceLastActiveTrkLayer(int lr)
{
  printf("Attention: overriding last active layer from %d to %d\n",fLastActiveTrkLayer,lr);
  fLastActiveTrkLayer = lr;
}

//=======================================================================================


