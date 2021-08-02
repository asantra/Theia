#include "Cluster.h"
#include <TString.h>

ClassImp(Cluster)


Cluster::~Cluster() {}

//_________________________________________________________________________
Cluster::Cluster(Cluster &src) 
:  TObject(src)
  ,fX(src.fX)
  ,fY(src.fY)
  ,fZ(src.fZ)
{}

//__________________________________________________________________________
Cluster& Cluster::operator=(const Cluster& src) 
{
  if (this!=&src) {
    TObject::operator=(src);
    fX = src.fX;
    fY = src.fY;
    fZ = src.fZ;
  }
  return *this;
}

//_________________________________________________________________________
void Cluster::Print(Option_t *opt) const 
{
  TString opts = opt;
  opts.ToLower();
  if (opts.Contains("lc")) 
    printf("Tr#%4d Loc (%+.4e,%+.4e %+.4e) %s",GetTrID(),fX,fY,fZ,IsKilled()?"Killed":""); 
  else 
    printf("Tr#%4d Lab (%+.4e,%+.4e %+.4e) %s",GetTrID(),GetXLab(),GetYLab(),GetZLab(),IsKilled()?"Killed":""); 
  if (opts.Contains("nl")) return;
  printf("\n");
}

//_________________________________________________________________________
Bool_t Cluster::IsEqual(const TObject* obj) const 
{
  // check if clusters are equal
  Cluster* cl = (Cluster*)obj;
  const double kTiny = 1e-12;
  if (TMath::Abs(GetXLab()-cl->GetXLab())>kTiny || 
      TMath::Abs(GetYLab()-cl->GetYLab())>kTiny) return kFALSE;
  return kTRUE;
}

//_________________________________________________________________________
Int_t Cluster::Compare(const TObject* obj) const
{
  // compare 1st labx, then laby
  Cluster* cl = (Cluster*)obj;
  if (GetXLab() > cl->GetXLab()) return -1; // tracking -Z = labX
  if (GetXLab() < cl->GetXLab()) return  1;
  if (GetYLab() < cl->GetYLab()) return -1; // tracking Y = labY
  if (GetYLab() > cl->GetYLab()) return  1;
  return 0;
}
