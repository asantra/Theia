#include "Detector.h"
#include "Probe.h"
#include "DrawProbe.h"

using namespace std;

typedef map<int,TString>             TMapiTS;
TMapiTS         islayers = { {1,"L1"}, {3,"L2"}, {5,"L3"}, {7,"L4"} };

Detector* det = 0;
vector<double>* zlayer = new vector<double>;

bool islayer(double z)
{
	for(int j=0 ; j<(int)zlayer->size() ; ++j)
	{
		double dz = abs(zlayer->at(j)-z);
		if(dz<1.e-6) return true;
	}
	return false;
}

TPolyLine3D* DrawProbe::TrackLine3d(const Probe* source, Double_t zMax, Double_t step=1, Color_t col=kBlack)
{
	double xyz[3];
	source->GetXYZ(xyz);
	double zCurr = xyz[2]; //source->GetZ();
	int nZ = (zMax - zCurr)/step + 1;
	if (nZ<2) {
		printf("bad limits\n");
		return 0;
	}
	Probe tmp(*source);
	double xp[nZ],yp[nZ],zp[nZ];
	xp[0] = xyz[0];
	yp[0] = xyz[1];
	zp[0] = xyz[2];
	int nz = 0;
	for (int iz=1;iz<nZ;iz++) {
		if(!det->PropagateToZBxByBz(&tmp, TMath::Min(tmp.GetZ()+step, zMax), step)) break; //propagation may fail..
		tmp.GetXYZ(xyz);
		xp[iz] = xyz[0];
		yp[iz] = xyz[1];
		zp[iz] = xyz[2];
		nz++;
	}
	TPolyLine3D *polyline = new TPolyLine3D(nz+1);
	polyline->SetLineColor(col);
	for (int i=0;i<nz+1;i++) {
		polyline->SetPoint(i,xp[i],yp[i],zp[i]);
	}
	return polyline;
}


TPolyMarker3D* DrawProbe::TrackMarker3d(const Probe* source, double zmin, double zmax, double zstep, Color_t col=kBlack)
{
	Probe tmp(*source);
	int nZ = (int)(zmax-zmin)/zstep;
	double xp[nZ],yp[nZ],zp[nZ];
	double xyz[3];
	tmp.GetXYZ(xyz);
	xp[0] = xyz[0];
	yp[0] = xyz[1];
	zp[0] = xyz[2];
	int nz = 0;
	for(int iz=1;iz<nZ;iz++) {
		if(!det->PropagateToZBxByBz(&tmp, tmp.GetZ()+zstep, zstep)) break; //propagation may fail...
		tmp.GetXYZ(xyz);
		xp[iz] = xyz[0];
		yp[iz] = xyz[1];
		zp[iz] = xyz[2];
		nz++;
	}
	TPolyMarker3D *polymarker = new TPolyMarker3D(zlayer->size());
	polymarker->SetMarkerColor(col);
	int n = 0;
	for(int i=0;i<nz+1;i++) {
		if(!islayer(zp[i])) continue;
		polymarker->SetPoint(n,xp[i],yp[i],zp[i]);
		n++;
	}
	return polymarker;
}


Color_t DrawProbe::trkcol(double E)
{
	if     (E>=14)           return kBlack;
	else if(E<14. and E>=12) return kRed;
	else if(E<12. and E>=10) return 95;
	else if(E<10. and E>=8.) return 91;
	else if(E<8.  and E>=7.) return 80;
	else if(E<7.  and E>=6.) return 71;
	else if(E<6.  and E>=5.) return 65;
	else if(E<5.  and E>=4.) return 60;
	else if(E<4.  and E>=3.) return 53;
	else if(E<3.  and E>=2.) return 51;
	else                     return 6;
	return kRed;
}
