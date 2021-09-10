#ifndef LUXEUTILITY_h
#define LUXEUTILITY_h
#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cmath>
#include "TVector2.h"
#include "TMath.h"
#include "TString.h"

using namespace std;
TVector2 rUnit2(TVector2 r1, TVector2 r2)
{
	TVector2 r = (r2-r1).Unit();
	return r;
}

float xofz(float* r1,float* r2, float z)
{
	float dz = r2[2]-r1[2];
	float dx = r2[0]-r1[0];
	if(dz==0)
	{
		cout << "ERROR in xofz: dz=0" << endl;
		exit(-1);
	}
	float a = dx/dz;
	float b = r1[0]-a*r1[2];
	float x = a*z+b;
	// cout << "in xofz: x=" << x << ", dz=" << dz << endl;
	return x;
}
	
float yofz(float* r1,float* r2, float z)
{
	float dz = r2[2]-r1[2];
	float dy = r2[1]-r1[1];
	if(dz==0)
	{
		cout << "ERROR in yofz: dz=0" << endl;
		exit(-1);
	}
	float a = dy/dz;
	float b = r1[1]-a*r1[2];
	float y = a*z+b;
	return y;
}

float zofx(float* r1,float* r2, float x)
{
	float dz = r2[2]-r1[2];
	float dx = r2[0]-r1[0];
	if(dx==0)
	{
		cout << "ERROR in zofx: dx=0" << endl;
		exit(-1);
	}
	float a = dz/dx;
	float b = r1[2]-a*r1[0];
	float z = a*x+b;
	return z;
}


void SetLogBins(Int_t nbins, Double_t min, Double_t max, Double_t* xpoints)
{
	Double_t logmin  = log10(min);
	Double_t logmax  = log10(max);
	Double_t logbinwidth = (Double_t)( (logmax-logmin)/(Double_t)nbins );
	xpoints[0] = min;
	for(Int_t i=1 ; i<=nbins ; i++) xpoints[i] = TMath::Power( 10,(logmin + i*logbinwidth) );
}


bool foundinvec(int x, vector<int>& v)
{
	vector<int>::iterator it = find(v.begin(),v.end(),x);
	return (it!=v.end());
}

int getvecindex(int x, vector<int>& v)
{
	vector<int>::iterator it = find(v.begin(),v.end(),x);
	return (it!=v.end()) ? distance(v.begin(),it) : -1;
}

int toint(TString str)
{
	stringstream strm;
	int x;
	strm << str;
	strm >> x;
	return x;
}

TString FormatEventID(int evnt)
{
	TString sevnt = "";
	if(evnt<10)                          sevnt = Form("000000%d", evnt);
	if(evnt>=10 && evnt<100)             sevnt = Form("00000%d", evnt);
	if(evnt>=100 && evnt<1000)           sevnt = Form("0000%d", evnt);
	if(evnt>=1000 && evnt<10000)         sevnt = Form("000%d", evnt);
	if(evnt>=10000 && evnt<100000)       sevnt = Form("00%d", evnt);
	if(evnt>=100000 && evnt<1000000)     sevnt = Form("0%d", evnt);
	if(evnt>=1000000 && evnt<10000000)   sevnt = Form("%d", evnt); 
	if(evnt>=10000000 && evnt<100000000) sevnt = Form("%d", evnt); // assume no more than 9,999,999,999 events...
	return sevnt;
}	

#endif
