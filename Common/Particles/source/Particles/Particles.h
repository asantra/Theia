#ifndef PARTICLES_H
#define PARTICLES_H
#include <iostream>
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"

class Particles{
public:  
  // Constructor
    Particles();
    Particles(float mass, float charge, float pdgid, float vtxx=0, float vtxy=0, float vtxz=0);
    
    // Destructor
    ~Particles();
    
    // Setter
    void SetMass(float mass);
    void SetPdgId(float pdgid);
    void SetCharge(float charge);
    void SetVtxPosition(float vtxx, float vtxy, float vtxz);
    void SetLocalPosition(float x, float y, float z);
    void SetMomentum(float px, float py, float pz, float Ez);
    // Getter
    TVector3 GetPosition();
    TVector3 GetVtxPosition();
    TLorentzVector GetMomentum();
  
    
private:
    float particleMass;
    float particleCharge;
    float particleId;
    float particleVtxx;
    float particleVtxy;
    float particleVtxz;
    float particleLocalX;
    float particleLocalY;
    float particleLocalZ;
    TVector3 vtxPosition;
    TVector3 localPosition;
    TLorentzVector momentum;
};

#endif
