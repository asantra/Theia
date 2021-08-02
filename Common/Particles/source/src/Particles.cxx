#include "Particles.h"


Particles::Particles() {}
Particles::~Particles() {}

Particles::Particles(float mass, float charge, float pdgid, float vtxx, float vtxy, float vtxz){
    particleMass   = mass;
    particleCharge = charge;
    particleId     = pdgid;
    particleVtxx   = vtxx;
    particleVtxy   = vtxy;
    particleVtxz   = vtxz;
}

// Setter
void Particles::SetMass(float mass){
    particleMass = mass;
}
void Particles::SetCharge(float charge){
    particleCharge = charge;
}
void Particles::SetPdgId(float pdgid){
    particleId = pdgid;
}
void Particles::SetVtxPosition(float vtxx, float vtxy, float vtxz){
    vtxPosition.SetXYZ(vtxx, vtxy, vtxz);
}
void Particles::SetLocalPosition(float x, float y, float z){
    localPosition.SetXYZ(x,y,z);
}
void Particles::SetMomentum(float px, float py, float pz, float E){
    momentum.SetPxPyPzE(px, py, pz, E);
}


// Getter
TVector3 Particles::GetPosition(){
    return localPosition;
}
TVector3 Particles::GetVtxPosition(){
    return vtxPosition;
}
TLorentzVector Particles::GetMomentum(){
    return momentum;
}
    
    
    
    
    
