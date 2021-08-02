#include "Magnets.h"

Magnets::Magnets() {}
Magnets::~Magnets() {}

/// get the diagram of the dipole
TPolyLine3D* Magnets::GetDipole(float xWdipole, float yHdipole, float z1dipole, float z2dipole, Color_t col)
{
    TPolyLine3D* polyline = new TPolyLine3D();
    polyline->SetPoint(0,-xWdipole/2,-yHdipole/2,z1dipole);
    polyline->SetPoint(1,-xWdipole/2,+yHdipole/2,z1dipole);
    polyline->SetPoint(2,+xWdipole/2,+yHdipole/2,z1dipole);
    polyline->SetPoint(3,+xWdipole/2,-yHdipole/2,z1dipole);
    polyline->SetPoint(4,-xWdipole/2,-yHdipole/2,z1dipole);

    polyline->SetPoint(5,-xWdipole/2,-yHdipole/2,z2dipole); // go up
    polyline->SetPoint(6,-xWdipole/2,+yHdipole/2,z2dipole); // move
    polyline->SetPoint(7,-xWdipole/2,+yHdipole/2,z1dipole); // go down
    polyline->SetPoint(8,-xWdipole/2,+yHdipole/2,z2dipole); // up again

    polyline->SetPoint(9,+xWdipole/2,+yHdipole/2,z2dipole); // move
    polyline->SetPoint(10,+xWdipole/2,+yHdipole/2,z1dipole); // go down
    polyline->SetPoint(11,+xWdipole/2,+yHdipole/2,z2dipole); // up again

    polyline->SetPoint(12,+xWdipole/2,-yHdipole/2,z2dipole); // move
    polyline->SetPoint(13,+xWdipole/2,-yHdipole/2,z1dipole); // go down
    polyline->SetPoint(14,+xWdipole/2,-yHdipole/2,z2dipole); // up again

    polyline->SetPoint(15,-xWdipole/2,-yHdipole/2,z2dipole); // move
    polyline->SetPoint(16,-xWdipole/2,-yHdipole/2,z1dipole); // go down
    polyline->SetPoint(17,-xWdipole/2,-yHdipole/2,z2dipole); // up again
    polyline->SetLineColor(col);
    return polyline;
    
}

// setter, set volume, xwidth, yheight and magnetic field.
void Magnets::SetVolume(float xWdipole, float yHdipole, float z1dipole, float z2dipole) {
    magnetVolume = xWdipole*yHdipole*std::abs(z1dipole-z2dipole);
}
void Magnets::SetXWidth(float xWdipole){
    xWidth = xWdipole;
}
void Magnets::SetYHeight(float yHdipole){
    yHeight = yHdipole;
}
/// magnetic field is hardcoded for now
void Magnets::SetMagneticField(TString process="bppp"){
    magneticField = (process="bppp") ? 1.7 : 1.0;
}
void Magnets::SetMagnetCenter(float x, float y, float z){
    magnetCenter.SetXYZ(x,y,z);
}

// Getter
float Magnets::GetVolume() {
    return magnetVolume;
}
float Magnets::GetXWidth(){
    return xWidth;
}
float Magnets::GetXHalfWidth(){
    return xWidth/2.0;
}
float Magnets::GetYHeight(){
    return yHeight;
}
float Magnets::GetYHalfHeight(){
    return yHeight/2.0;
}

float Magnets::GetMagneticField(){
    return magneticField;
}
void Magnets::SetMagneticFieldIrregularity(TString process="bppp"){
    float constmagneticField = (process="bppp") ? 1.7 : 1.0;
    float smearing = 0.5; // smearing in percentage, but this can be a position dependent value. Need to think later
}


