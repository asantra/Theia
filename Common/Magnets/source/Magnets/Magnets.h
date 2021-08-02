//////////////////////////////////////////////////////
//// this is for the creation of Magnet dipole   /////
//// prepared by arka.santra@cern.ch             /////
//////////////////////////////////////////////////////
#ifndef MAGNETS_H
#define MAGNETS_H

#include <iostream>
#include <cmath>
#include "TPolyLine3D.h"
#include "TColor.h"
#include "TString.h"
#include "TVector3.h"


using namespace std;

class Magnets {       // The class
public:             
    // Constructor
    Magnets();
    Magnets(float xWdipole, float yHdipole){
        xWidth = xWdipole;
        yHeight = yHdipole;
    }
    Magnets(float xWdipole, float yHdipole, float z1dipole, float z2dipole){     // Constructor
       std::cout << "The magnet dimension " << xWdipole << " X " << yHdipole << " X " << std::abs(z1dipole-z2dipole) << std::endl;
    }
    // Destructor
    ~Magnets();
    static TPolyLine3D* GetDipole(float xWdipole, float yHdipole, float z1dipole, float z2dipole, Color_t col);
    
    // Setter
    void SetVolume(float xWdipole, float yHdipole, float z1dipole, float z2dipole);
    void SetXWidth(float xWdipole);
    void SetYHeight(float yHdipole);
    void SetMagneticField(TString process);
    void SetMagnetCenter(float x, float y, float z);
    void SetMagneticFieldIrregularity(TString process);
    // Getter
    float GetVolume();
    float GetXWidth();
    float GetXHalfWidth();
    float GetYHeight(); 
    float GetYHalfHeight();
    float GetMagneticField();

private:
    /// private variables
    float magnetVolume;
    float xWidth;
    float yHeight;
    float magneticField;
    TVector3 magnetCenter;
};

#endif
