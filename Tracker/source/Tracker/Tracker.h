//////////////////////////////////////////////////////
//// this is for the creation of Magnet dipole   /////
//// prepared by arka.santra@cern.ch             /////
//////////////////////////////////////////////////////
#ifndef TRACKER_H
#define TRACKER_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <cctype>
#include <cmath>
#include "TString.h"


using namespace std;

class Tracker {       // The class
public:             
    // Constructor
    Tracker();
    // Destructor
    ~Tracker();
    
    float buildGeometry(float height, float width, float length);
    
    
private:
    /// private variables
    float volume;
};

#endif
