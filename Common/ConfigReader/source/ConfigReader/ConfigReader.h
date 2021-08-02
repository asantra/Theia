//////////////////////////////////////////////////////
//// this is for the creation of Magnet dipole   /////
//// prepared by arka.santra@cern.ch             /////
//////////////////////////////////////////////////////
#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <cctype>
#include <cmath>
#include "TString.h"


using namespace std;

class ConfigReader {       // The class
public:             
    // Constructor
    ConfigReader();
    // Destructor
    ~ConfigReader();
    
    void readConfigFile(string configFileName);
    
    
private:
    /// private variables
    string fileName;
};

#endif
