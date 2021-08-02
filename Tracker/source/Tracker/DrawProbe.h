#ifndef DRAWPROBE_H
#define DRAWPROBE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include "TPolyLine3D.h"
#include "TColor.h"
#include "TString.h"
#include "Probe.h"

class DrawProbe{
public:
    static TPolyLine3D* TrackLine3d(const Probe* source, Double_t zMax, Double_t step, Color_t col);
    static TPolyMarker3D* TrackMarker3d(const Probe* source, double zmin, double zmax, double zstep, Color_t col);
    static Color_t trkcol(double E);
};
#endif
