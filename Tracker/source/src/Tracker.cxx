#include "Tracker.h"

using namespace std;

Tracker::Tracker() {}
Tracker::~Tracker() {}

float Tracker::buildGeometry(float height, float width, float length){
    volume = height*width*length;
    return volume;
}
