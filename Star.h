#ifndef Star_h
#define Star_h


#include <string>

using namespace std;


struct Star
{
    Star();
    
    string mSpectralClass;
    string mLuminosityClass;
    double mLuminosity;
    double mTemperature;
    double mMass;
    unsigned mAge;
};


#endif
