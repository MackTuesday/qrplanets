#include <string>
#include <math.h>

#include "BrentsRand.h"
#include "GlobalConstants.h"
#include "Star.h"


Star::Star()
{
    double specClassDieRoll = BrentsUnitRand();
    double lumClassDieRoll = BrentsUnitRand();
        
    //var probs[0.01, 0.03, 0.0313, 0.0373, 0.0673, 0.1433, 0.2643];
    //var probs = [0.0, 0.0, 0.0, 0.0, 0.33, 0.67, 1.0];
    double probs[7] = {0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875};
    if (BrentsUnitRand() < probs[0]) {
        mSpectralClass = "M";
        mLuminosityClass = "III";
        double lumRoll = BrentsUnitRand();
        mLuminosity = pow(10.0, lumRoll * 1.04 + 3.66);
        mTemperature = pow(10.0, lumRoll * 0.188 + 3.38);
        mMass = pow(10.0, BrentsUnitRand() * 0.456 - 0.155);
    }
    else if (specClassDieRoll < probs[1]) {
        mSpectralClass = "D";
        mLuminosityClass = "V";
        double lumRoll = BrentsUnitRand();
        mLuminosity = pow(10.0, lumRoll * 2.7 - 4.25);
        mTemperature = pow(10.0, lumRoll * 1.57 + 3.6);
        mMass = pow(10.0, BrentsUnitRand() * 0.368 - 0.222);
    }
    else if (specClassDieRoll < probs[2]) {
        mSpectralClass = "B";
        mLuminosityClass = (lumClassDieRoll < 0.15) ? "IV" : "V";
        double lumRoll = BrentsUnitRand();
        mLuminosity = pow(10.0, lumRoll * 4.5 + 1.4);
        mLuminosity *= (lumClassDieRoll < 0.15) ? (BrentsUnitRand() * 0.25 + 1.0) : 1.0;
        mTemperature = pow(10.0, lumRoll * 0.3 + 4.0);
        mMass = pow(10.0, BrentsUnitRand() * 0.882 + 0.322);
    }
    else if (specClassDieRoll < probs[3]) {
        mSpectralClass = "A";
        mLuminosityClass = (lumClassDieRoll < 0.02) ? "IV" : "V";
        double lumRoll = BrentsUnitRand();
        mLuminosity = pow(10.0, lumRoll * 0.7 + 0.7);
        mLuminosity *= (lumClassDieRoll < 0.02) ? (BrentsUnitRand() * 0.7 + 0.8) : 1.0;
        mTemperature = pow(10.0, lumRoll * 0.125 + 3.875);
        mMass = pow(10.0, BrentsUnitRand() * 0.176 + 0.146);
    }
    else if (specClassDieRoll < probs[4]) {
        mSpectralClass = "F";
        mLuminosityClass = (lumClassDieRoll < 0.08) ? "IV" : "V";
        double lumRoll = BrentsUnitRand();
        mLuminosity = pow(10.0, lumRoll * 0.52 + 0.18);
        mLuminosity *= (lumClassDieRoll < 0.08) ? (BrentsUnitRand() * 0.7 + 1.1) : 1.0;
        mTemperature = pow(10.0, lumRoll * 0.0969 + 3.778);
        mMass = pow(10.0, BrentsUnitRand() * 0.129 + 0.017);
    }
    else if (specClassDieRoll < probs[5]) {
        mSpectralClass = "G";
        mLuminosityClass = (lumClassDieRoll < 0.3) ? "IV" : "V";
        double lumRoll = BrentsUnitRand();
        mLuminosity = pow(10.0, lumRoll * 0.4 - 0.22);
        mLuminosity *= (lumClassDieRoll < 0.3) ? (BrentsUnitRand() * 0.7 + 1.5) : 1.0;
        mTemperature = pow(10.0, lumRoll * 0.0621 + 3.716);
        mMass = pow(10.0, BrentsUnitRand() * 0.114 - 0.097);
    }
    else if (specClassDieRoll < probs[6]) {
        mSpectralClass = "K";
        mLuminosityClass = (lumClassDieRoll < 0.09) ? "IV" : "V";
        double lumRoll = BrentsUnitRand();
        mLuminosity = pow(10.0, lumRoll * 0.88 - 1.1);
        mLuminosity *= (lumClassDieRoll < 0.09) ? (BrentsUnitRand() * 0.3 + 0.6) : 1.0;
        mTemperature = pow(10.0, lumRoll * 0.148 + 3.568);
        mMass = pow(10.0, BrentsUnitRand() * 0.192 - 0.347);
    }
    else {
        mSpectralClass = "M";
        mLuminosityClass = "V";
        double lumRoll = BrentsUnitRand();
        mLuminosity = pow(10.0, lumRoll * 1.9 - 3.0);
        mTemperature = pow(10.0, lumRoll * 0.188  + 3.38);
        mMass = pow(10.0, BrentsUnitRand() * 0.751 - 1.097);
    }
    
    mMass *= kMassOfSol;
    
    if (mLuminosityClass == "III" || mLuminosityClass == "IV" ||
        mSpectralClass == "M" || mSpectralClass == "D") {
        mAge = kSystemAgesMature;
    }
    else if (mSpectralClass == "B" || mSpectralClass == "A") {
        mAge = kSystemAgesYoung;
    }
    else {
        double ageDieRoll = BrentsUnitRand();
        if (mSpectralClass == "F") {
            mAge = (ageDieRoll < 0.33) ? kSystemAgesYoung : kSystemAgesMature;
        }
        else if (mSpectralClass == "G") {
            mAge = (ageDieRoll < 0.1) ? kSystemAgesYoung : kSystemAgesMature;
        }
        else if (mSpectralClass == "K") {
            mAge = (ageDieRoll < 0.02) ? kSystemAgesYoung : kSystemAgesMature;
        }
        else {
            mAge = kSystemAgesMature;
        }
    }
}
