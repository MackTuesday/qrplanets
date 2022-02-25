#include <math.h>

#include "BrentsRand.h"
#include "GlobalConstants.h"
#include "ChooseMass.h"


double kMassBoundsSubjovianMin = 0.8;
double kMassBoundsSubjovianRange = 1.0;
double kMassBoundsJovianMin = 1.8;
double kMassBoundsJovianRange = 1.4;
double kMassBoundsTransjovianMin = 3.2;
double kMassBoundsTransjovianRange = 0.41;
double kMassBoundsSelenianMin = -1.5;
double kMassBoundsSelenianRange = 0.8;
double kMassBoundsGlacianMin = -1.5;
double kMassBoundsGlacianRange = 1.3;
double kMassBoundsTerranMin = -0.7;
double kMassBoundsTerranRange = 0.9;
double kMassBoundsNereanMin = -0.7;
double kMassBoundsNereanRange = 0.9;
double kMassBoundsAquarianMin = 0.3;
double kMassBoundsAquarianRange = 0.3;
double kMassBoundsHadeanMin = 0.3;
double kMassBoundsHadeanRange = 0.3;
double kMassBoundsVestanMin = 0.85;
double kMassBoundsVestanRange = 0.3;

double kDensityBoundsLesserSubjovianMin = 2.80;
double kDensityBoundsLesserSubjovianRange = 0.40;
double kDensityBoundsSubjovianMin = 2.90;
double kDensityBoundsSubjovianRange = 0.30;
double kDensityBoundsJovianMin = 2.95;
double kDensityBoundsJovianRange = 0.30;
// We don't use density bounds for 
// transjovians -- see code
double kDensityBoundsSelenianMin = 3.52;
double kDensityBoundsSelenianRange = 0.24;
double kDensityBoundsGlacianMin = 3.11;
double kDensityBoundsGlacianRange = 0.46;
double kDensityBoundsTerranMin = 3.52;
double kDensityBoundsTerranRange = 0.24;
double kDensityBoundsNereanMin = 3.52;
double kDensityBoundsNereanRange = 0.24;
double kDensityBoundsAquarianMin = 3.3;
double kDensityBoundsAquarianRange = 0.3;
double kDensityBoundsHadeanMin = 3.85;
double kDensityBoundsHadeanRange = 0.3;
double kDensityBoundsVestanMin = 3.75;
double kDensityBoundsVestanRange = 0.3;


//double jovianProb = 0.25;

// Assumes values in planetProbabilities are organized so the 
// probabilities of jovian classes are all at the end of the list.
/*void AdjustJovianProbs(double* planetProbabilities, unsigned numProbabilities,
                       unsigned indexOfFirstJovian, double newProbability,
                       double* adjustedProbabilities)
{
    double baseProbability = planetProbabilities[indexOfFirstJovian];
    // Not-jovian is the first part of the list, so increasing the
    // jovian probability means decreasing newProbability.
    double backwardProbability = 1.0 - newProbability;
    double probabilityRatio = backwardProbability / baseProbability;
    // p <== p * n / b
    // p <== (p - b) * (1 - n) / (1 - b) + n
    //       (p - 1 + 1 - b) * (1 - n) / (1 - b) + n
    //       (p - 1) * (1 - n) / (1 - b) + (1 - n) + n
    for (unsigned i = 0; i < indexOfFirstJovian; i++) {
        adjustedProbabilities[i] = planetProbabilities[i] * probabilityRatio;
    }
    for (unsigned i = indexOfFirstJovian; i < numProbabilities; i++) {
        adjustedProbabilities[i] = (planetProbabilities[i] - baseProbability) *
                                   (1.0 - backwardProbability) / 
                                   (1.0 - baseProbability) + backwardProbability;
    }
} */


void ChooseJovianNormalEpistellarMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[4] = {0.27, 0.46, 0.62, 0.92};
    //double newProbs[4];
    //AdjustJovianProbs(probs, 4, 1, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, massDieRoll * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTerranRange + kMassBoundsTerranMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsTerranRange + 
                            kDensityBoundsTerranMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + kMassBoundsSubjovianMin);
        if (mass < 17.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsLesserSubjovianRange + 
                                kDensityBoundsLesserSubjovianMin);
        }
        else {
            density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                                kDensityBoundsSubjovianMin);
        }
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsJovianRange + 
                         kMassBoundsJovianMin);
        if (mass < 160.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsJovianRange + 
                                kDensityBoundsJovianMin);
        }
        else {
            density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
        }
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTransjovianRange + 
                         kMassBoundsTransjovianMin);
        density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
    }
}



void ChooseJovianNormalSupertorridMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[7] = {0.27, 0.46, 0.49, 0.52, 0.61, 0.62, 0.92};
    //double newProbs[7];
    //AdjustJovianProbs(probs, 7, 4, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTerranRange + kMassBoundsTerranMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsTerranRange + 
                            kDensityBoundsTerranMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsHadeanRange + kMassBoundsHadeanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsHadeanRange + 
                            kDensityBoundsHadeanMin);
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsVestanRange + kMassBoundsVestanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsVestanRange + 
                            kDensityBoundsVestanMin);
    }
    else if (planetTypeDieRoll < probs[4]) {
        mass = 0.0;
        density = 0.0;
    }
    else if (planetTypeDieRoll < probs[5]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + 
                         kMassBoundsSubjovianMin);
        if (mass < 17.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsLesserSubjovianRange + 
                                kDensityBoundsLesserSubjovianMin);
        }
        else {
            density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                                kDensityBoundsSubjovianMin);
        }
    }
    else if (planetTypeDieRoll < probs[6]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsJovianRange + 
                         kMassBoundsJovianMin);
        if (mass < 160.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsJovianRange + 
                                kDensityBoundsJovianMin);
        }
        else {
            density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
        }
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTransjovianRange + 
                         kMassBoundsTransjovianMin);
        density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
    }
}


void ChooseJovianNormalTorridMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[7] = {0.27, 0.46, 0.49, 0.50, 0.59, 0.62, 0.92};
    //double newProbs[7];
    //AdjustJovianProbs(probs, 7, 4, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTerranRange + kMassBoundsTerranMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsTerranRange + 
                            kDensityBoundsTerranMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsHadeanRange + kMassBoundsHadeanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsHadeanRange + 
                            kDensityBoundsHadeanMin);
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsVestanRange + kMassBoundsVestanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsVestanRange + 
                            kDensityBoundsVestanMin);
    }
    else if (planetTypeDieRoll < probs[4]) {
        mass = 0.0;
        density = 0.0;
    }
    else if (planetTypeDieRoll < probs[5]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + kMassBoundsSubjovianMin);
        if (mass < 17.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsLesserSubjovianRange + 
                                kDensityBoundsLesserSubjovianMin);
        }
        else {
            density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                                kDensityBoundsSubjovianMin);
        }
    }
    else if (planetTypeDieRoll < probs[6]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsJovianRange + kMassBoundsJovianMin);
        if (mass < 160.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsJovianRange + 
                                 kDensityBoundsJovianMin);
        }
        else {
            density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
        }
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTransjovianRange + 
                              kMassBoundsTransjovianMin);
        density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
    }
}

// Get probs from if statement
// Note index of jovian threshold
// Make array from probs
// Adjust probs
// Put newProbs in if statements
void ChooseJovianNormalTemperateMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    
    double probs[8] = {0.27, 0.36, 0.46, 0.49, 0.52, 0.53, 0.62, 0.92};
    if (age == kSystemAgesMature) {
        probs[1] -= 0.04;
    }
    
    //double newProbs[8];
    //AdjustJovianProbs(probs, 8, 5, jovianProb, newProbs);

    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTerranRange + kMassBoundsTerranMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsTerranRange + 
                            kDensityBoundsTerranMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsNereanRange + kMassBoundsNereanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsNereanRange + 
                            kDensityBoundsNereanMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsHadeanRange + kMassBoundsHadeanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsHadeanRange + 
                            kDensityBoundsHadeanMin);
    }
    else if (planetTypeDieRoll < probs[4]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsVestanRange + kMassBoundsVestanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsVestanRange + 
                            kDensityBoundsVestanMin);
    }
    else if (planetTypeDieRoll < probs[5]) {
        mass = 0.0;
        density = 0.0;
    }
    else if (planetTypeDieRoll < probs[6]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + kMassBoundsSubjovianMin);
        if (mass < 17.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsLesserSubjovianRange + 
                                kDensityBoundsLesserSubjovianMin);
        }
        else {
            density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                                kDensityBoundsSubjovianMin);
        }
    }
    else if (planetTypeDieRoll < probs[7]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsJovianRange + kMassBoundsJovianMin);
        if (mass < 160.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsJovianRange + 
                                kDensityBoundsJovianMin);
        }
        else {
            density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
        }
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTransjovianRange + 
                         kMassBoundsTransjovianMin);
        density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
    }
}

void ChooseJovianNormalFrigidMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[4] = {0.27, 0.59, 0.62, 0.92};
    //double newProbs[4];
    //AdjustJovianProbs(probs, 4, 1, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsGlacianRange + kMassBoundsGlacianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsGlacianRange + 
                            kDensityBoundsGlacianMin);
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = 0.0;
        density = 0.0;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + kMassBoundsSubjovianMin);
        if (mass < 17.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsLesserSubjovianRange + 
                                kDensityBoundsLesserSubjovianMin);
        }
        else {
            density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                                kDensityBoundsSubjovianMin);
        }
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsJovianRange + kMassBoundsJovianMin);
        if (mass < 160.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsJovianRange + 
                                kDensityBoundsJovianMin);
        }
        else {
            density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
        }
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTransjovianRange + 
                         kMassBoundsTransjovianMin);
        density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
    }
}


void ChooseJovianNormalSuperfrigidMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[4] = {0.27, 0.59, 0.62, 0.92};
    //double newProbs[4];
    //AdjustJovianProbs(probs, 4, 1, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsGlacianRange + kMassBoundsGlacianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsGlacianRange + 
                            kDensityBoundsGlacianMin);
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = 0.0;
        density = 0.0;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + 
                         kMassBoundsSubjovianMin);
        if (mass < 17.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsLesserSubjovianRange + 
                                kDensityBoundsLesserSubjovianMin);
        }
        else {
            density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                                kDensityBoundsSubjovianMin);
        }
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsJovianRange + kMassBoundsJovianMin);
        if (mass < 160.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsJovianRange + 
                                kDensityBoundsJovianMin);
        }
        else {
            density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
        }
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTransjovianRange + 
                         kMassBoundsTransjovianMin);
        density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
    }
}

// Get probs from if statement
// Note index of jovian threshold
// Make array from probs
// Adjust probs
// Put newProbs in if statements
void ChooseJovianPostStellarInnerMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[6] = {0.27, 0.46, 0.49, 0.52, 0.62, 0.92};
    //double newProbs[6];
    //AdjustJovianProbs(probs, 6, 4, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTerranRange + kMassBoundsTerranMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsTerranRange + 
                            kDensityBoundsTerranMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsHadeanRange + kMassBoundsHadeanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsHadeanRange + 
                            kDensityBoundsHadeanMin);
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsVestanRange + kMassBoundsVestanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsVestanRange + 
                            kDensityBoundsVestanMin);
    }
    else if (planetTypeDieRoll < probs[4]) {
        mass = 0.0;
        density = 0.0;
    }
    else if (planetTypeDieRoll < probs[5]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsJovianRange + kMassBoundsJovianMin);
        if (mass < 160.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsJovianRange + 
                                kDensityBoundsJovianMin);
        }
        else {
            density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
        }
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTransjovianRange + 
                              kMassBoundsTransjovianMin);
        density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
    }
}

void ChooseJovianPostStellarMidOuterMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[7] = {0.27, 0.46, 0.49, 0.50, 0.58, 0.61, 0.94};
    //double newProbs[7];
    //AdjustJovianProbs(probs, 7, 4, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTerranRange + kMassBoundsTerranMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsTerranRange + 
                            kDensityBoundsTerranMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsHadeanRange + kMassBoundsHadeanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsHadeanRange + 
                            kDensityBoundsHadeanMin);
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsVestanRange + kMassBoundsVestanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsVestanRange + 
                            kDensityBoundsVestanMin);
    }
    else if (planetTypeDieRoll < probs[4]) {
        mass = 0.0;
        density = 0.0;
    }
    else if (planetTypeDieRoll < probs[5]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + kMassBoundsSubjovianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                            kDensityBoundsSubjovianMin);
    }
    else if (planetTypeDieRoll < probs[6]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsJovianRange + kMassBoundsJovianMin);
        if (mass < 160.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsJovianRange + 
                                     kDensityBoundsJovianMin);
        }
        else {
            density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
        }
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTransjovianRange + 
                              kMassBoundsTransjovianMin);
        density = mass / 160.0 * (densityDieRoll * 0.1 + 0.95) * 1000.0;
    }
}

void ChooseJovianPostStellarDeepMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[2] = {0.56, 0.72};
    //double newProbs[2];
    //AdjustJovianProbs(probs, 2, 1, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = 0.0;
        density = 0.0;
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + kMassBoundsSubjovianMin);
        if (mass < 17.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsLesserSubjovianRange + 
                                kDensityBoundsLesserSubjovianMin);
        }
        else {
            density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                                kDensityBoundsSubjovianMin);
        }
    }
}

void ChooseNonjovianNormalSupToTempMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[4] = {0.27, 0.49, 0.62, 0.83};
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTerranRange + kMassBoundsTerranMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsTerranRange + 
                            kDensityBoundsTerranMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsHadeanRange + kMassBoundsHadeanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsHadeanRange + 
                            kDensityBoundsHadeanMin);
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsVestanRange + kMassBoundsVestanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsVestanRange + 
                            kDensityBoundsVestanMin);
    }
    else {
        mass = 0.0;
        density = 0.0;
    }
}

void ChooseNonjovianNormalFrigidMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[2] = {0.56, 0.72};
    //double newProbs[2];
    //AdjustJovianProbs(probs, 2, 1, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsGlacianRange + kMassBoundsGlacianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsGlacianRange + 
                            kDensityBoundsGlacianMin);
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = 0.0;
        density = 0.0;
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + kMassBoundsSubjovianMin);
        if (mass < 17.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsLesserSubjovianRange + 
                                kDensityBoundsLesserSubjovianMin);
        }
        else {
            density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                                kDensityBoundsSubjovianMin);
        }
    }
}

void ChooseNonjovianNormalSuperfrigidMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[2] = {0.56, 0.72};
    //double newProbs[2];
    //AdjustJovianProbs(probs, 2, 1, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsGlacianRange + kMassBoundsGlacianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsGlacianRange + 
                            kDensityBoundsGlacianMin);
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = 0.0;
        density = 0.0;
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + kMassBoundsSubjovianMin);
        if (mass < 17.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsLesserSubjovianRange + 
                                kDensityBoundsLesserSubjovianMin);
        }
        else {
            density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                                kDensityBoundsSubjovianMin);
        }
    }
}

void ChooseNonjovianPostStellarInnerMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[4] = {0.27, 0.49, 0.62, 0.83};
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTerranRange + kMassBoundsTerranMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsTerranRange + 
                            kDensityBoundsTerranMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsHadeanRange + kMassBoundsHadeanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsHadeanRange + 
                            kDensityBoundsHadeanMin);
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsVestanRange + kMassBoundsVestanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsVestanRange + 
                            kDensityBoundsVestanMin);
    }
    else {
        mass = 0.0;
        density = 0.0;
    }
}

void ChooseNonjovianPostStellarMidOuterMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[4] = {0.27, 0.49, 0.62, 0.78};
    //double newProbs[4];
    //AdjustJovianProbs(probs, 4, 3, jovianProb, newProbs);
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTerranRange + kMassBoundsTerranMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsTerranRange + 
                            kDensityBoundsTerranMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[2]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsHadeanRange + kMassBoundsHadeanMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsHadeanRange + 
                            kDensityBoundsHadeanMin);
    }
    else if (planetTypeDieRoll < probs[3]) {
        mass = 0.0;
        density = 0.0;
    }
    else {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSubjovianRange + 
                         kMassBoundsSubjovianMin);
        if (mass < 17.0) {
            density = pow(10.0, densityDieRoll * kDensityBoundsLesserSubjovianRange + 
                                kDensityBoundsLesserSubjovianMin);
        }
        else {
            density = pow(10.0, densityDieRoll * kDensityBoundsSubjovianRange + 
                                kDensityBoundsSubjovianMin);
        }
    }
}

void ChooseNonjovianPostStellarDeepMass(unsigned age, double& mass, double& density)
{
    double planetTypeDieRoll = BrentsUnitRand();
    double massDieRoll = BrentsUnitRand();
    double densityDieRoll = BrentsUnitRand();
    double probs[2] = {0.66, 0.84};
    if (planetTypeDieRoll < probs[0]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsSelenianRange + kMassBoundsSelenianMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsSelenianRange + 
                            kDensityBoundsSelenianMin) + 100.0 * mass;
    }
    else if (planetTypeDieRoll < probs[1]) {
        mass = pow(10.0, BrentsUnitRand() * kMassBoundsTerranRange + kMassBoundsTerranMin);
        density = pow(10.0, densityDieRoll * kDensityBoundsTerranRange + 
                            kDensityBoundsTerranMin) + 100.0 * mass;
    }
    else {
        mass = 0.0;
        density = 0.0;
    }
}
