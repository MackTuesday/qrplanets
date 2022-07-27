// TODO: Do the whole life thing.

// TODO: Consider melting and boiling points at tropopause as well as
// surface. Use tropopausePressure constant unless surface pressure is less.
// Consider entire range when deciding on presence of clouds.
// If component is clouds up high but gas below, it precipitates without 
// particles reaching the surface. Use blackbody equilibrium temperature
// for tropopause temperature.

// TODO: Model daily and seasonal changes in temperature.

// TODO: It doesn't get a frost layer or hydrosphere if it doesn't have a solid
// surface. Maybe we say this is true when the pressure is greater than 1000 bars
// or so.

#include <math.h>

#include <vector>

#include "BrentsRand.h"
#include "ChooseMass.h"
#include "GlobalConstants.h"
#include "StarSystem.h"


StarSystem::StarSystem(unsigned seed, bool justTheStar)
{
    BrentsRandSet(seed);
    mStar = new Star;
    if (!justTheStar) {
        mSystemHasPrimeJovian = (BrentsUnitRand() <= 0.25);
        mPrimeJovianOrbitNumber = 0xffffffff;
        CreatePlanets();
        for (unsigned i = 0; i < mPlanets.size(); i++) {
            mPlanets[i]->ComputeChemistry();
        }
    }
}


void StarSystem::ChooseBaseAU(double rangeFactor, double& saInAUs, unsigned& theZone)
{
    unsigned zones[8] = {kPlZonesEpi, kPlZonesEpi, kPlZonesInn, kPlZonesInn,
                         kPlZonesMid, kPlZonesMid, kPlZonesMid, kPlZonesOut};
    unsigned zoneDieRoll = BrentsRand(8);
    theZone = zones[zoneDieRoll];
    double distDieRoll = BrentsUnitRand();
    saInAUs = 0.0;
    switch (theZone) {
        case kPlZonesEpi:
        if (BrentsUnitRand() < 0.2) {
            saInAUs = 0.01 * ((1.5*distDieRoll + 1.0)*(1.5*distDieRoll + 2.0) * 0.5 + 1.5);
        }
        else {
            // r = 1  defines the value of the kEpistellarBound variable
            saInAUs = 0.01 * ((3.5*distDieRoll + 2.0)*(3.5*distDieRoll + 3.0) * 0.5 + 0.5) -
                      0.00075;
        }
        break;
        
        case kPlZonesInn:
            saInAUs = kEpistellarBound + 2.178 * distDieRoll*distDieRoll;
        break;
            
        case kPlZonesMid:
            saInAUs = kInnerBound + 4.0 * distDieRoll;
        break;
            
        default:
            saInAUs = kMiddleBound + kOuterBoundCoefficient * distDieRoll * rangeFactor;
    }
}


void StarSystem::MakePrimeJovian(double& jovianMass, double& density)
{
    jovianMass = -1.0;
    if (mStar->mLuminosityClass == "III") {
        jovianMass = pow(10.0, BrentsUnitRand() * 0.5 + BrentsUnitRand() * 1.5 - 0.75) * 
                     kMassOfJupiter / kMassOfEarth;
    }
    else {
        char spectralClass = mStar->mSpectralClass[0];
        switch (spectralClass) {
            case 'B':
            jovianMass = pow(10.0, BrentsUnitRand() * 1.5 - 0.25) * kMassOfJupiter / kMassOfEarth;
            break;
                
            case 'A':
            case 'F':
            case 'G':
            case 'K':
            jovianMass = pow(10.0, BrentsUnitRand() * 0.5 + BrentsUnitRand() * 1.5 - 0.75) * 
                         kMassOfJupiter / kMassOfEarth;
            break;
                
            case 'M':
            jovianMass = pow(10.0, BrentsUnitRand() * 0.5 + BrentsUnitRand() - 0.75) * 
                         kMassOfJupiter / kMassOfEarth;
            break;
                
            default:
            jovianMass = pow(10.0, BrentsUnitRand() * 0.5 + BrentsUnitRand() * 1.5 - 0.75) * 
                         kMassOfJupiter / kMassOfEarth;
        }
    }

    density = 0.0;
    if (jovianMass < 17.0) {
        density = pow(10.0, BrentsUnitRand() * 2.73 + 0.48);
    }
    else if (jovianMass < 160.0) {
        density = pow(10.0, BrentsUnitRand() * 2.40 + 0.48);
    }
    else {
        density = pow(jovianMass / 318.0 * (BrentsUnitRand() * 0.2 + 0.9), 1.1) * 1330.0;
    }
}


void MarkForDeletion(unsigned age, double& mass, double& density)
{
    mass = -1.0;
    density = -1.0;
}


bool StarSystem::PlanetsTooClose(double lesserDist, double greaterDist,
                                 double closerMass, double furtherMass, double primaryMass)
{
    // From a paper by one Gladman
    // For sufficient separation we require (where a1 > a0)
    // mu1 = m1/m3 and mu2 = m2/m3; m1 is the closer mass, m2 is the further mass,
    // and m3 is the primary mass.
    // q = 3^1/6 (mu1 + mu2)^1/3
    // (a1 - a0)/a0 > 2q + 2qq - (11 mu1 + 7 mu2) / (3^(5/3) q)
    if (lesserDist == 0.0 || primaryMass == 0.0) {
        return false;
    }

    double mu1 = closerMass / primaryMass;
    double mu2 = furtherMass / primaryMass;
    
    double q = 1.2009370 * pow(mu1 + mu2, 1.0/3.0); 
    double relativeDifference = (greaterDist - lesserDist) / lesserDist;
    return relativeDifference <= 2.0*q * (1.0 + q) - (11.0*mu1 + 7.0*mu2) / (6.2402515*q);
}


void StarSystem::CreatePlanets()
{
    vector<double> semimajorAxes;
    double plRangeFactor = pow(mStar->mMass/kMassOfSol, 1.0/3.0);
    double outerBound = kMiddleBound + kOuterBoundCoefficient * plRangeFactor;
    double maxPlanetDistanceAUs = 30.0 * plRangeFactor;

    // Compute distances
    double a = 0.0;
    if (mSystemHasPrimeJovian) {
        double jovianSAInAUs = 0.0;
        unsigned jovianZone = 0xffffffff;
        ChooseBaseAU(plRangeFactor, jovianSAInAUs, jovianZone);
        unsigned numBeforeJovian = 0xffffffff;
        unsigned jovianRoll = BrentsRand(4);
        switch (jovianZone) {
            case kPlZonesEpi:
            case kPlZonesInn:
            numBeforeJovian = (jovianRoll > 2) ? (jovianRoll-2) : 0;
            break;
                
            case kPlZonesMid:
            case kPlZonesOut:
            default:
            numBeforeJovian = jovianRoll + 1;
            break;
        }

        double lastPosition = 0.0;
        // This is a way of choosing n values of uniform distribution
        // in guaranteed increasing order.
        // We might throw away some of these later in the function
        // depending on their proximity to each other. If so we don't do
        // it here because the decision depends on their relative masses
        // and we aren't doing that part here.
        for (unsigned i = 0; i < numBeforeJovian; i++) {
            double position = 
               jovianSAInAUs - (jovianSAInAUs - lastPosition) *
                               pow(BrentsUnitRand(), 1.0/(numBeforeJovian-i));
            semimajorAxes.push_back(position);
            lastPosition = position;
        }
        mPrimeJovianOrbitNumber = numBeforeJovian;
        semimajorAxes.push_back(jovianSAInAUs);

        a = jovianSAInAUs * (1.6 + 0.8 * BrentsUnitRand());
    }
    else {
        unsigned dummy = 0;
        ChooseBaseAU(plRangeFactor, a, dummy);
    }

    while (a < maxPlanetDistanceAUs) {
        semimajorAxes.push_back(a);
        a *= 1.6 + 0.8 * BrentsUnitRand();
    }

    // Figure out temperature
    vector<double> temperatures;
    for (unsigned i = 0; i < semimajorAxes.size(); i++) {
        temperatures.push_back(Planet::ComputeTemp(semimajorAxes[i], 1.0, mStar->mLuminosity, 
                                                   0.0, 0.0, 0.0));
    }

    // Now compute masses
    vector<double> masses;
    vector<double> densities;
    double tempZones[4] = {430, 280, 160, 80};
    double planetTypeDieRoll = BrentsUnitRand();
    void (*lastMassFunc)(unsigned, double&, double&) = NULL; // We use this later
    double m = 0.0;
    double d = 0.0;
    if (mSystemHasPrimeJovian) {
        MakePrimeJovian(m, d);
        for (unsigned i = 0; i < mPrimeJovianOrbitNumber; i++) {
            masses.push_back(0.0);
            densities.push_back(0.0);
        }
        masses.push_back(m);
        densities.push_back(d);

        if (mStar->mSpectralClass == "D" || mStar->mLuminosityClass == "III") {
            for (unsigned i = 0; i < semimajorAxes.size(); i++) {
/*              if (i == mPrimeJovianOrbitNumber) {
                    MarkForDeletion(mStar->mAge, m, d);
                    masses[i] = m;
                    densities[i] = d;
                }
*/              if (semimajorAxes[i] < kInnerBound) {
                    lastMassFunc = MarkForDeletion;
                }
                else if (semimajorAxes[i] < outerBound) {
                    lastMassFunc = ChooseJovianPostStellarMidOuterMass;
                }
                else {
                    lastMassFunc = ChooseJovianPostStellarDeepMass;
                }
                lastMassFunc(mStar->mAge, m, d);
                if (i < masses.size()) {
                    masses[i] = m;
                    densities[i] = d;
                }
                else {
                    masses.push_back(m);
                    densities.push_back(d);
                }
            }
        }
        else {
            for (unsigned i = 0; i < semimajorAxes.size(); i++) {
                if (i == mPrimeJovianOrbitNumber) {
                    continue;
                }
                if (temperatures[i] > 430) {
                    if (semimajorAxes[i] < kEpistellarBound) {
                        lastMassFunc = ChooseJovianNormalEpistellarMass;
                    }
                    else {
                        lastMassFunc = ChooseJovianNormalSupertorridMass;
                    }
                }
                else if (temperatures[i] > 280) {
                    lastMassFunc = ChooseJovianNormalTorridMass;
                }
                else if (temperatures[i] > 160) {
                    lastMassFunc = ChooseJovianNormalTemperateMass;
                }
                else if (temperatures[i] > 80) {
                    lastMassFunc = ChooseJovianNormalFrigidMass;
                }
                else {
                    lastMassFunc = ChooseJovianNormalSuperfrigidMass;
                }
                lastMassFunc(mStar->mAge, m, d);
                if (i < masses.size()) {
                    masses[i] = m;
                    densities[i] = d;
                }
                else {
                    masses.push_back(m);
                    densities.push_back(d);
                }
            }
        }
    }
    else {
        if (mStar->mSpectralClass == "D" || mStar->mLuminosityClass == "III") {
            for (unsigned i = 0; i < semimajorAxes.size(); i++) {
                if (semimajorAxes[i] < kInnerBound) {
                    lastMassFunc = MarkForDeletion;
                }
                else if (semimajorAxes[i] < outerBound) {
                    lastMassFunc = ChooseNonjovianPostStellarMidOuterMass;
                }
                else {
                    lastMassFunc = ChooseNonjovianPostStellarDeepMass;
                }
                lastMassFunc(mStar->mAge, m, d);
                if (i < masses.size()) {
                    masses[i] = m;
                    densities[i] = d;
                }
                else {
                    masses.push_back(m);
                    densities.push_back(d);
                }
            }
        }
        else {
            for (unsigned i = 0; i < semimajorAxes.size(); i++) {
                if (temperatures[i] > 160) {
                    lastMassFunc = ChooseNonjovianNormalSupToTempMass;
                }
                else if (temperatures[i] > 80) {
                    lastMassFunc = ChooseNonjovianNormalFrigidMass;
                }
                else {
                    lastMassFunc = ChooseNonjovianNormalSuperfrigidMass;
                }
                lastMassFunc(mStar->mAge, m, d);
                if (i < masses.size()) {
                    masses[i] = m;
                    densities[i] = d;
                }
                else {
                    masses.push_back(m);
                    densities.push_back(d);
                }
            }
        }
    }

    // Last planet can't be an asteroid belt
    while (masses[masses.size()-1] == 0.0) {
        lastMassFunc(mStar->mAge, m, d);
        masses[masses.size()-1] = m;
        densities[masses.size()-1] = d;
    }
        
    // Clear out planets marked for deletion
    unsigned i = 0;
    while (i < masses.size()) {
        if (masses[i] < 0.0) {
            semimajorAxes.erase(semimajorAxes.begin()+i);
            masses.erase(masses.begin()+i);
            densities.erase(densities.begin()+i);
            temperatures.erase(temperatures.begin()+i);
        }
        else {
            i++;
        }
    }
        
    // Throw away planets too close to more massive ones. This will
    // happen only within the orbit of a prime jovian because of the
    // way we're generating planets in that case.
    if (mSystemHasPrimeJovian) {
        i = 1;
        while (i <= mPrimeJovianOrbitNumber) {
            // If it's too close to its neighbors, throw it away.
            // Notice that either i increments or primeJovianOrbitNumber
            // decrements, so this loop is guaranteed to terminate.
            if (PlanetsTooClose(semimajorAxes[i-1], semimajorAxes[i],
                                masses[i-1]*kMassOfEarth, masses[i]*kMassOfEarth, mStar->mMass)) {
                unsigned lesserMassIdx = (masses[i-1] < masses[i]) ? (i-1) : i;
                semimajorAxes.erase(semimajorAxes.begin()+lesserMassIdx);
                masses.erase(masses.begin()+lesserMassIdx);
                temperatures.erase(temperatures.begin()+lesserMassIdx);
                densities.erase(densities.begin()+lesserMassIdx);
                mPrimeJovianOrbitNumber--;
            }
            else {
                i++;
            }
        }
    }

    vector<double> magRetentions;
    for (unsigned i = 0; i < masses.size(); i++) {
        // Totally ad hoc function chosen for its properties
        magRetentions.push_back(tanh(log(masses[i] * 1.22 + 1.0e-300)) * 0.5 + 0.5);
    }
        
    vector<double> planets;
    for (unsigned i = 0; i < semimajorAxes.size(); i++) {
        Planet* planet = new Planet(mStar, i, semimajorAxes[i], masses[i] * kMassOfEarth, 
                                    densities[i], magRetentions[i], temperatures[i]);
        mPlanets.push_back(planet);
    }
}
