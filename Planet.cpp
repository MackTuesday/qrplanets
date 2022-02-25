#include <math.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "GlobalConstants.h"
#include "Planet.h"

using namespace std;


void DumpPlanetChem(Planet* planet, vector<double>& atm, vector<double>& hyd, 
                                    vector<double>& cld, vector<double>& fst)
{
    stringstream ss;
    ss << setprecision(3);
    ss << "Albedo:     " << planet->mAlbedo << '\n';
    ss << "Greenhouse: " << planet->mGreenhouse << '\n';
    ss << "Temp:       " << planet->mTemperature << '\n';
    ss << "Pressure:   " << planet->mPressure << '\n';
    ss << "MinMolMass: " << planet->mMinMolMass << '\n';
    ss << "Life:       " << planet->mLife << '\n';
    for (unsigned i = 0; i < atm.size(); i++) {
        if (atm[i] > 0.0 || hyd[i] > 0.0 ||
            cld[i] > 0.0 || fst[i] > 0.0) {
            ss << i << ' ' << kTheComponentNames[i] << ":   ";
            if (atm[i] > 0.0) {
                ss << " atm " << atm[i] << "  ";
            }
            if (hyd[i] > 0.0) {
                ss << " hyd " << hyd[i] << "  ";
            }
            if (cld[i] > 0.0) {
                ss << " cld " << cld[i] << "  ";
            }
            if (fst[i] > 0.0) {
                ss << " fst " << fst[i] << "  ";
            }
            ss << '\n';
        }
    }
    ss << '\n';
    string dump = ss.str();
    cout << dump;
/*
    double mAtmAbundance;
    double mHydAbundance;
    double mCldAbundance;
    double mFstAbundance;
    unsigned mComponentIndex;

    vector<PlanetComponent*> mVolComponents;
    double mDensity;
    double mAlbedo;
    double mGreenhouse;
    double mRadius;
    double mGravity;
    double mEscapeVelocity;
    double mSurfaceArea;
    double mMinMolMass;
    double mOrigVolMass;
    double mMagRetention;
    double mSemimajorAxis;
    double mMass;
    double mIntrinsicLum;
    double mTemperature;
    double mPressure;
    double mLife;
    unsigned mNumber;
*/
}


PlanetComponent::PlanetComponent(unsigned i, double a)
{
    mComponentIndex = i;
    mAtmAbundance = a;
    mHydAbundance = 0.0;
    mCldAbundance = 0.0;
    mFstAbundance = 0.0;
}


double Planet::ComputeTemp(double sa, double r, double stellarLum, double intrinsicLum,
                           double albedo, double greenhouse)
{
    double temp = sqrt(sqrt(6.0783e9 * stellarLum * (1.0 - albedo + greenhouse) / (sa*sa) + 
                           55.14 +   // CMB term
                            1.328e19 * intrinsicLum / (r*r)));
    return temp;
}


double Planet::ComputeMinMolMass(double temp, double escapeVelocity)
{
    // m = ln (t/2.8e8) 3RT / v^2
    // t = retention time, seconds
    // R = gas constant
    // v = escape velocity
        
    return 200.0 * 3.0 * kGasConstant * temp / (escapeVelocity * escapeVelocity);
}


Planet::Planet(Star* star, unsigned n, double sa, double m, double d, double mr, double temp)
{
    mStar = star;
    mNumber = n;
    mSemimajorAxis = sa;
    mMass = m;
    mDensity = d;
    mMagRetention = mr;
    mTemperature = temp;
    mAlbedo = 0.3;
    mGreenhouse = 0.0;
    mRadius = (m > 0.0 && d > 0.0) ? pow(m / d * 0.75 / kPI, 1.0/3.0) : 100.0;
    mGravity = kGravConstant * m / (mRadius * mRadius);
    mEscapeVelocity = sqrt(2.0 * kGravConstant * m / mRadius);
    mSurfaceArea = 4.0 * 3.14159 * mRadius * mRadius;
    mIntrinsicLum = ComputeIntrinsicLum();
    
//printf("%.3le %.3le %.3le %.3le %.3le %.3le %.3le %.3le\n", sa, m, d, mr, temp, mRadius, mGravity, mEscapeVelocity);
        
    mMinMolMass = ComputeMinMolMass(temp, mEscapeVelocity);
        
    double sumOfAbundances = 0.0;
    for (unsigned i = 0; i < kNumComponents; i++) {
        PlanetComponent* pc = new PlanetComponent(i, kTheComponents[i][kCompFieldAbundance]);
        // We don't play with the abundance of hydrogen or helium
        if (i > 1) {
            double abundanceAdjustment = pow(2.0, BrentsUnitRand() * 2.0 - 1.0);
            
            double adjustedAbundance = pc->mAtmAbundance * abundanceAdjustment;
            if (adjustedAbundance > 0.2) {
                adjustedAbundance = 0.25 * pow(0.8, BrentsUnitRand());
            }
            
            pc->mAtmAbundance = adjustedAbundance;
        }
        mVolComponents.push_back(pc);
        sumOfAbundances += pc->mAtmAbundance;
    }
        
    // Clear out bio components because we're about to set them properly
    for (unsigned i = 0; i < kNumBioComponents; i++) {
        unsigned bioIdx = kBioComponents[i];
        PlanetComponent* bioComp = mVolComponents[bioIdx];
        bioComp->mAtmAbundance = 0.0;
    }

    // Is there life?
    mLife = 0.0;
    // The probability we use here is misleading because life often gets
    // killed off due to atmospheric changes.
    if (BrentsUnitRand() < 0.33) {
        unsigned bioCompRoll = BrentsRand(kNumBioComponents);
        double abundance = BrentsUnitRand() * 0.2 + 0.05;
        mVolComponents[kBioComponents[bioCompRoll]]->mAtmAbundance = abundance;
        // Select a common component
        unsigned nonBioCompRoll = BrentsRand(kNumNonBioComponents);
        while (mVolComponents[kNonBioComponents[nonBioCompRoll]]->mAtmAbundance < 1e-4) {
            nonBioCompRoll = BrentsRand(kNumNonBioComponents);
        }
        mVolComponents[kNonBioComponents[nonBioCompRoll]]->mAtmAbundance *= 
          (BrentsUnitRand() + 0.618) * 0.001;
    }
        
    if (sumOfAbundances == 0.0) {
        sumOfAbundances = 1e-15;
    }
        
    // 1000 is the density of water
    mOrigVolMass = mMass * sumOfAbundances * mMagRetention * kVolatilesMassScalar;
    mPressure = mOrigVolMass * mGravity / mSurfaceArea;
}


void Planet::ComputeChemistry()
{
    for (unsigned i = 0; i < 10; i++) {
        ReactionStep();
    }

    HandleDecompositions();
    
    if (mPressure >= 1000.0) {
        mGreenhouse = 0.0;
        mTemperature = ComputeTemp(mSemimajorAxis, mRadius, mStar->mLuminosity, 
                                   mIntrinsicLum, mAlbedo, mGreenhouse);
        for (unsigned i = 0; i < mVolComponents.size(); i++) {
            mVolComponents[i]->mHydAbundance = 0.0;
            mVolComponents[i]->mFstAbundance = 0.0;
        }
    }

    mLife = ComputeLife();
}


double Planet::ComputeLife()
{
    double bioAbundance = 0.0;
/*    if (atmAbundances) {
        for (i = 0; i < bioComponents.length; i++) {
            var bioIdx = bioComponents[i];
            bioAbundance += atmAbundances[bioIdx] + hydAbundances[bioIdx] + 
                            cldAbundances[bioIdx] + fstAbundances[bioIdx];
        }
    }
    else { */
        for (unsigned i = 0; i < kNumBioComponents; i++) {
            unsigned bioIdx = kBioComponents[i];
            PlanetComponent* comp = mVolComponents[bioIdx];
            bioAbundance += comp->mAtmAbundance + comp->mHydAbundance + 
                            comp->mCldAbundance + comp->mFstAbundance;
        }
//    }

    // kg m/s^2/m^2 / m/s^2 = kg
    // bars * paPerBar / (g * m/s^2/g)
    // bars * paPerBar / (g * m/s^2/g)
    // bars * paPerBar / (g * accPerG)
    // bars * paPerBar / g * gPerAcc
    //var paPerBar = 100000.0;
    //var gPerAcc = 10.0;
    //var tier = Math.log(bioAbundance * (this.pressure * paPerBar / this.gravity * gPerAcc) * 
    //                    Math.E + 1.0) * 2.9 - 1; 
    //tier = (bioAbundance <= 0.0) ? 0.0 : Math.max(0, Math.min(5.9, tier));

    double tier = 0.0;
    if (bioAbundance > 0.0) {
        // Subtract one so minimum is 0
        // Divide by (Math.exp(3.0) - 1.0) to normalize to [0, 1)
        // Multiply by (5.9 - 1.0) so the interval has the right width
        // Add 1 so the minimum is 1
        tier = (exp(BrentsUnitRand() * 3.0) - 1.0) / (exp(3.0) - 1.0) * (5.9 - 1.0) + 1.0;
    }
    return tier;
}


double Planet::ComputeMeltingPoint(double* s, double p)
{
    double tabMP = s[kCompFieldMeltPt];
    double tpTemp = s[kCompFieldTripPtTemp];
    double tpPres = s[kCompFieldTripPtPress];
    double cpTemp = s[kCompFieldCritPtTemp];
    double cpPres = s[kCompFieldCritPtPress];
    
    if (p < tpPres) {
        return ComputeBoilingPoint(s, p);
    }
    else {
        double logRelPres = log(p / tpPres);
        // This curve is linear as x-->inf and constant as x--> -inf
        return (log(exp(logRelPres) + 1.0) - log(2.0)) * 0.005 + tpTemp;
    }
}


double Planet::ComputeBoilingPoint(double* s, double p)
{
    double tabBP = s[kCompFieldBoilPt];
    double tpTemp = s[kCompFieldTripPtTemp];
    double tpPres = s[kCompFieldTripPtPress];
    double cpTemp = s[kCompFieldCritPtTemp];
    double cpPres = s[kCompFieldCritPtPress];
    
    // if the pressure is greater than critical pressure,
    // return critical temperature
    if (p > cpPres) {
        return cpTemp;
    }
    
    if (p < tpPres) {
        if (p <= 1.0) {
            return tabBP * p;
        }
        else if (tpPres > 1.0) {
            return (tpTemp - tabBP) / (tpPres - 1.0) * (p - 1.0) + tabBP;
        }
        else {
            return (tabBP - tpTemp) / (1.0 - tpPres) * (1.0 - p) + tpTemp;
        }
    }
    
    // if the triple point pressure is greater than 1 bar, we don't
    // have separate MP and BP in our table. So interpolate between tp and cp.
    if (tpPres > 1.0) {
        double slope = log(cpPres / tpPres) / (cpTemp - tpTemp);
        double temp = log(p / tpPres) / slope + tpTemp;
        return temp;
    }
    
    // If pressure is greater than one bar, interpolate between table value and cp
    if (p > 1.0) {
        double slope = (log(cpPres) - log(1.0)) / (cpTemp - tabBP);
        double temp = (log(p) - log(1.0)) / slope + tabBP;
        return temp;
    }
    // Else interpolate between table and (0,0)
    else if (p < 1.0) {
        double slope = (log(1.0) - log(tpPres)) / (tabBP - tpTemp);
        double temp = -(log(1.0) - log(p)) / slope + tabBP;
        return temp;
    }
    else {
        return tabBP;
    }
}


void Planet::HandleDecompositions()
{
    double sumOfAbundances = 0.0;
    for (unsigned i = 0; i < mVolComponents.size(); i++) {
        sumOfAbundances += mVolComponents[i]->mAtmAbundance;
    }
    
    if (sumOfAbundances == 0.0) {
        return;
    }

    double decompAbundance = 0.0;
    for (unsigned i = 0; i < kNumDecompositions; i++) {
        double* decomp = kDecompositions[i];
        if (mTemperature > decomp[1]) {
            decompAbundance += mVolComponents[decomp[0]]->mAtmAbundance;
            mVolComponents[decomp[0]]->mAtmAbundance = 0.0;
            mVolComponents[decomp[0]]->mHydAbundance = 0.0;
            mVolComponents[decomp[0]]->mCldAbundance = 0.0;
            mVolComponents[decomp[0]]->mFstAbundance = 0.0;
        }
    }
    
    mPressure *= (sumOfAbundances - decompAbundance) / sumOfAbundances;
}


double Planet::ComputeIntrinsicLum()
{
    double intrinsicLum = (mMass <= 2.56968e28) ? 1.717e-10*mMass : 9.422e-57*pow(mMass, 2.63);
    intrinsicLum *= 0.8 + 0.45 * BrentsUnitRand();
    return intrinsicLum/kLuminosityOfSol;
}


void Planet::ReactionStep()
{
    // Compute pressure based on mass of atmosphere
    vector<double> atmComponents;
    vector<double> cldComponents;
    vector<double> hydComponents;
    vector<double> fstComponents;
    vector<double> atmHydAbundances;
    double sumAtmAbundances = 0.0;
    double atmMass = 0.0;
    double totalLiq = 0.0;
    for (unsigned i = 0; i < mVolComponents.size(); i++) {
        atmComponents.push_back(0.0);
        cldComponents.push_back(0.0);
        hydComponents.push_back(0.0);
        fstComponents.push_back(0.0);
        
        double atmHydAbundance = mVolComponents[i]->mAtmAbundance + 
                                 mVolComponents[i]->mHydAbundance;
        atmHydAbundances.push_back(atmHydAbundance);
        
        if (atmHydAbundance <= 0.0) {
            continue;
        }
        
        double meltingPoint = ComputeMeltingPoint(kTheComponents[i], mPressure);
        double boilingPoint = ComputeBoilingPoint(kTheComponents[i], mPressure);
        
        // Components with low enough BP are atmosphere
        double gasFuzzinessRange = tanh(BrentsUnitRand() * 2.0 + 0.2) * 0.5;
        double minGasTemp = boilingPoint * (1 - gasFuzzinessRange*0.5);
        // The gasFuzziness correction proportion goes into atmosphere, and if the
        // temperature falls into both gas and liquid ranges, the rest goes into
        // hydrosphere.
        double atmHydProportion = (mTemperature - minGasTemp) / (boilingPoint - minGasTemp);
        atmHydProportion = min(atmHydProportion, 1.0);
        if (minGasTemp < mTemperature) {
            atmComponents[i] = atmHydAbundances[i] * atmHydProportion;
            sumAtmAbundances += atmComponents[i];
            atmMass += atmComponents[i] * mOrigVolMass;
        }
        
        // Components with borderline BP are clouds
        double liqFuzzinessRange = tanh(BrentsUnitRand() * 2.0 + 0.2) * 0.5;
        double maxLiqTemp = boilingPoint * (1 + liqFuzzinessRange*0.5);
        cldComponents[i] = 0.0;
        if (boilingPoint < mTemperature && mTemperature < maxLiqTemp) {
            cldComponents[i] = (mTemperature - boilingPoint) / (maxLiqTemp - boilingPoint);
        }
        else if (minGasTemp < mTemperature && mTemperature <= boilingPoint) {
            cldComponents[i] = (boilingPoint - mTemperature) / (boilingPoint - minGasTemp);
        }
        
        // Components with low enough MP are hydrosphere
        double minLiqTemp = meltingPoint * (1 - liqFuzzinessRange*0.5);
        if (minLiqTemp < mTemperature && mTemperature < maxLiqTemp) {
            hydComponents[i] = atmHydAbundances[i] * (1.0 - atmHydProportion);
            totalLiq += hydComponents[i] * mOrigVolMass;
        }
        
        // Components with borderline MP are frost
        double maxSolTemp = meltingPoint * (1 + liqFuzzinessRange*0.5);
        fstComponents[i] = 0.0;
        if (meltingPoint < mTemperature && mTemperature < maxSolTemp) {
            fstComponents[i] = (mTemperature - meltingPoint) / (maxSolTemp - meltingPoint);
        }
        else if (mTemperature <= meltingPoint && atmHydAbundances[i] > 0.0) {
            fstComponents[i] = 1.0;
        }
    }
    
    double hydPerSqMeter = totalLiq / mSurfaceArea * 0.001; // 0.001 m^3/kg, density of water
    // A smattering of seas is easy to make. Conversely, a smattering of islands
    // is hard to cover up.
    mNormedHydDepth = hydPerSqMeter / kHydDepthNorm;
    double relToShelf = 2.0 * mNormedHydDepth - kShelfParam;
    double twoMinusShelf = 2.0 - kShelfParam;
    // Desmos equation
    // y=\left(\left(\operatorname{abs}(2x-a)\right)^{b}\cdot\operatorname{sgn}(2x-a)+a^{b}+cx\ \right)/\ \left(\left(\operatorname{abs}(2-a)\right)^{b}\cdot\operatorname{sgn}(2-a)+a^{b}+c\right)
    mHydCover = (pow(fabs(relToShelf), kTerrainPointiness) * kFsgn(relToShelf) +
                 pow(kShelfParam, kTerrainPointiness) + kShelfSlope * mNormedHydDepth) /
                (pow(fabs(twoMinusShelf), kTerrainPointiness) * kFsgn(twoMinusShelf) +
                 pow(kShelfParam, kTerrainPointiness) + kShelfSlope);
//    mHydCover = pow(0.5 + 0.5*tanh(log(mNormedHydDepth+1e-12)/log(10.0)), 6.0);
//    mHydCover = (mHydCover < 0.1) ? max(mNormedHydDepth, mHydCover) : mHydCover;
    mHydCover = (mHydCover < 0.0) ? 0.0 : mHydCover;
    mHydCover = (mHydCover > 1.0) ? 1.0 : mHydCover;

    // Estimate albedo based on clouds and frost
    double sqrSumCloud = 0.0;
    double sqrSumFrost = 0.0;
    double coverage = 0.0;
    double atmMassPerArea = (mSurfaceArea > 0.0) ? atmMass/mSurfaceArea : 0.0;
    for (unsigned i = 0; i < mVolComponents.size(); i++) {
        // The following is totally hacked
        coverage = pow(cldComponents[i] * atmMassPerArea * 0.003 * 
                       (kTheComponents[i])[kCompFieldAbundance], 3.0);
        coverage = std::min(coverage, 1.0);
        sqrSumCloud += coverage * coverage;
        coverage = pow(fstComponents[i] * atmMassPerArea * 0.003 * 
                       (kTheComponents[i])[kCompFieldAbundance], 3.0);
        coverage = std::min(coverage, 1.0);
        sqrSumFrost += coverage * coverage;
    }
    sqrSumCloud = std::min(sqrSumCloud, 1.0);
    sqrSumFrost = std::min(sqrSumFrost, 1.0);
    mCloudCover = sqrt(sqrSumCloud);
    mFrostCover = sqrt(sqrSumFrost);
    double sqrAtmReflect = pow(1.0 - exp(atmMassPerArea*kAtmReflectCoeff), 2.0) * 0.09;
    mAlbedo = tanh((0.1 + BrentsUnitRand() * BrentsUnitRand() * 0.1 +
                   sqrt(sqrAtmReflect + sqrSumCloud + sqrSumFrost) * 0.75) / 0.85) * 0.85;

    // Estimate greenhouse effect
    double ghSqrSum = 0.0;
//printf("Greenhouse: ");
    if (sumAtmAbundances > 0.0) {
        for (unsigned i = 0; i < kNumComponents; i++) {
            double* comp = kTheComponents[i];
            
            // If is a greenhouse gas
            if (comp[kCompFieldGreenhouse] > 0.0) {
                double partialMass = atmMassPerArea * atmComponents[i] / sumAtmAbundances;
//printf("%u %.3lf ", i, partialMass);
                ghSqrSum += partialMass * partialMass;
            }
        }
    }
//printf("\n");

    // Adjust temperature based on albedo and greenhouse gases
    mGreenhouse = sqrt(ghSqrSum) * kGreenhouseEffectScalar;
    double surfaceTemp = ComputeTemp(mSemimajorAxis, mRadius, mStar->mLuminosity, 
                                     mIntrinsicLum, mAlbedo, mGreenhouse);
    
    // Cull lightest atmospheric components based on escape velocity
    mMinMolMass = ComputeMinMolMass(surfaceTemp, mEscapeVelocity);
    for (unsigned i = 0; i < kNumComponents; i++) {
        double* comp = kTheComponents[i];
        // table molar masses are in g/mol
        if (comp[kCompFieldMolWeight] < mMinMolMass * 1000.0) {
            //atmComponents[i] *= 1.0 - gasFuzzinessCorrections[i];
            atmComponents[i] = 0.0;
        }
    }

    // Store the numbers we computed
    for (unsigned i = 0; i < kNumComponents; i++) {
        mVolComponents[i]->mAtmAbundance = atmComponents[i];
        mVolComponents[i]->mHydAbundance = hydComponents[i];
        mVolComponents[i]->mCldAbundance = cldComponents[i];
        mVolComponents[i]->mFstAbundance = fstComponents[i];
    }
    mPressure = atmMass > 0.0 ? (atmMass * mGravity / mSurfaceArea * 1e-5) : 0.0;
    mTemperature = surfaceTemp;

    //DumpPlanetChem(this, atmComponents, hydComponents, cldComponents, fstComponents);
}
