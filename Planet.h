#include <vector>

#include "BrentsRand.h"
#include "Star.h"


#define kHydDepthNorm         3.0e+04  // meters
#define kShelfParam           1.0      // offset (dimensionless)
#define kShelfSlope           0.75     // slope (dimensionless)
#define kTerrainPointiness    4.0      // exponent


struct PlanetComponent
{
    double mAtmAbundance;
    double mHydAbundance;
    double mCldAbundance;
    double mFstAbundance;
    unsigned mComponentIndex;
    PlanetComponent(unsigned i, double a);
};


struct Planet
{
    static double ComputeTemp(double sa, double r, double stellarLum, double intrinsicLum, 
                              double albedo, double greenhouse);
    static double ComputeMinMolMass(double temp, double escapeVelocity);
    
    Planet(Star* star, unsigned number, double semimajorAxis, double mass, double density,
           double magRetention, double temp);
    void ComputeChemistry();
    void ReactionStep();
    double ComputeMeltingPoint(double* s, double p);
    double ComputeBoilingPoint(double* s, double p);
    void HandleDecompositions();
    double ComputeLife();
    Star* mStar;
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
    double mHydCover;
    double mCloudCover;
    double mFrostCover;
    double mLife;
    double mNormedHydDepth;
    unsigned mNumber;

private:
    double ComputeIntrinsicLum();
};
