#include <vector>

#include "Planet.h"
#include "Star.h"


struct StarSystem
{
    StarSystem(unsigned seed=423732843, bool justTheStar=false);
    void CreatePlanets();
    void MakePrimeJovian(double& mass, double& density);
    double ComputeTemp(double sa, double stellarLum);
    void ChooseBaseAU(double rangeFactor, double& saInAUs, unsigned& theZone);
    bool PlanetsTooClose(double lesserDist, double greaterDist,
                         double closerMass, double furtherMass, double primaryMass);

    Star* mStar;
    vector<Planet*> mPlanets;

    bool mSystemHasPrimeJovian;
    unsigned mPrimeJovianOrbitNumber;
};
