#ifndef GlobalConstants_h
#define GlobalConstants_h


#include <string>

using namespace std;


#define kNumComponents       53
#define kNumBioComponents    17
#define kNumNonBioComponents 36
#define kNumDecompositions    5


enum {
    kCompFieldMolWeight,
    kCompFieldMolFormHeat,
    kCompFieldMassFormHeat,
    kCompFieldMeltPt,
    kCompFieldBoilPt,
    kCompFieldTripPtTemp,
    kCompFieldTripPtPress,
    kCompFieldCritPtTemp,
    kCompFieldCritPtPress,
    kCompFieldAbundance,
    kCompFieldGreenhouse,
    kNumComponentFields
};


extern double kTheComponents[kNumComponents][kNumComponentFields];
extern string kTheComponentNames[kNumComponents];
extern unsigned kNonBioComponents[kNumNonBioComponents];
extern unsigned kBioComponents[kNumBioComponents];
extern double kDecompositions[kNumDecompositions][2];


/*
var componentHTMLNames = [
    "H<sub>2</sub>",                   "He",           
    "H<sub>2</sub>O",                  "N<sub>2</sub>",         
    "CH<sub>4</sub>",                  "NH<sub>3</sub>",          
    "CO<sub>2</sub>",                  "C<sub>2</sub>H<sub>6</sub>",           
    "Ne",                              "CH<sub>3</sub>OH",         
    "Ar",                              "H<sub>2</sub>S", 
    "CO",                              "SO<sub>x</sub>",    
    "CHOOH",                           "CH<sub>2</sub>O",         
    "SiH<sub>4</sub>",                 "CH<sub>2</sub>(OH)<sub>2</sub>",      
    "C<sub>2</sub>H<sub>2</sub>O<sub>4</sub>",      
    "CH<sub>3</sub>COOH",              
    "CH<sub>3</sub>NH<sub>2</sub>",    "HF",                              
    "Si(OH)<sub>x</sub>",              "(Na,K)OH",                        
    "CH<sub>3</sub>SH",                "HCl",                             
    "H<sub>2</sub>SO<sub>4</sub>",     "CH<sub>3</sub>CH<sub>2</sub>SH",      
    "(Fe,Ni)(CO)<sub>x</sub>",         "CH<sub>x</sub>Cl<sub>4-x</sub>",   
    "HNO<sub>3</sub>",                 "H<sub>3</sub>PO<sub>4</sub>",  
    "SiH<sub>x</sub>Cl<sub>4-x</sub>", "C<sub>2</sub>H<sub>x</sub>Cl<sub>6-x</sub>",    
    "S<sub>x</sub>",                   "P<sub>x</sub>",       
    "P<sub>2</sub>H<sub>4</sub>",      "O<sub>2</sub>",               
    "Cl<sub>2</sub>",                  "F<sub>2</sub>",               
    "H<sub>2</sub>O<sub>2</sub>",      "C<sub>2</sub>H<sub>4</sub>",            
    "C<sub>2</sub>H<sub>2</sub>",      "N<sub>2</sub>H<sub>4</sub>",         
    "NO<sub>2</sub>",                  "NO",      
    "N<sub>2</sub>O",                  "HCN",  
    "HNC",                             "COS",  
    "COCl<sub>2</sub>",                "Cl<sub>2</sub>O",
    "OF<sub>2</sub>"
];
*/


// 1 - exp ma = 0.3
// exp ma = 0.7
// ma = ln 0.7
// a = ln 0.7 / 1e7
const double kPI = 3.14159265358979;
const double kSBConstant = 5.67e-8;       // W/mmKKKK
const double kGravConstant = 6.67e-11;    // mmm/kgss
const double kGasConstant = 8.314;        // J/molK
/*const kRetentionScalar = 200.0;         // ad hoc */
const double kVolatilesMassScalar = 0.002; // ad hoc
/*const kPrecCoverageScalar = 0.03;       // ad hoc
const kCloudAlbedoScalar = 0.22;        // ad hoc
const kFrostAlbedoScalar = 0.35;        // ad hoc */
const double kGreenhouseEffectScalar = 1.0e-6; // ad hoc 
const double kAtmReflectCoeff = -1e-5;    // ad hoc
const double kMassOfEarth = 5.976e+24;    // kg
/*const radiusOfEarth = 6.378e+6;       // m   */
const double kMassOfJupiter = 1.9e+27;    // kg
const double kMassOfSol = 1.99e+30;       // kg
const double kAU = 1.496e+11;             // m
const double kLuminosityOfSol = 3.8e+26;  // W
/*const tropopausePressure = 0.1;        // bar -- amazingly, we find roughly this pressure at
                                       // at the tropopauses of Venus, Earth, Jupiter,
                                       // Saturn, Uranus, and Neptune.
*/
#define kEpistellarBound       0.1555  // AUs
#define kInnerBound            2.361
#define kMiddleBound           6.361
#define kOuterBoundCoefficient 6.0

enum {
    kPlZonesEpi = 0,
    kPlZonesInn = 1,
    kPlZonesMid = 2,
    kPlZonesOut = 3,
    kNumPlZones
};
/*
//Mass Ranges
//===========
//Subjovian:   0.02-  0.2 Jupiters
//Jovian:      0.2 -  5   Jupiters
//Transjovian: 5   - 13   Jupiters
//Selenian:    0.01-  0.1  Earths
//Glacian:     0.01-  0.3  Earths
//Terran:      0.3 -  1.0  Earths
//Nerean:      0.3 -  1.0  Earths
//Aquarian:    2   -  4    Earths
//Hadean:      2   -  4    Earths
//Vestan:      7   - 14    Earths
*/

enum {
    kSystemAgesNewborn = 0,
    kSystemAgesYoung = 1,
    kSystemAgesMature = 2,
    kNumSystemAges
};

#define kFsgn(x) ((0.0 < (x)) - ((x) < 0.0))


#endif
