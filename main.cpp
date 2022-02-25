// Yes, I know I haven't done things in the most elegant way. I got it all to
// work and have moved on to more important things.

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "BrentsRand.h"
#include "GlobalConstants.h"
#include "StarSystem.h"

using namespace std;


void init();


/*
function consoleReport(starSystem) {
    var htmlString = "";
    htmlString += "Star:\n";
    htmlString += "Type " + starSystem.star.spectralClass + starSystem.star.luminosityClass + "\n";
    htmlString += "Luminosity " + starSystem.star.luminosity.toPrecision(3) + 
                  " times the sun\n";
    htmlString += "Planets: semimajor axis, mass, temperature\n";
    for (var i = 0; i < starSystem.planets.length; i++) {
        htmlString += starSystem.planets[i].semimajorAxis.toPrecision(3) + " " + 
                      starSystem.planets[i].mass.toPrecision(3) + " " +
                      starSystem.planets[i].temperature.toPrecision(3) + " ";
        if (i == starSystem.primeJovianOrbitNumber) {
            htmlString += "* ";
        }
        htmlString += "\n";
    }
    return htmlString;
}
*/

string ConsoleReport(StarSystem* starSystem) {
    stringstream ss;
    ss << setprecision(3);
    ss << "Star:\n";
    ss << "Type " << starSystem->mStar->mSpectralClass << starSystem->mStar->mLuminosityClass <<
          "\n";
    ss << "Luminosity " << starSystem->mStar->mLuminosity << " times Sol\n";
    ss << "Planets: semimajor axis, mass, temperature\n";
    for (unsigned i = 0; i < starSystem->mPlanets.size(); i++) {
        ss << starSystem->mPlanets[i]->mSemimajorAxis << " ";
        ss << starSystem->mPlanets[i]->mMass/kMassOfEarth << " ";
        ss << starSystem->mPlanets[i]->mTemperature << " ";
        if (i == starSystem->mPrimeJovianOrbitNumber) {
            ss << "* ";
        }
        ss << "\n";
    }
    string report = ss.str();
    return report;
}

/*
function getPrincipalAtmComponents(components) {
    var sortedComponents = components.sort( 
        function(a,b) { return b.atmAbundance - a.atmAbundance; } 
                                          );
    var sumOfAbundances = 0.0;
    for (var i = 0; i < sortedComponents.length; i++) {
        sumOfAbundances += sortedComponents[i].atmAbundance;
    }
    
    if (sumOfAbundances == 0.0) {
        return [];
    }
    
    var indices = [];
    var sumNormalized = 0.0;
    var idx = 0;
    while (sumNormalized < 0.99 && idx < sortedComponents.length && idx < 4) {
        sumNormalized += sortedComponents[idx].atmAbundance / sumOfAbundances;
        indices.push(sortedComponents[idx].componentIndex);
        idx++;
    }
    
    return indices;
}

function getPrincipalHydComponents(components) {
    var sortedComponents = components.sort( 
        function(a,b) { return b.hydAbundance - a.hydAbundance; } 
                                          );
    var sumOfAbundances = 0.0;
    for (var i = 0; i < sortedComponents.length; i++) {
        sumOfAbundances += sortedComponents[i].hydAbundance;
    }
    
    if (sumOfAbundances == 0.0) {
        return [];
    }
    
    var indices = [];
    var sumNormalized = 0.0;
    var idx = 0;
    while (sumNormalized < 0.99 && idx < sortedComponents.length && idx < 4) {
        sumNormalized += sortedComponents[idx].hydAbundance / sumOfAbundances;
        indices.push(sortedComponents[idx].componentIndex);
        idx++;
    }
    
    return indices;
}

function getPrincipalCldComponents(components) {
    var sortedComponents = components.sort( 
        function(a,b) { return b.cldAbundance - a.cldAbundance; } 
                                          );
    var indices = [];
    for (var i = 0; i < 3 && i < sortedComponents.length; i++) {
        if (sortedComponents[i].cldAbundance > 1e-6) {
            indices.push(sortedComponents[i].componentIndex);
        }
    }
    
    return indices;
}

function getPrincipalFstComponents(components) {
    var sortedComponents = components.sort( 
        function(a,b) { return b.fstAbundance - a.fstAbundance; } 
                                          );
    var indices = [];
    for (var i = 0; i < 3 && i < sortedComponents.length; i++) {
        if (sortedComponents[i].fstAbundance > 1e-6) {
            indices.push(sortedComponents[i].componentIndex);
        }
    }
    
    return indices;
}

function planetsTable(starSystem) {
    var star = starSystem.star;
    var planets = starSystem.planets;
    
    var htmlString = "";
//    htmlString += "<p>";
    htmlString += "Star:<br />";
    htmlString += "Type " + starSystem.star.spectralClass + starSystem.star.luminosityClass + "<br />";
    htmlString += "Luminosity " + starSystem.star.luminosity.toPrecision(3);
//    htmlString += "</p>";
    htmlString += "<br />";

    // semimajorAxis, mass, radius, temperature, pressure, atm, clouds, hyd, frosts
    htmlString += "<table>";
    
    htmlString += "<tr>";
    htmlString += "<td><b>" + "Dist"  + "</b></td>";
    htmlString += "<td><b>" + "Mass"  + "</b></td>";
    htmlString += "<td><b>" + "Rad"   + "</b></td>";
    htmlString += "<td><b>" + "Dens"  + "</b></td>";
    htmlString += "<td><b>" + "Temp"  + "</b></td>";
    htmlString += "<td><b>" + "Press" + "</b></td>";
    htmlString += "<td><b>" + "Alb"   + "</b></td>";
    htmlString += "<td><b>" + "Atm"   + "</b></td>";
    htmlString += "<td><b>" + "Cloud" + "</b></td>";
    htmlString += "<td><b>" + "Hyd"   + "</b></td>";
    htmlString += "<td><b>" + "Frost" + "</b></td>";
    htmlString += "<td><b>" + "Life" + "</b></td>";
    htmlString += "</tr>";
        
    for (var i = 0; i < planets.length; i++) {
        var atmIndicesList = getPrincipalAtmComponents(starSystem.planets[i].volComponents);
        var cldIndicesList = getPrincipalCldComponents(starSystem.planets[i].volComponents);
        var hydIndicesList = getPrincipalHydComponents(starSystem.planets[i].volComponents);
        var fstIndicesList = getPrincipalFstComponents(starSystem.planets[i].volComponents);
        
        var atmHTML = "";
        if (atmIndicesList.length > 0) {
            atmHTML += componentHTMLNames[atmIndicesList[0]];
            for (var j = 1; j < atmIndicesList.length; j++) {
                atmHTML += ", " + componentHTMLNames[atmIndicesList[j]];
            }
        }
        
        var cldHTML = "";
        if (cldIndicesList.length > 0) {
            cldHTML += componentHTMLNames[cldIndicesList[0]];
            for (var j = 1; j < cldIndicesList.length; j++) {
                cldHTML += ", " + componentHTMLNames[cldIndicesList[j]];
            }
        }
        
        var hydHTML = "";
        if (hydIndicesList.length > 0) {
            hydHTML += componentHTMLNames[hydIndicesList[0]];
            for (var j = 1; j < hydIndicesList.length; j++) {
                hydHTML += ", " + componentHTMLNames[hydIndicesList[j]];
            }
        }
        
        var fstHTML = "";
        if (fstIndicesList.length > 0) {
            fstHTML += componentHTMLNames[fstIndicesList[0]];
            for (var j = 1; j < fstIndicesList.length; j++) {
                fstHTML += ", " + componentHTMLNames[fstIndicesList[j]];
            }
        }
        
        htmlString += "<tr>";
        htmlString += "<td>" + planets[i].semimajorAxis.toPrecision(3) + "</td>";
        htmlString += "<td>" + (planets[i].mass/massOfEarth).toPrecision(3) + "</td>";
        if (planets[i].mass > 0.0) {
            htmlString += "<td>" + (planets[i].radius/radiusOfEarth).toPrecision(3) + "</td>";
            htmlString += "<td>" + (planets[i].mass/(Math.pow(planets[i].radius,3.0)*
                                                    4/3*Math.PI)/1000).toPrecision(3) + "</td>";
        }
        else {
            htmlString += "<td>" + "n/a" + "</td>";
            htmlString += "<td>" + "n/a" + "</td>";
        }
        htmlString += "<td>" + planets[i].temperature.toPrecision(3) + "</td>";
        if (planets[i].pressure == 0.0) {
            htmlString += "<td>" + planets[i].pressure.toPrecision(2) + "</td>";
        } else if (planets[i].pressure < 0.001) {
            var pressureValue = "<td>" + planets[i].pressure.toExponential(2) + "</td>";
            var unbreakable = pressureValue.replace("-", "\u2011");
            htmlString += unbreakable;
        } else if (planets[i].pressure < 1000.0) {
            htmlString += "<td>" + planets[i].pressure.toPrecision(3) + "</td>";
        }
        else {
            htmlString += "<td>" + "&gt; 1000" + "</td>";
        }
        htmlString += "<td>" + planets[i].albedo.toPrecision(3)        + "</td>";
        htmlString += "<td>" + atmHTML + "</td>";
        htmlString += "<td>" + cldHTML + "</td>";
        htmlString += "<td>" + hydHTML + "</td>";
        htmlString += "<td>" + fstHTML + "</td>";
        if (planets[i].life <= 0) {
            htmlString += "<td>" + "&nbsp;" + "</td>";
        }
        else {
            htmlString += "<td>" + planets[i].life.toPrecision(2) + "</td>";
        }
        htmlString += "</tr>";
    }
    
    htmlString += "</table>";
    document.getElementById("planetsTableDiv").innerHTML += htmlString;
    return htmlString;
}

function init() {
    var starSystem = new StarSystem();
    console.log(consoleReport(starSystem));
    document.getElementById("planetsTableDiv").innerHTML = planetsTable(starSystem);
}

function makeNewSystemClicked() {
    init(false);
}
*/

void init()
{
    StarSystem* starSystem = new StarSystem();
    cout << ConsoleReport(starSystem);
}

/*
function theMainFunction() {
    document.getElementById("makeNewSystemButton").innerHTML += 
        "<button onclick=\"makeNewSystemClicked()\">Make New System</button></br>";
    init(true);
}
*/

int main(int argc, char** argv)
{
    init();
    return 0;
}
