//===============================================================================
//
// PolEM  -  Polarized Electron Mott Scattering Model for Geant4
//
// File MottCrossSection.hh
//
// Authors: Michal Dragowski, Gunter Weber and Marta Wlodarczyk
//
// Copyright (c) 2020 Michal Dragowski, Gunter Weber and Marta Wlodarczyk
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software. Publication of the results
// obtained with the Software shall contain a reference to
// M. Dragowski et al., Nucl. Instr. Meth. B 488, 37 (2021)
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//===============================================================================

#ifndef MottCrossSection_h
#define MottCrossSection_h 1

#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cmath>

class Result{
private:
    G4double theta;
    G4double S;
    G4double T;
    G4double U;
    G4ThreeVector newMomentum;
    G4ThreeVector newPol;
public:
    Result();
    Result(G4double, G4double, G4double, G4double);
    ~Result();
    inline G4ThreeVector GetMomentum(){ return newMomentum; };
    inline G4ThreeVector GetPol(){ return newPol; };
    void Scattering(G4ThreeVector, G4ThreeVector);
};

class MottTables{
private:
    static const G4int nAngularSteps = 606;
    static const G4int nEnergySteps = 500;
    static const G4int nJumperSteps = 100;
    G4int nEnergyLossSteps;
    G4double minEnergy;
    G4double maxEnergy;
    G4double energyStep;
    std::vector<G4double> energyV;
    std::vector<G4double> totalV;
    std::vector<G4double> angleV;
    std::vector<G4double> csda1V;
    std::vector<G4double> csda2V;
    std::vector< std::vector<G4double> > diffV;
    std::vector< std::vector<G4double> > sV;
    std::vector< std::vector<G4double> > tV;
    std::vector< std::vector<G4double> > uV;
    std::vector< std::vector<G4int> > jumperV;
    G4double FindValue(G4double* xArray, G4double* yArray, G4double x, G4int nMax, G4int start);
    G4double Interpolate(G4double x1, G4double x2, G4double y1, G4double y2, G4double x);
public:
    MottTables();
    ~MottTables();
    void Initialise(G4int Z);
    G4double CrossSection(G4double kineticEnergy);
    G4double ELoss(G4double kineticEnergy);
    Result GetTheta(G4double kineticEnergy);
    inline G4double GetMin(){ return minEnergy; };
    inline G4double GetMax(){ return maxEnergy; };
};

class MottCrossSection{
private:
    std::unordered_map<G4int, MottTables*> mottTablesMap;
public:
    MottCrossSection();
    ~MottCrossSection();
    void Initialise(G4int Z);
    inline G4double CrossSection(G4double kineticEnergy, G4int Z){ return mottTablesMap.at(Z)->CrossSection(kineticEnergy); };
    inline G4double ELoss(G4double kineticEnergy, G4int Z){ return mottTablesMap.at(Z)->ELoss(kineticEnergy); };
    inline Result GetTheta(G4double kineticEnergy, G4int Z){ return mottTablesMap.at(Z)->GetTheta(kineticEnergy); };
    inline G4double GetMin(G4int Z){ return mottTablesMap.at(Z)->GetMin(); };
    inline G4double GetMax(G4int Z){ return mottTablesMap.at(Z)->GetMax(); };
};

#endif
