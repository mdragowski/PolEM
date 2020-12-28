//=======================================================================================================
//
// File MottScatteringModel.hh created by Michal Dragowski based on the Geant4 software
//
// This file includes software developed by Members of the Geant4 Collaboration ( http://cern.ch/geant4 )
//
// The Geant4 software is copyright of the Copyright Holders of the Geant4 Collaboration.
// It is provided under the terms and conditions of the Geant4 Software License,
// available at http://cern.ch/geant4/license .
// See the license for a list of copyright holders, disclaimer and the limitation of liability.
//
//=======================================================================================================

#ifndef MottScatteringModel_h
#define MottScatteringModel_h 1

#include "MottCrossSection.hh"
#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4DataVector.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"
#include "G4DynamicParticle.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include <vector>

class MottScatteringModel : public G4VEmModel{
private:
    G4ParticleChangeForGamma* fParticleChange;
    MottCrossSection* mottCrossSection;
public:
    MottScatteringModel();
    virtual ~MottScatteringModel();
    virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);
    virtual void InitialiseLocal(const G4ParticleDefinition*, G4VEmModel*);
    virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                G4double kinEnergy,
                                                G4double Z,
                                                G4double A,
                                                G4double cut,
                                                G4double emax);
    virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                   const G4MaterialCutsCouple*,
                                   const G4DynamicParticle*,
                                   G4double tmin,
                                   G4double maxEnergy);
    inline MottCrossSection* GetPointer(){ return mottCrossSection; }
};

#endif
