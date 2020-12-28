//=======================================================================================================
//
// File MottScatteringModel.cpp created by Michal Dragowski based on the Geant4 software
//
// This file includes software developed by Members of the Geant4 Collaboration ( http://cern.ch/geant4 )
//
// The Geant4 software is copyright of the Copyright Holders of the Geant4 Collaboration.
// It is provided under the terms and conditions of the Geant4 Software License,
// available at http://cern.ch/geant4/license .
// See the license for a list of copyright holders, disclaimer and the limitation of liability.
//
//=======================================================================================================

#include "MottScatteringModel.hh"

MottScatteringModel::MottScatteringModel() : G4VEmModel("MottScattering"){
    fParticleChange = 0;
    mottCrossSection = 0;
}

MottScatteringModel::~MottScatteringModel(){
    if(IsMaster())
        delete mottCrossSection;
}

void MottScatteringModel::Initialise(const G4ParticleDefinition* part, const G4DataVector&){
    if(part->GetParticleName() != "e-"){
        G4cerr << "EXCEPTION Trying to initialise MottScatteringModel for particle other than e-\n";
    }

    if(!fParticleChange)
        fParticleChange = GetParticleChangeForGamma();

    if(IsMaster() && !mottCrossSection){
        mottCrossSection = new MottCrossSection();
        G4double minEnergy=0, maxEnergy=1e+10;
        
        // initialize Mott cross section for each element used in detector construction
        
        for(const G4Element* element : *G4Element::GetElementTable()){
            G4int Z = G4lrint(element->GetZ());
            mottCrossSection->Initialise(Z);
            G4double minTmp = mottCrossSection->GetMin(Z) * keV;
            if(minTmp > minEnergy)
                minEnergy = minTmp;
            G4double maxTmp = mottCrossSection->GetMax(Z) * keV;
            if(maxTmp < maxEnergy)
                maxEnergy = maxTmp;
        }
        
        // set the energy range

        if(minEnergy > LowEnergyLimit())
            SetLowEnergyLimit(minEnergy);
        if(maxEnergy < HighEnergyLimit())
            SetHighEnergyLimit(maxEnergy);
        if(minEnergy > LowEnergyActivationLimit())
            SetActivationLowEnergyLimit(minEnergy);
        if(maxEnergy < HighEnergyActivationLimit())
            SetActivationHighEnergyLimit(maxEnergy);
    }
}

void MottScatteringModel::InitialiseLocal(const G4ParticleDefinition*, G4VEmModel* masterModel){
    if(!mottCrossSection){
        MottScatteringModel* masterMott = (MottScatteringModel*) masterModel;
        mottCrossSection = masterMott->GetPointer();
    }
    SetLowEnergyLimit(masterModel->LowEnergyLimit());
    SetHighEnergyLimit(masterModel->HighEnergyLimit());
    SetActivationLowEnergyLimit(masterModel->LowEnergyActivationLimit());
    SetActivationHighEnergyLimit(masterModel->HighEnergyActivationLimit());
}

G4double MottScatteringModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                         G4double kinEnergy,
                                                         G4double Z, G4double,
                                                         G4double, G4double){
    return cm2*mottCrossSection->CrossSection(kinEnergy/keV, Z);
}

void MottScatteringModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                            const G4MaterialCutsCouple* couple,
                                            const G4DynamicParticle* dp,
                                            G4double cutEnergy,
                                            G4double){
    // elastic scattering
    
    G4double kinEnergy = dp->GetKineticEnergy();
    G4int Z = SelectRandomAtom(couple,dp->GetDefinition(),kinEnergy,cutEnergy,kinEnergy)->GetZ();

    Result result = mottCrossSection->GetTheta(kinEnergy/keV, Z);
    result.Scattering(dp->GetMomentumDirection(), dp->GetPolarization());
    G4ThreeVector newMomentum = result.GetMomentum();
    G4ThreeVector newPol = result.GetPol();

    // energy loss
    
    G4double eDep = mottCrossSection->ELoss(kinEnergy/keV, Z) * keV;
    if(eDep){
        G4double stepSize = fParticleChange->GetCurrentTrack()->GetStepLength() / nm;
        const G4double density = couple->GetMaterial()->GetDensity() * cm3/g;
        eDep *= stepSize * density;
    }

    fParticleChange->ProposeMomentumDirection(newMomentum);
    fParticleChange->ProposePolarization(newPol);
    fParticleChange->SetProposedKineticEnergy(kinEnergy-eDep);
    fParticleChange->ProposeNonIonizingEnergyDeposit(eDep);
    fParticleChange->ProposeLocalEnergyDeposit(eDep);
}
