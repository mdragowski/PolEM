//===============================================================================
//
// PolEM  -  Polarized Electron Mott Scattering Model for Geant4
//
// File MottCrossSection.cpp
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

// correct errors from ELSEPA
//#define PolEM_CORRECT

// implement average energy loss according to csda
//#define PolEM_CSDA

#include "MottCrossSection.hh"

Result::Result() : theta(0.), S(0.), T(0.), U(0.) {}

Result::Result(G4double thetaV, G4double SV, G4double TV, G4double UV) : theta(thetaV), S(SV), T(TV), U(UV) {}

Result::~Result(){}

void Result::Scattering(G4ThreeVector oldMomentum, G4ThreeVector oldPol){
    const G4ThreeVector zAxis(0., 0., 1.);
    // change reference frame
    G4double alpha = oldMomentum.theta(); // angle between momentum and Z axis
    G4ThreeVector rotAxis = (oldMomentum.cross(zAxis)).unit(); // normalized cross product of momentum and Z axis
    // sample phi    
    G4double phi = 2.*M_PI*G4UniformRand();
    // if momentum and polarization are not parallel, azimutal anisotropy arises
    const G4double smallValue = std::numeric_limits< G4double >::min();
    if(std::abs( oldMomentum.angle(oldPol) ) >= smallValue){
        G4ThreeVector rotPol = oldPol;
        rotPol.rotate(alpha, rotAxis);
        G4double ex = 1. + std::abs(S) * oldPol.mag();
        G4double rnd = ex*G4UniformRand();
        while(rnd > (1. + S*(rotPol[1]*cos(phi) - rotPol[0]*sin(phi)))){
            phi = 2.*M_PI*G4UniformRand();
            rnd = ex*G4UniformRand();
        }
    }
    // execute scattering
    newMomentum.setRThetaPhi(1., theta, phi);
    newMomentum.rotate(-alpha, rotAxis);
    // obtain the normal vector of the scattering plane
    G4ThreeVector normal = (oldMomentum.cross(newMomentum)).unit();
    // modification of the electron polarization vector
    G4ThreeVector vector1 = normal.cross( oldPol.cross(normal) );
    G4ThreeVector vector2 = normal.cross(oldPol);
    G4double fact = oldPol * normal;
    G4double denom = 1. + S * fact;
    fact += S;
    newPol = (fact*normal + T*vector1 + U*vector2)/denom;
}


MottTables::MottTables() : nEnergyLossSteps(100) {
    angleV = std::vector<G4double>(nAngularSteps);
    energyV = std::vector<G4double>(nEnergySteps);
    totalV = std::vector<G4double>(nEnergySteps);
    csda1V = std::vector<G4double>(nEnergyLossSteps);
    csda2V = std::vector<G4double>(nEnergyLossSteps);
    diffV = std::vector< std::vector<G4double> >(nEnergySteps,std::vector<G4double>(nAngularSteps));
    sV = std::vector< std::vector<G4double> >(nEnergySteps,std::vector<G4double>(nAngularSteps));
    tV = std::vector< std::vector<G4double> >(nEnergySteps,std::vector<G4double>(nAngularSteps));
    uV = std::vector< std::vector<G4double> >(nEnergySteps,std::vector<G4double>(nAngularSteps));
    jumperV = std::vector< std::vector<G4int> >(nJumperSteps,std::vector<G4int>(nEnergySteps));
}

MottTables::~MottTables(){}

void MottTables::Initialise(G4int Z){
    // import total cross section
    std::string fileString = "data/t" + std::to_string(Z) + ".dat";
    FILE *f = fopen(fileString.c_str(), "rb");
    if(f==nullptr)
        G4cerr << "ERROR: Cannot open binary file!\n";
    if(!fread(&energyV[0], nEnergySteps*sizeof(energyV[0]), 1, f))
        G4cerr << "ERROR: Incomplete file (energy)!\n";
    if(!fread(&totalV[0], nEnergySteps*sizeof(totalV[0]), 1, f))
        G4cerr << "ERROR: Incomplete file (total)!\n";
    if(fgetc(f) != EOF){
        if(!feof(f))
            G4cerr << "ERROR: File too large!\n";
        else
            G4cerr << "ERROR: Corrupted file!\n";
    }
    fclose(f);
    
    // import differential cross section
    fileString = "data/d" + std::to_string(Z) + ".dat";
    f = fopen(fileString.c_str(), "rb");
    if(f==nullptr)
        G4cerr << "ERROR: Cannot open binary file!\n";
    if(!fread(&angleV[0], nAngularSteps*sizeof(angleV[0]), 1, f))
        G4cerr << "ERROR: Incomplete file (angle)!\n";
    for(G4int i=0; i<nEnergySteps; i++){
        if(!fread(&diffV[i][0], nAngularSteps*sizeof(diffV[i][0]), 1, f))
            G4cerr << "ERROR: Incomplete file (diff)!\n";
    }
    for(G4int i=0; i<nEnergySteps; i++){
        if(!fread(&sV[i][0], nAngularSteps*sizeof(sV[i][0]), 1, f))
            G4cerr << "ERROR: Incomplete file (S)!\n";
    }
    for(G4int i=0; i<nEnergySteps; i++){
        if(!fread(&tV[i][0], nAngularSteps*sizeof(tV[i][0]), 1, f))
            G4cerr << "ERROR: Incomplete file (T)!\n";
    }
    for(G4int i=0; i<nEnergySteps; i++){
        if(!fread(&uV[i][0], nAngularSteps*sizeof(uV[i][0]), 1, f))
            G4cerr << "ERROR: Incomplete file (U)!\n";
    }
    if(fgetc(f) != EOF){
        if(!feof(f))
            G4cerr << "ERROR: File too large!\n";
        else
            G4cerr << "ERROR: Corrupted file!\n";
    }
    fclose(f);
    
    // check minimum and maximum energy vaules in the cross section table
    minEnergy=energyV[0];
    maxEnergy=energyV[nEnergySteps-1];
    // compute cross section table energy step
    energyStep=log(maxEnergy+1.-minEnergy)/(nEnergySteps-1);
    
#ifdef PolEM_CSDA
    // import csda energy loss
    G4double num;
    G4int stepN=0;
    std::string line;
    fileString = "data/csda" + std::to_string(Z) + ".dat";
    std::ifstream in(fileString.c_str());
    while(getline(in, line) && line.size()!=0 && stepN<nEnergyLossSteps){
        std::istringstream str(line);
        if(str >> num)
            csda1V[stepN] = num*1000.; // [keV]
        else
            G4cerr << "ERROR: Incomplete file (CSDA)!\n";
        if(str >> num)
            csda2V[stepN] = num/10000.; // [keV / nm / (g/cm3)]
        else
            G4cerr << "ERROR: Incomplete file (CSDA)!\n";
        stepN++;
    }
    in.close();
    if(stepN<nEnergyLossSteps)
        nEnergyLossSteps=stepN;

    if(csda1V[0]>minEnergy)
        G4cerr << "WARNING: Lowest energy of CSDA data set is higher than energy threshold!\n";
    if(csda1V[nEnergyLossSteps-1]<maxEnergy)
        G4cerr << "WARNING: Highest energy of CSDA data set is lower than electron energy!\n";
#endif

#ifdef PolEM_CORRECT
    // correct errors from ELSEPA
    for(G4int j=0; j<nAngularSteps; j++){
        for(G4int i=2; i<nEnergySteps-1; i++){
            if( energyV[i]>10 && std::abs(sV[i][j]-sV[i-1][j]) > 2*std::abs(sV[i-1][j]-sV[i-2][j]) ){
                G4double sTmp = sV[i-1][j] + (energyV[i]-energyV[i-1]) * (sV[i+1][j]-sV[i-1][j]) / (energyV[i+1]-energyV[i-1]);
                G4double tTmp = tV[i-1][j] + (energyV[i]-energyV[i-1]) * (tV[i+1][j]-tV[i-1][j]) / (energyV[i+1]-energyV[i-1]);
                G4double uTmp = uV[i-1][j] + (energyV[i]-energyV[i-1]) * (uV[i+1][j]-uV[i-1][j]) / (energyV[i+1]-energyV[i-1]);
                G4double dTmp = diffV[i-1][j] + (energyV[i]-energyV[i-1]) * (diffV[i+1][j]-diffV[i-1][j]) / (energyV[i+1]-energyV[i-1]);
                for(G4int n=1; n<5; n++){
                    if( std::abs(sTmp-sV[i-1][j]) < std::abs(sV[i][j]-sV[i-1][j]) || i>nEnergySteps-2-n )
                        break;
                    G4double sNew = sV[i-1][j] + (energyV[i]-energyV[i-1]) * (sV[i+n+1][j]-sV[i-1][j]) / (energyV[i+n+1]-energyV[i-1]);
                    G4double tNew = tV[i-1][j] + (energyV[i]-energyV[i-1]) * (tV[i+n+1][j]-tV[i-1][j]) / (energyV[i+n+1]-energyV[i-1]);
                    G4double uNew = uV[i-1][j] + (energyV[i]-energyV[i-1]) * (uV[i+n+1][j]-uV[i-1][j]) / (energyV[i+n+1]-energyV[i-1]);
                    G4double dNew = diffV[i-1][j] + (energyV[i]-energyV[i-1]) * (diffV[i+n+1][j]-diffV[i-1][j]) / (energyV[i+n+1]-energyV[i-1]);
                    if( std::abs(sNew-sV[i-1][j]) < std::abs(sTmp-sV[i-1][j]) ){
                        sTmp=sNew;
                        tTmp=tNew;
                        uTmp=uNew;
                        dTmp=dNew;
                    }
                }
                if( std::abs(sTmp-sV[i-1][j]) < std::abs(sV[i][j]-sV[i-1][j]) ){
                    sV[i][j]=sTmp;
                    tV[i][j]=tTmp;
                    uV[i][j]=uTmp;
                    diffV[i][j]=dTmp;
                }
                else{
                    sTmp = sV[i-1][j] + (energyV[i]-energyV[i-1]) * (sV[i-1][j]-sV[i-2][j]) / (energyV[i-1]-energyV[i-2]);
                    tTmp = tV[i-1][j] + (energyV[i]-energyV[i-1]) * (tV[i-1][j]-tV[i-2][j]) / (energyV[i-1]-energyV[i-2]);
                    uTmp = uV[i-1][j] + (energyV[i]-energyV[i-1]) * (uV[i-1][j]-uV[i-2][j]) / (energyV[i-1]-energyV[i-2]);
                    dTmp = diffV[i-1][j] + (energyV[i]-energyV[i-1]) * (diffV[i-1][j]-diffV[i-2][j]) / (energyV[i-1]-energyV[i-2]);
                    if( std::abs(sTmp-sV[i-1][j]) < std::abs(sV[i][j]-sV[i-1][j]) ){
                        sV[i][j]=sTmp;
                        tV[i][j]=tTmp;
                        uV[i][j]=uTmp;
                        diffV[i][j]=dTmp;
                    }
                }
            }
        }
    }
    for(G4int i=0; i<nEnergySteps; i++){
        for(G4int j=1; j<nAngularSteps; j++){
            if( diffV[i][j] < diffV[i][j-1] ){
                diffV[i][j] = (diffV[i][j-1]+diffV[i][j+1])/2;
                if( diffV[i][j] < diffV[i][j-1] )
                    G4cerr << "WARNING: Differential cross section table inconsistent!\n";
                sV[i][j] = (sV[i][j-1]+sV[i][j+1])/2;
                tV[i][j] = (tV[i][j-1]+tV[i][j+1])/2;
                uV[i][j] = (uV[i][j-1]+uV[i][j+1])/2;
            }
        }
    }
#endif

    // define starting points for the interpolation during simulation
    for(G4int i=0; i<nAngularSteps; i++){
        for(G4int v=0; v<nEnergySteps; v++){
            jumperV[0][v]=0;
            for(G4int j=1; j<nJumperSteps; j++){
                if(diffV[v][i] < (G4double) j/nJumperSteps)
                    jumperV[j][v]=i;
            }
        }
    }
}

// elastic cross section
G4double MottTables::CrossSection(G4double kineticEnergy){
    if(kineticEnergy < minEnergy)
        return 0;
    else if(kineticEnergy > maxEnergy){
        G4cerr << "WARNING Trying to initialise MottScatteringModel for energies outside the cross section table range\n"
               << "Cross section at " << kineticEnergy << " keV will be set to the value at " << maxEnergy << " keV\n";
        return totalV[nEnergySteps-1];
    }
    else
        return FindValue(&energyV[0],&totalV[0],kineticEnergy,nEnergySteps,0);
}

G4double MottTables::ELoss(G4double kineticEnergy){
#ifdef PolEM_CSDA
    return FindValue(&csda1V[0],&csda2V[0],kineticEnergy,nEnergyLossSteps,0);
#else
    return 0;
#endif
}

// obtain polar scattering angle theta and electron polarization transfer parameters
Result MottTables::GetTheta(G4double kineticEnergy){
    G4double theta, S, T, U;
    G4double* bestV;
    // find the cross section table for the nearest electron energy
    G4int best1 = std::max(0,(G4int)ceil(log(1.+kineticEnergy-minEnergy)/energyStep)-1);
    G4int best2 = std::min(nEnergySteps-1,best1+1);
    // sample the scattering angle theta and obtain S, T, U functions
    G4double rnd = G4UniformRand();
    if(best1 != best2){ // linear interpolation
        // theta
        G4double W1 = std::abs(kineticEnergy-energyV[best1]);
        G4double W2 = std::abs(kineticEnergy-energyV[best2]);
        G4double W3 = std::abs(energyV[best2]-energyV[best1]);
        W1 = std::max(0.,(1.-W1/W3));
        W2 = std::max(0.,(1.-W2/W3));
        bestV = &diffV[best1][0];
        G4double theta1 = FindValue(bestV,&angleV[0],rnd,nAngularSteps,
                                    jumperV[std::min(nJumperSteps-1,G4int(rnd*nJumperSteps))][best1]);
        bestV = &diffV[best2][0];
        G4double theta2 = FindValue(bestV,&angleV[0],rnd,nAngularSteps,
                                    jumperV[std::min(nJumperSteps-1,G4int(rnd*nJumperSteps))][best2]);
        theta = (theta1*W1+theta2*W2)/(W1+W2);
        // polarization transfer
        bestV = &sV[best1][0];
        G4double S1 = FindValue(&angleV[0],bestV,theta,nAngularSteps,0);
        bestV = &tV[best1][0];
        G4double T1 = FindValue(&angleV[0],bestV,theta,nAngularSteps,0);
        bestV = &uV[best1][0];
        G4double U1 = FindValue(&angleV[0],bestV,theta,nAngularSteps,0);
        bestV = &sV[best2][0];
        G4double S2 = FindValue(&angleV[0],bestV,theta,nAngularSteps,0);
        bestV = &tV[best2][0];
        G4double T2 = FindValue(&angleV[0],bestV,theta,nAngularSteps,0);
        bestV = &uV[best2][0];
        G4double U2 = FindValue(&angleV[0],bestV,theta,nAngularSteps,0);
        S = S1*W1+S2*W2;
        T = T1*W1+T2*W2;
        U = U1*W1+U2*W2;
    }
    else{ // data table for the nearest electron energy
        // theta
        bestV = &diffV[best1][0];
        theta = FindValue(bestV,&angleV[0],rnd,nAngularSteps,
                          jumperV[std::min(nJumperSteps-1,G4int(rnd*nJumperSteps))][best1]);
        // polarization transfer
        bestV = &sV[best1][0];
        S = FindValue(&angleV[0],bestV,theta,nAngularSteps,0);
        bestV = &tV[best1][0];
        T = FindValue(&angleV[0],bestV,theta,nAngularSteps,0);
        bestV = &uV[best1][0];
        U = FindValue(&angleV[0],bestV,theta,nAngularSteps,0);
    }
    // check normalization
    G4double norm = sqrt(S*S+T*T+U*U);
    S = S/norm;
    T = T/norm;
    U = U/norm;
    // until now theta is in deg
    theta = theta * M_PI / 180.;
    return Result(theta, S, T, U);
}

G4double MottTables::FindValue(G4double* xArray, G4double* yArray, G4double x, G4int nMax, G4int start){
    G4int i, mode;
    // check the order of the x-array
    if(xArray[start] < xArray[nMax-1]){
        mode=+1;
        if(start > 1)
            i=start;
        else
            i=1;
        if(x < xArray[start] || x > xArray[nMax-1])
            G4cerr << "ERROR in FindValue (+)\n";
    }
    else{
        mode=-1;
        i=nMax-1;
        if(x > xArray[0] || x < xArray[nMax-1])
            G4cerr << "ERROR in FindValue (-)\n";
    }
    // look for the right position within the x-array
    while(x > xArray[i] && i < nMax){
        i=i+mode;
        if(xArray[i] < xArray[i-mode])
            G4cerr << "ERROR in FindValue\n";
    }
    return Interpolate(xArray[i],xArray[i-mode],yArray[i],yArray[i-mode],x);
}

G4double MottTables::Interpolate(G4double x1, G4double x2, G4double y1, G4double y2, G4double x){
    if((std::abs(x1-x) <= std::abs(x1-x2)) && (std::abs(x2-x) <= std::abs(x1-x2))){
        if(x1 != x2)
            return y1+(y2-y1)/(x2-x1)*(x-x1);
        else
            return (y1+y2)/2.;
    }
    else{
        G4cerr << "ERROR in interpolate\n";
        return NAN;
    }
}


MottCrossSection::MottCrossSection(){}

MottCrossSection::~MottCrossSection(){}

void MottCrossSection::Initialise(G4int Z){
    MottTables* mottTable = new MottTables();
    mottTable->Initialise(Z);
    mottTablesMap[Z]=mottTable;
}
