const fvMesh& mesh = solidRegions[i];

// solidThermo& thermo = thermoSolid[i];

// tmp<volScalarField> trho = thermo.rho();
// const volScalarField& rho = trho();

//volScalarField& rhorho = RhoRhoSolid[i]; //TO CHECK
//volScalarField& rho2 = rhoFluid[i];


// volScalarField& e = thermo.he();

// const Foam::fvModels& fvModels = fvModelsSolid[i];
// Foam::fvConstraints& fvConstraints = fvConstraintsSolid[i];

// volScalarField& volPower = volPowerSolid[i];

#include "checkRadiationModel.H"

solidNoLoopControl& pimple = pimples.solid(i);
