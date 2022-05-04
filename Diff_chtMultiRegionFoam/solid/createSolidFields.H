// Initialise solid field pointer lists
PtrList<solidThermo> thermoSolid(solidRegions.size());
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size());
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size());
PtrList<volScalarField> volPowerSolid(solidRegions.size()); //new piece, WORKS! It was added to take into account the volumetric power deposition into the solid due to nuclear reactions. The same is done for the fluid regions
PtrList<volScalarField> T(solidRegions.size());
PtrList<volScalarField> rho(solidRegions.size());

//PtrList<volScalarField> rhoSolid(solidRegions.size()); //new piece. It is added to print the rho fields also for region, in order to automatize the coupling process with Serpent. NEED TO BE CHECKED IF WORKS!!

// Populate solid field pointer lists
forAll(solidRegions, i)
{
/*    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoSolid\n" << endl;
    thermoSolid.set(i, solidThermo::New(solidRegions[i]));
*/
    Info<< "    Adding to fvModelsSolid\n" << endl;
    fvModelsSolid.set
    (
        i,
        new Foam::fvModels(solidRegions[i])
    );

    Info<< "    Adding fvConstraintsSolid\n" << endl;
    fvConstraintsSolid.set
    (
        i,
        new fvConstraints(solidRegions[i])
    );

    T.set				 // start new piece, WORKS!
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    rho.set				 // start new piece, WORKS!
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );


/*    Info<< "    Adding to volPowerSolid\n" << endl;
    volPowerSolid.set				 // start new piece, WORKS!
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "volPower",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );  					// ends new piece, WORKS!
*/

/*    Info<< "    Adding to rhoSolid\n" << endl;   // start new piece, CHECK IT


    (
        i,
        volScalarField thermoSolid[i].rho 
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoSolid[i].rho()
        )
    );						// end new piece, CHECK IT
*/

}


#include "createSolidNuclearFields.H"
