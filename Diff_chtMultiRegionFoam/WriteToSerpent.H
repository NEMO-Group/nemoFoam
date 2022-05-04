forAll(solidRegions, i)
{

    if (runTime.writeTime())
    {
        volScalarField rho
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
        );

        runTime.write();
    }
}
