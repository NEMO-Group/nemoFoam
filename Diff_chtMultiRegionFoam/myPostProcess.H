Switch DoPostProcess(runTime.controlDict().lookup("DoPostProcess"));

if (DoPostProcess)
{
	forAll(solidRegions, jj)
	{
		Info<< "\n PostProcess for solid region: " << solidRegions[jj].name() << endl;

		for(label iiii = 0; iiii<nEnergyGroupsInit; iiii++)
		{
			const label fluxIndex = iiii+jj*nEnergyGroupsInit;

			Info<< nl << "\n     Integral of:" << fluxSolid[fluxIndex].name() << " = " << fvc::domainIntegrate(fluxSolid[fluxIndex]).value() << endl;
		}
	}


Info<< "\n ---- POWER DEPOSITION " << endl;

	forAll(solidRegions, jj)
	{
		Info<< "\n Total power deposition in solid region: " << solidRegions[jj].name() << " = " << fvc::domainIntegrate(qTot[jj]).value() << " W."<< endl;
	}

/*
Info<< "\n ---- POWER DEPOSITION " << endl;

		Info<< "\n Total power deposition in solid = " << fvc::domainIntegrate(qTot) << " W."<< endl;

*/

/*
Info<< "\n ---- POWER DEPOSITION bis" << endl;

	forAll(solidRegions, jj)
	{
		Info<< "\n Total power deposition in solid region: " << solidRegions[jj].name() << " = " << gSum(qTot[jj]*solidRegions[jj].V()) << " W."<< endl;
	}
*/
}


