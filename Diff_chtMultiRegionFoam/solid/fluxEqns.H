{ // Remember: i is the region index, ii is the group index

    for(label ii = i*nEnergyGroupsInit; ii<i*nEnergyGroupsInit+nEnergyGroupsInit; ii++)
    {
        volScalarField& fluxi = fluxSolid[ii];

        fluxi.correctBoundaryConditions();

        fvScalarMatrix fluxEqni
        (
	    fvm::ddt(invVel[ii],fluxi)
          - fvm::laplacian(D[ii], fluxi)
          + fvm::Sp(Sink[ii], fluxi)
         ==
	    neutronSource[ii] +
	    ExternalNeutronSource[ii]
            //fvModels.source(fluxi) +

        );

        fluxEqni.relax();
        fluxEqni.solve();
    }

    // Update total flux
    fluxTotSolid[i] *= 0;

    for(label ii = i*nEnergyGroupsInit; ii<i*nEnergyGroupsInit+nEnergyGroupsInit; ii++)
    {
	   fluxTotSolid[i] += fluxSolid[ii];
    }
}
