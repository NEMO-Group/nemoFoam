{
	for(label ii = i*nPrecGroupsInit; ii<i*nPrecGroupsInit+nPrecGroupsInit; ii++)
	{
		volScalarField& preci = precSolid[ii];

		preci.correctBoundaryConditions();

		fvScalarMatrix precEqni
		(
			fvm::ddt(preci)
		+	fvm::Sp(precLambda[ii], preci)
		==
			precBeta[ii]*NuRfTot[i]
		//	fvModels.source(preci)
		);

		precEqni.relax();
		precEqni.solve();
	}
}
