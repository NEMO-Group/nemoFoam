{
	for(label ii = i*nDecGroupsInit; ii<i*nDecGroupsInit+nDecGroupsInit; ii++)
	{
		volScalarField& deci = decSolid[ii];

		deci.correctBoundaryConditions();

		fvScalarMatrix decEqni
		(
			fvm::ddt(deci)
		+	fvm::Sp(decLambda[ii], deci)
		==
			decBeta[ii]*EfRfTot[i]
		//	fvModels.source(deci)
		);

		decEqni.relax();
		decEqni.solve();
	}
}
