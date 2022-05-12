{
if (Debug)  {Info << "\n\nupdatePowerSource: START " << endl;}
if (Debug)  {Info << "\n - i =  "<< i << ", region is: " << solidRegions[i].name() << "." << endl ; }


	// Reset total fields
	qTot[i] *= 0;
		if (Debug)  {Info << "\n --- qTot["<< i << "] = 0" << endl ; }

	QTot[i] *= 0;
		if (Debug)  {Info << "\n --- QTot["<< i << "] = 0" << endl ; }


	for(label ii = 0; ii<nEnergyGroupsInit; ii++)
	{
	  	label posEqn = ii + nEnergyGroupsInit*i;
			if (Debug)  {Info << "\n --- posEqn = "<< posEqn << endl ; }

		// Reset field
		q_PromptFission[posEqn] *= 0;
			if (Debug)  {Info << "\n --- q_PromptFission["<< posEqn << "] = 0" << endl ; }

		q_Kerma[posEqn] *= 0;
			if (Debug)  {Info << "\n --- q_Kerma["<< posEqn << "] = 0" << endl ; }

		// Evaluate power deposition based on current fluxes and cross-sections
		q_PromptFission[posEqn] = Ef[posEqn] * (1-decBetaTot[i]) * Rf[posEqn];
			if (Debug)  {Info << "\n --- q_PromptFission["<<posEqn<<"] = Ef["<<posEqn<<"] * (1-decBetaTot["<<i<<"]) * Rf["<<posEqn<<"]"<< endl ; }

		q_Kerma[posEqn] = Kerma[posEqn] * fluxSolid[posEqn];//(1-decBetaTot[i]) * Ri[posEqn];
  			if (Debug)  {Info << "\n --- q_Kerma["<<posEqn<<"] = kerma["<<posEqn<<"] * (1-decBetaTot["<<i<<"]) * Ri["<<posEqn<<"]" << endl ; }

		qTot[i] += (q_PromptFission[posEqn] + q_Kerma[posEqn]) ;
  			if (Debug)  {Info << "\n --- qTot["<<i<<"] += (q_PromptFission["<<posEqn<<"] + q_Kerma["<<posEqn<<"])" << endl ; }
	}

	for(label ii = 0; ii<nDecGroupsInit; ii++)
	{
	  	label posEqn = ii + nDecGroupsInit*i;
			if (Debug)  {Info << "\n --- posEqn = "<< posEqn << endl ; }

		// Reset field
		q_Decay[posEqn] *= 0;
			if (Debug)  {Info << "\n --- q_Decay["<< posEqn << "] = 0" << endl ; }

		// Evaluate power deposition based on current decay heat precursors
		q_Decay[posEqn] = decLambda[posEqn]*decSolid[posEqn];
			if (Debug)  {Info << "\n --- q_Decay["<<posEqn<<"] = decLambda["<<posEqn<<"]*decSolid["<<posEqn<<"]"<< endl ; }

		qTot[i] += q_Decay[posEqn] ;
  			if (Debug)  {Info << "\n --- qTot["<<i<<"] += (q_PromptFission["<<posEqn<<"] + q_Kerma["<<posEqn<<"])" << endl ; }
	}

//	QTot[i] = unitField[i]*fvc::domainIntegrate(qTot[i]).value()*dimensionedScalar("", dimPower, 1);
//	QTot[i] = qTot[i]*solidRegions[i].V();
//  			if (Debug)  {Info << "\n --- QTot["<<i<<"] = qTot["<<i<<"]*solidRegions["<<i<<"].V()" << endl ; }
//  			if (Debug)  {Info << "\n --- QTot["<<i<<"] = unitField["<<i<<"]*fvc::domainIntegrate(qTot["<<i<<"]).value()*dimensionedScalar("", dimPower, 1)" << endl ; }
}






















