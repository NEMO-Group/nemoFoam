{
if (Debug)  {Info << "\n\nupdateFissionRate: START " << endl;}
if (Debug)  {Info << "\n - i =  "<< i << ", region is: " << solidRegions[i].name() << "." << endl ; }

	// Reset total fields
	RfTot[i] *= 0;
		if (Debug)  {Info << "\n --- RfTot["<< i << "] = 0" << endl ; }

	RiTot[i] *= 0;
		if (Debug)  {Info << "\n --- RiTot["<< i << "] = 0" << endl ; }

	NuRfTot[i] *= 0;
		if (Debug)  {Info << "\n --- NuRfTot["<< i << "] = 0" << endl ; }

	EfRfTot[i] *= 0;
		if (Debug)  {Info << "\n --- EfRfTot["<< i << "] = 0" << endl ; }



	for(label ii = 0; ii<nEnergyGroupsInit; ii++)
	{
	  	label posEqn = ii + nEnergyGroupsInit*i;
  			if (Debug)  {Info << "\n --- posEqn = "<< posEqn << endl ; }

	// Reset field
		Rf[posEqn] *= 0;
			if (Debug)  {Info << "\n --- Rf["<< posEqn << "] = 0" << endl ; }

		Ri[posEqn] *= 0;
			if (Debug)  {Info << "\n --- Ri["<< posEqn << "] = 0" << endl ; }

	// Evaluate new fission rates based on current fluxes and cross-sections
		Rf[posEqn] = Sf[posEqn]*fluxSolid[posEqn];
  			if (Debug)  {Info << "\n --- Rf["<<posEqn<<"] = Sf["<<posEqn<<"]*fluxSolid["<<posEqn<<"]" << endl ; }

	// Evaluate new interaction rates based on current fluxes and cross-sections
		Ri[posEqn] = (Stot[posEqn] - Sf[posEqn])*fluxSolid[posEqn];
  			if (Debug)  {Info << "\n --- Ri["<<posEqn<<"] = (Stot["<<posEqn<<"] - Sf["<<posEqn<<"]) *fluxSolid["<<posEqn<<"]" << endl ; }


	// Compute total fields

		RfTot[i] += Rf[posEqn];
  			if (Debug)  {Info << "\n --- RfTot["<<i<<"] += Rf["<<posEqn<<"]" << endl ; }
		RiTot[i] += Ri[posEqn];
  			if (Debug)  {Info << "\n --- RiTot["<<i<<"] += Ri["<<posEqn<<"]" << endl ; }
		NuRfTot[i] += Nu[posEqn]*Rf[posEqn];
  			if (Debug)  {Info << "\n --- NuRfTot["<<i<<"] += Nu["<<posEqn<<"]*Rf["<<posEqn<<"]" << endl ; }
		EfRfTot[i] += Ef[posEqn]*Rf[posEqn];;
  			if (Debug)  {Info << "\n --- EfRfTot["<<i<<"] += Ef["<<posEqn<<"]*Rf["<<posEqn<<"]" << endl ; }
	}
if (Debug)  {Info << "\n\nupdateFissionRate: END " << endl;}
}
