
#include "./custom_main.h"

/**
 *	\main drug_AGS custom main file
 *	\brief Custom module file for drug_AGS example
 * 
 *	\details Modules needed for the drug_AGS example. This custom module can be used to study the inhibition of AGS cell lines with AKT, beta-catenin and TAK inhibitors.
 *
 *
 *	\date 19/10/2020
 *	\author Arnau Montagud, BSC-CNS, with code previously developed by Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */

void inject_density_sphere(int density_index, double concentration, double membrane_lenght) 
{
	// Inject given concentration on the extremities only
	#pragma omp parallel for
	for( int n=0; n < microenvironment.number_of_voxels() ; n++ )
	{
		auto current_voxel = microenvironment.voxels(n);
		std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};

		if ((membrane_lenght - norm(cent)) <= 0)
			microenvironment.density_vector(n)[density_index] = concentration; 	
	}
}

void remove_density( int density_index )
{	
	for( int n=0; n < microenvironment.number_of_voxels() ; n++ )
		microenvironment.density_vector(n)[density_index] = 0; 	
	std::cout << "Removal done" << std::endl;
}

			// TNF 
			// if ( PhysiCell_globals.current_time >= time_put_tnf )
			// {
			// 	time_tnf_next = PhysiCell_globals.current_time + duration_add_tnf;
			// 	time_put_tnf += time_add_tnf;
			// }

			// if ( PhysiCell_globals.current_time >= time_remove_tnf )
			// {
			// 	int k = microenvironment.find_density_index("tnf");
			// 	if ( k >= 0 )
			// 		remove_density(k);
			// 	time_remove_tnf += PhysiCell_settings.max_time;
			// }

			// if ( PhysiCell_globals.current_time <= time_tnf_next )
			// {
			// 	int k = microenvironment.find_density_index("tnf");
			// 	if ( k >= 0 ) 
			// 		inject_density_sphere(k, concentration_tnf, membrane_lenght);
			// }

// void set_density_for_current_time (std::string dens, double time_add_dens, double time_put_dens, double duration_add_dens, double time_remove_dens, double time_dens_next, double concentration_dens, double membrane_length )
// {
// 	if (PhysiCell_globals.current_time >= time_put_dens)
// 	{
// 		time_dens_next = PhysiCell_globals.current_time + duration_add_dens;
// 		time_put_dens += time_add_dens;
// 	}

// 	if (PhysiCell_globals.current_time >= time_remove_dens)
// 	{
// 		int k = microenvironment.find_density_index(dens);
// 		if (k >= 0)
// 			remove_density(k);
// 		time_remove_dens += PhysiCell_settings.max_time;
// 	}

// 	if (PhysiCell_globals.current_time <= time_dens_next)
// 	{
// 		int k = microenvironment.find_density_index(dens);
// 		if (k >= 0)
// 			inject_density_sphere(k, concentration_dens, membrane_length);
// 	}
// }


/* Change the current value of the input coefficient, increase or decrease according to up value */
void evolve_coef( int up, double* coef, double dt )
/**{ 
	// increase exponentially
	if ( up )
	{
		if ( (*coef) < EPSILON ) 
			(*coef) = EPSILON; 	
		(*coef) = std::sqrt( (*coef) );
		(*coef) = (*coef) > 1 ? (1-EPSILON) : (*coef);
	}
	else
	{
		// decrease exponentially
		if ( (*coef) >= 1 )
			(*coef) = 1 - EPSILON;
		(*coef) *= (*coef);	
		(*coef) = (*coef) < 0 ? EPSILON : (*coef);
	}
}*/
{ 
	// if up, increase, else decrease
	if ( !up )
		dt = -dt;

	(*coef) +=  (*coef) * (1 - (*coef)) * dt/10.0 ; //what is the dt/10.0?

	(*coef) = (*coef) > 1 ? (1-EPSILON) : (*coef);
	(*coef) = (*coef) < 0 ? (EPSILON) : (*coef);
}


void evolve_linear( int up, double* var, double dt )
{
		(*var) = 1.5 * (*var) ;
}