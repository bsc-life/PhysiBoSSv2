#include "./boolean_model_interface.h"


void update_custom_variables( Cell* pCell )
{
	// first density is oxygen - shouldn't be changed: index from 1
	for (int i = 0; i < microenvironment.number_of_densities(); i++) 
	{
		//TODO: check that update_custom_variables of drugs work
		std::string drug_name = microenvironment.density_names[i];
		if (drug_name != "oxygen")
		{
			int drug_index = microenvironment.find_density_index(drug_name);
			int index_drug_conc = pCell->custom_data.find_variable_index("concentration_reporter_"+drug_name);
			// int index_drug_conc = pCell->custom_data.find_variable_index(drug_name + "_concentration");
			int index_drug_node = pCell->custom_data.find_variable_index(drug_name + "_node");
			string drug_target = get_value(drug_targets, drug_name);
			// pCell->custom_data.variables.at(index_drug_conc).value = pCell->nearest_density_vector()[drug_index];
			// pCell->custom_data.variables.at(index_drug_node).value = pCell->boolean_network.get_node_value("anti_" + drug_target);
			pCell->custom_data[index_drug_conc] = pCell->nearest_density_vector()[drug_index];
			pCell->custom_data[index_drug_node] = pCell->boolean_network.get_node_value("anti_" + drug_target);
		}	
	}
	//TODO: should multiplier_reporter update be here?
}


void set_boolean_node (Cell* pCell, std::string drug_name, int index, double threshold) {
	if (index != -1)
		{
			string drug_target = get_value(drug_targets, drug_name);
			std::string node_name = "anti_" + drug_target;
			 // get internalized substrate concentration. not really, you are getting the nearest voxel values.
    		double drug_conc = pCell->nearest_density_vector()[index];
			double cell_viability = get_cell_viability_for_drug_conc(pCell, drug_conc, parameters.strings("cell_line"), drug_name, index);
			double cell_inhibition = 1 - cell_viability;
			double random_num = (double) rand()/RAND_MAX;
			if (random_num <= cell_inhibition) 
			{
			
				pCell->boolean_network.set_node_value(node_name, 1);
			}
			else 
			{
				pCell->boolean_network.set_node_value(node_name, 0);
			}
			
		}
		// old:
		// 	if (index != -1)
		// {
		// 	std::string node_name = "anti_" + drug_name;
		// 	 double drug_conc = pCell->phenotype.molecular.internalized_total_substrates[index];
		// 	if (drug_conc > threshold) 
		// 	{
			
		// 		pCell->boolean_network.set_node_value(node_name, 1);
		// 	}
		// 	else 
		// 	{
		// 		pCell->boolean_network.set_node_value(node_name, 0);
		// 	}
			
		// }

}

void set_input_nodes(Cell* pCell) {

	if (PhysiCell::parameters.ints("simulation_mode") == 0)
	{	
		// single inhibition - just one drug is present 

		std::string drug_name = microenvironment.density_names[1];
		int drug_index = microenvironment.find_density_index(drug_name);
		int drug_threshold = PhysiCell::parameters.doubles( "threshold_" + drug_name);

		if (pCell->type_name == drug_name + "_sensitive")
		{
			// cell is sensitive to the drug -> set boolean node
			set_boolean_node(pCell, drug_name, drug_index, drug_threshold);
		}
	}
	else if (PhysiCell::parameters.ints("simulation_mode") == 1)
	{	
		// double inhibition - two drugs are present - we have 4 cell strains 

		std::string drug1_name = microenvironment.density_names[1];
		int drug1_index = microenvironment.find_density_index(drug1_name);
		int drug1_threshold = PhysiCell::parameters.doubles( "threshold_" + drug1_name);

		std::string drug2_name = microenvironment.density_names[2];
		int drug2_index = microenvironment.find_density_index(drug2_name);
		int drug2_threshold = PhysiCell::parameters.doubles( "threshold_" + drug2_name);
	
		if (pCell->type_name == drug1_name + "_sensitive")
		{	
			// cell is sensitive to both drugs
			set_boolean_node(pCell, drug1_name, drug1_index, drug1_threshold);
			set_boolean_node(pCell, drug2_name, drug2_index, drug2_threshold);
		}
		else if (pCell->type_name == drug2_name + "_resistant")
		{
			// cell is only sensitive to the first drug
			set_boolean_node(pCell, drug1_name, drug1_index, drug1_threshold);
		}
		else if (pCell->type_name == drug2_name +  "_sensitive")
		{
			// cell is only sensitive to the second drug
			set_boolean_node(pCell, drug2_name, drug2_index, drug2_threshold);
		}

		// else: cell is sensitive to no drug --> no boolean node is set
	}
	//old:
	// 	if (PhysiCell::parameters.ints("simulation_mode") == 0)
	// {	
	// 	// single inhibition - just one drug is present 

	// 	std::string drug_name = microenvironment.density_names[1];
	// 	int drug_index = microenvironment.find_density_index(drug_name);
	// 	int drug_threshold = PhysiCell::parameters.doubles( "threshold_" + drug_name);

	// 	if (pCell->type_name == drug_name + "_sensitive")
	// 	{
	// 		// cell is sensitive to the drug -> set boolean node
	// 		set_boolean_node(pCell, drug_name, drug_index, drug_threshold);
	// 	}
	// }
	// else if (PhysiCell::parameters.ints("simulation_mode") == 1)
	// {	
	// 	// double inhibition - two drugs are present - we have 4 cell strains 

	// 	std::string drug1_name = microenvironment.density_names[1];
	// 	int drug1_index = microenvironment.find_density_index(drug1_name);
	// 	int drug1_threshold = PhysiCell::parameters.doubles( "threshold_" + drug1_name);

	// 	std::string drug2_name = microenvironment.density_names[2];
	// 	int drug2_index = microenvironment.find_density_index(drug2_name);
	// 	int drug2_threshold = PhysiCell::parameters.doubles( "threshold_" + drug2_name);
	
	// 	if (pCell->type_name == drug1_name + "_sensitive")
	// 	{	
	// 		// cell is sensitive to both drugs
	// 		set_boolean_node(pCell, drug1_name, drug1_index, drug1_threshold);
	// 		set_boolean_node(pCell, drug2_name, drug2_index, drug2_threshold);
	// 	}
	// 	else if (pCell->type_name == drug2_name + "_resistant")
	// 	{
	// 		// cell is only sensitive to the first drug
	// 		set_boolean_node(pCell, drug1_name, drug1_index, drug1_threshold);
	// 	}
	// 	else if (pCell->type_name == drug2_name +  "_sensitive")
	// 	{
	// 		// cell is only sensitive to the second drug
	// 		set_boolean_node(pCell, drug2_name, drug2_index, drug2_threshold);
	// 	}
	// 	// else: cell is sensitive to no drug --> no boolean node is set
	// }
}

void from_nodes_to_cell (Cell* pCell, Phenotype& phenotype, double dt)
{
	std::vector<bool>* nodes = pCell->boolean_network.get_nodes();
	int bn_index;

	// For AGS model

	std::string prosurvival_basename = "Prosurvival_b";
	std::string antisurvival_basename = "Antisurvival_b";
	double prosurvival_value = 0.0;
	double antisurvival_value = 0.0;

	for(int i=1; i<=3; i++)
	{
		bn_index = pCell->boolean_network.get_node_index( prosurvival_basename + std::to_string(i) );
		if ( (*nodes)[bn_index] > 0)
		{
			prosurvival_value += (*nodes)[bn_index];
		}
		bn_index = pCell->boolean_network.get_node_index( antisurvival_basename + std::to_string(i) );
		if ( (*nodes)[bn_index] > 0)
		{
			antisurvival_value += (*nodes)[bn_index];
		}
	}
	int start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live ); // Q_phase_index; 
	int end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live ); // K_phase_index;
	int apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );  

	// live model
	if( pCell->phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
	{
		// multiplier implementation, old AGS
		// phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier * pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index);
		// pCell->phenotype.death.rates[apoptosis_index] = ( 1 / multiplier ) * pCell->phenotype.death.rates[apoptosis_index];
		//alternative: // phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier *	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
		//alternative: // phenotype.death.rates[apoptosis_index] = ( 1 / multiplier ) * pCell->parameters.pReference_live_phenotype->.death.rates[apoptosis_index]; // for completeness sake, will it work?

		// switch implementation:

		// double multiplier = 1.0;
		// multiplier = ( prosurvival_value + 1 ) / ( antisurvival_value + 1 ) ; //[0.25, 0.33, 0.5, 0.58, 0.66, 0.75, 1, 1.5, 2, 3, 4]
		// // TODO: connect multiplier_reporter with multiplier
		// static int multiplier_reporter_index = pCell->custom_data.find_variable_index("multiplier_reporter");
	    // pCell->custom_data[multiplier_reporter_index] = multiplier;

		// if ( multiplier == 1 ) // same prosurvival_value and antisurvival_value
		// {
		// 	pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) = PhysiCell::parameters.doubles("base_transition_rate");
		// }
		// else if ( multiplier > 1 ) // bigger prosurvival_value than antisurvival_value
		// {
		// 	double high_transition_rate = PhysiCell::parameters.doubles("base_transition_rate") * PhysiCell::parameters.doubles("transition_rate_multiplier");
		// 	pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) = high_transition_rate;
		// }
		// else if ( multiplier < 1 ) // smaller prosurvival_value than antisurvival_value
		// {
		// 	pCell->phenotype.death.rates[apoptosis_index] = PhysiCell::parameters.doubles("apoptosis_rate_multiplier") * pCell->phenotype.death.rates[apoptosis_index];
		// }
		//TODO: substitute transition_rate_multiplier and apoptosis_rate_multiplier with the multiplier defined by prosurvival_value than antisurvival_value

		// TODO: fer figura de saturació i pujar i baixar amb tokens
		// https://web.expasy.org/cellosaurus/CVCL_0139
		// normalitzar dades Asmund i eixides de simulació

		// void simple_volume_function( Cell* pCell, Phenotype& phenotype, double dt )
		// {
		// 	double V_target = phenotype.volume.target_solid_nuclear * (1.0 + phenotype.volume.target_cytoplasmic_to_nuclear_ratio) / (1 - phenotype.volume.fluid_fraction);
		// 	double rate = phenotype.cycle.model().transition_rate(0,0) * log(0.1);
		// 	double addme = V_target;
		// 	addme -= phenotype.volume.total;
		// 	addme *= rate;
		// 	addme *= dt; // dt*rate*( V_target - V )
		// 	phenotype.volume.total += addme;
		// 	return;
		// }

		// hyperbole implementation:
		// no base_transition_rate, no transition_rate_multiplier, no apoptosis_rate_multiplier
		// fit target_prolif, pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate, multiplier_divisor
		// AGS doubling time: 20 - 24 hours depending on source. rate is "1/min". thus, 20 hours = 0.00083

		// prosurvival_value i antisurvival_value = [0, 1, 2, 3], so we add 3 to the base multiplier_exp to have it between 0 and 6, instead of -3 and +3.
		// double multiplier_exp = 3.0; // [0, 6]
		// multiplier_exp += prosurvival_value;
		// multiplier_exp -= antisurvival_value;
		// double target_prolif = 0.00083;
		// logistic curve: max_value/(1 + e^(-(x- xmid/scale))
		// el rate multiplicador ha de ser <1 per a tenir el comportament exponencial, si és 1, quan satura i ja no baixa 
		// double multiplier_divisor = 60;
		// double rhs = multiplier_exp/multiplier_divisor * ( target_prolif - pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) );
		// pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) += rhs;
		// double rhs = multiplier_exp/6 * ( target_prolif - pCell->phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) );
		// pCell->phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) += rhs;
		
		// rhs ara és positiu i negatiu per a que puge i baixe:
		// fit target_prolif, pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate
		// multiplier_divisor
		double multiplier_exp = 0.0; // [-3, +3], but then transition_rate can be negative!
		multiplier_exp += prosurvival_value;
		multiplier_exp -= antisurvival_value;
		double target_prolif = 0.00083;
		// logistic curve: max_value/(1 + e^(-(x- xmid/scale))
		// rate * (xmax - x)
		// el rate multiplicador ha de ser <1 per a tenir el comportament exponencial, si és 1,  satura a la primea i ja no baixa 
		// double divisor_prolif = 30;
		double divisor_prolif = PhysiCell::parameters.doubles("divisor_prolif");
		// double rhs = multiplier_exp/divisor_prolif * ( target_prolif - pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) );
		double rhs = prosurvival_value/divisor_prolif * ( target_prolif - pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) );
		pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) += rhs;
		if ( pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) < 0)
		{
			pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) = 0;
		}

		double multiplier_apop = PhysiCell::parameters.doubles("multiplier_apop");
		pCell->phenotype.death.rates[apoptosis_index] = (antisurvival_value * -1 * multiplier_apop) * pCell->phenotype.death.rates[apoptosis_index];
		// if ( multiplier_exp < 1 )
		// {
		// 	pCell->phenotype.death.rates[apoptosis_index] = (multiplier_exp * -1 * multiplier_apop) * pCell->phenotype.death.rates[apoptosis_index];
		// }

		//  connect multiplier_reporter with multiplier
		static int multiplier_reporter_index = pCell->custom_data.find_variable_index("multiplier_reporter");
		pCell->custom_data[multiplier_reporter_index] = multiplier_exp;
		static int multiplier_prolif_index = pCell->custom_data.find_variable_index("divisor_prolif_reporter");
		pCell->custom_data[multiplier_prolif_index] = divisor_prolif;
		static int multiplier_apop_index = pCell->custom_data.find_variable_index("multiplier_apop_reporter");
		pCell->custom_data[multiplier_apop_index] = multiplier_apop;

		// TODO: overwrite pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) if the sims do not grow enough
	}

	if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
	{
		std::cout << pCell->phenotype.cycle.current_phase().name << "\n" << std::endl;
	}
	pCell->set_internal_uptake_constants(dt);
}

// Following is NOT used
void from_nodes_to_cell_prostate (Cell* pCell, Phenotype& phenotype, double dt)
{
	std::vector<bool>* nodes = pCell->boolean_network.get_nodes();
	int bn_index;

	// Prostate live model
	// map apoptosis, proliferation and invasion values to agent-based model

	if( pCell->phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
	{
		int start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		int end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		int apoptosis_index = phenotype.death.find_death_model_index(PhysiCell_constants::apoptosis_death_model);
		
		// Update Apoptosis 
		if(pCell->boolean_network.get_node_value("Apoptosis"))
		{
			// simple implementation, just lead immediately to death
			// pCell->start_death(apoptosis_index);

			// increase death rate whenever the node is ON 
			pCell->phenotype.death.rates[apoptosis_index] = PhysiCell::parameters.doubles("apoptosis_rate_multiplier") * pCell->phenotype.death.rates[apoptosis_index];
		}

		// Update Adhesion
		if( pCell->boolean_network.get_node_value("EMT"))
		{
			// reduce cell-cell adhesion 
			// pCell->evolve_coef_sigmoid( 
			// 	pCell->boolean_network.get_node_value("EMT"), phenotype.mechanics.cell_cell_adhesion_strength, dt 
			// );

			//phenotype.mechanics.cell_cell_adhesion_strength = PhysiCell::parameters.doubles("homotypic_adhesion_max") * theta 
		}


		// Update pReference_live_phenotype for proliferation node 

		if (pCell->boolean_network.get_node_value("Proliferation")) 
		{
			// multiplier implementation
			//pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) *= 2.5;
			//std::cout << "Rate up! " << std::endl;

			//switch implementation
			double high_transition_rate = PhysiCell::parameters.doubles("base_transition_rate") * PhysiCell::parameters.doubles("transition_rate_multiplier");
			pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) = high_transition_rate;
		}
		else 
		{
			//multiplier implementation 
			//pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) *= 0.4;
			//std::cout << "Rate down! " << std::endl;


			//switch implementation 
			pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index) = PhysiCell::parameters.doubles("base_transition_rate");
		}

		// Update Migration
		if(pCell->boolean_network.get_node_value("Migration"))
		{ 
			// pCell->phenotype.motility.is_motile = true;
		 	// pCell->phenotype.motility.migration_speed = PhysiCell::parameters.doubles("migration_speed");
			// pCell->phenotype.motility.migration_bias = PhysiCell::parameters.doubles("migration_bias");
			// pCell->phenotype.motility.persistence_time = PhysiCell::parameters.doubles("persistence");

			// pCell->evolve_coef(pCell->boolean_network.get_node_value("Migration"),	phenotype.motility.migration_speed, dt 
			// );
			// pCell->phenotype.motility.migration_speed = PhysiCell::parameters.doubles("max_motility_speed") * migration_coeff
		}
		else 
		{
			pCell->phenotype.motility.is_motile = false;
		}


		// Update Invasion
		if(pCell->boolean_network.get_node_value("Invasion"))
		{
			// nothing happens for now 
		}

	}

	pCell->set_internal_uptake_constants(dt);
}


void boolean_model_interface_main (Cell* pCell, Phenotype& phenotype, double dt){
    static int index_next_physiboss_run = pCell->custom_data.find_variable_index("next_physiboss_run");

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	if (PhysiCell_globals.current_time >= pCell->custom_data.variables.at(index_next_physiboss_run).value)
	{
		set_input_nodes(pCell);

		pCell->boolean_network.run_maboss();
		// Get noisy step size
		double next_run_in = pCell->boolean_network.get_time_to_update();
		pCell->custom_data.variables.at(index_next_physiboss_run).value = PhysiCell_globals.current_time + next_run_in;
		
		update_custom_variables(pCell);

		from_nodes_to_cell(pCell, phenotype, dt);
	}
}