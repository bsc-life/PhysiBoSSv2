
#include "./custom.h"

/**
 *	\main gastric custom
 *	\brief Custom module file for gastric example
 * 
 *	\details Modules needed for the gastric example. This custom module can be used to study the inhibition of AGS cell line with AKT, beta-catenin and TAK inhibitors.
 *
 *
 *	\date 17/03/21
 *	\author Arnau Montagud, BSC-CNS, with code previously developed by Annika Meert, Gerard Pradas and Miguel Ponce de Leon, BSC-CNS
 */

// declare cell definitions here 

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling;
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.functions.set_orientation = NULL;

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	initialize_cell_definitions_from_pugixml();

	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// set molecular properties for each drug defined in PhysiCell_settings.xml
	// starts with the second density because the first density is oxygen
	for (int i = 0; i < microenvironment.number_of_densities(); i++) 
	{
		if (microenvironment.density_names[i] != "oxygen") 
		{
			cell_defaults.phenotype.molecular.fraction_released_at_death[i] = 0.0;
		}
	}

	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 

	return; 
}

double get_decay_rate(double half_life)
{
	// natural logarithm of 2 / half-life = k (reaction rate coefficient)
	double decay_rate = log(2)/(half_life);
	return decay_rate;
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}	

	// set intial conditions and dirichlet boundary conditions for the drugs; vector already contains the condition for oxygen
	// TODO add oxygen_condition to XML if possible
	double oxygen_condition = 160.0;
	vector<double> condition_vector = {oxygen_condition};
	vector<bool> activation_vector {1};
	vector<double> decay_vector {0.1};
	for (int i = 0; i < microenvironment.number_of_densities(); i++)
	{
		std::string drug_name = microenvironment.density_names[i];
		if (drug_name != "oxygen") 
		{
			int total_drug_levels = parameters.ints("total_concentration_levels");
			if ( total_drug_levels == 0 ) 
			{
				double drug_concentration= parameters.doubles("concentration_" + drug_name);
				condition_vector.push_back(drug_concentration);
				activation_vector.push_back(1);
			}
			else 
			{
				int current_drug_level= parameters.ints("current_concentration_level_" + drug_name);
				// int total_drug_levels = parameters.ints("total_concentration_levels");
				string cell_line = parameters.strings("cell_line");
				int simulation_mode = parameters.ints("simulation_mode");
				double drug_concentration = get_drug_concentration_from_level_initial(cell_line, drug_name, current_drug_level, total_drug_levels, simulation_mode);
				condition_vector.push_back(drug_concentration);
				activation_vector.push_back(1);
			}
			// double half_life = get_value(half_lives, drug_name);
			// double decay_rate = get_decay_rate(half_life);
			// decay_vector.push_back(decay_rate);
		}
	}
	default_microenvironment_options.Dirichlet_activation_vector = activation_vector;
	default_microenvironment_options.Dirichlet_condition_vector = condition_vector;
	default_microenvironment_options.Dirichlet_zmax_values = condition_vector;
	default_microenvironment_options.Dirichlet_zmin_values = condition_vector;
	default_microenvironment_options.Dirichlet_ymax_values = condition_vector;
	default_microenvironment_options.Dirichlet_ymin_values = condition_vector;
	default_microenvironment_options.Dirichlet_xmax_values = condition_vector;
	default_microenvironment_options.Dirichlet_xmin_values = condition_vector;

	// initialize BioFVM 
	initialize_microenvironment();

	// microenvironment.decay_rates = decay_vector;
	
	return; 
}

void setup_tissue( void )
{
	Cell* pC;

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	std::string bnd_file = PhysiCell::parameters.strings("bnd_file");
	std::string cfg_file = PhysiCell::parameters.strings("cfg_file");
	BooleanNetwork prostate_network;
	double maboss_time_step = PhysiCell::parameters.doubles("maboss_time_step");
	prostate_network.initialize_boolean_network(bnd_file, cfg_file, maboss_time_step);

	for (int i = 0; i < cells.size(); i++)
	{
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		float radius = cells[i].radius;
		int phase = cells[i].phase;
		double elapsed_time = cells[i].elapsed_time;

		double random_num_1 = (double) rand()/RAND_MAX;
		double random_num_2 = (double) rand()/RAND_MAX;

		// double drug_sim = microenvironment.number_of_densities();
		// if ( drug_sim == 1 )
		if ( microenvironment.number_of_densities() == 1 )
		{
			pC = create_cell(get_cell_definition("default"));
		}
		else
		{
			if (PhysiCell::parameters.ints("simulation_mode") == 0)
			{
				// single inhibition - just one drug is present 
				std::string drug_name1 = microenvironment.density_names[1];
				double eval1 = PhysiCell::parameters.doubles( "prop_drug_resistant_" + drug_name1);
				if (random_num_1 < eval1)
				// if (random_num_1 < PhysiCell::parameters.doubles("prop_drug_resistant_" + microenvironment.density_names[1]))
				{
					// cell is NOT sensitive to the drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_resistant"));
				}
				else 
				{
					// cell is sensitive to the drug
					pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_sensitive"));
				}
			}
			else if (PhysiCell::parameters.ints("simulation_mode") == 1)
			{
				// double inhibition - two drugs are present - we have 4 cell strains 
				std::string drug_name1 = microenvironment.density_names[1];
				double eval1 = PhysiCell::parameters.doubles( "prop_drug_resistant_" + drug_name1);
				std::string drug_name2 = microenvironment.density_names[2];
				double eval2 = PhysiCell::parameters.doubles( "prop_drug_resistant_" + drug_name2);
				// if (random_num_1 < PhysiCell::parameters.doubles("prop_drug_resistant_" + microenvironment.density_names[1]))
				if (random_num_1 < eval1)
				{
					// TODO: check that these gates makes sense. Right now drug1 rules: if a cell is drug1_sensitive it will also be drug2_sensitive. There is no way that a cell is drug1_sensitive and drug2_resistant. (see boolean_model_interface.cpp, line 91). Also, in this line, check that the cell definitions make sense.
					// Maybe change names to double_resistant, drug1 or 2_resistant, and double sensitive
					if (random_num_2 < eval2)
					{
						// cell is resistant to both drugs
						pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_resistant"));
					}
					else 
					{
						// cell is only resistant to the first drug
						pC = create_cell(get_cell_definition(microenvironment.density_names[2] + "_sensitive"));
					}
				}
				else
				{
					if (random_num_2 < eval2)
					{
						// cell is only resistant to the second drug
						pC = create_cell(get_cell_definition(microenvironment.density_names[2] + "_resistant"));
					}
					else
					{
						// cell is resistant to no drug
						pC = create_cell(get_cell_definition(microenvironment.density_names[1] + "_sensitive"));
					}
					
				}
				
			}
			else
			{
				pC = create_cell(get_cell_definition("default"));
			}
		}
		pC->assign_position( x, y, z );
		// pC->set_total_volume(sphere_volume_from_radius(radius));
		
		// pC->phenotype.cycle.data.current_phase_index = phase;
		pC->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;

		pC->boolean_network = prostate_network;
		pC->boolean_network.restart_nodes();
		static int index_next_physiboss_run = pC->custom_data.find_variable_index("next_physiboss_run");
		pC->custom_data.variables.at(index_next_physiboss_run).value = pC->boolean_network.get_time_to_update();
		update_custom_variables(pC);
	}

	return; 
}

// custom cell phenotype function to run PhysiBoSS when is needed
void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if ( pCell->phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model) 
	// {
		// std::cout << pCell->phenotype.cycle.current_phase().name << " 0,0: " << pCell->phenotype.cycle.data.transition_rate(0, 0) << "\n" << std::endl;
	// }
	// if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
	// {
	// 	std::cout << pCell->phenotype.cycle.current_phase().name << " 0,1: " << pCell->phenotype.cycle.data.transition_rate(0, 1) << " / " << pCell->phenotype.cycle.current_phase().name << " 1,0: " << pCell->phenotype.cycle.data.transition_rate(1, 0)  << "\n" << std::endl;
	// }

	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
	boolean_model_interface_main (pCell, phenotype, dt);

	// static int index_next_physiboss_run = pCell->custom_data.find_variable_index("next_physiboss_run");

	// if( phenotype.death.dead == true )
	// {
	// 	pCell->functions.update_phenotype = NULL;
	// 	return;
	// }

	// if (PhysiCell_globals.current_time >= pCell->custom_data.variables.at(index_next_physiboss_run).value)
	// {
	// 	set_input_nodes(pCell);

	// 	pCell->boolean_network.run_maboss();
	// 	// Get noisy step size
	// 	double next_run_in = pCell->boolean_network.get_time_to_update();
	// 	pCell->custom_data.variables.at(index_next_physiboss_run).value = PhysiCell_globals.current_time + next_run_in;
		
	// 	update_custom_variables(pCell);

	// 	from_nodes_to_cell(pCell, phenotype, dt);
	// }
}

std::vector<std::string> prolif_apoptosis_coloring( Cell* pCell )
{
	// start with live coloring 
	// std::vector<std::string> output = false_cell_coloring_live_dead(pCell); 
	// return output; 
	std::vector<std::string> output;
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptosis_death_model)
	{
		//apoptotic cells are colored red
		output = {"crimson", "black","darkred", "darkred"};
	}

	else if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrosis_death_model)
	{
		//necrotic cells are colored brown
		output = {"peru", "black","saddlebrown", "saddlebrown"};
	}

	else if (PhysiCell::parameters.ints("simulation_mode") == 0) 
	{
		std::string drug_name = microenvironment.density_names[1];
		if (pCell->type_name == drug_name + "_sensitive")
		{
			//drug sensitive living cells are colored blue
			output = {"deepskyblue", "black", "darkblue", "darkblue"};
		} 
		else 
		{
			//drug resistant living cells are colored green
			output = {"limegreen", "black", "darkgreen", "darkgreen"};
		}
	}
	else if (PhysiCell::parameters.ints("simulation_mode") == 1) 
	{
		// // color living cells just in one color 
		// output = {"limegreen", "black", "darkgreen", "darkgreen"};

		// In case we want to color all 4 strains differently:
		// double inhibitions --> 4 cell strains 
		std::string drug1_name = microenvironment.density_names[1];
		std::string drug2_name = microenvironment.density_names[2]; 
		if (pCell->type_name == drug1_name + "_sensitive")
		{
			//cells that are sensitive to both drugs are colored blue
			output = {"deepskyblue", "black", "darkblue", "darkblue"};
		}
		else if (pCell->type_name == drug2_name + "_resistant")
		{
			// cells that are just sensitive to the first drug 
			output = {"limegreen", "black", "darkgreen", "darkgreen"};

		}
		else if (pCell->type_name == drug2_name + "_sensitive")
		{
			// cells are just sensitive to the second drug
			output = {"gold", "black", "orange", "orange"};
		}
		else 
		{
			// cells aren't sensitive to any drug
			output = {"mediumorchid", "black", "mediumpurple", "mediumpurple"};
		}
	}
	else 
	{
		// no drug simulation - living cells are colored green
		output = {"limegreen", "black", "darkgreen", "darkgreen"};
	}
	return output;

}


// void set_boolean_node (Cell* pCell, std::string drug_name, int index, double threshold) {
// 	if (index != -1)
// 		{
// 			std::string node_name = "anti_" + drug_name;
// 			 double drug_conc = pCell->phenotype.molecular.internalized_total_substrates[index];
// 			if (drug_conc > threshold) 
// 			{
			
// 				pCell->boolean_network.set_node_value(node_name, 1);
// 			}
// 			else 
// 			{
// 				pCell->boolean_network.set_node_value(node_name, 0);
// 			}
			
// 		}
// }

// void set_input_nodes(Cell* pCell) {
// 	if (PhysiCell::parameters.ints("simulation_mode") == 0)
// 	{	
// 		// single inhibition - just one drug is present 

// 		std::string drug_name = microenvironment.density_names[1];
// 		int drug_index = microenvironment.find_density_index(drug_name);
// 		int drug_threshold = PhysiCell::parameters.doubles( "threshold_" + drug_name);

// 		if (pCell->type_name == drug_name + "_sensitive")
// 		{
// 			// cell is sensitive to the drug -> set boolean node
// 			set_boolean_node(pCell, drug_name, drug_index, drug_threshold);
// 		}
// 	}
// 	else if (PhysiCell::parameters.ints("simulation_mode") == 1)
// 	{	
// 		// double inhibition - two drugs are present - we have 4 cell strains 

// 		std::string drug1_name = microenvironment.density_names[1];
// 		int drug1_index = microenvironment.find_density_index(drug1_name);
// 		int drug1_threshold = PhysiCell::parameters.doubles( "threshold_" + drug1_name);

// 		std::string drug2_name = microenvironment.density_names[2];
// 		int drug2_index = microenvironment.find_density_index(drug2_name);
// 		int drug2_threshold = PhysiCell::parameters.doubles( "threshold_" + drug2_name);
	
// 		if (pCell->type_name == drug1_name + "_sensitive")
// 		{	
// 			// cell is sensitive to both drugs
// 			set_boolean_node(pCell, drug1_name, drug1_index, drug1_threshold);
// 			set_boolean_node(pCell, drug2_name, drug2_index, drug2_threshold);
// 		}
// 		else if (pCell->type_name == drug2_name + "_resistant")
// 		{
// 			// cell is only sensitive to the first drug
// 			set_boolean_node(pCell, drug1_name, drug1_index, drug1_threshold);
// 		}
// 		else if (pCell->type_name == drug2_name +  "_sensitive")
// 		{
// 			// cell is only sensitive to the second drug
// 			set_boolean_node(pCell, drug2_name, drug2_index, drug2_threshold);
// 		}

// 		// else: cell is sensitive to no drug --> no boolean node is set
// 	}

// }

// void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt)
// {
// 	std::vector<bool>* nodes = pCell->boolean_network.get_nodes();
// 	int bn_index;

// 	// For AGS model

// 	std::string prosurvival_basename = "Prosurvival_b";
// 	std::string antisurvival_basename = "Antisurvival_b";
// 	double prosurvival_value = 0.0;
// 	double antisurvival_value = 0.0;

// 	for(int i=1; i<=3; i++)
// 	{
// 		bn_index = pCell->boolean_network.get_node_index( prosurvival_basename + std::to_string(i) );
// 		if ( (*nodes)[bn_index] > 0)
// 		{
// 			prosurvival_value += (*nodes)[bn_index];
// 		}
// 		bn_index = pCell->boolean_network.get_node_index( antisurvival_basename + std::to_string(i) );
// 		if ( (*nodes)[bn_index] > 0)
// 		{
// 			antisurvival_value += (*nodes)[bn_index];
// 		}
// 	}

// 	static int start_phase_index; // Q_phase_index; 
// 	static int end_phase_index; // K_phase_index;
// 	static int apoptosis_index; 
// 	double multiplier = 1.0;

// 	// live model 
			
// 	if( pCell->phenotype.cycle.model().code == PhysiCell_constants::live_cells_cycle_model )
// 	{
// 		start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
// 		apoptosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
// 		end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );

// 		multiplier = ( prosurvival_value + 1 ) / ( antisurvival_value + 1 ) ; //[0.25, 0.33, 0.5, 0.58, 0.66, 0.75, 1, 1.5, 2, 3, 4]
// 		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier *	phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
// 		pCell->phenotype.death.rates[apoptosis_index] = ( 1 / multiplier ) * pCell->phenotype.death.rates[apoptosis_index];

// 	}

// 	if( phenotype.cycle.model().code == PhysiCell_constants::advanced_Ki67_cycle_model || phenotype.cycle.model().code == PhysiCell_constants::basic_Ki67_cycle_model )
// 	{
// 		// std::cout << pCell->phenotype.cycle.current_phase().name << " 0,1: " << pCell->phenotype.cycle.data.transition_rate(0, 1) <<  ". " << prosurvival_value << ", " << antisurvival_value << " / " << pCell->phenotype.cycle.current_phase().name << " 1,0: " << pCell->phenotype.cycle.data.transition_rate(1, 0) <<  ". " << prosurvival_value << ", " << antisurvival_value << "\n" << std::endl;

// 	}

// 	// int tnf_substrate_index = microenvironment.find_density_index( "tnf" );
// 	// static double tnf_secretion = parameters.doubles("tnf_secretion_rate");
// 	// double tnf_secretion_rate = 0;
// 	// // produce some TNF
// 	// if ( pCell->boolean_network.get_node_value( "NFkB" ) )
// 	// {
// 	// 	tnf_secretion_rate = (tnf_secretion / microenvironment.voxels(pCell->get_current_voxel_index()).volume);

// 	// }
// 	// pCell->phenotype.secretion.secretion_rates[tnf_substrate_index] = tnf_secretion_rate;
// 	pCell->set_internal_uptake_constants(dt);
// }


// ***********************************************************
// * NOTE: Funtion to read init files created with PhysiBoSS *
// ***********************************************************
std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header) 
{ 
	// File pointer 
	std::fstream fin; 
	std::vector<init_record> result;

	// Open an existing file 
	fin.open(filename, std::ios::in); 

	// Read the Data from the file 
	// as String Vector 
	std::vector<std::string> row; 
	std::string line, word;

	if(header)
		getline(fin, line);

	do 
	{
		row.clear(); 

		// read an entire row and 
		// store it in a string variable 'line' 
		getline(fin, line);

		// used for breaking words 
		std::stringstream s(line); 

		// read every column data of a row and 
		// store it in a string variable, 'word' 
		while (getline(s, word, delimiter)) { 

			// add all the column data 
			// of a row to a vector 
			row.push_back(word); 
		}

		init_record record;
		record.x = std::stof(row[2]);
		record.y = std::stof(row[3]);
		record.z = std::stof(row[4]);
		record.radius = std::stof(row[5]);
		record.phase = std::stoi(row[13]);
		record.elapsed_time = std::stod(row[14]);

		result.push_back(record);
	} while (!fin.eof());
	
	return result;
}