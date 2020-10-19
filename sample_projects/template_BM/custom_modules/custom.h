
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

#include "../addons/PhysiBoSS/src/boolean_network.h"

using namespace BioFVM; 
using namespace PhysiCell;

struct init_record
{
	float x;
	float y;
	float z;
	float radius;
	int phase;
	double elapsed_time;
};

// setup functions to help us along 
void create_cell_types( void );
void setup_tissue( void ); 

// set up the BioFVM microenvironment 
void setup_microenvironment( void ); 

// custom pathology coloring function 
std::vector<std::string> my_coloring_function( Cell* );

// custom cell phenotype functions could go here 
void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt );

void update_custom_variables( Cell* pCell );

void set_input_nodes(Cell* pCell); 
void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt);
void do_proliferation( Cell* pCell, Phenotype& phenotype, double dt );

std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header);

inline float sphere_volume_from_radius(float radius) {return 4/3 * PhysiCell_constants::pi * std::pow(radius, 3);}
