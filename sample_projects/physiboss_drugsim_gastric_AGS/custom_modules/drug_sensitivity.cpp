
#include "./drug_sensitivity.h"

/**
 *	\drug sensitivity
 *	\brief drug sensitivity module for prostate example
 * 
 *	\details Module needed for the prostate example. 
 *
 *
 *	\date 19/03/2021
 *	\author Annika Meert, BSC-CNS, with some additions from Arnau Montagud for the AGS-specific cell line and drugs
 */

using namespace std;

vector<pair<string, vector<double>>> read_csv(string filename){
    // Reads a CSV file into a vector of <string, vector<int>> pairs where
    // each pair represents <column name, column values>

    // Create a vector of <string, int vector> pairs to store the result
    vector<pair<string, vector<double>>> result;

    // Create an input filestream
    ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw runtime_error("Could not open file");

    // Helper vars
    string line, colname;
    double val;

    // Read the column names
    if(myFile.good())
    {
        // Extract the first line in the file
        getline(myFile, line);

        // Create a stringstream from line
        stringstream ss(line);

        // Extract each column name
        while(getline(ss, colname, ',')){
            
            // Initialize and add <colname, double vector> pairs to result
            result.push_back({colname, vector<double> {}});
        }
    }

    // Read data, line by line
    while(getline(myFile, line))
    {
        // Create a stringstream of the current line
        stringstream ss(line);
        
        // Keep track of the current column index
        int colIdx = 0;
        string token;
        
        // Extract each double
        while(std::getline(ss, token, ',')){
            
            // Add the current double to the 'colIdx' column's values vector
            double value;
            try
            {
                value = std::stod(token);
            }
            catch(std::exception& e)
            {
                colIdx++;
                continue;
            }
            result.at(colIdx).second.push_back(value);
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
            
            // Increment the column index
            colIdx++;
           
        }
    }

    // Close file
    myFile.close();

    return result;
}

const vector<pair<string, int>> cell_line_ids = {
    { "LNCaP", 907788},
    { "BPH1", 924105},
    { "DU145", 905935},
    { "22Rv1", 924100},
    { "VCaP", 1299075},
    { "PC3", 905934},
    { "AGS",906790}
};

const vector<pair<string, int>> drug_ids = {
    { "Ipatasertib", 1924},
    { "Afuresertib", 1912},
    { "Afatinib", 1032},
    { "Erlotinib", 1168},
    { "Ulixertinib", 2047},
    { "Luminespib", 1559},
    { "Trametinib", 1372},
    { "Selumetinib", 1736},
    { "Pictilisib", 1058},
    { "Alpelisib", 1560},
    { "BIBR1532", 2043},  //No growth data
    { "PI103", 302},
    { "PD0325901", 1060}, 
    { "AKT_inhibitor_VIII", 171},
    { "BIRB0796", 1042},
    { "CT99021", 1241},
    { "Oxozeaenol", 1242},
    { "PKF118", 999999}
};

const vector<pair<string, string>> drug_targets = {
    { "Ipatasertib", "AKT"},
    { "Afuresertib", "AKT"},
    { "Afatinib", "EGFR"},
    { "Erlotinib", "EGFR"},
    { "Ulixertinib", "ERK"},
    { "Luminespib", "HSPs"},
    { "Trametinib", "MEK1_2"},
    { "Selumetinib", "MEK1_2"},
    { "Pictilisib", "PI3K"},
    { "Alpelisib", "PI3K"},
    { "BIBR1532", "TERT"},
    { "PI103", "PI3K"},
    { "PD0325901", "MEK"},
    { "AKT_inhibitor_VIII", "AKT"},
    { "BIRB0796", "p38alpha"},
    { "CT99021", "GSK3"},
    { "Oxozeaenol", "TAK1"},
    { "PKF118", "betacatenin"}
};
// not used anymore
const vector<pair<string, int>> half_lives = {
    { "Ipatasertib", 2748},
    { "Afuresertib", 2448},
    { "Afatinib", 2220},
    { "Erlotinib", 2172},
    { "Ulixertinib", 105},
    { "Luminespib", 7200},
    { "Trametinib", 5760},
    { "Selumetinib", 822},
    { "Pictilisib", 1062},
    { "Alpelisib", 822},
    { "PI103", 2000},
    { "PD0325901", 2000},
    { "AKT_inhibitor_VIII", 2000},
    { "BIRB0796", 2000},
    { "CT99021", 2000},
    { "Oxozeaenol", 2000},
    { "PKF118", 2000}
};

// TODO: change to reading it from XML
const vector<pair<string, vector<double>>> csv_file = read_csv( "./gastric_drug_sensitivity2.csv");

string get_value (const vector<pair<string, string>> dict, string key) {
    vector< pair<string, string>>::const_iterator dict_iterator = find_if( dict.begin(), dict.end(),[&key](const pair < string, string>& element){ return element.first  == key;} );
    return (*dict_iterator).second;
}

int get_value (const vector<pair<string, int>> dict, string key) {
    vector< pair<string, int>>::const_iterator dict_iterator = find_if( dict.begin(), dict.end(),[&key](const pair < string, int>& element){ return element.first  == key;} );
    return (*dict_iterator).second;
}

vector<double> get_value (const vector<pair<string, vector<double>>> dict, string key) {
    vector< pair<string, vector<double>>>::const_iterator dict_iterator = find_if( dict.begin(), dict.end(),[&key](const pair < string, vector<double>>& element){ return element.first  == key;} );
    return (*dict_iterator).second;
}

// TODO get this from XML
int get_index(string drug_name, string cell_line_name) { 
    // retrieve the id for the drug and cell line
    int drug_identifier = get_value(drug_ids, drug_name);
    int cell_identifier = get_value(cell_line_ids, cell_line_name);
   
    // retrieve the index for the csv datastructure where both identifiers are met
    int index = 0;
    //static int index_CL = cell_defaults.custom_data.find_vector_variable_index("\"CL\"");
	//vector <double > cell_line_vector = cell_defaults.custom_data.vector_variables.at(index_CL).value;
    vector <double> cell_line_vector = get_value(csv_file, "\"CL\"");
    // static int index_drug = cell_defaults.custom_data.find_vector_variable_index("\"DRUG_ID_lib\"");
	// vector <double > drug_vector = cell_defaults.custom_data.vector_variables.at(index_drug).value;
    vector <double> drug_vector = get_value(csv_file, "\"DRUG_ID_lib\"");
    for (int i = 0; i < drug_vector.size(); i++) {
        if (cell_line_vector[i] == cell_identifier && drug_vector[i] == drug_identifier)
        {
            index = i;
            // cout << "Index found: drug was found in the csv file!" << endl;
            break;
        }
    }
  
    return index;
}

// TODO get this from XML if possible
vector<double> get_drug_sensitivity_values_initial (string drug_name, string cell_line_name) {
    int index = get_index(drug_name, cell_line_name);
    // static int index_max_conc = cell_defaults.custom_data.find_vector_variable_index("\"maxc\"");
    // static int index_xmid = cell_defaults.custom_data.find_vector_variable_index("\"xmid\"");
    // static int index_scale = cell_defaults.custom_data.find_vector_variable_index("\"scal\"");
    // vector <double > max_conc_vector = cell_defaults.custom_data.vector_variables.at(index_max_conc).value;
    // vector <double > xmid_vector = cell_defaults.custom_data.vector_variables.at(index_xmid).value;
    // vector <double > scale_vector = cell_defaults.custom_data.vector_variables.at(index_scale).value;
    vector <double> max_conc_vector = get_value(csv_file, "\"maxc\"");
    vector <double> xmid_vector = get_value(csv_file, "\"xmid\"");
    vector <double> scale_vector = get_value(csv_file, "\"scal\"");
    
    return {max_conc_vector[index], xmid_vector[index], scale_vector[index]};
}

// TODO get this from XML
vector<double> get_drug_sensitivity_values (Cell* pCell, string drug_name, string cell_line_name) {
    static int max_conc_vector_index = pCell->custom_data.find_variable_index(drug_name+"_maxc");
    static double  max_conc_vector = pCell->custom_data[max_conc_vector_index];
    static int xmid_vector_index = pCell->custom_data.find_variable_index(drug_name+"_xmid");
    static double  xmid_vector = pCell->custom_data[xmid_vector_index];
    static int scale_vector_index = pCell->custom_data.find_variable_index(drug_name+"_scal");
    static double  scale_vector = pCell->custom_data[scale_vector_index];
    return {max_conc_vector, xmid_vector, scale_vector};
}

// TODO get this from xml if possible
double get_drug_concentration_from_level_initial (string cell_line, string drug_name, int conc_level, int num_of_conc_levels, int simulation_mode) {
    // IC10 --> cell viability = 0.9, lowest drug concentration
    double highest_limit = 0.9;
    // IC90 --> cell viability = 0.1, highest drug concentration
    double lowest_limit = 0.1;

    // divide the range into the total number of levels
    double range_size = (highest_limit - lowest_limit) / (num_of_conc_levels - 1);
    double final_viability = highest_limit - (conc_level - 1) * range_size;
    
   // call functions to retrieve data from datastructure cell line and drug
    vector<double> drug_sens_vals = get_drug_sensitivity_values_initial(drug_name, cell_line);
    // call linear_mixed_model_function
    double max_conc = drug_sens_vals[0];
    double xmid = drug_sens_vals[1];
    double scale = drug_sens_vals[2];

    // get the drug concentration for the cell viability
    double x = get_x_for_cell_viability(xmid, scale, final_viability);
    double drug_conc = get_conc_from_x(x, max_conc);

    // if simulation mode is on double drugs half the drug concentration
    // if (simulation_mode == 1) {
    //     drug_conc = drug_conc / 2;
    // }
    return drug_conc;
}

// returns x: the concentration scaled for 9 different concentrations
double get_x_from_conc(double x_conc, double max_conc) {
    double x = (log (x_conc / max_conc) / log (2) ) + 9;
    return(x);
}

// returns conc: concentration in micromolar 
double get_conc_from_x (double x, double max_conc) {
    double x_conc = max_conc * pow(2 , (x - 9));
    return(x_conc);
}

// returns the natural logarithm for x
double get_lx_from_x (double x, double max_conc) {
    double lx = log (get_conc_from_x(x, max_conc));
    return(lx);
} 

// two parameter logistic function used to model GDSC dose response data 
// this model is described further in published in Vis, D.J. et al Pharmacogenomics 2016, 17(7):691-700) 
// double logist3 (double x, double xmid, double scale) {
//     //https://github.com/CancerRxGene/gdscIC50/blob/master/R/nlme_fit.R

//     return;
// }


// this function then has to be called whenever we translate concentration into maboss
// check in the datastructure for cell line and drug for the important values and enter concentration
// then we receive a y_hat that we can use to calculate the maboss activity
// TODO: remove max_conc from this function as it is not needed
double get_cell_viability_for_x(double x, double max_conc, double xmid, double scale) 
{
    // double y_hat = 1 - logist3(x, xmid, scale);
    double y_hat = 1/(1 + pow( exp(1), (x- xmid) /scale));
    return(y_hat);
}

double derivative_linear_mixed_model(double lx, double max_conc, double xmid, double scale){
    double x = get_x_from_conc( exp(lx), max_conc);
    double derivative_y_hat = pow(exp(1), (x-xmid) /scale) /  (scale * pow((pow (exp(1), (x-xmid) /scale) + 1), 2));
    return derivative_y_hat;
}

double get_x_for_cell_viability (double xmid, double scale, double cell_viability) {
    double x = (log((1/cell_viability) - 1) * scale) + xmid;
    return x;
}

// TODO remove index as it is not used
double get_cell_viability_for_drug_conc (Cell* pCell, double drug_conc, string cell_line, string drug_name, int index) 
{
    //std::cout << "Concentration of " << drug_name << ": " << drug_conc << std::endl;
    // call functions to retrieve data from datastructure cell line and drug
    vector<double> drug_sens_vals = get_drug_sensitivity_values(pCell, drug_name, cell_line);
    // call linear_mixed_model_function
    double max_conc = drug_sens_vals[0];
    double xmid = drug_sens_vals[1];
    double scale = drug_sens_vals[2];
    // transfor drug concentration into lx
    double x = get_x_from_conc(drug_conc, max_conc);
    double y_hat = get_cell_viability_for_x(x, max_conc, xmid, scale);
    // return cell viability value y_hat
    //std::cout << "Cell viability for " << drug_name << ": " << y_hat << std::endl;
    return y_hat;
}



