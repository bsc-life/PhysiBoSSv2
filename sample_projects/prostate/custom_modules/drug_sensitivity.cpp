
#include "./drug_sensitivity.h"

/**
 *	\drug sensitivity
 *	\brief drug sensitivity module for prostate example
 * 
 *	\details Module needed for the prostate example. 
 *
 *
 *	\date 15/02/2020
 *	\author Annika Meert, BSC-CNS
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
    { "PC3", 905934}
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
    { "BIBR1532", 2043}
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
    { "BIBR1532", "TERT"}
};

int get_index(Cell* pCell, string drug_name, string cell_line_name) { 
    // retrieve the identifier for the drug and cell line
    vector< pair<string, int>>::const_iterator drug_identifier_it = find_if( drug_ids.begin(), drug_ids.end(),[&drug_name](const pair < string, int>& element){ return element.first  == drug_name;} );
    vector< pair<string, int>>::const_iterator cell_identifier_it = find_if( cell_line_ids.begin(), cell_line_ids.end(),[&cell_line_name](const pair<string, int>& element){ return element.first == cell_line_name;} );

    // retrieve the index for the custom vector data where both identifiers are met
    int index = 0;
    static int index_CL = pCell->custom_data.find_variable_index("\"CL\"");
	vector <double > cell_line_vector = pCell->custom_data.vector_variables.at(index_CL).value;
    static int index_drug = pCell->custom_data.find_variable_index("\"DRUG_ID_lib\"");
	vector <double > drug_vector = pCell->custom_data.vector_variables.at(index_drug).value;
    for (int i = 0; i < drug_vector.size(); i++) {
        if (cell_line_vector[i] == (*cell_identifier_it).second && drug_vector[i] == (*drug_identifier_it).second)
        {
            index = i;
            break;
        }
    }

    return index;
}


vector<double> get_drug_sensitivity_values (Cell* pCell, string drug_name, string cell_line_name) {
    int index = get_index(pCell, drug_name, cell_line_name);
    static int index_max_conc = pCell->custom_data.find_variable_index("\"maxc\"");
    static int index_xmid = pCell->custom_data.find_variable_index("\"xmid\"");
    static int index_scale = pCell->custom_data.find_variable_index("\"scal\"");

	vector <double > max_conc_vector = pCell->custom_data.vector_variables.at(index_max_conc).value;
    vector <double > xmid_vector = pCell->custom_data.vector_variables.at(index_xmid).value;
    vector <double > scale_vector = pCell->custom_data.vector_variables.at(index_scale).value;
    
    return {max_conc_vector[index], xmid_vector[index], scale_vector[index]};
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
double linear_mixed_model_function(double lx, double max_conc, double xmid, double scale) 
{
    double x = get_x_from_conc( exp(lx), max_conc);
    // double y_hat = 1 - logist3(x, xmid, scale);
    double y_hat = 1/(1 + pow( exp(1), (x- xmid) /scale));
    return(y_hat);
}

// drug concentration is entered in micromolar
double get_cell_viability_for_drug_conc (Cell* pCell, string cell_line, string drug_name, int index) 
{
    // get internalized substrate concentration
    double drug_conc = pCell->phenotype.molecular.internalized_total_substrates[index];
    // call functions to retrieve data from datastructure cell line and drug
    vector<double> drug_sens_vals = get_drug_sensitivity_values(pCell, drug_name, cell_line);
    // call linear_mixed_model_function
    double max_conc = drug_sens_vals[0];
    double xmid = drug_sens_vals[1];
    double scale = drug_sens_vals[2];
    // transfor drug concentration into lx
    double x = get_x_from_conc(drug_conc, max_conc);
    double lx = get_lx_from_x(x, max_conc);
    double y_hat = linear_mixed_model_function(lx, max_conc, xmid, scale);
    // return cell viability value y_hat
    return 1 - y_hat;
}



