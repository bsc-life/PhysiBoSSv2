
#include "./drug_sensitivity.h"

/**
 *	\main prostate custom
 *	\brief Custom module file for prostate example
 * 
 *	\details Modules needed for the prostate example. 
 *
 *
 *	\date 14/12/2020
 *	\author Annika Meert, BSC-CNS
 */

using namespace std;

std::vector<std::pair<std::string, std::vector<int>>> read_csv(std::string filename){
    // Reads a CSV file into a vector of <string, vector<int>> pairs where
    // each pair represents <column name, column values>

    // Create a vector of <string, int vector> pairs to store the result
    std::vector<std::pair<std::string, std::vector<int>>> result;

    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line, colname;
    int val;

    // Read the column names
    if(myFile.good())
    {
        // Extract the first line in the file
        std::getline(myFile, line);

        // Create a stringstream from line
        std::stringstream ss(line);

        // Extract each column name
        while(std::getline(ss, colname, ',')){
            
            // Initialize and add <colname, int vector> pairs to result
            result.push_back({colname, std::vector<int> {}});
        }
    }

    // Read data, line by line
    while(std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);
        
        // Keep track of the current column index
        int colIdx = 0;
        
        // Extract each integer
        while(ss >> val){
            
            // Add the current integer to the 'colIdx' column's values vector
            result.at(colIdx).second.push_back(val);
            
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

std::vector<std::pair<std::string, int>> cell_line_identifiers = {
    { "LNCaP", 907788},
    { "BPH1", 924105},
    { "DU145", 905935},
    { "22RV1", 924100},
    { "VCaP", 1299075},
    { "PC3", 905934}
};

std::vector<std::pair<std::string, int>> drug_identifiers = {
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

// int get_drug_index(int drug_identifier, std::vector<std::pair<std::string, std::vector<int>>> csv_file) {

// }

// int get_cell_line_index (int cell_line_identifier, csv_file) {

// }

double get_x_from_conc(double x_conc, double max_conc) {
    double x = (log (x_conc / max_conc) / log (2) ) + 9;
    return(x);
}

double get_conc_from_x (double x, double max_conc) {
    double x_conc = pow(max_conc * 2 , (x - 9));
    return(x_conc);
}

double get_lx_from_x (double x, double max_conc) {
    double lx = log (get_conc_from_x (x, max_conc));
    return(lx);
} 

// two parameter logistic function used to model GDSC dose response data 
// this model is described further in published in Vis, D.J. et al Pharmacogenomics 2016, 17(7):691-700) 
double logist3 (double x, double xmid, double scale) {
    //https://github.com/CancerRxGene/gdscIC50/blob/master/R/nlme_fit.R
    return;
}

double linear_mixed_model_function(double lx, double max_conc, double xmid, double scale) 
{
    double x = get_x_from_conc( exp(lx), max_conc);
    double y_hat = 1 - logist3(x, xmid, scale);
    return(y_hat);
}



