#!/usr/bin/env python3
# coding: utf-8

# Tool: physiboss_drugsim.py
# Date: January 2021
# Author: Annika Meert
# E-mail: annikameert@protonmail.com

#%% Imports 
import os, sys
import argparse
import re
import shutil
from itertools import combinations, product
from lxml import etree

current_wd = os.getcwd()
arg = sys.argv
parser = argparse.ArgumentParser()

parser.add_argument("-p", "--project", required=True, help="Name of project that drug simulations are based on (ex. 'prostate').")
parser.add_argument("-d", "--drugs", required=True, help="Names of nodes affected by drug, comma separated (ex. 'MYC_MAX, ERK').")
parser.add_argument("-c", "--drug_conc", default="0, 0.2, 0.4, 0.6, 0.8, 1.0", help="Levels of drug inhibition, between 0 and 1, comma separated (ex: '0, 0.2, 0.4, 0.6, 0.8, 1.0').")
parser.add_argument("-m", "--mode", default="both", choices=['single', 'double', 'both'], help="Mode of simulation for drug inhibition: single, double or both.")
parser.add_argument("-i", "--input_cond", default='00', nargs='?', choices=['00', 'AR', 'AR_EGF', 'EGF'], help="Initial condition for drug simulation.")
parser.add_argument("-cl", "--cluster", default=False ,type=bool, help="Use of cluster or not.") 
# example: physiboss_drugsim.py -p prostate -d "MYC_MAX, ERK" -c "0.2, 0.8" -m "single" -i "00" -cl yes

args = parser.parse_args()
#%% Process Arguments

####################################################################
# Process arguments
####################################################################

print("Arguments:\n")
print(args)
print("\n")

project=args.project
# check if the project is in sample projects 
if not os.path.isdir("sample_projects/" + project):
    print(project + " folder is not found")
    sys.exit(1)

drugs=args.drugs
node_list = str(drugs.replace (" ", "").replace(",", ", ")).split(", ")
print("Nodes to be inhibited: "+str(node_list).replace("[", "").replace("]","").replace("'",""))

# specify boolean model path
input_cond = args.input_cond
bool_model_path_dir = "{}/{}/{}/{}".format("sample_projects", project, "config", "boolean_network")
bool_model_filename = "LNCaP_mut_RNA_00"
bool_model = "{}/{}".format(bool_model_path_dir, bool_model_filename)


node_list_1 = [x.split(", ") for x in node_list]

# create two list containing the value and the name of the chosen drug concentrations 
drug_conc = args.drug_conc
drug_conc_value_list = [float(i) for i in drug_conc.replace(" ","").replace(",",", ").split(", ")]
drug_conc_name_list = [str(i).replace(".", "_") for i in drug_conc_value_list]
print("Inhibited nodes' levels: "+str(drug_conc_value_list).replace("[","").replace("]",""))

# set base paths for output and project folders 
sample_project_path = "sample_projects"
prostate_path = "{}/{}".format(sample_project_path, project)
project_path = "{}/{}_{}_{}".format(sample_project_path, "physiboss_drugsim", project, "LNCaP") 
single_run_path = "{}/{}".format(project_path, "single_runs")
double_run_path = "{}/{}".format(project_path, "double_runs")


####################################################################
# Function definitions
####################################################################

# adds a drug to a network by modifying and adding corresponding nodes in .bnd and .cfg files 
def add_drugs_to_network(bool_model, druglist):
    # creating new .bnd file containing modified nodes and anti-nodes for all drugs
    bnd_file = "{}.{}".format(bool_model,"bnd")
    bnd_file_content = open(bnd_file).read()

    bnd_file_content = bnd_file_content.replace(" {", "{")
    bnd_file_content = bnd_file_content.replace("\n{", "{")
    bnd_file_content = bnd_file_content.replace("{", "\n{")

    bnd_nodes = bnd_file_content.split("Node")[1:]
    for item_index, item in enumerate(bnd_nodes):
        bnd_nodes[item_index] = "Node" + item 

    # filter nodes that are selected in the node_list (drugs)
    bnd_nodes_list = [re.findall(r'Node.*',line) for line in open(str(bnd_file))]
    bnd_nodes_list = list(filter(None, bnd_nodes_list))
    bnd_nodes_list = [[w.replace("Node ", "").replace(" {", "") for w in line] for line in bnd_nodes_list]

    # check if all selected drugs are present in the .bnd file 
    drugs_missing_bnd = [item for item in node_list_1 if item not in bnd_nodes_list]
    if len(drugs_missing_bnd) > 0:
        print("Nodes NOT in bnd file: " + str(drugs_missing_bnd).replace("[", "").replace("]", "").replace("'",""))

    new_nodes = str()
    for drug in druglist:
        print("Modifying bnd file: Adding inhibitor node: " + str(drug))
        string = "Node {}\n".format(drug)
        item_count = 0
        for item_index, item in enumerate(bnd_nodes):
            if string in (bnd_nodes[item_index]):
                print("{} found".format(drug))
                item_count += 1
                for line_index, line in enumerate(bnd_nodes[item_index].split("\n")):
                    if " logic" in line:
                        new_nodes += line.split(";")[0].split(" = ")[0] + " = (" + line.split(";")[0].split(" = ")[1] + ") AND NOT anti_{};\n".format(drug)
                    else:
                        new_nodes += line + "\n"
                bnd_nodes.pop(item_index)
        if item_count == 0:
            print("{} not found".format(drug))
    
    for last_item in bnd_nodes:
        new_nodes += last_item 
    
    for inhibitor_to_add in druglist:
        new_nodes += "Node anti_{}".format(inhibitor_to_add) + " \n{\n" + "\tlogic = (anti_{});\n".format(inhibitor_to_add) + "\trate_up = @logic ? 0 : 0;\n\trate_down = @logic ? 0 : 0;\n}\n\n\n"

    new_bnd_name = "{}_{}.{}".format(bool_model, "all_drugs", "bnd")
    fw1 = open(new_bnd_name,"w")
    fw1.write(new_nodes)
    fw1.close()

    # create new .cfg file setting anti-nodes on 0 and drug-nodes on random 
    cfg_file = "{}.{}".format(bool_model,"cfg")
    cfg_file_content = open(cfg_file).read()

    cfg_file_content = cfg_file_content.replace("[1]", " [1]").replace("[0]", " [0]").replace("] ,", "],").replace("=", " = ").replace("  ", " ")
    new_cfg = str()

    for line_index, line in enumerate(cfg_file_content.split("\n")):
        if ".istate = TRUE" in line:
            line = "[" + line.replace(".istate = TRUE", "].istate = 1 [1], 0 [0]")
        elif ".istate = FALSE" in line:
            line = "[" + line.replace(".istate = FALSE", "].istate = 0 [1], 1 [0]")
        elif "].istate = 0 [0], 1 [1]" in line:
            line = line.replace("].istate = 0 [0], 1 [1]", "].istate = 1 [1], 0 [0]")
        elif "].istate = 1 [0], 0 [1]" in line:
            line = line.replace("].istate = 1 [0], 0 [1]", "].istate = 0 [1], 1 [0]")
        else:
            line = line
        new_cfg += line + "\n"  

    b1 = [re.findall(r'.*istate.*',line) for line in open(str(bool_model+".cfg"))]
    b2 = list(filter(None, b1))
    cfg_nodes = [[re.sub(r"].istate.*","",w) for w in line] for line in b2]
    cfg_nodes = [[re.sub(r"\[","",w) for w in line] for line in cfg_nodes]
    missing_nodes_cfg = [item for item in bnd_nodes_list if item not in cfg_nodes]
    if (len(missing_nodes_cfg) > 0):
        missing_line = str([re.sub(r"$",".istate = 0.5 [0] , 0.5 [1];",str(line)) for line in missing_nodes_cfg]).replace("'","").replace('"[',"").replace('"]',"").replace(";",";\n").replace('", ','[')
        new_cfg += missing_line
    drugs_in_cfg = [item for item in druglist if [item] in cfg_nodes]

            
    name_cfg = "{}_{}.cfg".format(bool_model,"all_drugs")
    fw2 = open(name_cfg, "w")
    fw2.write(new_cfg)
    fw2.close()

    for drug in druglist:
        f = open(name_cfg, "r")
        lines = f.readlines()
        f.close()
        nf = open(name_cfg, "w")
        for line in lines:
            if "[{}].istate = ".format(drug) in line:
                new_line = line
                new_line = new_line = "[{}].istate = {} [1], {} [0];\n[anti_{}].istate = {} [1], {} [0];\n".format(drug, "0.5", "0.5", drug, "0", "1")
                nf.write(str(new_line))
            else:
                nf.write(str(line))
        nf.close()



# adds a drug to a physicell xml file
def add_drug_to_xml(drug, conc, path_to_xml, config_path, model_name, mode, output_path):
    parser = etree.XMLParser(remove_blank_text=True)
    root = etree.parse(path_to_xml, parser).getroot()

    # set the output directory for the current run
    save = root.find("save")
    output_folder = save.find("folder")
    output_folder.text(output_path)

    # insert drug as a density in the microenvironment 
    # 1. childs 
    microenv_setup = root.find('microenvironment_setup')

    # 2. childs
    variables = microenv_setup.findall('variable')
    variable_count = 0
    for variable in variables:
        variable_count += 1
    new_variable = etree.SubElement(microenv_setup, "variable")
    new_variable.set('name', drug)
    new_variable.set('units', 'mmol')
    new_variable.set("ID", str(variable_count))
    physical_parameter_set = etree.SubElement(new_variable, "physical_parameter_set")
    initial_cond = etree.SubElement(new_variable, "initial_condition")
    initial_cond.set('units', 'mmHg')
    initial_cond.text = str(0.0)
    dich = etree.SubElement(new_variable, "Dirichlet_boundary_condition")
    dich.set('units', 'mmHg')
    dich.set('enabled', 'true')
    dich.text = str(0.0)

    # 3. childs 
    diffusion_coeff = etree.SubElement(physical_parameter_set, "diffusion_coefficient")
    diffusion_coeff.set('units', 'micron^2/min')
    diffusion_coeff.text = str(1200.0)
    decay_rate = etree.SubElement(physical_parameter_set, "decay_rate")
    decay_rate.set('units', '1/min')
    decay_rate.text = str(0.0275)
    
    cell_definitions = root.find('cell_definitions')
    # insert drug as a secreting substance for the default strain 
    cell_definition = cell_definitions.find('cell_definition')
    phenotype = cell_definition.find('phenotype')
    if (phenotype != None):
        secretion = phenotype.find('secretion')
    if (secretion != None):
        substrate = etree.SubElement(secretion, "substrate")
        substrate.set("name", drug)
        secretion_rate = etree.SubElement(substrate, "secretion_rate")
        secretion_rate.set("units", "1/min")
        secretion_rate.text = str(0.0)
        secretion_target = etree.SubElement(substrate, "secretion_target")
        secretion_target.set("units", "substrate density")
        secretion_target.text = str(1.0)
        uptake_rate = etree.SubElement(substrate, "uptake_rate")
        uptake_rate.set("units", "1/min")
        uptake_rate.text = str(0.0025)
        net_export_rate = etree.SubElement(substrate, "net_export_rate")
        net_export_rate.set("units", "total substrate/min")
        net_export_rate.text = str(0.0)

    # insert drug in custom data 
    custom_data = cell_definition.find('custom_data')
    if (custom_data != None):
        new_drug_conc = etree.SubElement(custom_data, drug + "_concentration")
        new_drug_conc.set("units", "dimensionless")
        new_drug_conc.text = str(0.0)
        new_drug_node = etree.SubElement(custom_data, drug + "_node")
        new_drug_node.set("units", "dimensionless")
        new_drug_node.text = str(0.0)

    # insert two new cell strains for the drug 
    new_cell_def_1 = etree.SubElement(cell_definitions, "cell_definition")
    new_cell_def_1.set("name", drug + "_sensitive")
    new_cell_def_1.set("ID", str(len(cell_definitions.getchildren())-1))
    new_cell_def_1.set("parent_type", "default")

    new_cell_def_2 = etree.SubElement(cell_definitions, "cell_definition")
    new_cell_def_2.set("name", drug + "_insensitive")
    new_cell_def_2.set("ID", str(len(cell_definitions.getchildren())-1))
    new_cell_def_2.set("parent_type", "default")

    # insert drug in user parameters
    user_parameters = root.find('user_parameters')

    uptake_rate = etree.SubElement(user_parameters, drug + "_uptake_rate")
    uptake_rate.set("type", "double")
    uptake_rate.set("units", "")
    uptake_rate.text = str(0.0025)

    secretion_rate = etree.SubElement(user_parameters, drug + "_secretion_rate")
    secretion_rate.set("type", "double")
    secretion_rate.set("units", "fg/cell/min")
    secretion_rate.text = str(0.1)

    duration_add = etree.SubElement(user_parameters, "duration_add_" + drug)
    duration_add.set("type", "int")
    duration_add.set("units", "min")
    duration_add.text = str(8000)

    time_remove = etree.SubElement(user_parameters, "time_remove_" + drug)
    time_remove.set("type", "int")
    time_remove.set("units", "min")
    time_remove.text = str(8000)

    time_add = etree.SubElement(user_parameters, "time_add_" + drug)
    time_add.set("type", "int")
    time_add.set("units", "min")
    time_add.text = str(0)

    threshold = etree.SubElement(user_parameters, "threshold_" + drug)
    threshold.set("type", "double")
    threshold.set("units", "dimensionless")
    threshold.text = str(0.14)

    drug_conc = etree.SubElement(user_parameters, "concentration_" + drug)
    drug_conc.set("type", "double")
    drug_conc.set("units", "ng/mL")
    drug_conc.text = str(0.5)

    # set the proportion of drug inhibition for the drug 

    inhibition_level = etree.SubElement(user_parameters, "prop_drug_sensitive_" + drug)
    inhibition_level.set("type", "double")
    inhibition_level.set("units", "dimensionless")
    inhibition_level.text = str(conc.replace("_", "."))

    # set the new bnd and cfg files 
    bnd_file = user_parameters.find('bnd_file') 
    # this path is for later when i have in the makefile saved where the files are 
    # bnd_file.text = "{}/{}/{}/{}_{}.{}".format(".", config_path, "boolean_network", model_name, "all_drugs", "bnd")
    bnd_file.text = "{}/{}/{}/{}_{}.{}".format(".", "config", "boolean_network", model_name, "all_drugs", "bnd")
    cfg_file = user_parameters.find('cfg_file') 
    # cfg_file.text = "{}/{}/{}/{}_{}.{}".format(".", config_path, "boolean_network", model_name, "all_drugs", "cfg")
    cfg_file.text = "{}/{}/{}/{}_{}.{}".format(".", "config", "boolean_network", model_name, "all_drugs", "cfg")

    # set the chosen simulation mode 
    simulation_mode = user_parameters.find("simulation_mode")
    if (simulation_mode == None):
        simulation_mode = etree.SubElement(user_parameters, "simulation_mode")
    if (mode == "single"):
        simulation_mode.text = "0"
    else:
        simulation_mode.text = "1"

    et = etree.ElementTree(root)
    et.write(path_to_xml, pretty_print=True)   


def setup_drug_simulations(druglist, bool_model_name, bool_model, base_project_path, conc_list, mode):

    # add drugs to network files
    add_drugs_to_network(bool_model, druglist)
    
    translation_table = dict.fromkeys(map(ord, "[()'[] ]"), None)

    # modify druglist and concentrationlist if mode is double
    if (mode == "double"):
        drug_combinations = combinations(druglist, 2)
        conc_combinations = product(conc_list, repeat=2)
        druglist = drug_combinations
        conc_combinations_list = []
        for elem in conc_combinations:
            conc_combinations_list.append(elem)
        conc_list = conc_combinations_list

    # set the output and config path 
    if (mode == "single"):
        output_base_path = "output/single"
        config_base_path = "config/single"
    else:
        output_base_path = "output/double"
        config_base_path = "config/double"
    
    for drug in druglist:

        for conc in conc_list:
            filtered_drugname = str(drug).translate(translation_table)
            filtered_conc = str(conc).translate(translation_table)
            drugsim_path = "{}/{}_{}_{}".format(base_project_path, bool_model_name,filtered_drugname.replace(",","_"), filtered_conc.replace(",","_"))
            output_path = "{}/{}_{}_{}".format(output_base_path, bool_model_name,filtered_drugname.replace(",","_"), filtered_conc.replace(",","_"))
            config_path = "{}/{}_{}_{}".format(config_base_path, bool_model_name, filtered_drugname.replace(",","_"), filtered_conc.replace(",","_"))

            # create the corresponding output folder
            if os.path.exists(output_path):
                shutil.rmtree(output_path)
            os.makedirs(output_path)  

            # create the corresponding config folder
            if os.path.exists(config_path):
                shutil.rmtree(config_path)
            os.makedirs(config_path)

            # create the project folder and copy the prostate project files in it
            shutil.copytree(prostate_path, drugsim_path)

            # modify the .xml file for the current run
            xml_path = "{}/{}/{}".format(drugsim_path, "config", "PhysiCell_settings.xml")
            if (type(drug) is tuple):
                # for the tuples the first two elements of drug and conc belong together
                add_drug_to_xml(drug[0], conc[0],  xml_path, config_path, bool_model_filename, mode, output_path)
                add_drug_to_xml(drug[1], conc[1], xml_path, config_path, bool_model_filename, mode, output_path)
            else: 
                add_drug_to_xml(drug, conc, xml_path, config_path, bool_model_filename, mode, output_path)

            # add the new sample project to the Makefile in the main Physiboss folder

            # modify the Makefile for the current run // add also that the files have to be stored in the new config and new output folders 

            # run physiboss

    # delete created folders again 
    # shutil.rmtree(project_path)
    # delete new .cfg and .bnd in the original prostate folder 
    os.remove("{}_{}.cfg".format(bool_model,"all_drugs"))
    os.remove("{}_{}.bnd".format(bool_model,"all_drugs"))

         

####################################################################
# Run drug simulations 
####################################################################

mode = args.mode
if (mode == "single"):
    setup_drug_simulations(node_list, bool_model_filename, bool_model, single_run_path, drug_conc_name_list,mode)
elif (mode == "double"):
    setup_drug_simulations(node_list, bool_model_filename, bool_model, double_run_path, drug_conc_name_list, mode)
elif (mode == "both"):
    setup_drug_simulations(node_list, bool_model_filename, bool_model, single_run_path, drug_conc_name_list, "single")
    setup_drug_simulations(node_list, bool_model_filename, bool_model, double_run_path, drug_conc_name_list, "double")

# for each run: each drug and each concentration (and each initial condition) depending on the mode create the different run folders in sample_projects
# named for example "MYC_MAX_00_0_8" or "MYC_MAX_ERK_00_0_4_0_8"
# rule: first initial condition, then used drugs then all concentrations in the same order than the drugs
# copy the files of the prostate project into every created folder

# create in the same loop the output folder structure (depending on which mode was chosen)
# output
#     single
#         --MYC_MAX_00_0_8
#     double
#         --MYC_MAX_ERK_00_0_4_0_8


# modify copied prostate files accordingly:
# i just need to create one config file right? - a config file that contains all specified drug nodes and anti nodes
# anti nodes are set to 0[1] in the beginning and the drug itself to random : it's like a simulation without the drug, 0.0 inhibition
# to simulate different concentrations of the drug .cpp file is modified when cells are initialized 

# then add the used densities in the files and all the other things that have to be modified 


# in .xml file specify the corresponding output folder 

# when everything is set: sample project contains all folders for each run and output folders are set too then run physicell

# then call another script "create_all_sims_txt.py" that creates txt file with:
    # make run0
    # ./run0
    # make run1
    # ./run1

# then call another script "run_drug_sim.sh"/"run_drug_sim_mn4.sh" - that can be either for cluster or not that runs the created .txt file
# the mn4 one calls it with greasy to parallelize 
