# COMPILES (COMPutIng anaLytic positivE Steady states)

This is a MATLAB code to determine the analytic positive steady state solution of a system of ordinary differential equations which is expressed using its dynamically equivalent chemical reaction network. The main function used is steadyState.m (see examples for sample codes on how to input the chemical reaction network and how to use the function). The ZIP file contains all the functions, classes, and examples used in this package.

## Code Description

**steadyState.m**

This is the main function for this package. Users need to run this function once the chemical reaction network has been inputted using the function addReaction. Running steadyState automatically displays the results of the package but output variables are also available if users need to access them. steadyState.m uses several functions which are written after the main funtion in the same MATLAB file.

The function steadyState returns the steady state solution of a chemical reaction network parametrized by rate constants (and/or sigma's). The free parameters and conservation laws are also listed after the solution. If there are subnetworks that could not be solved because translation is taking too long, solving the subnetwork is skipped and a message saying it could not be solved is displayed. In this case, the parametrization of the steady state of the entire network is not completed but the lists of solved and unsolved subnetworks are displayed. In the case where all subnetworks are solved but the solution of the entire network cannot be parametrized in terms of the free parameters (due to two or more species dependent on each other), the solution returned is in terms of both free parameters and other "nonfree" species. If there are subnetworks containing only 1 reaction, it means that there are species with 0 steady steady; hence, the network has not positive steady state and a message appears saying so. The output variables 'equation', 'species', 'free_parameter', 'conservation_law', and 'model' allow the user to view the following, respectively:

  * List of parametrization of the steady state of the system
  * List of steady state species of the network
  * List of free parameters of the steady state
  * List of conservation laws of the system
  * Complete network with all the species listed in the 'species' field of the structure 'model'

**addReaction.m**

To add the reactions to the network, use the function addReaction where the output is 'model'. See examples for sample codes.

This function returns a structure called 'model' with added field 'reaction' with subfields 'id', 'reactant', 'product', 'reversible', and 'kinetic'. The output variable 'model' allows the user to view the network with the added reaction. The following are the inputs to the function:

  * model: a structure, representing the CRN
  * id: visual representation of the reaction, e.g., reactant -> product (string)
  * reactant_species: species of the reactant complex (cell)
  * reactant_stoichiometry: stoichiometry of the species of the reactant complex (cell)
  * reactant_kinetic: kinetic orders of the species of the reactant complex (array)
  * product_species: species of the product complex (cell)
  * product_stoichiometry: stoichiometry of the species of the product complex (cell)
  * product_kinetic: "kinetic orders" of the species of the product complex, if the reaction is reversible (array); if the reaction in NOT reversible, leave blank
  * reversible: logical; whether the reaction is reversible or not (true or false)

**Classes**

The following are included as separate MATLAB files because they are classes that need to be saved individually: **edge.m**, **graph_.m**, **vertex.m**. These classes are used in the function analyticSolution.

*Make sure that these 5 MATLAB files are in the same working directory when using the function steadyState.m.*

## Examples

15 examples are included in the package. Examples 0 to 6, 10b, and 11b include checker files to show that the solutions returned by the function steadyState are indeed steady states of their respective ordinary differential equations.

## Limitations

1. The package assumes that the chemical reaction network inputted has mass action kinetics.

2. In some instances, all subnetworks are solved but the solution of the entire network could not be parametrized in terms of the free parameters (due to two or more species dependent on each other). Renaming the variables can sometimes solve the problem: upon doing this, some subnetworks may be solved for different species depending on the variable assigned to them since the selection of species to solve is based on alphabetical order (see Examples 10 and 10b, and Examples 11 and 11b).
