#!/bin/sh
 
# Using the different read and write state functionality to efficiently pace a large tissue model
# to stable state.

# There are multiple ways to read and write state files using tissue models.
# Setting Write_state to On will write the state of the whole tissue, which can then be read in with 
# Read_state On. Setting Write_state to ave will write the state of the centre cell in the tissue, 
# which can then be read into the whole tissue with Read_state ave (you will need to create a state 
# file for every celltype and modulation condition used in your tissue model). Finally, you can also 
# read in state files created by the single cell models using Read_state single cell.

# Read/Write_state On can be used to save whole tissue states for more efficient simulations 
# (e.g. pace to stable state, then test S2 intervals and locations, or pharma action, or similar). 
# The single cell and ave functionality can be used sequentially to to get to stable state on a large 
# tissue model more efficiently.

# Let’s explore that, using Tissue_model functional_model_test at a rapid rate, with heterogeneity turned on, 
# and ISO applied according to a map (reproducing a complex situation).

# First, we’ll create state files using the single cell models paced to stable states. We need a state file 
# for each celltye in the model (EPI, M and RA) and with and without ISO.
# Then we’ll use these state files as a starting point, and pace a 1D homogeneous tissue for each condition 
# to output the coupled cell properties. We’ll use fewer beats as we already paced the single cell models to 
# stable state over 200 beats. Note: we need to set some stimulus properties, as the default for this tissue 
# model (centre and sphere) does not create a stimulus region in 1D (due to the sphere; we also move region 
# to edge so that we are not writing the state of a cell in the stimulus region).
for celltype in ENDO M RA
do
	for i in 0 1
    do
        iso=$(printf "%.2f" $i)        
        ./model_single_native Model minimal BCL 400 Celltype ${celltype} ISO ${iso} Beats 200 Write_state On Reference single_cell_pre_pace Results_Reference ${celltype}_${iso}
        ./model_tissue_native Tissue_model functional_model_test Tissue_order 1D Model minimal BCL 400 Stimulus_location_type edge S1_shape cuboid Tissue_type homogeneous Celltype ${celltype} ISO ${iso} Beats 20 Read_state single_cell Write_state ave Reference Tissue_state_read Results_Reference ave_coupled_write_${celltype}_${iso}
    done
done

# And finally, we’ll read the average coupled state files into our full 2D or 3D tissue model, pace for five beats, 
# and write the whole tissue state:
./model_tissue_native Tissue_model functional_model_test Tissue_order 2D Tissue_type heterogeneous Model minimal BCL 400 ISO 1 ISO_map On Beats 5 Read_state ave Write_state On Reference Tissue_state_read Results_Reference ave_coupled_read_whole_tissue_write

# Which can now be read in to start a new simulation from that point:
./model_tissue_native Tissue_model functional_model_test Tissue_order 2D Tissue_type heterogeneous Model minimal BCL 400 ISO 1 ISO_map On Beats 2 Read_state On Reference Tissue_state_read Results_Reference whole_tissue_read

# We can see the effect this has by inspecting the properties of the cells. "Cell 2" corresponds to the ENDO cell in
# a region without ISO, so comparing:
# "Outputs_tissue_native_Tissue_state_read/Results_ave_coupled_write_ENDO_0.00/Properties_cell2.dat" time and APD (1 vs 14)
# with "Outputs_tissue_native_Tissue_state_read/Results_ave_coupled_read_whole_tissue_write/Properties_cell2.dat"
# and "Outputs_tissue_native_Tissue_state_read/Results_whole_tissue_read/Properties_cell2.dat", we can see that
# by the time we read in the whole state, the APD is stable

# This process is limited if we are using complex tissue models with continuous maps, as creating state files 
# for every individual condition becomes infeasible or impossible. The ISO and ACh concentrations are written 
# to two decimal places in the state files, so if a continuous ISO or ACh map is used, we would need to create 
# a state file for every celltype and every concentration in two decimal place intervals between 0 and max 
# concentration; other modulation proportions or maps (Remodelling_prop, Direct_modulation map) are not encoded 
# into state files and so tissue models which contain these cannot read from single_cell or ave coupled cell 
# state files; only full state reading and writing will work. For heterogeneity and discrete ISO, ACh, mutation 
# and Remodelling maps, we can feasibly use the above process to create state files for every celltype and 
# modulation condition in the tissue model.

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents
