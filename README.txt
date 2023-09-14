******** Release version 07/09/2023 up to here ***********

Summary of whats new in terms of models:
    Network model main files added and all underlying functions
    Tissue moodels:
        Vector_field_and_fibrosis_2D_sims
        2D_human_atria_fibrosis_300x300_OY
        2D_human_atria_fibrosis_300x300_field_control
        2D_human_atria_fibrosis_300x300_field_remodelled
        Rat_vent_three_eigenvectors_control_base
    Windows compilation is now much cleaner (please see full documentation for instructions)

Please read documentation for installation and use instructions:
    Qucik_use_docs:
        AVAILABLE_MODEL_LIST.txt            -   Lists all of the cell models that can be used within this code
        BASIC_INSTRUCTIONS_USE.txt          -   Simple instructions to get setup and running
        BASIC_INSTRUCTIONS_MODIFICATION.txt -   Simple instructions to modify the code for your own usage
        NETWORK_MODEL_OPTIONS.txt           -   Options that are specific for use of the new network model and the tool to create heterogeneous maps of network connections

    Full_documentation:
        Documentation.pdf                   -   More detailed and extensive documentation for use and modifications of the code

Legacy versions of the code (V1.1 and older) are still contained in the folder "Legacy"


Version tracking and updates:
=====Version 1.4 SRF read for reproducibility added
    30/06/2023  - new rat vent DTI geo Benson added for network model sims
    
    27/06/2023  - Tidy of tissue IF selection statements to avoid else if
                    - t->Tissue_settings_set  flag added
                - Maximum duration of SRF updated from 1000ms to 1500ms (this would be slowest wave allowed)
                    - if (srf->duration > 1500)   srf->duration = 1500; in Spontaneous_release_functions.cpp
                - Also minimum updated to assign to max, so that we don't get over-representation of short waves, as would happen if assigned to smallest
                    -  if (srf->duration < 20)     srf->duration = 1500

    20/06/2023  - Functionality to uncouple (electrically) physically-connected regions
                    - Add settings if wanted:
                        t->disconnect_regions_flag         = true;
                        t->Ndisconnected_regions           = 2; // region junctions really, as it is N pairs of regions to disconnect
                        t->disconnect_regions[0][0]        = 2;
                        t->disconnect_regions[0][1]        = 3;
                    - Function "Modify_neighbours_region_disconnect" added to lib/Tissue.cpp and called in both FDM tissue mains
                    - Additional clause added to network model create junctions function to make false if disconnected regions
                    - This functionality has been added to canine atria models to ensure PV and RA do not directly couple electrically

    19/06/2023  - AF remodelling for minimal RA model added
                - minimal RA model also updated to have a higher resting potential

    12/06/2023  - Functilnality to read in SRF from file added -> "SRF_Mode Read" now controls this
                    - Argument to set filename to read in "SRF_read_filename"
                    - Note: SRF properties output when read in may be slighlty different due to t_init, but none of those that affect the read will be different (all actual waveform parameters are the same)
                    - Default ARG value for SRF mode set, as main files use Arg rather that SRF mode to determine whether to read in
                - Tissue model "basic_small" added; tissue dimensions added to screen and file settings outputs
                - Default spatial outputs changed: now 1 ms for binary and 5 ms for vtk


=====Version 1.3 Update of network model handling no fibres
    01/04/2023  - Change to the way network model handles nodes without fibres
                    - now assigns all directions as fibre, so that it does not intoduce an artificial region of less coupling (assuming less coupling more arrhythmic than over coupling)
                    - fibreflag and output to screen that this has been done


=====Version 1.2 Network model included, conduction velocity calculation improved, OS Compatability issues addressed
    18/10/2022 - Network model of tissue excitation included, associated with Tissue_native_network and Tissue_integrated_network
                    - New tissue state read/write functions added for network
                    - CV calculation has been improved based on curvature arguments
                    - New tissue models have been added - fibre vector field in 2D and 3D sheets 
                                                         - rat DTI three eigenvector 3D ventricle
                    - New arguments associated with the network model
                        - network connection maps
                        - symmetry factor
                    - Model types for network have been added to the bin to vtk tools (to identify correct output folders)
                - New tool to create heterogeneous network maps has been included
                - Output of three different fibre files added
                
    30/11/2022  - output files renamed from ".dat" to ".txt" for better portability with Windows
                - Main files now determine which OS is being used, and create folders/files differently for the different OS for compatability -> no need for different main files for Windows
                - Compile.bat file added to compile in Windows
                - Should now be fully compatible with Windows, Mac and Linux (other OS have not been considered)
                - *Intention to also optimise memory useage, but not been done yet* -> tried but it seemed to make no worthwhile difference


=====Version 1.1 Sub-cellular heterogeneity included, associated with Methods paper
    03/02/2021 - Extracellular concentrations added as run-time controllable parameters
	05/04/2020 - Functionality for sub-cellular heterogeneity and reading in maps added; documentation updated; "Sub_cellular_het_geometries" directory and single example map file added to "MSCSF_state_and_geometry_files"; example of using this added to Example_scripts/3D_cell; no 7 
    01/04/2020 - Arguments.c file updated to Windows format as default (works in all OS).


=====Version 1.0 Release.
    07/10/2019 - Fix:   In function "output_settings()", the order of variables Gleak and Gup was the incorrect way round relative to the text
    fprintf(so,"\t\t\t\tJup_scale  = %.02f || Jleak_scale = %.02f || || Jrel_scale = %0.2f\n\n", p.Gleak, p.Gup, p.Grel); 
    >>  fprintf(so,"\t\t\t\tJup_scale  = %.02f || Jleak_scale = %.02f || || Jrel_scale = %0.2f\n\n", p.Gup, p.Gleak, p.Grel);

    24/06/2019 - Documentation for installing the package on Windows computers has been included (thank you to Jakub Tomek for this implementation and associated documentation). This is provided ahead of a full version of the framework with included Visual Studio project files, which I intend to finish and upload soon.

    23/06/2019 - Fix:   Small issue with writing the settings file for all implementations other than single native, so it writes now to correct folder.
                        [E.g. output_settings_tissue(Sim, Tissue, directory) -> output_settings_tissue(Sim, Tissue, res_dir_full) in Tissue_native_main.cc and similar equivalents for all main files other than Single_cell_native.cc]

    20/06/2019 - Release version uploaded, with accompanying author revised version of the manuscript (currently under review at PLOS Comp. Biol).
