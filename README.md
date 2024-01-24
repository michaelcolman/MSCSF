____________________________________________________________________
Source code for the "Multi-scale cardiac simulation framework", including multiple models and associated with multiple publications.
    
    Please cite using the Zenodo resource: 10.5281/zenodo.10204624 (all versions doi)

This version is released with the publications:
    "Patchy fibrosis promotes trigger-substrate interactions that both generate and maintain atrial fibrillation"
    Michael A Colman, Roshan Sharma, Oleg V Aslanidi, Jichao Zhao
    Interface Focus. 2023 Dec 15;13(6):20230041. doi: 10.1098/rsfs.2023.0041

    "Interactions between calcium-induced arrhythmia triggers and the electrophysiological-anatomical substrate underlying induction of atrial fibrillation"
    Colman, Michael A., Varela, Marta, McLeod, Rob, Hancox, Jules C and Aslanidi, Oleg V
    The Journal of Physiology (in press. doi: 10.1113/JP285740)

Original source code presented with:
    "Arrhythmia Mechanisms and Spontaneous Calcium Release: Bi-directional coupling between re-entry and focal excitation”
    Michael A Colman, PLOS Comp Biol 2019 Aug 8;15(8):e107260


NOT included in this package:
    The folder containing all geometry files (tissue and cell structures) and state files. These can be found here:
    https://github.com/michaelcolman/MSCSF_state_and_geometry_files
    and is kept up-to-date with additional geometry files associated with new publications
____________________________________________________________________

____________________________________________________________________
Basic instructions for compiling and running the code under different conditions.
Please see also "BASIC_INSTRUCTIONS_MODDING.txt" for how to update the code with new models and functionality
Please see also example script files for examples and illustrations of performing different functionality
See "Full_documentation.pdf" for the full disclaimer including citation list, and for more detailed breakdown
code structure for more complex modding.
See "Single_cell_model_list.txt" and "Tissue_model_list.txt"" for all available cell and tissue models, thier identifiers and filenames.
____________________________________________________________________

____________________________________________________________________
0) Package contents:
    BASIC_USE.txt (this file)
    BASIC_INSTRUCTIONS_MODIFICATION.txt
    Full_documentation.pdf
    CODE                            - source code for the framework
    MSCSF_state_and_geometry_files  - directory containing all state and geometry files (see instruction 2).
    Simulations                     - simulations folder
    Example_scripts                 - further directories with example scripts for most functionality
    Single_cell_model_list.txt
    Tissue_model_list.txt
____________________________________________________________________

____________________________________________________________________
1) Compile the code (mac and Linux)
    
    Navigate to "CODE"
    Copy the Makefile_linux or Makefile_mac to a file called "Makefile"
    Type "make"
    This will compile everything:

    •	“model_single_native”   – Single cell implementations for traditional cell models
    •	“model_tissue_native”   – Tissue implementations using traditional cell models
    •	“model_tissue_network”  – Tissue implementations using traditional cell models, network model of coupling
    •	“model_single_3D”       – 3D stochastic spatial calcium handling single cell models
    •	“model_single_0D”       – Non-spatial equivalent to the 3D stochastic model
    •	“model_tissue_0D”       – Tissue models using the 0D derived calcium system
    •	“model_tissue_0D_network” – Tissue models using the 0D derived calcium system, network model of coupling
    •	“model_Ca_clamp_3D”     – Ca2+ clamp protocol for 3D single cell
    •	“model_Ca_clamp_OD”     – Ca2+ clamp protocol for 0D single cell + Spontaneous Release Functions
    •   "bin_to_vtk_tissue"     - converts binary data to plain text and/or vtk data files; tissue models
    •   "bin_to_vtk_3Dcell"     - converts binary data to plain text and/or vtk data files; 3D single cell models

    You can also type   "make x" where x is the executable name without the "model_" prefix to compile just that implementation.
                        "make bin_to_vtk_{tissue/3Dcell}" to compile just these post-processing tools

1b) Compile the code (Windows)

    If you do not have a compiler, download MinGW (a version has been provided for you). Place this in a sensible location on your computer
    Ensure the Compile.bat points to the path of your compiler (right click and select "edit", then put the path of your compiler)
    Click Compile.bat to compile all versions of the code, or Compile_single_cell_only.bat to compile just the native single cell models
    NOTE: no parallelisation included in these versions - add OMP as you see fit, and the code should work fine
____________________________________________________________________


____________________________________________________________________
2) Set your path - this is the location your tissue geometry and all state files will be stored/read/written
    
    Move the folder "MSCSF_state_and_geometry_files" to a sensible location
    Add the full path to PATH.txt e.g. "/Users/username/MSCSF_state_and_geometry_files"
____________________________________________________________________

____________________________________________________________________
3) Run code
    
    Can run within the CODE folder, but it is suggested to run elsewhere (e.g. simulations folder).
    Copy the executable desired and PATH.txt to your simulations folder
    Execute with default setup: ./model_X
____________________________________________________________________

____________________________________________________________________
4) Running the code with non-default options

    Multiple example scripts have been provided in Example_scripts, to demonstrate full functionality.
    Use the below instructions and the exmaples as a guide to perform any simulation you like.
    
    Executable can be run with any number of options for settings, passed as command line arguments:

    • "./model_X ARG1 VALUE ARG2 VALUE ..."
    
    Any number of arguments  may be passed and in any order - as long as the value follows the associated argument
    run "./model_X Display_args" to print to screen all argument options for that implementation
    Full options and description is in the documentation; below (point 7) are given just the basics.
    Also refer to the example scripts, or copy desired scripts into your simulations folder and modify them as desired.
____________________________________________________________________

____________________________________________________________________
5) Main outputs (full list in documentation):
    
    The simulations will produce multiple outputs.
    • In the parent directory, "Log.dat" contains a log of all simulations performed, the time set, and settings
    • In the parent directory will be created an "Outputs_X" directory where X is determined by 1) executable run and 2) any reference passed in.
    This folder contains all the outputs for simulations associated with X.
    Within this folder are:

    • All models:
        • "Outputs_X/Results_Y/         - Directory (containing time traces for cellular data - ("Y" is determined by a results reference if passed)
            • Currents.dat              - time courses of voltage, membrane currents, gating variables and ionic concentrations
            • Properties.dat            - measured properties (dvdt, APD, CaT_min/max etc) for the final two beats of a simulation
            • Vm_linescan_x.dat         - Idealised tissue models only; contains rows of voltage for an x-linescan, one row per ms
            • CRU.dat                   - 3Dcell/integrated models only; detailed time courses of calcium concentration and calcium handling dynamics
            • Ca_linescan_{x/y/z}.dat"  - 3Dcell models only; linescans crossing through the centre of the cell in x, y and z 
            • Note:                     - For tissue models, data for three different cells within the tissue will be written, with "_cell_{1-3}.dat"
                                          appended to the filenames for Currents, Properties and CRU.
    
        • "Outputs_X/Results_Y/Settings.dat" - text file - record of all of the settings for the most recent simulation X+Y; equivilent to screen outputs
        • "Outputs_X/Properties_log.dat"     - text file - contains a single line entry of final measured properties for every simulation within Outputs_X - appended.

    • Single-cell models only:
        • "Outputs_X/Parameters_Y/
            • Ix.dat                        - All voltage dependent parameters/variables for current Ix
            • Magnitude_parameters.dat      - A list of the final values of all ion channel conductances and calcium handling maximal flux rates (all het and mod included)
            • Modifier_parameters.dat"      - A list of the final values of all modifiers (scale factors, voltage shifts, time-constant scaling etc)
            • Voltage clamp files (if Vclamp is set to On):
                • I{CaL/K/K1}_IV.dat       - Current (peak) voltage relationship for ICaL, Ito, IKur, IKr, IKs and IK1
                • I{CaL/K}_Vclamp_trace_Z.dat  - time courses of currents for voltage step Z
                • I{CaL/K}_Vclamp_traces.dat   - time courses for all voltage steps, in sequence 
       
    • 3D cell and tissue models: 
        • "Outputs_X/Spatial_Results_Y" - Directory (containing spatial output data where relevant - tissue and 3D cell models)
            • Note:                     - By default, outputs binary files. These can be converted to plain text or vtk file data 
                                          using the supplied tools (see Section 6 below) - binary data of course saves space.
            • Tissue models:
                • Vm_output_t.{bin/vtk}      - where t is the time in ms; contains the voltage for all tissue (.bin) and all space (.vtk). 
                • Ca/CaSR_output_t.{bin/vtk} - integrated (0D) models only; spatial data of Cai or CaSR

            • 3D cell models:
                • {Ca/CaSR/CaDS}_output_t.{bin/vtk} - spatial calcium data within the cell
 
            • Note:                      - you can output any spatial variable in native or integrated models if you like:
                                           in the relevant main file (Tissue_native_main.cc/Tissue_integrated_main.cc/Single_cell_3D_main.cc)
                                           navigate to the spatial data outputs (search function "vtk_3D_output" or "array_1D_output") and
                                           follow the format of the variables already output, passing in an appropriate string reference.

            • All maps (tissue and 3D cell models), geometries and fibres and final activation pattern (tissue models only) 
                will also be output in "Outputs_X" for inspection, verification, and figures.
 
• Key output file contents (complete list in Full_documentation.pdf):

    • Outputs_X/Results_Y/Currents.dat
    
        column      variable    unit
        1           time        ms
        2           voltage     mV
        3           INa         A/F
        8           Ito         A/F
        12          ICaL        A/F
        19          IKr         A/F
        22          IKs         A/F
        24          IK1         A/F
        34          Cai         mM
    
    • Outputs_X/Results_Y/Properties.dat
    
        1           time        ms
        2           voltage     mV
        5           dv/dt_max   mV/ms
        6           Vmin        mV
        7           Vmax        mV
        8           Vamp        mV
        11          APD_30      ms
        12          APD_50      ms
        13          APD_70      ms
        14          APD_90      ms

     • Outputs_X/Properties_log.dat

        1           BCL                         ms
        2           S2                          ms
        5           APD_30, final beat          ms
        7           APD_50, final beat          ms
        9           APD_70, final beat          ms
        11          APD_90, final beat          ms
        13          dv/dt_max, final beat       mV/ms
        15          Vmin, final                 mV
        17          Vmax, final                 mV
        19          Vamp, final                 mV
        21          CaT min, final              microM
        23          CaT_max, final              microM

    • Outputs_X/Results_Y/CRU.dat  (spatial/integrated cell models only)

        1           time        ms
        2           voltage     mV
        3           [Ca]ds      microM
        4           [Ca]ss      microM
        5           [Ca]cyto    microM
        6           [Ca]jSR     microM
        7           [Ca]nSR     microM
        8           J_rel       microM/ms
        9           NRyR_Open   -
        16          J_CaL       microM/ms
        18          J_SERCA     microM/ms

    • Outputs_X/Parameters_Y/ICaL_Vclamp_trace{s/_Y}.dat

        1           time        ms
        2           voltage     mV
        3           ICaL        A/F
        4           Cai         mM 
    
    • Outputs_X/Parameters_Y/IK_Vclamp_trace{s/_Y}.dat

        1           time        ms
        2           voltage     mV
        3           Ito         A/F
        4           IKur        A/F
        5           Ito+IKur    A/F
        6           IKr         A/F
        7           IKs         A/F

    • Outputs_X/Parameters_Y/Ix.dat: (Not all variables are computed for all models, so check which are relevant to your model)
        Note: these voltage dependent functions are exaclty as produced under simulation conditions, with all modulation incorporated  
                Use to both check/verify the implementation of modulation, and for convenience for plotting and figures etc
        
        Current/file        1       2       3       4       5       6       7        8       9        10      11      12      13
        INa.dat             Vm      va_a    va_b    va_ss   va_tau  vi_1_a  vi_1_b   vi_1_ss vi_1_tau vi_2_a  vi_2_b  vi_2_ss vi_2_tau
        INaL.dat            Vm      va_a    va_b    va_ss   va_tau  vi_a    vi_b     vi_ss   vi_tau
        Ito.dat             Vm      va_ss   va_tau  vi_ss   vi_tau  vi_3_ss vi_s_tau vi_Fs
        ICaL.dat            Vm      va_ss   va_tau  vi_ss   vi_tau  vi_s_tau
        IKur.dat            Vm      va_ss   va_tau  vi_ss   vi_tau  dynamic_g
        IKr.dat             Vm      va_ss   va_tau  v_ti
        IKs.dat             Vm      va_ss   va_tau  va_2_tau
        IKACh.dat           Vm      va_ss   va_tau  v_ti 

    • Outputs_X/Parameters_Y/Magnitide_parameters.dat and Modifier_parameters.dat contain strings with the variables and values in the files

    • When using 0D single cell models with the SRF turned on, the SRF model distributions will also be output to Parameters_[results_reference]/SRF_distributions.
        “Parameters_[results_reference ]/SRF_distributions/ti_distribution.dat” – histogram of ti
        “Parameters_[results_reference ]/SRF_distributions/duration_distribution.dat” – histogram of duration
        “Parameters_[results_reference ]/SRF_distributions/NRyRo_peak_distribution.dat” – histogram of NRyRo
        If a dynamic mode is selected, these above files will be appended with various CaSR concentrations, and contain the distribution at that CaSR. 
            The additional file “../CaSR_dependence.dat” will also be created which tracks how the SRF parameters vary with CaSR (1 – CaSR; 2 – PSCRE; 3 – ti_sep; 
            4 – ti_width_1; 5 – ti_width_2; 6 – MD; 7 – duration_width_1; 8 – duration width 2).
____________________________________________________________________

____________________________________________________________________
6) Visualsing spatial data (3D-cell models and tissue models):

    If running these models, with default settings, the simulation will write only binary data for the spatial outputs 
        (vtk files can be output directly by the model - see arguments below)

    Use the provided tools to convert the binary data to plain text or vtk data files:
    From within the same directory where the simulation was performed, or in a parent directory containing the Outputs directory you want to analyse, 
    run: "./bin_to_vtk_{tissue/3Dcell}” with arguments identifying the simulation and conversion settings:

        •	Any Reference and Results_Reference passed (so it knows which directory to look in to read binary and write data/vtk files)
        •	The variable you want to convert (i.e. Vm, Ca, CaSR etc): Variable [V]
        •	The Tissue_order and Tissue_model (tissue) or Cell_size and Sim_cell_size (3Dcell) used to perform the simulation (so it knows geometry sizes and files etc).
        •	Which type of data to write: Write_data [On/Off] (plain text) and/or Write_vtk [On/Off] 
        •	The time range and time interval over which to convert data: start_time [n1] end_time [n2] interval [n3]
        •	For tissue models, you need to also specify the model type (Model_type [native/integrated])
        •	For 3Dcell models, you can also write plain text data for specified 2D slices:
            o	Write_slices [On/Off]
            o	XZ_slice_y [x] to define the y coordinate of the XZ slice (and similar for all directions).
        
    Examples: 
        "./bin_to_vtk_tissue start_time 100 end_time 200 interval 10 Model_type integrated Variable Cai Write_vtk On Write_data On. Tissue_order 2D Tissue_model basic Reference X Results_Reference Y" 
        "./bin_to_vtk_3Dcell start_time 100 end_time 200 interval 10 Variable CaSR Write_vtk On Write_data On Write_slices On YZ_slice_x 5 Cell_size standard Sim_cell_size full Reference X Results_Reference Y"
____________________________________________________________________

____________________________________________________________________
7) Basic arguments

    • All implementations:
        Model                   [string model identifier]
        BCL                     [x ms]
        Beats                   [n]
        Total_time              [x ms]
        S2                      [x ms]
        Reference               [any string]  -> names the "Outputs_X" directory
        Results_Reference       [any string]  -> names the "Results_Y" and "Spatial_Results_Y" directory within Outputs_X
        Read_state              [On/Off + other options for tissue - see example scripts for tissue read/write state]
        Write_state             [On/Off + other options for tissue]
        State_Reference_write   [string]    -> adds a reference to the state-file written
        State_Reference_read    [string]    -> reads in the state file with reference
        Settings_file           [filename]  -> read options from a settings file (see below for writing and using settings files)
        Vclamp                  [On/Off]    -> performs simple voltage clamp for ICaL and potassium currents (will need to update if more complex protocol is required)
     
    • All tissue models:
        Tissue_order                    [1D/2D/3D/geo]  -> idealised 1-3D models, or geo where a geoemtry file is read in
        Tissue_model                    [string]        -> selects settings for specific tissue models, includes idealised: "basic", "conduction_velocity", "re-entry"
        Tissue_type                     [homogeneous/heterogeneous] -> electrical heterogeneity
        Orientation_type                [isotropic/anisotropic] -> whether to run the isotropic FDM (without orientation), or anistropic (with orientation)
        D_uniformity                    [uniform/regional/map/regional_map]  -> Whether D1, D2, D_AR vary in space for the baseline model (not associated with remodelling);
                                            implemented as D1 and D_AR scale-factors, relative to the baseline values set. regional = scales by celltype;
                                             map = values set by map; regional_map = both
        D1                              [x] sets D1 explicitly
        D_AR                            [x] sets anisotropy ratio (D2 = D1/D_AR) explicitly
        Dscale_base_map_file            [filename] -> define explicitly filename which sets Dscale from map (if "map" or "regional_map" is set by D_uniformity)
        D_AR_scale_base_map_file        [filename] -> define explicitly filename which sets D_AR_scale from map (if "map" or "regional_map" is set by D_uniformity)
        Stimulus_location_type          [edge/centre/cross_field/{other specific string}]
        S2_Stimulus_location_type       [S1/cross_field/edge/{other specific string}]
        {S1/S2}_shape                   [cuboid/sphere]
        {stim/S2_stim}_map_file         [filename] -> define explicitly map filename to be read in for S1 or S2 stimuli
        Spatial_output_interval_data    [n ms]     -> interval to output binary spatial data (default is 5 ms)
        Spatial_output_interval_vtk     [n ms]     -> interval to output vtk data directly (default is 0, which is off)
        Read_state                      [Off/On/phase/single_cell/ave]  -> phase = read state files for phase-distribution re-entry; 
                                                                           single_cell = read in from single_cell written file; 
                                                                           ave = read in from single coupled cell; 
                                                                           On = whole tissue read
        Wirte_state                     [Off/On/phase/ave]              -> same as read, except no single_cell as cannot write single_cell state using a tissue model!

        And for idealised tissue models (Tissue_order = 1D, 2D or 3D but not geo):
            OX                          [x]     -> value of x-component of orientation (globally applied)
            OY                          [x]     -> value of y-component of orientation (globally applied)
            OZ                          [x]     -> value of z-component of orientation (globally applied)  ; only 2 need be specified, as third will be defined from normalisation constraint
            Global_orientation_direction    X/Y/Z/{XY/XZ/YZ}_plus/{XY/XZ/YZ}_minus/XYZ_{ppp/ppm/pmp/mpp}  -> applies predefined OX, OY and OZ for specific global orientation directions
            
    • 3D single cell models:
        Cell_size                       [string e.g. standard/thin..]   -> determines cell dimensions
        Sim_cell_size                   [full/portion/testing]          -> determines whether full or just a portion of cell is simulated 
        Spatial_output_interval_data    [n ms]  -> interval to output binary spatial data (default is 5 ms)
        Spatial_output_interval_vtk     [n ms]  -> interval to output vtk data directly (default is 0, which is off)
        {volds/RyR/LTCC}_het            [Off/random]    -> to apply volds, NRyR and LTCC homogeneously in tissue, or with random variation around a mean
        Detub                           [On/Off]        -> apply variable TT denisty
        {SERCA/NCX}_het                 [Off/On]        -> apply a sub-cellular heterogneous SERCA or NCX scale map 
        {RyR/LTCC}_het                  [map]           -> (additional to Off/random as above); applies NRyR or LTCC heterogeneity from map (cannot be map and random simultaneous) 
        {RyR_het/LTCC/SERCA/NCX}_map_file   [filename]  -> define explicitly map file to scale NRyR/NLTCC or Gup/GNCX
        Read_state                      [Off/On/ave]    -> On = reads in spatial state of cell; ave = reads in average values (and applies homogeneously within cell)
        Write_state                     [OFf/On/ave]    -> same as above, but writing equivilents.

    • 0D single cell and tissue models - spontaneous release functions
        SRF_mode                [Off/Direct_Control/Dynamic]
        SRF_model               [3D_cell/General]   -> applies to "Dymamic" mode only
        SRF_Pset                [string]            -> selects pre-defined settings or User_Control to set all parameters from arguments

    • Modulation/modification (can have all set at once; only one of each type at a time; modifications applied in sequence):
        Celltype                [string]    -> selects regional cell model
        ISO                     [0-1]       -> represents 0 to 1microM or saturating solution
        ISO_model               [string]    -> which ISO model to apply, if multiple are available
        ACh                     [0-1]       -> represents 0 to 1microM or saturating solution
        ACh_model               [string]    -> which ACh model to apply, if multiple are available
        Agent                   [string]    -> apply pharmacologiocal agent 
        Remdelling              [string]    -> apply remodelling 
        Mutation                [string]    -> apply mutation 
        Spatial_gradient        [string]    -> apply a spatial gradient model (e.g apico-basal het); in tissue, requires a map; in single cell, requires also:
            Spatial_gradient_proportion [0-1] to determine where the single cell lies within the gradient

        Direct_control modulation:
            X_scale                 [x]         -> scale current or flux X (e.g. INa_scale, Jup_scale)
            Ix_va_ss_shift          [x mV]      -> shifts the voltage dependence of the steady-state of activation gate for current x
            Ix_vi_ss_shift          [x mV]      -> shifts the voltage dependence of the steady-state of inactivation gate for current x
            Ix_va_tau_scale         [x]         -> scales time constant of activation gate for current x
            Ix_vi_tau_scale         [x]         -> scales time constant of inactivation gate for current x
            RyR_Po, LTCC_po         [x]         -> integrated cell models only; scales the open transition rate. 

    • Tissue modulation (Tissue model only!):
        Dscale                  [x]         -> scales D1 (homogeneously if no map is set); D_AR preserved, so D2 calculated from scaled D1
        D_AR_scale              [x]         -> scales D_AR (D1 = D/D_AR) (homogeneously if no mpa is set)
        Dscale_mod_map          [On/Off]    -> Apply Dscale according to map; map is floats, where map = 0, Dscale_local = 1(=D1) (i.e. control), 
                                                where map = 1, Dscale_local = Dscale(=D1*Dscale) (i.e. fully scaled according to value in Dscale)
        Dscale_mod_map_file     [string]    -> Define explicitly the filename containing the Dscale modulation map
        D_AR_scale_mod_map      [On/Off]    -> Apply D_AR_scale according to map; map is floats, where map = 0, D_AR_scale_local = 1(=D_AR) (i.e. control), 
                                                where map = 1, D_AR_scale_local = D_AR_scale(=D_AR*D_AR_scale) (i.e. fully scaled according to value in Dscale)
        D_AR_scale_mod_map_file [string]    -> Define explicitly the filename containing the D_AR_scale modulation map
        {ISO/ACh/Remodelling/Dscale_mod/D_AR_scale_mod}_map  [On/Off]  -> Whether to apply single-cell modulations homogeneously to tissue, or by map 
                                                                        (exactly as with D_scale: map determines where zero-maximum effect of modulation is applied.
        Direct_modulation_map   [On/Off]    -> Whether to apply direct modulation (X_scale;Ix_va_ss_shift;Ix_vi_ss_shift;Ix_va_tau_scale;Ix_vi_tau_scale) 
                                                homogeneously or by map. MAP MUST BE INTEGER for direct_mod on or off, not continuos (as with above maps)
        {X}_map_file            [filename]  -> explicitly set the file for the modulation map X (X=ISO/ACh/Remodelling .....)
____________________________________________________________________

____________________________________________________________________
8) Using a settings file
    
    You can use settings files to store the options passed in, for quick reproduction etc
    These files can be anywhere, where the full path of the file must be passed in after "Settings_file" argument.
    You can use a settings file AND command line arguments:
        "Settings_file" argument must be passed FIRST (this is the only restriction on order of arguments)
        Any arguments passed after this will supplement or overwrite the options in the file 
        -> command line args supersede settings files if the same argument is in both.

    Settings files must have the following format as a plain text document,
    new line for each argument, spaces or tabs to separate arg and value:
    
        Number of args
        ARG 1   Value 1
        ARG 2   Value 2 
    
    Example: An example settings_file could have contents:
        3
        BCL 345
        Model hAM_CRN
        Results_Reference file_test

    And performing:
        ./model_single_native Settings_file [filename] Remodelling AF_GB
        will apply BCL = 345; Model = hAM_CRN; Results_Reference = file_test; Remodelling = AF_GB
____________________________________________________________________
        
