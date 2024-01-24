#!/bin/bash

# It is certainly not suggested to simply run this whole script and produce all data (that would be thousands of simulations). 
# Rather, this is intended to make it as easy as possible for you to run whatever simulations you chose. 
# Please explore and use for whatever purposes you may need!

# Create CaSR SRF files================
./model_single_0D Total_time 1 SRF_mode Dynamic SRF_model General SRF_Pset User_control SRF_Dyn_PSCR_threshold 1.05 SRF_Dyn_CaSR_max 1.8 SRF_Dyn_CaSR_Prange 0.175 SRF_Dyn_ti_sep_max 1800 SRF_Dyn_ti_sep_min 950 SRF_Dyn_ti_width_max 1250 SRF_Dyn_ti_width_min 250 SRF_Dyn_MD_max 800 SRF_Dyn_MD_min 100 SRF_Dyn_duration_width_max 500 SRF_Dyn_duration_width_min 25 SRF_Dyn_H 1.5 Write_SRF_settings On reference SRF_write

mkdir Settings_files
cp -r Outputs_0Dcell_SRF_write/Parameters/SRF_distributions/SRF_Settings_file_*.txt Settings_files

# Create independent SRF files
mkdir SRF_settings_files_DC
clang -o Create_SRF_settings_linear_variation Create_SRF_settings_linear_variation.c
./Create_SRF_settings_linear_variation

# Create state files ==================
# Single cell ========
# loop over cycle length
for i in 1000 900 800 700 600 500 400 300 250 200 150
do
    BCL=$(printf "%04d" $i)
 
    #TT=$(( $i * 250 )) # i not BCL as leading 0 causes problems
    TT=$(( 1 )) #just for testing the script. Obviously state files cannot be created with a total time of 1

    ./model_single_0D Model minimal Celltype RA Beats 250 BCL $BCL Total_time $TT Write_state On Remodelling none ISO 0.0 Reference Control_control Results_Reference $BCL
    ./model_single_0D Model minimal Celltype RA Beats 250 BCL $BCL Total_time $TT Write_state On Remodelling AF ISO 0.0 Reference AF_control Results_Reference $BCL
    ./model_single_0D Model minimal Celltype RA Beats 250 BCL $BCL Total_time $TT Write_state On Remodelling none ISO 1.0 Reference Control_ISO Results_Reference $BCL
    ./model_single_0D Model minimal Celltype RA Beats 250 BCL $BCL Total_time $TT Write_state On Remodelling AF ISO 1.0 Reference AF_ISO Results_Reference $BCL
done

# Tissue ===========
for k in 400 200 150
do
    BCL=$(printf "%04d" $k)
    #TT=$(( $k * 10 )) # i not BCL as leading 0 causes problems
    TT=$(( 1 )) #just for testing the script. Obviously state files cannot be created with a total time of 1

    # ave
    for rem in none AF
    do
        for i in 0.0 1.0
        do
            ISO=$(printf "%2f" $i)

            ./model_tissue_0D_network Tissue_order geo Tissue_model 2D_human_atria_fibrosis_300x300_OY Model minimal Celltype RA BCL $BCL Remodelling $rem ISO $ISO Read_state single_cell Write_state ave Tissue_type homogeneous Total_time $TT Beats 10 Spatial_output_interval_data 0 Spatial_output_interval_vtk 0 D1 0.4 D_AR 7 reference Tissue_ave_state_write_${rem}_${ISO}
        done
    done

    #full tissue
    for tissue in OY field_control field_remodelled
    do
        for rem in none AF
        do
            for i in 0.0 1.0
            do
                ISO=$(printf "%2f" $i)

                # control geometry
                ./model_tissue_0D_network Tissue_order geo Tissue_model 2D_human_atria_fibrosis_300x300_$tissue Model minimal Celltype RA BCL $BCL Remodelling $rem ISO $ISO Read_state ave Write_state On Tissue_type homogeneous Total_time $TT Beats 10 Spatial_output_interval_data 0 Spatial_output_interval_vtk 0 D1 0.4 D_AR 7 reference Tissue_state_write_Control State_Reference_write Control

                #fibrosis geometries
                for map in 10_25 10_50 25_25 25_50 37_25 37_50 50_25 50_50
                do
                    ./model_tissue_0D_network Tissue_order geo Tissue_model 2D_human_atria_fibrosis_300x300_$tissue Model minimal Celltype RA BCL $BCL Remodelling $rem ISO $ISO Read_state ave Write_state On Tissue_type homogeneous Total_time $TT Beats 10 Spatial_output_interval_data 0 Spatial_output_interval_vtk 0 D1 0.4 D_AR 7 Junction_modulation_map On Junction_modulation_map_file 2D_human_atria_fibrosis_300x300_${tissue}_JNmap_${map}.txt   reference Tissue_state_write_Fib_${map} State_Reference_write fib_${map}
                done
            done
        done
    done
done


# Arrhythmia substrate ================
TT=$(( 1 ))
#TT=$(( 3000 ))

for k in 200 150
do
    BCL=$(printf "%04d" $k)
    for tissue in OY field_control field_remodelled
    do
        for rem in none AF
        do
            for i in 0.0 1.0
            do
                ISO=$(printf "%2fd" $i)

                # control geometry
                ./model_tissue_0D_network Tissue_order geo Tissue_model 2D_human_atria_fibrosis_300x300_$tissue Tissue_type homogeneous Model minimal Celltype RA BCL $BCL NBeats 5 Total_time $TT Remodelling $rem ISO $ISO Read_state On state_reference_read Control reference Rapid_pace_Control_${tissue} results_reference rem_${rem}_IS0_${ISO}_BCL_${BCL}

                 #fibrosis geometries
                 for map in 10_25 10_50 25_25 25_50 37_25 37_50 50_25 50_50
                 do
                     ./model_tissue_0D_network Tissue_order geo Tissue_model 2D_human_atria_fibrosis_300x300_$tissue Tissue_type homogeneous Model minimal Celltype RA BCL $BCL NBeats 5 Total_time $TT Remodelling $rem ISO $ISO Read_state On state_reference_read fib_${map} Junction_modulation_map On Junction_modulation_map_file 2D_human_atria_fibrosis_300x300_${tissue}_JNmap_${map}.txt reference Rapid_pace_Fib_${map}_${tissue} results_reference rem_${rem}_IS0_${ISO}_BCL_${BCL}
                 done
             done
         done
     done
 done

#CaSR dependence of SCRE ===============
TT=$(( 1 ))
#TT=$(( 2000 ))

for i in `seq 1 1 1`
do
    run=$(printf "%04d" $i)
    for j in `seq 900 100 1300`
    do
        CaSR=$(printf "%04d" $j)
        for tissue in field_remodelled #OY field_control field_remodelled
        do
            for rem in none AF
            do
                for i in 0.0 #1.0
                do
                    ISO=$(printf "%2fd" $i)
                    # control geometry
                    ./model_tissue_0D_network Settings_file Settings_files/SRF_Settings_file_${CaSR}.txt Tissue_order geo Tissue_model 2D_human_atria_fibrosis_300x300_${tissue} Tissue_type homogeneous Model minimal Celltype RA BCL 400 NBeats 1 Total_time ${TT} Remodelling $rem ISO $ISO Read_state On state_reference_read Control Delayed_CaSR_IC On CaSR_IC_delay 500 CaSR 1000 Cai 0.1 reference PTA_Control_${tissue} results_reference rem_${rem}_IS0_${ISO}_CaSR_${CaSR}_run_${run}

                     #fibrosis geometries
                     for map in 10_25 10_50 25_25 25_50 37_25 37_50 50_25 50_50
                     do
                         ./model_tissue_0D_network Settings_file Settings_files/SRF_Settings_file_${CaSR}.txt Tissue_order geo Tissue_model 2D_human_atria_fibrosis_300x300_${tissue} Tissue_type homogeneous Model minimal Celltype RA BCL 400 NBeats 1 Total_time ${TT} Remodelling $rem ISO $ISO Read_state On state_reference_read Control Delayed_CaSR_IC On CaSR_IC_delay 500 CaSR 1000 Cai 0.1 Junction_modulation_map On Junction_modulation_map_file 2D_human_atria_fibrosis_300x300_${tissue}_JNmap_${map}.txt reference PTA_Fib_${map}_${tissue} results_reference rem_${rem}_IS0_${ISO}_CaSR_${CaSR}_run_${run}
                     done
                 done
             done
         done
     done
 done

# Fibroblast coupling =================
# Create state files 
for k in 400
do
    BCL=$(printf "%04d" $k)
    #TT=$(( $k * 10 )) # i not BCL as leading 0 causes problems
    TT=$(( 1 )) #just for testing the script. Obviously state files cannot be created with a total time of 1

    for tissue in OY field_control field_remodelled
    do
        for rem in none AF
        do
            for i in 0.0 1.0
            do
                ISO=$(printf "%2f" $i)

                for l in 10 25 37 50
                do
                    ln=$(printf "%02d" $l)
                    for m in 25 50
                    do
                        area=$(printf "%02d" $m)
                        map=$(printf "%s_%s" $ln $area);
                        echo $map

                        ./model_tissue_0D_network Tissue_order geo Tissue_model 2D_human_atria_fibrosis_300x300_$tissue Model minimal Celltype RA BCL $BCL Remodelling $rem ISO $ISO Read_state ave Write_state On Tissue_type homogeneous Total_time $TT Beats 10 Spatial_output_interval_data 0 Spatial_output_interval_vtk 0 D1 0.4 D_AR 7 Junction_modulation_map On Junction_modulation_map_file 2D_human_atria_fibrosis_300x300_${tissue}_JNmap_${map}.txt Direct_modulation_map   On Direct_modulation_map_file 2D_human_atria_fibrosis_map_ln_${ln}_perc_${area}.txt IK1_Erev_shift 10 reference Tissue_state_write_Fib_${map} State_Reference_write fib_FB_${map}
                    done
                done
            done
        done
    done
done

# And then run whatever data simulations you like using these maps and state files.


# 3D atrial models and creating connection maps =========
# The 3D atrial models have been provided in a few different formats. 
#   The original
#   Cleaned to remove isolated nodes, as used in this study
#   Reflected into the correct R-L. This was performed after data in this study, and therefore you will have to create your own fibrosis connection maps.
#   We will do that here

#Ensure you compile the tool for creating connection maps. 
#   "make create_connection_map"

# Now we will run it, tell it to take in the fibrosis map and remove a set proportion of axial and transverse connections within that map area
# If you want to remove connections globally, simply do not pass in the "Apply_map_file" argument

./create_connection_map Tissue_order geo Tissue_model Human_atria_AKL_heart_1_corrected Apply_map_file Human_atria_AKL_heart_1_fibrosis_map_reflected.dat axial_perc 0.2 trans_perc 0.8

# NOTE: this can take a while, as it is no parallel and is using random numbers. Only has to be done once to create the map. 

# And now we could perform simulations in 3D, using the above tissue models and created connection maps in place of the 2D tissue models used throughout the remainder of this script.
