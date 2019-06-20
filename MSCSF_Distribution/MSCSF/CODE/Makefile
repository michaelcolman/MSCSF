# Compiler and flags
CC = clang++#g++
CFLAGS = -O3 -w #-std=c++11
CFLAGS2 = -Xpreprocessor -fopenmp -lomp

# build options
all: single_native tissue_native single_3D single_0D tissue_0D Ca_clamp_0D Ca_clamp_3D bin_to_vtk_dat_tissue bin_to_vtk_dat_3Dcell

# compilation file lists
common = lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp
SC = lib/Spatial_coupling.cpp
tissue = lib/Tissue.cpp
spatial_Ca = lib/CRU.cpp lib/myofilament.cpp
SRF = lib/Spontaneous_release_functions.cpp
dyad = lib/Single_dyad.cpp

# Compile
single_native: $(common) Single_cell_native_main.cc
	$(CC) $(CFLAGS) -o model_single_native $(common) Single_cell_native_main.cc

tissue_native: $(common) $(SC) $(tissue) Tissue_native_main.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_tissue_native $(common) $(SC) $(tissue) Tissue_native_main.cc

single_3D: $(common) $(SC) $(spatial_Ca) Single_cell_3D_main.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_single_3D $(common) $(SC) $(spatial_Ca) Single_cell_3D_main.cc

single_0D: $(common) $(spatial_Ca) Single_cell_0D_main.cc
	$(CC) $(CFLAGS) -o model_single_0D $(common) $(SC) $(spatial_Ca) $(SRF) Single_cell_0D_main.cc

tissue_0D: $(common) $(SC) $(tissue) Tissue_integrated_main.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_tissue_0D $(common) $(SC) $(tissue) $(spatial_Ca) $(SRF) Tissue_integrated_main.cc

Ca_clamp_0D: $(common) $(spatial_Ca) Single_cell_Ca_clamp_0D.cc
	$(CC) $(CFLAGS) -o model_Ca_clamp_0D $(common) $(SC) $(spatial_Ca) $(SRF) Single_cell_Ca_clamp_0D.cc

Ca_clamp_3D: $(common) $(SC) $(spatial_Ca) Single_cell_Ca_clamp_3D.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o model_Ca_clamp_3D $(common) $(SC) $(spatial_Ca) Single_cell_Ca_clamp_3D.cc

bin_to_vtk_dat_tissue: $(common) $(SC) $(tissue) Data_convert_binary_to_vtk_text_tissue.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o bin_to_vtk_tissue $(common) $(SC) $(tissue) Data_convert_binary_to_vtk_text_tissue.cc

bin_to_vtk_dat_3Dcell: $(common) $(SC) $(spatial_Ca) Data_convert_binary_to_vtk_text_3Dcell.cc
	$(CC) $(CFLAGS) $(CFLAGS2) -o bin_to_vtk_3Dcell $(common) $(SC) $(spatial_Ca) Data_convert_binary_to_vtk_text_3Dcell.cc

clean:
	rm model_*
