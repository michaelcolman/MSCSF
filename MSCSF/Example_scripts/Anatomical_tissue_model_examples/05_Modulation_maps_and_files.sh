#!/bin/sh

# Applying modulation maps by coordinates or file, and specifying the filename

# As with idealised tissue models, we can apply modulation (ISO, ACh, Remodelling, Mutation, Direct_modulation, D_scale_mod, 
# D_AR_scale_mod) either homogeneous or by a map. We can set the map using coordinates, exactly as with the idealised model 
# example 9, using map_in_type coords (at least one map has to be turned on for it to do anything!):
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 ISO 1 ISO_map On map_in_type coords Reference ISO_map_coords

# And we can explicitly set the map coordinates and sizes as in that previous example. However, this is limiting in two major ways:
#   1) we are limited to discrete maps (0 or 1 for off or on, multiplied by total amount of ISO);
#   2) all modulation must be subject to the same map (if map is set to on for that modulation – we can still combine some by maps 
#   and others homogeneously).
# To overcome these limitations, we can use map files – which we can create to contain any shape and continuity we like 
# (note: Direct_modulation maps must still be discrete, as continuous imposition of these modulations is not implemented). We can 
# either apply all modulations by coordinates or files – we cannot combine coordinates and files for different modulation maps.

# Map files can be named to indicate the modulation they are designed for, but, as with stimulus maps, the actual name of the map 
# file does not limit where it can it can be used. The functional_model_test has example maps files for all modulations, and are 
# defaulted in the tissue mode settings. So let first apply all modulations, set the modulation map to “file” and set all modulations 
# to be applied by a map. We won’t pass any map files (so it reads the default settings) and will set D_uniformity to uniform (so we 
# can easily interpret the D1 and D2 maps from the modulation). We need to pass some direct modulation in order for that map to have 
# an effect. Note: the maps simply define the regions and proportion of maximum effect to which the modulation applies, so for ISO, 
# ACh and Remodelling the maximum value of the map should be = 1, which will be the region where ISO, ACh or remodelling are applied 
# according to their concentration/proportion variable (i.e. if ISO = 0.5, then the map region of 1 will apply ISO at a concentration 
# of 0.5, and the map region 0.5 will apply concentration 0.25). So passing different values of ISO, ACh and Remodelling_prop will 
# affect your simulation, but will not be reflected in the maps. Similarly, the D_scale and D_AR_scale maps must also be between 0 and 1, 
# where map=1 applies the full value of Dscale or D_AR_scale (which may be larger or smaller than 1). This is different to the baseline 
# D1 and D_AR maps (example 4), which are absolute scale factors to multiply the parameters D1 and D_AR.
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 D_uniformity uniform ISO 1 ISO_map On ACh 0.5 ACh_map On Remodelling test Remodelling_map On Remodelling_proportion 0.75 ICaL_scale 0.5 IK1_scale 2 Direct_modulation_map On Dscale 2 Dscale_mod_map On D_AR_scale 3 D_AR_scale_mod_map On map_in_type file Reference All_mod_maps_by_file

# We can see on the run-time outputs that different map files are loaded for each modulation. We can now check all the map files, and 
# that they are all different (note the Direct_modulation map is not continuous) and that the D1 and D2 maps correlate with 
# “Map_Duniformity.vtk” and “Map_D_AR_uniformity.vtk”, as in the previous example (again, remember D2 depends on D1 map and D_AR map).

# Sometimes, we may want the same map to apply to multiple modulations – for example a remodelled region could be associated with 
# electrical and diffusion dysfunction. We will use the functionality to set the map file at run-time in order to impose the same 
# map on remodelling and Dscale:
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 D_uniformity uniform Remodelling test Remodelling_map On remod_map_file functional_model_test_Dscale_mod_map.dat Dscale 0.5 Dscale_mod_map On Dscale_mod_map_file functional_model_test_Dscale_mod_map.dat Reference Remodelling_and_Dscale_same_map_file

# We can of course combine D_uniformity (by region or maps) with Dscale and D_AR_scale maps, in order to apply baseline variation in D,
# with a remodelled variation imposed on top.

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


