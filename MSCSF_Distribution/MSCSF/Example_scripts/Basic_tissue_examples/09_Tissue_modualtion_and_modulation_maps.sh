#!/bin/sh

# Cellular modulation (ISO, Remodelling, ICaL_scale etc) can be applied to tissue models exactly as with single-cell models. 
# By default, they apply homogeneously to the tissue, on top of any local heterogeneity due to celltype. 
# We can also apply tissue modulation – i.e. Dscale and/or D_AR_scale, to multiply the values of D1 and D2 or the anisotropy ratio, relative to the set values.

# Let’s create and control and modulation simulation: 
./model_tissue_native Tissue_order 2D Tissue_model functional_model_test Total_time 200 Reference Tissue_control
./model_tissue_native Tissue_order 2D Tissue_model functional_model_test Total_time 200 ISO 1 ICaL_scale 0.5 Dscale 2 D_AR_scale 4 Reference Tissue_modulation


# Wec an also define whether to apply the modulation homogeneously, or by a spatial map. Whereas we can individually turn the modulation map on or off for all 
# different modulations (Dscale_mod_map, D_AR_mod_map, Remodelling_map, ISO_map, ACh_map, Direct_modulation_map [Off/On]), we can only use a single map itself 
# for idealised models. If you require different maps for different modulations, you must use an anatomically detailed model and provide map files. 
# Note: “Direction_modulation_map” will control where any direct_modulations are applied (ICaL_scale, INa_va_shift, IKur_vi_tau_scale etc). 

# Let’s apply Dscale and ISO by a spatial map, but leave D_AR_scale and ICaL_scale to apply homogeneously. 
# Maps are defaulted to off, but we will include explicit settings for each map here just for clarity. 
./model_tissue_native Tissue_order 2D Tissue_model functional_model_test Total_time 200 ISO 1 ICaL_scale 0.5 Dscale 2 D_AR_scale 4 Dscale_mod_map On D_AR_scale_mod_map Off ISO_map On Direct_modulation_map Off Reference Tissue_modulation_map

# Note that in “Outputs_tissue_native_Tissue_modulation_map” there are now two additional files: Map_Dmod .vtk and Map_ISO.vtk, each with a region of 1s and a 
# region of 0s. This number is what determines the magnitude of the modulation – so where it is 1, ISO will be set to the value of ISO and Dscale_mod will be 
# set to the value of Dscale. Where it is 0, ISO = 0 and Dscale = 0 and thus parameters = control. You can also visualise D1 (D1_map.vtk) to see if the 
# modulation map has the intended effect. 

# We can control the shape, size and location of the map at run-time in an exactly analogous way to setting the stimulus spatial properties, by using map_shape, 
# map_x_loc, map_x_size etc. As with stimulus, this sets the centre and size from centre to edge; spherical maps can only have an x_size which determines radius.
./model_tissue_native Tissue_order 2D Tissue_model functional_model_test Total_time 200 ISO 1 ICaL_scale 0.5 Dscale 2 D_AR_scale 4 Dscale_mod_map On D_AR_scale_mod_map Off ISO_map On Direct_modulation_map Off Reference Tissue_modulation_map_controlled map_shape sphere map_x_loc 90 map_x_size 5 map_y_loc 7

# Finally, let’s also turn D_uniformity to regional to demonstrate how regional D and Dscale with a map can combine:
./model_tissue_native Tissue_order 2D Tissue_model functional_model_test Total_time 0 ISO 1 ICaL_scale 0.5 Dscale 2 D_AR_scale 4 Dscale_mod_map On D_AR_scale_mod_map Off ISO_map On Direct_modulation_map Off Reference Tissue_modulation_map_and_regional_D map_shape sphere map_x_loc 90 map_x_size 35 map_y_loc 7 D_uniformity regional 

# And you may notice now that D1_map.vtk has clear regional variation, and an additional area where the modulation map is applied (we made the map larger so we can see it affecting different regions).

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


