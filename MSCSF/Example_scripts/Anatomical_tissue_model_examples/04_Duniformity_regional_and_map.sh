#!/bin/sh

# Scaling the baseline D1 and D_AR by region, map and both

# The default setting for the functional_mode_test is for D to vary according to region/celltype – the same functionality as 
# is available with idealised cell models (so long as the heterogeneity settings exist). A difference here, of course, is that 
# our heterogeneity settings (tissue segmentation) is absorbed into the tissue geometry file – it either contains 1s and 0s only, 
# or further numbers for segmentation – rather than having to set the segmentation properties as with an idealised model. So 
# long as our tissue geometry file has multiple regions, we can impose electrical heterogeneity and/or D_uniformity by region. 
# However, with anatomically detailed tissue models, we can also impose D_uniformity by a map, or by both a region and map 
# (applied sequentially). The map multiples D1 and/or D_AR by the map value for every cell (so can be used to impose a D spatial 
# variability which does not depend on celltype and/or which varies continuously or checkerboard etc – anything you can make with a map!). 
# Remember, this is the baseline D1 and D_AR associated with the model – not a remodelled region which we want to change the diffusion 
# parameters (which is controlled by Dscale, D_AR_scale, and Dscale_mod_map, D_AR_scale_mod_map – see next example) – and so is 
# to correspond different diffusion parameters around the tissue model as part of its control setup. When using a map file, we must 
# provide a map file for both D1 and D_AR – if we want only one to vary spatially, just create a map of all 1s for the one we don’t want varied. 

# So let’s compare our D1 and D2 maps under uniform, regional, map and regional_map D_uniformity. 
# We’ll again set Total_time to 0 as we just need the maps for this illustration:
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 D_uniformity uniform Reference D_uniformity_uniform
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 D_uniformity regional Reference D_uniformity_regional
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 D_uniformity map Reference D_uniformity_map
./model_tissue_native Tissue_order geo Tissue_model functional_model_test Total_time 0 D_uniformity regional_map Reference D_uniformity_regional_map

# We notice that “../Map_{D1/D2}.vtk” for the uniform simulation are homogeneous (as we would hope!), whereas for the regional case, 
# they clearly vary according to celltype (compare with the segmentation in “Geometry_anatomy.vtk”). In the map simulation, we can see 
# there are now two map files: “Map_Duniformity.vtk” and “Map_D_AR_uniformity.vtk”, and we can easily see how these maps correlate with 
# the actual values of D1 and D2: remember, D2 depends on both D1 map and D_AR map, as D_AR is relative to D1. Finally, we can see how applying 
# the map and regional variation combine sequentially in the regional_map simulation example (seeing the smooth gradients in D1 in this 
# last example may require rescaling between the ranges for each celltype, as the smooth gradient is imposed within discrete celltype variation.) 

# Check "BASIC_INSTRUCTIONS_USE.txt" and Full_documentation.pdf for output file contents


