# MSCSF
The full package for the "multi-scale cardiac simulation framework" in C/C++

23/06/2019 - Fix: Small issue with writing the settings file for all implementations other than single native, so it writes now to correct folder.
             [E.g. output_settings_tissue(Sim, Tissue, directory) -> output_settings_tissue(Sim, Tissue, res_dir_full) in Tissue_native_main.cc and similar equivalents for all main files other than Single_cell_native.cc]

24/06/2019 - Documentation for installing the package on Windows computers has been included (thank you to Jakub Tomek for this implementation and associated documentation). This is provided ahead of a full version of the framework with included Visual Studio project files, which I intend to finish and upload soon.  

20/06/2019 - Release version uploaded, with accompanying author revised version of the manuscript (currently under review at PLOS Comp. Biol).


