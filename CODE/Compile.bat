echo.
PATH C:\Users\fbsmac\Documents\MinGW\bin
:: Single cell: native (standard non-spatial)
g++ Single_cell_native_main.cc lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp -o model_single_cell_native.exe

:: Tissue native: Note: no parallelisation here -> add open MP yourself to this compile line if you have it installed (it is suggested you do install it)
g++ Tissue_native_main.cc lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp lib/Spatial_coupling.cpp lib/Tissue.cpp -o model_tissue_native.exe

:: Tissue network: Note: no parallelisation here -> add open MP yourself to this compile line if you have it installed (it is suggested you do install it)
g++ Tissue_native_network_main.cc lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp lib/Spatial_coupling.cpp lib/Tissue.cpp -o model_tissue_network.exe

:: Single cell: spatial cell
g++ Single_cell_3D_main.cc lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp lib/Spatial_coupling.cpp lib/CRU.cpp lib/myofilament.cpp -o model_single_cell_3D.exe

:: Single cell: non-spatial reduction of spatial cell (for spontaneous release functions)
g++ Single_cell_0D_main.cc lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp lib/Spatial_coupling.cpp lib/CRU.cpp lib/myofilament.cpp lib/Spontaneous_release_functions.cpp -o model_single_cell_0D.exe

g:: Single cell: spatial cell -> Ca clamp
g++ Single_cell_Ca_clamp_3D.cc lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp lib/Spatial_coupling.cpp lib/CRU.cpp lib/myofilament.cpp -o model_Ca_clamp_3D.exe

:: Single cell: non-spatial reduction of spatial cell (for spontaneous release functions) -> Ca clamp
g++ Single_cell_Ca_clamp_0D.cc lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp lib/Spatial_coupling.cpp lib/CRU.cpp lib/myofilament.cpp lib/Spontaneous_release_functions.cpp -o model_Ca_clamp_0D.exe

:: Tissue integrated for spontanoeus release
g++ Tissue_integrated_main.cc lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp lib/Spatial_coupling.cpp lib/Tissue.cpp lib/CRU.cpp lib/myofilament.cpp ib/Spontaneous_release_functions.cpp -o model_tissue_0D.exe

:: Tissue integrated for spontanoeus release - network model
g++ Tissue_integrated_network.cc lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp lib/Spatial_coupling.cpp lib/Tissue.cpp lib/CRU.cpp lib/myofilament.cpp ib/Spontaneous_release_functions.cpp -o model_tissue_0D_network.exe
