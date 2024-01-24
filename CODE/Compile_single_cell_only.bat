echo.
PATH C:\Users\fbsmac\Documents\MinGW\bin
:: Single cell: native (standard non-spatial)
g++ Single_cell_native_main.cc lib/Arguments.c lib/Initialisation.c  lib/Model.c lib/Model*.cpp lib/Read_write_state.c lib/Outputs.cpp -o model_single_cell_native.exe
