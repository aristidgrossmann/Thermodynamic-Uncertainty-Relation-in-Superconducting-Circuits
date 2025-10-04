
@echo off

SET filename=Components

SET cpp_filename=%filename%.cpp
SET exe_filename=%filename%.exe
SET compile_prompt=g++ -I../../../SDE_Solver/include ../../../SDE_Solver/src/SDE_Solver.cpp %cpp_filename% -o %exe_filename%

%compile_prompt%
%exe_filename%