#! /bin/csh -f

cp ../../../UPSY_component_test_program_mesh_creation .
mpiexec -n 2 UPSY_component_test_program_mesh_creation

mv  test_meshes_and_grids  results