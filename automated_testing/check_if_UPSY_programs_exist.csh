#! /bin/csh -f

set test_dir = "build/src/UPSY"
set expected_programs = ( \
    UPSY_unit_test_program \
    UPSY_multinode_unit_test_program \
    UPSY_component_test_program_mesh_creation \
    UPSY_component_test_program_mesh_discretisation \
    UPSY_component_test_program_mesh_focussing \
    UPSY_component_test_program_mesh_remapping_mesh_grid \
    UPSY_component_test_program_mesh_remapping_mesh_mesh \
    HENKIEDEPENKIE \
)

set all_exist = 1
set missing_programs = ()

foreach program ($expected_programs)
    if ( -x "$test_dir/$program" ) then
        echo "Found: $program"
    else
        echo "Missing: $program"
        set all_exist = 0
        set missing_programs = ($missing_programs $program)
    endif
end

if ( $all_exist == 1 ) then
    exit 0
else
    exit 1
endif