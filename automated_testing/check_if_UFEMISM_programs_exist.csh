#! /bin/csh -f

set test_dir = "build/src/UFEMISM"
set expected_programs = ( \
    UFEMISM_program \
    UFEMISM_unit_test_program
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