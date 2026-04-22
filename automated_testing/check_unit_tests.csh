#!/bin/csh

# Script to check unit test results and fail if any tests failed

set output_file = "automated_testing/UPSY/unit_tests/results/unit_tests_output.txt"

# Check if the file exists
if (! -e $output_file) then
    echo "Error: Unit test output file '$output_file' not found."
    exit 1
endif

# Find failed tests
set failed_count = `grep -c "^ Unit test failed:" $output_file`

if ($failed_count > 0) then
    echo "The following unit tests failed:"
    grep "^ Unit test failed:" $output_file
    echo ""
    echo "Unit tests failed. Exiting with error code 1."
    exit 1
else
    echo "All unit tests passed."
    exit 0
endif