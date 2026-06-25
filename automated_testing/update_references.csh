#!/bin/csh

# Define colors
set RED = "\033[0;31m"
set GREEN = "\033[0;32m"
set MAG = "\033[0;35m"
set RESET = "\033[0m"

set os_name = `uname -s`
if ("$os_name" == "Darwin") then
    set reference_dir = "reference_macos"
else if ("$os_name" == "Linux") then
    set reference_dir = "reference_linux"
else
    printf "${RED}Error: unsupported OS \"${os_name}\".${RESET}\n"
    exit 1
endif

# Run UPSY tests

printf "${GREEN}Starting UPSY tests (target: ${reference_dir})... ${RESET}\n"

csh automated_testing/UPSY/update_references.csh

if ($status != 0) then
    printf "${RED}Error: UPSY tests failed with status $status ${RESET}\n"
    exit 1
else
    printf "${GREEN}UPSY tests completed successfully.${RESET}\n"
endif

# Run UFEMISM tests

printf "${MAG}Starting UFEMISM tests (target: ${reference_dir})...${RESET}\n"

csh automated_testing/UFEMISM/update_references.csh

if ($status != 0) then
    printf "${RED}Error: UFEMISM tests failed with status $status ${RESET}\n"
    exit 1
else
    printf "${GREEN}UFEMISM tests completed successfully.${RESET}\n"
endif

