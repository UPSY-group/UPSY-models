#!/bin/csh

# Define colors
set RED = "\033[0;31m"
set GREEN = "\033[0;32m"
set MAG = "\033[0;35m"
set RESET = "\033[0m"

# Run UPSY tests

printf "${GREEN}Starting UPSY tests... ${RESET}\n"

csh automated_testing/UPSY/update_references.csh

if ($status != 0) then
    printf "${RED}Error: UPSY tests failed with status $status ${RESET}\n"
    exit 1
else
    printf "${GREEN}UPSY tests completed successfully.${RESET}\n"
endif

# Run UFEMISM tests

printf "${MAG}Starting UFEMISM tests...${RESET}\n"

csh automated_testing/UFEMISM/update_references.csh

if ($status != 0) then
    printf "${RED}Error: UFEMISM tests failed with status $status ${RESET}\n"
    exit 1
else
    printf "${GREEN}UFEMISM tests completed successfully.${RESET}\n"
endif

