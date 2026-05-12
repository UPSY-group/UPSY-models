#!/usr/bin/env bash

set -o pipefail

version="$1"
if [[ -z "$version" ]]; then
  echo "Usage: $0 <dev|perf>" >&2
  exit 2
fi

mkdir -p .vscode/task-logs
logfile=".vscode/task-logs/ubuntu-upsy-unit-tests-${version}-$(date +%Y%m%d-%H%M%S).log"

act workflow_call \
  -W .github/workflows/UPSY_unit_tests.yml \
  --matrix os:ubuntu-latest \
  --matrix build_type:"$version" \
  --env OMPI_ALLOW_RUN_AS_ROOT=1 \
  --env OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  --container-architecture linux/amd64 \
  -P ubuntu-latest=catthehacker/ubuntu:act-latest \
  2>&1 | tee "$logfile"
