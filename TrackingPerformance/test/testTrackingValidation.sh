#!/bin/bash
##
## Copyright (c) 2020-2024 Key4hep-Project.
##
## This file is part of Key4hep.
## See https://key4hep.github.io/key4hep-doc/ for further info.
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##

##
## Test for TrackingValidation.
##
## This test intentionally does not run DDSim, digitization, track finding,
## fitting, or perfect tracking. It runs only TrackingValidation on a
## pre-produced CLD reconstruction file provided through CMake ExternalData.
##

set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

RUN_FILE="${SCRIPT_DIR}/runTrackingValidation.py"

INPUT_FILE="${1:-}"

CLD_XML="${TRACKINGPERF_CLD_COMPACT_FILE:-${K4GEO}/FCCee/CLD/compact/CLD_o3_v01/CLD_o3_v01.xml}"

VALIDATION_FILE="${VALIDATION_FILE:-validation_output_test.root}"


TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/tracking_validation_cld.XXXXXX")"
trap 'rm -rf "${TMPDIR}"' EXIT

OUTPUT_FILE="${TMPDIR}/out_validation_only_cld.root"

echo "=== TrackingValidation test ==="
echo "Working directory: $(pwd)"
echo "Run file:         ${RUN_FILE}"
echo "Input file:       ${INPUT_FILE}"
echo "Compact file:     ${CLD_XML}"
echo "Validation file:  ${VALIDATION_FILE}"
echo "Temporary EDM output file: ${OUTPUT_FILE}"

if [ -z "${INPUT_FILE}" ]; then
  echo "ERROR: no input file provided"
  echo "Usage: $0 /path/to/MuGuns_CLD_o3_v01_2026_07_01.root"
  exit 1
fi

if [ ! -f "${RUN_FILE}" ]; then
  echo "ERROR: steering file not found: ${RUN_FILE}"
  exit 1
fi

if [ ! -f "${INPUT_FILE}" ]; then
  echo "ERROR: input file not found: ${INPUT_FILE}"
  exit 1
fi

if [ ! -f "${CLD_XML}" ]; then
  echo "ERROR: CLD compact file not found: ${CLD_XML}"
  exit 1
fi

if ! command -v k4run >/dev/null 2>&1; then
  echo "ERROR: k4run not found in PATH"
  exit 1
fi

rm -f "${VALIDATION_FILE}"

echo "=== Running TrackingValidation ==="

k4run "${RUN_FILE}" \
  --inputFile "${INPUT_FILE}" \
  --outputFile "${OUTPUT_FILE}" \
  --validationFile "${VALIDATION_FILE}" \
  --compactFile "${CLD_XML}" \
  --runDigi 0 \
  --runFinder 0 \
  --runFitter 0 \
  --runPerfectTracking 0 \
  --runValidation 1 \
  --useDCH 0 \
  --mode 0 \
  --doPerfectFit 0 \
  --mcParticles "MCPhysicsParticles" \
  --hitSimLinks "VXDTrackerHitRelations" \
  --finderTracks "SiTracks" \
  --fittedTracks "FittedTracks" \
  --finderEfficiencyDefinition 1 \
  --finderPurityThreshold 0.75

echo "=== Checking validation output ==="

if [ ! -f "${VALIDATION_FILE}" ]; then
  echo "ERROR: validation output was not created: ${VALIDATION_FILE}"
  exit 1
fi

if [ ! -s "${VALIDATION_FILE}" ]; then
  echo "ERROR: validation output is empty: ${VALIDATION_FILE}"
  exit 1
fi

echo "Validation test completed successfully."
echo "Validation file: ${VALIDATION_FILE}"