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
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

MODEL_FILE="${1:-}"
ENERGY_GEV="${2:-15}"
THETA_DEG="${3:-20}"
N_EVENTS="${4:-10}"
OUT_TAG="${5:-${ENERGY_GEV}GeV}"

TRACKINGPERF_RUN_SIM="${TRACKINGPERF_RUN_SIM:-1}"
TRACKINGPERF_INPUT_FILE_OVERRIDE="${TRACKINGPERF_INPUT_FILE_OVERRIDE:-}"
TRACKINGPERF_KEEP_TMP="${TRACKINGPERF_KEEP_TMP:-0}"

if [ -z "${MODEL_FILE}" ]; then
  echo "ERROR: missing ONNX model path argument"
  echo "Usage: $0 /full/path/to/SimpleGatrIDEAv3o1.onnx ENERGY_GEV THETA_DEG N_EVENTS OUT_TAG"
  exit 1
fi

if [ ! -f "${MODEL_FILE}" ]; then
  echo "ERROR: ONNX model file not found: ${MODEL_FILE}"
  exit 1
fi

if [ -z "${K4GEO:-}" ]; then
  echo "ERROR: K4GEO is not set"
  exit 1
fi

if [ ! -d "${K4GEO}" ]; then
  echo "ERROR: K4GEO does not point to a valid directory: ${K4GEO}"
  exit 1
fi

if ! command -v k4run >/dev/null 2>&1; then
  echo "ERROR: k4run not found in PATH"
  exit 1
fi

if [ "${TRACKINGPERF_RUN_SIM}" -eq 1 ]; then
  if ! command -v ddsim >/dev/null 2>&1; then
    echo "ERROR: ddsim not found in PATH"
    exit 1
  fi

  if ! command -v curl >/dev/null 2>&1; then
    echo "ERROR: curl not found in PATH"
    exit 1
  fi
fi

VALIDATION_FILE="${VALIDATION_FILE:-validation.root}"

XML_FILE="${K4GEO}/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"
RUN_FILE="${SCRIPT_DIR}/runTrackingValidation.py"

VAL_FILE="${VALIDATION_FILE}"
LOG_FILE="${LOG_FILE:-run_${OUT_TAG}.log}"

if [ ! -f "${XML_FILE}" ]; then
  echo "ERROR: geometry XML file not found: ${XML_FILE}"
  exit 1
fi

if [ ! -f "${RUN_FILE}" ]; then
  echo "ERROR: tracking validation run file not found: ${RUN_FILE}"
  exit 1
fi

TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/tracking_validation_${OUT_TAG}.XXXXXX")"
echo "TMPDIR = ${TMPDIR}"

if [ "${TRACKINGPERF_KEEP_TMP}" -eq 1 ]; then
  echo "Keeping temporary directory after exit"
else
  trap 'rm -rf "${TMPDIR}"' EXIT
fi

STEERING_FILE="${TMPDIR}/SteeringFile_IDEA_o1_v03.py"
SIM_FILE="${TMPDIR}/out_sim_edm4hep.root"
RECO_FILE="${TMPDIR}/out_reco.root"

rm -f "${VAL_FILE}"

SEED=42

if [ "${TRACKINGPERF_RUN_SIM}" -eq 1 ]; then
  INPUT_FILE="${SIM_FILE}"
else
  if [ -z "${TRACKINGPERF_INPUT_FILE_OVERRIDE}" ]; then
    echo "ERROR: TRACKINGPERF_RUN_SIM=0 but TRACKINGPERF_INPUT_FILE_OVERRIDE is empty"
    echo "Please provide an existing EDM4hep file, e.g."
    echo "  TRACKINGPERF_RUN_SIM=0 TRACKINGPERF_INPUT_FILE_OVERRIDE=/path/to/input.root $0 /path/to/model.onnx ${ENERGY_GEV} ${THETA_DEG} ${N_EVENTS} ${OUT_TAG}"
    exit 1
  fi

  if [ ! -f "${TRACKINGPERF_INPUT_FILE_OVERRIDE}" ]; then
    echo "ERROR: TRACKINGPERF_INPUT_FILE_OVERRIDE does not exist: ${TRACKINGPERF_INPUT_FILE_OVERRIDE}"
    exit 1
  fi

  INPUT_FILE="${TRACKINGPERF_INPUT_FILE_OVERRIDE}"
fi

echo "=== Test configuration ==="
echo "Script dir:                 ${SCRIPT_DIR}"
echo "Temporary dir:              ${TMPDIR}"
echo "Geometry XML:               ${XML_FILE}"
echo "Run script:                 ${RUN_FILE}"
echo "ONNX model:                 ${MODEL_FILE}"
echo "TRACKINGPERF_RUN_SIM:       ${TRACKINGPERF_RUN_SIM}"
echo "TRACKINGPERF_KEEP_TMP:      ${TRACKINGPERF_KEEP_TMP}"
echo "Energy [GeV]:               ${ENERGY_GEV}"
echo "Theta [deg]:                ${THETA_DEG}"
echo "Events:                     ${N_EVENTS}"
echo "Seed:                       ${SEED}"
echo "Input file:                 ${INPUT_FILE}"
echo "Reco file:                  ${RECO_FILE}"
echo "Validation file:            ${VAL_FILE}"
echo "Log file:                   ${LOG_FILE}"
echo "runDigi:                    1"
echo "runFinder:                  1"
echo "runFitter:                  1"
echo "runPerfectTracking:         0"
echo "runValidation:              1"
echo "useDCH:                     1"
echo "mode:                       0"
echo "doPerfectFit:               0"
echo "finderEfficiencyDefinition: 1"
echo "finderPurityThreshold:      0.75"

if [ "${TRACKINGPERF_RUN_SIM}" -eq 1 ]; then
  echo "=== Downloading DDSim steering file ==="
  curl -L \
    -o "${STEERING_FILE}" \
    https://raw.githubusercontent.com/key4hep/k4geo/master/example/SteeringFile_IDEA_o1_v03.py

  if [ ! -f "${STEERING_FILE}" ]; then
    echo "ERROR: failed to download DDSim steering file"
    exit 1
  fi

  echo "=== Step 1: DDSim ==="
  ddsim \
    --steeringFile "${STEERING_FILE}" \
    --compactFile "${XML_FILE}" \
    -G \
    --gun.particle mu- \
    --gun.energy "${ENERGY_GEV}*GeV" \
    --gun.distribution uniform \
    --gun.thetaMin "${THETA_DEG}*deg" \
    --gun.thetaMax "${THETA_DEG}*deg" \
    --gun.phiMin "0*deg" \
    --gun.phiMax "0*deg" \
    --random.seed "${SEED}" \
    --numberOfEvents "${N_EVENTS}" \
    --outputFile "${SIM_FILE}"

  if [ ! -f "${SIM_FILE}" ]; then
    echo "ERROR: simulation output was not created: ${SIM_FILE}"
    exit 1
  fi
else
  echo "=== Step 1: DDSim skipped (TRACKINGPERF_RUN_SIM=0) ==="
fi

echo "=== Step 2: digitization + tracking + validation ==="
k4run "${RUN_FILE}" \
  --inputFile "${INPUT_FILE}" \
  --modelPath "${MODEL_FILE}" \
  --outputFile "${RECO_FILE}" \
  --validationFile "${VAL_FILE}" \
  --compactFile "${XML_FILE}" \
  --runDigi 1 \
  --runFinder 1 \
  --runFitter 1 \
  --runPerfectTracking 0 \
  --runValidation 1 \
  --useDCH 1 \
  --mode 0 \
  --doPerfectFit 0 \
  --finderEfficiencyDefinition 1 \
  --finderPurityThreshold 0.75
if [ ! -f "${RECO_FILE}" ]; then
  echo "ERROR: reconstruction output was not created: ${RECO_FILE}"
  exit 1
fi

if [ ! -f "${VAL_FILE}" ]; then
  echo "ERROR: validation output was not created: ${VAL_FILE}"
  exit 1
fi

echo "=== Step 3: check outputs ==="
if [ "${TRACKINGPERF_RUN_SIM}" -eq 1 ]; then
  test -f "${SIM_FILE}"
fi
test -f "${RECO_FILE}"
test -f "${VAL_FILE}"

echo "Test completed successfully."
echo "Validation file: ${VAL_FILE}"
