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
## Test for VertexValidation.
##
## Runs VertexValidation on a small hand-made EDM4hep file, with and without
## the vertex -> particle track weights, and compares the output with the
## expected classification, matching and residuals.
##

set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/vertex_validation.XXXXXX")"
trap 'rm -rf "${TMPDIR}"' EXIT

INPUT_FILE="${TMPDIR}/vertex_validation_input.root"

echo "=== Writing VertexValidation test input ==="
python3 "${SCRIPT_DIR}/makeVertexValidationInput.py" "${INPUT_FILE}"

for mode in weighted unweighted; do
  VALIDATION_FILE="${TMPDIR}/vertex_validation_${mode}.root"
  LINK_ARGS=()
  if [ "${mode}" = weighted ]; then
    LINK_ARGS=(--vertexParticleLinks PrimaryVertices_ParticleLinks)
  fi

  echo "=== Running VertexValidation (${mode}) ==="
  k4run "${SCRIPT_DIR}/runVertexValidation.py" \
    --inputFile "${INPUT_FILE}" \
    --vertices PrimaryVertices \
    "${LINK_ARGS[@]}" \
    --validationFile "${VALIDATION_FILE}"

  echo "=== Checking VertexValidation output (${mode}) ==="
  python3 "${SCRIPT_DIR}/checkVertexValidation.py" "${VALIDATION_FILE}" "${mode}"
done

echo "VertexValidation test completed successfully."
