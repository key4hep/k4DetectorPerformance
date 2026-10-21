#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from pyLCIO import EVENT, UTIL, IOIMPL, IMPL

import podio
from podio.root_io import Writer
import edm4hep

import math
import random
from array import array
import argparse

import numpy as np

### Author: P.McKeown, CERN, Aug 2024
### Adapted File to create MC Particles for passing to ddsim
### Currently designed to create two particles (e.g. photons) orthogonal to the face of the ECAL Barrel, in this case of ILD.
### The separation between the two particles is varied uniformly in the range [0, maxSep], with the midpoint between the two particles along the global Z axis being set by centralPosZ
### The energy of each particle is varied uniformly (and independently) between Emax and Emin

def get_parser():
    parser = argparse.ArgumentParser(
        description='Generation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--pdg', action='store',
                        type=int, default=22,
                        help='Particle Type')

    parser.add_argument('--mass', action='store',
                        type=int, default=0,
                        help='Mass')

    parser.add_argument('--charge', action='store',
                        type=int, default=0,
                        help='Charge')

    '''
    parser.add_argument('--angleMin_theta', action='store',
                        type=int, default=30,
                        help='Minimum theta angle of incident particle input: deg, converted to [Rad]')

    parser.add_argument('--angleMax_theta', action='store',
                        type=int, default=90,
                        help='Maximum theta angle of incident particle input:deg, converted to [Rad]')

    parser.add_argument('--angleMin_phi', action='store',
                        type=int, default=30,
                        help='Minimum phi angle of incident particle input: deg, converted to [Rad]')

    parser.add_argument('--angleMax_phi', action='store',
                        type=int, default=90,
                        help='Maximum phi angle of incident particle input:deg, converted to [Rad]')
    '''

    parser.add_argument('--maxSep', action='store',
                        type=int, default=90,
                        help='Maximum separation between the two photons created [mm]')
    parser.add_argument('--centralPosX', action='store',
                        type=int, default=0,
                        help='Central position of the patch at the global x axis [mm]')
    parser.add_argument('--centralPosY', action='store',
                        type=int, default=1804.7,
                        help='Central position of the patch at the global y axis [mm]')
    parser.add_argument('--centralPosZ', action='store',
                        type=int, default=150,
                        help='Central position (i.e. mid-way point) between two incident photons along the global z axis [mm]')
    parser.add_argument('--Emax', action='store',
                        type=int, default=5,
                        help='Maximum energy of incident particle [GeV]')
    parser.add_argument('--Emin', action='store',
                        type=int, default=5,
                        help='Minimum energy of incident particle [GeV]')
    parser.add_argument('--output', action='store',
                        type=str, default='Di-particle_MC_Particles',
                        help='output edm4hep file name')

    return parser


def initialize_writer(outfile):
    wrt = IOIMPL.LCFactory.getInstance().createLCWriter()
    wrt.open(outfile, EVENT.LCIO.WRITE_NEW)
    print("Opened outfile:", outfile)
    return wrt

def write_run_header(wrt, pdg, charge, mass):
    run = IMPL.LCRunHeaderImpl()
    run.setRunNumber(0)
    run.parameters().setValue("Generator", "${lcgeo}_DIR/examples/lcio_particle_gun.py")
    run.parameters().setValue("PDG", pdg)
    run.parameters().setValue("Charge", charge)
    run.parameters().setValue("Mass", mass)
    wrt.writeRunHeader(run)

def initialize_event(j):
    col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE)
    evt = IMPL.LCEventImpl()
    evt.setEventNumber(j)
    evt.addCollection(col, "MCParticle")
    return evt, col

def calculate_energy_and_momentum(mass, p, position, y_dir):
    x, z = position
    # Compute direction from the local (x,z) coordinates.
    direction = np.array([x, y_dir, z], dtype=float)
    direction = direction / np.linalg.norm(direction)
    energy = math.sqrt(mass * mass + p * p)
    momentum = array('f', [p * direction[0], p * direction[1], p * direction[2]])
    return energy, momentum

def sample_points(radius, n_points, sep):
    points = []
    for _ in range(n_points):
        # Sample a random center within the radius
        angle = np.random.uniform(0, 2 * np.pi)
        center_dist = np.sqrt(np.random.uniform(0, 1)) * radius
        center_x = center_dist * np.cos(angle)
        center_y = center_dist * np.sin(angle)
        center = np.array([center_x, center_y])

        # Sample a random distance between 0 and sep
        dist = np.random.uniform(0, sep)

        # Calculate the two points equidistant from the center along a line
        dx = dist / 2
        point1 = np.array([center_x + dx, center_y])
        point2 = np.array([center_x - dx, center_y])

        # Rotate the points around the center by a random angle
        rotation_angle = np.random.uniform(0, 2 * np.pi)
        rotation_matrix = np.array([
            [np.cos(rotation_angle), -np.sin(rotation_angle)],
            [np.sin(rotation_angle),  np.cos(rotation_angle)]
        ])
        point1 = center + np.dot(rotation_matrix, point1 - center)
        point2 = center + np.dot(rotation_matrix, point2 - center)

        points.append((point1, point2))

    return np.array(points)  # shape (n_points, 2, 2)

def create_vertex_and_endpoint(positions, centralPosY, centralPosZ):
    vy = centralPosY

    vx_1, vz_1 = positions[0]
    vx_2, vz_2 = positions[1]

    vertex_1 = array('d', [vx_1, vy, vz_1])
    endpoint_1 = array('d', [vx_1, vy, vz_1])

    vertex_2 = array('d', [vx_2, vy, vz_2])
    endpoint_2 = array('d', [vx_2, vy, vz_2])

    return vertex_1, endpoint_1, vertex_2, endpoint_2

def create_mcparticle(genstat, mass, pdg, momentum, charge, vertex, endpoint):
    mcp = IMPL.MCParticleImpl()
    mcp.setGeneratorStatus(genstat)
    mcp.setMass(mass)
    mcp.setPDG(pdg)
    mcp.setMomentum(momentum)
    mcp.setCharge(charge)
    mcp.setVertex(vertex)
    mcp.setEndpoint(endpoint)
    return mcp

def write_to_lcio(outfile, nevt, pdg, mass, charge, maxSep, centralPosX, centralPosY, centralPosZ, Emax, Emin):
    wrt = initialize_writer(outfile)
    random.seed()

    genstat = 1
    write_run_header(wrt, pdg, charge, mass)

    positions = sample_points(10, nevt, maxSep)
    positions[:, :, 0] = positions[:, :, 0] + centralPosX
    positions[:, :, 1] = positions[:, :, 1] + centralPosZ

    for j in range(nevt):
        evt, col = initialize_event(j)

        # Sample energies (or momentum magnitudes) for each particle
        p1, p2 = random.uniform(Emin, Emax), random.uniform(Emin, Emax)

        energy1, momentum1 = calculate_energy_and_momentum(mass, p1, positions[j][0])
        energy2, momentum2 = calculate_energy_and_momentum(mass, p2, positions[j][1])

        vertex1, endpoint1, vertex2, endpoint2 = create_vertex_and_endpoint(positions[j], centralPosY, centralPosZ)

        mcp1 = create_mcparticle(genstat, mass, pdg, momentum1, charge, vertex1, endpoint1)
        mcp2 = create_mcparticle(genstat, mass, pdg, momentum2, charge, vertex2, endpoint2)

        col.addElement(mcp1)
        col.addElement(mcp2)

        wrt.writeEvent(evt)

        # Compute normalized direction vectors from the momentum vectors.
        norm1 = np.linalg.norm(momentum1)
        norm2 = np.linalg.norm(momentum2)
        if norm1 != 0:
            dir1 = np.array(momentum1, dtype=float) / norm1
        else:
            dir1 = np.array([0.0, 0.0, 0.0])
        if norm2 != 0:
            dir2 = np.array(momentum2, dtype=float) / norm2
        else:
            dir2 = np.array([0.0, 0.0, 0.0])

    wrt.close()
    return 0

def initialize_writer_edm4hep(outfile):
    writer = Writer(outfile)
    print("Created writer for outfile: ", outfile)
    return writer

def configure_frame(writer):
    frame = podio.Frame()

    frame.put_parameter(key="Run", value=0, as_type="int")
    frame.put_parameter(key="File", value=__file__, as_type="str")

    return  frame

def initialize_event_edm4hep(j):
    evth_col = edm4hep.EventHeaderCollection()
    evth = evth_col.create()
    evth.setEventNumber(j)
    return evth_col

def create_mcparticle_edm4hep(genstat, mass, pdg, momentum, charge, vertex, endpoint, MC_col):
    mcp = MC_col.create()
    mcp.setGeneratorStatus(genstat)
    mcp.setMass(mass)
    mcp.setPDG(pdg)
    mcp.setMomentum(edm4hep.Vector3d(momentum[0], momentum[1], momentum[2]))
    mcp.setCharge(charge)
    mcp.setVertex(edm4hep.Vector3d(vertex[0], vertex[1], vertex[2]))
    mcp.setEndpoint(edm4hep.Vector3d(endpoint[0], endpoint[1], endpoint[2]))
    return mcp, MC_col


def write_to_edm4hep(outfile, nevt, pdg, mass, charge, maxSep, centralPosX, centralPosY, centralPosZ, Emax, Emin):
    writer = initialize_writer_edm4hep(outfile)
    random.seed()

    genstat = 1

    positions = sample_points(10, nevt, maxSep) ## ToDo: make radius(10) a configurable parameter
    ## Place in global coordinate system
    positions[:, :, 0] = positions[:, :, 0] + centralPosX
    positions[:, :, 1] = positions[:, :, 1] + centralPosZ

    for j in range(nevt):

        frame= configure_frame(writer)
        evth_col = initialize_event_edm4hep(j)
        MC_col = edm4hep.MCParticleCollection()

        # Sample energies (or momentum magnitudes) for each particle
        p1, p2 = random.uniform(Emin, Emax), random.uniform(Emin, Emax)
        energy1, momentum1 = calculate_energy_and_momentum(mass, p1, positions[j][0], centralPosY)
        energy2, momentum2 = calculate_energy_and_momentum(mass, p2, positions[j][1], centralPosY)

        vertex1, endpoint1, vertex2, endpoint2 = create_vertex_and_endpoint(positions[j], centralPosY, centralPosZ)

        mcp1 = create_mcparticle_edm4hep(genstat, mass, pdg, momentum1, charge, vertex1, endpoint1, MC_col)
        mcp2 = create_mcparticle_edm4hep(genstat, mass, pdg, momentum2, charge, vertex2, endpoint2, MC_col)

        frame.put(evth_col, "EventHeader")
        frame.put(MC_col, "MCParticles")

        writer.write_frame(frame, "events")

    return 0


if __name__ == "__main__":

    parser = get_parser()
    parse_args = parser.parse_args()

    pdgid = parse_args.pdg
    mass = parse_args.mass
    charge = parse_args.charge

    mSep = parse_args.maxSep
    cPosX = parse_args.centralPosX
    cPosY = parse_args.centralPosY
    cPosZ = parse_args.centralPosZ

    engrMax = parse_args.Emax
    engrMin = parse_args.Emin

    output = parse_args.output

    np.random.seed(42)

    for i in range(1, 2):
        filename = f"{output}"
        write_to_edm4hep(filename+'.edm4hep.root', 100, pdgid, mass, charge, mSep, cPosX, cPosY, cPosZ, engrMax, engrMin)
