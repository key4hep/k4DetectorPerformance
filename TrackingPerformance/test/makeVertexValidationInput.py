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

# Write a small, hand-made EDM4hep file for testVertexValidation.sh. The
# expected VertexValidation output for each event is in
# checkVertexValidation.py.
#
#   event 0: vertex 0 flagged primary, 4 truth-PV tracks with weight 1, one
#            secondary track with weight 0.5 and one with 0.05 (below
#            MinTrackWeight); vertex 1 matches a secondary truth vertex;
#            vertex 2 has only a track without truth link
#   event 1: nothing flagged primary; vertex 0 is the clean truth-PV match,
#            vertex 1 (softer) splits from it, vertex 2 is merged
#   event 2: vertex 0 (soft) matches the truth PV first, vertex 1 (harder,
#            flagged primary) takes over and vertex 0 becomes split

import sys

import edm4hep
from podio import Frame, root_io

X, Y, Z = edm4hep.Cartesian.x, edm4hep.Cartesian.y, edm4hep.Cartesian.z


class Event:
    def __init__(self):
        self.mcs = edm4hep.MCParticleCollection()
        self.tracks = edm4hep.TrackCollection()
        self.mcLinks = edm4hep.TrackMCParticleLinkCollection()
        self.particles = edm4hep.ReconstructedParticleCollection()
        self.vertices = edm4hep.VertexCollection()
        self.vertexLinks = edm4hep.VertexRecoParticleLinkCollection()

    def mc(self, pos, pt, parent=None):
        mc = self.mcs.create()
        mc.setVertex(edm4hep.Vector3d(*pos))
        mc.setMomentum(edm4hep.Vector3d(pt, 0.0, 1.0))
        if parent is not None:
            mc.addToParents(parent)
        return mc

    def particle(self, mc, pt):
        """Reconstructed particle with one track, truth-linked to mc unless None"""
        track = self.tracks.create()
        if mc is not None:
            link = self.mcLinks.create()
            link.setFrom(track)
            link.setTo(mc)
            link.setWeight(1.0)
        particle = self.particles.create()
        particle.setMomentum(edm4hep.Vector3f(pt, 0.0, 1.0))
        particle.addToTracks(track)
        return particle

    def vertex(self, pos, sigma, particles, primary=False):
        """particles: list of (particle, track weight); sigma in mm"""
        vtx = self.vertices.create()
        vtx.setPosition(edm4hep.Vector3f(*pos))
        for dim in (X, Y, Z):
            vtx.setCovMatrix(sigma * sigma, dim, dim)
        vtx.setChi2(1.0)
        vtx.setNdf(2)
        vtx.setPrimary(primary)
        for particle, weight in particles:
            vtx.addToParticles(particle)
            link = self.vertexLinks.create()
            link.setFrom(vtx)
            link.setTo(particle)
            link.setWeight(weight)

    def write(self, writer):
        frame = Frame()
        frame.put(self.mcs, "MCParticles")
        frame.put(self.tracks, "Tracks")
        frame.put(self.mcLinks, "SiTracksMCTruthLink")
        frame.put(self.particles, "ReconstructedParticles")
        frame.put(self.vertices, "PrimaryVertices")
        frame.put(self.vertexLinks, "PrimaryVertices_ParticleLinks")
        writer.write_frame(frame, "events")


def event0(writer):
    ev = Event()
    pv = (0.010, -0.020, 0.500)
    sv = (1.0, 0.5, 2.0)
    pvMCs = [ev.mc(pv, 2.0) for _ in range(4)]
    svMCs = [ev.mc(sv, 1.0, parent=pvMCs[0]) for _ in range(4)]
    pvParts = [ev.particle(mc, 2.0) for mc in pvMCs]
    svParts = [ev.particle(mc, 1.0) for mc in svMCs]
    unlinked = ev.particle(None, 1.0)
    # residual (2, 0, 3) um with sigma 2 um
    ev.vertex(
        (0.012, -0.020, 0.503),
        0.002,
        [(p, 1.0) for p in pvParts] + [(svParts[0], 0.5), (svParts[1], 0.05)],
        primary=True,
    )
    ev.vertex((1.001, 0.5, 2.0), 0.01, [(svParts[2], 1.0), (svParts[3], 1.0)])
    ev.vertex((0.5, 0.5, 0.5), 0.01, [(unlinked, 1.0)])
    ev.write(writer)


def event1(writer):
    ev = Event()
    pv = (-0.005, 0.003, -1.2)
    sv = (2.0, -1.0, 3.0)
    pvMCs = [ev.mc(pv, 5.0) for _ in range(5)]
    svMCs = [ev.mc(sv, 1.0, parent=pvMCs[0]) for _ in range(2)]
    hard = [ev.particle(mc, 10.0) for mc in pvMCs[:3]]
    soft = ev.particle(pvMCs[3], 1.0)
    mixed = [ev.particle(pvMCs[4], 1.0)] + [ev.particle(mc, 1.0) for mc in svMCs]
    # residual (0, 0, 4) um with sigma 4 um
    ev.vertex((-0.005, 0.003, -1.196), 0.004, [(p, 1.0) for p in hard])
    ev.vertex((-0.005, 0.003, -1.1), 0.01, [(soft, 1.0)])
    ev.vertex((2.0, -1.0, 3.0), 0.01, [(p, 1.0) for p in mixed])
    ev.write(writer)


def event2(writer):
    ev = Event()
    pv = (0.0, 0.0, 0.0)
    pvMCs = [ev.mc(pv, 5.0) for _ in range(3)]
    soft = ev.particle(pvMCs[0], 1.0)
    hard = [ev.particle(mc, 10.0) for mc in pvMCs[1:]]
    ev.vertex((0.0, 0.0, 0.1), 0.01, [(soft, 1.0)])
    # residual (-1, 1, 0) um with sigma 1 um
    ev.vertex((-0.001, 0.001, 0.0), 0.001, [(p, 1.0) for p in hard], primary=True)
    ev.write(writer)


if __name__ == "__main__":
    out = root_io.Writer(sys.argv[1] if len(sys.argv) > 1 else "vertex_validation_input.root")
    for make in (event0, event1, event2):
        make(out)
    # podio finishes the file at exit
