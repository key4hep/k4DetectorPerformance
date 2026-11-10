# v00-01

* 2026-11-10 Juan Miguel Carceller ([PR#12](https://github.com/key4hep/k4DetectorPerformance/pull/12))
  - Avoid unnecessary C compiler detection when configuring the project.

* 2026-11-10 Joshua Beirer ([PR#10](https://github.com/key4hep/k4DetectorPerformance/pull/10))
  - TrackingValidation: fix the truth reference parameters and the track-particle matching, and align the matching with the ACTS definitions
  		  - Truth z0 is now taken at the closest approach to the beamline (was to the origin, wrong for vertices away from z = 0), and truth d0 is computed in double precision
  		  - Truth d0/z0/phi are taken on the helix of the particle's production vertex and momentum at its closest approach to the beamline, so they are also correct for displaced particles. The fitter trees gain `vertexR`/`vertexZ` to select prompt or displaced particles
  		  - MC particles are identified by ObjectID, so subset collections such as `MCPhysicsParticles` are handled correctly; reco hits with several contributing particles count for each of them
  		  - A track is matched only to its majority particle and only if it passes the `FinderEfficiencyDefinition` criterion, which now also selects the tracks in `fitter_vs_mc` and `fitter_vs_perfect`
  		  - `finder_particle_to_tracks` lists only tracks for which the particle is the majority particle; `finder_track_to_particles` writes -1 for particles outside the input MC collection
  		  - `g_efficiency_vs_p` is now a `TGraphAsymmErrors` with Clopper-Pearson uncertainties
  		  - The helper `truthFromMC_GenfitConvention` is renamed to `truthPerigeeFromMC`
  		  - The CI test uses all CLD tracker link collections and checks the output content

* 2026-10-21 Giovanni Marchiori ([PR#8](https://github.com/key4hep/k4DetectorPerformance/pull/8))
  - Use updated name of DCH digitizer

* 2026-08-03 Arina Ponomareva ([PR#4](https://github.com/key4hep/k4DetectorPerformance/pull/4))
  * Add TrackingValidation consumer for track-finder and track-fitter performance studies
  * Add finder-level validation through particle-to-track and track-to-particle association trees
  * Provide fitter-level validation with residuals computed with respect to MC truth
  * Support optional comparison to perfectly associated fitted tracks
  * Add summary plots written to the output ROOT file (tracking efficiency, d0 resolution, momentum resolution, transverse-momentum resolution)
  * Support configurable validation modes (full, finder-only, fitter-only)
  * Support configurable finder-efficiency definitions and purity threshold
  * Add helper and plotting utilities used by the validation workflow
  * Add test steering and end-to-end test for reconstruction and validation
  * Add documentation in the README

# v00-01

* This file is also automatically populated by the tagging script