# SU(3) plaquette-oracle circuits

Companion artifacts to *Block Encoding Non-Abelian Lattice Gauge Theory*
(arXiv:2608.17115). This package contains
the compiled magnetic-plaquette block-encoding oracle for SU(3) lattice
gauge theory at truncation B = 6 with b_rot = 38 rotation bits: six
explicit circuits in OpenQASM 3, the transition tables the circuits load,
the ground-truth data used to check them, and a verifier that runs the
check.

## Layout

```
qasm/       six OpenQASM 3 circuits and their manifests
tables/     per-plane transition tables and shared constants
data/       ground-truth matrix elements and near-miss targets
verify.py   verifies the circuits against data/
```

## qasm/

One circuit per plane and direction: `up_{12,13,23}.qasm.gz` implement the
oracle U_p on the three plaquette planes, and `up_dag_{12,13,23}.qasm.gz`
implement its adjoint. The adjoint is not a mechanical gate reversal of the
forward circuit, since the uncompute gadgets use mid-circuit measurement
with classical feedforward; the six files are separate compilations.

Each file is gzip-compressed OpenQASM 3.0 text using `stdgates.inc`, about
485,000 statements. Registers are one flat `qubit[1055] q`, a
`qubit[37] car` carry register for the phase-gradient adder, and a classical
register for the discarding measurements. Four subroutines carry the
recurring gadgets:

- `temp_and`: the AND-compute (one `ccx`; accounted at 4 T),
- `unand`: the measurement-based uncompute, an X-basis measurement and a
  classically controlled `cz`, 0 T,
- `add_phi`: the in-place ripple-carry adder that realizes every
  data-dependent rotation as an addition into the phase-gradient register,
- `gradient_prep`: the phase-gradient state preparation, run once per
  circuit; its 38 `rz(2*pi*2^j/2^38)` gates are the only rotation gates in
  the artifact.

`qasm_provenance.json` is a per-circuit manifest: the operation census, the
T count under the accounting stated in the paper, register widths, and file
checksums.

`qasm_readout_maps.json` holds, for each circuit, the register wire map
(which physical wires carry the active links, control links, and
multiplicities of the input and output basis states) and the injection
points where an input state is written before replay. `verify.py` uses
this file to build and read off basis states without touching the circuit
text.

## tables/

`b6_plane12.json`, `b6_plane13.json`, `b6_plane23.json` hold the per-plane
transition tables (proposal rows, corner data, class labels) and
`b6_shared.json` the shared constants. These are the data the compiler
loaded at build time to produce the circuits, and the values behind the
circuits' lookups.

## data/

`b6_compose_sectors_12.json`, `b6_compose_sectors_13.json`, and
`b6_compose_sectors_23.json` hold the ground-truth matrix elements h_st:
for each plane, a set of control sectors and, within each sector, the
edges s -> t with their matrix element and the block-encoding scale uhat.
These come from the HDF5 matrix-element pipeline, through a code path
independent of the tables in `tables/` and of the angle synthesis that
produced the circuits.

`zero_checks_b6.json` holds, for every edge, an additional readout target
that should come back at zero: a near-miss basis state the circuit must
not connect to. These targets apply to the three forward circuits.

## Verification

`python verify.py` parses each circuit's OpenQASM 3 text and replays it as
a sparse statevector simulation, restricted to the basis states one
fixture edge touches. It checks:

- the amplitude <0,t|U_p|s,0> at fixture edges against h_st/uhat, to
  1e-9, on all six circuits, every edge under `--full` (the adjoint
  circuits are checked the same way, with the prepared and read-out
  states exchanged);
- the near-miss zero amplitudes in `zero_checks_b6.json`, below 1e-9;
- a classical check that the `add_phi` subroutine body implements
  addition mod 2^38, the ripple-carry adder underlying the phase-gradient
  rotation.

The one input the verifier takes from the paper without checking it
directly is the phase-gradient kickback identity: that a controlled
addition into the phase-gradient register realizes the intended `rz`
rotation. Every other check runs against the fixture data.

By default, `python verify.py` runs `--sample 8`: the 8 cheapest edges per
circuit. This takes roughly 90 s for a cheap edge and roughly 7 minutes for
a typical edge, single core, both measured rates. `--full --workers N`
checks every edge, with checkpointing so a run can resume after an
interruption. At the measured rates the full census is of order 100
core-hours across all six circuits; the heaviest control sectors dominate,
so the figure is rough. A full run of all six circuits took 2 h 02 m wall
clock on a 16-core desktop at 30 workers, 61 core-hours total.

## License

MIT; see `LICENSE`.

## Citation

If you use these circuits or the verifier, cite the paper:

```bibtex
@article{draper2026blockencoding,
  author        = {Draper, Patrick},
  title         = {Block Encoding Non-Abelian Lattice Gauge Theory},
  eprint        = {2608.17115},
  archivePrefix = {arXiv},
  year          = {2026}
}
```
