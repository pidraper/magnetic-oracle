# SU(3) plaquette-oracle circuits

Companion to *Block Encoding Non-Abelian Lattice Gauge Theory* (arXiv:2608.17115). This package contains the compiled plaquette operator block-encoding oracle for SU(3) lattice gauge theory. The local truncation is a per-vertex electric casimir cut, 
C_2(tot) <= 6, and rotations are stored in b_rot = 38 bits. Six explicit circuits are provided in OpenQASM 3, corresponding to three planes, the oracle, and its adjoint. Transition tables used to generate numbers in the circuits' lookup tables are also included. A verification script verify.py can be used to check that the amplitudes produced by the circuit match the target values.

## Layout

```
qasm/       six OpenQASM 3 circuits + manifests
tables/     per-plane transition tables + shared constants
data/       target matrix elements
verify.py   verifies the circuits against data/
```

## qasm/

One circuit per plane and direction: `up_{12,13,23}.qasm.gz` implement the oracle U_p on the three plaquette planes, and `up_dag_{12,13,23}.qasm.gz` implement the adjoint. (The adjoint is not a mechanical gate reversal of the forward circuit because the uncompute gadgets use mid-circuit measurement with classical feedforward.)

Each file is gzip-compressed OpenQASM 3.0 text using `stdgates.inc`, about 485,000 statements. Registers are one flat `qubit[1055] q`, a `qubit[37] car` carry register for the phase-gradient adder, and a classical register for the discarding measurements. Four subroutines are used throughout the circuits:

- `temp_and`: toffoli implementation (4 T)
- `unand`: measurement-based uncompute, an X-basis measurement and a classically controlled `cz` (0 T)
- `add_phi`:  in-place ripple-carry adder used in the programmed rotations
- `gradient_prep`:  phase-gradient state preparation, run once per circuit. its 38 `rz(2*pi*2^j/2^38)` gates are the only rotation gates in the circuits.

`qasm_provenance.json` : the per-circuit manifest

`qasm_readout_maps.json` holds, for each circuit, the register wire map (which physical wires carry the active links, control links, and multiplicities of the input and output basis states) and the injection points where an input state is written before replay. `verify.py` uses this file to build and read off basis states without touching the circuit text.

## tables/

`b6_plane12.json`, `b6_plane13.json`, `b6_plane23.json` provide the per-plane transition tables (proposal rows, corner data, class labels). 
`b6_shared.json` contains the shared constants. These are the data the compiler
loaded at build time to produce the circuits.

## data/

`b6_compose_sectors_12.json`, `b6_compose_sectors_13.json`, and `b6_compose_sectors_23.json` contain the true matrix elements h_st. For each plane, transitions are specified by control sectors, and within each sector, core data transitions s -> t. These come from a pyclebsch computation that is independent from the circuits.

`zero_checks_b6.json` contains, for every edge, an additional readout target that should come back at zero: a near-miss basis state the circuit must not populate. These targets apply to the three forward circuits.

## Verification

`python verify.py` parses each circuit's OpenQASM 3 text and replays it as a sparse statevector simulation. It checks:

- the amplitude <0,t|U_p|s,0> for tabulated transitions s -> t in `data/` against h_st/uhat, on all six circuits. Use `--full` to check every transition. (The adjoint circuits are checked the same way, with the prepared and read-out states exchanged)
- the near-miss zero amplitudes in `zero_checks_b6.json`
- a classical check that the `add_phi` subroutine body implements   addition mod 2^38, the ripple-carry adder used in the programmed rotations.

A check that runs in a couple of seconds, the classical adder check alone:

```
$ python3 verify.py --circuits up_12 --adder-only
[up_12] adder check: 204 trials OK
RESULT: PASS
```

Two transitions on one circuit, a few minutes on a single core:

```
$ python3 verify.py --circuits up_12 --sample 2 --workers 1
[up_12] adder check: 204 trials OK
  up_12 sector=vacuum edge=0 delta=4.269e-12 zero=0.000e+00 PASS
  up_12 sector=vacuum edge=1 delta=6.393e-12 zero=0.000e+00 PASS
[up_12] summary: n_edges=2 max_delta=6.393e-12 zero=2/2 adder_trials=204 -- PASS
RESULT: PASS
```

In each transition line, `delta` is |replayed amplitude - h_st/uhat| and `zero` is the amplitude at the near-miss target (`-` where none applies, as for the adjoint circuits). The full census, with checkpointing so an interrupted run resumes where it left off:

```
$ python3 verify.py --full --workers 30 --checkpoint ckpt.json
...
[up_dag_23] summary: n_edges=566 max_delta=6.393e-12 zero=0/0 adder_trials=204 -- PASS
RESULT: PASS
```

By default, `python verify.py` runs `--sample 8`: the 8 cheapest edges per circuit. This took roughly 90 s for a cheap edge and roughly 7 minutes for a typical edge, single core, on an Apple M4. `--full --workers N` checks every edge, with checkpointing so a run can resume after an interruption. At the measured rates the full census is of order 100 core-hours across all six circuits. A full run of all six circuits took 2h wall clock at 30 workers on a 16-core AMD Ryzen 9 9950X3D, 61 core-hours total.

## License

MIT; see `LICENSE`.

## Citation

If you use these circuits, please cite the paper:

```bibtex
@article{draper2026blockencoding,
  author        = {Draper, Patrick},
  title         = {Block Encoding Non-Abelian Lattice Gauge Theory},
  eprint        = {2608.17115},
  archivePrefix = {arXiv},
  year          = {2026}
}
```
