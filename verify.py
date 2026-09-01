"""Standalone verifier for the exported OpenQASM 3 U_p streams.

Stdlib-only, runnable without ymcirc or any other dependency installed.
It does not use a general OpenQASM 3 parser; it accepts exactly the
narrow statement subset the exporter emits and hard-errors on anything
else, since a verifier that silently tolerates unrecognized syntax could
silently skip miscompiled statements.
"""
import argparse
import cmath
import concurrent.futures
import gzip
import json
import math
import os
import random
import re
import sys


class VerifyError(Exception):
    """Raised when the input doesn't match the accepted QASM subset, or
    when a replay or key check fails."""


B = 38  # phi/operand register width; fixed by the exporter, not derived.
PRUNE_TOL = 1e-16  # sparse-branch amplitude prune floor for replay() below.

# One shared token pattern for a gate argument: `q[123]`, `car[7]`, or a
# bare def-body parameter name like `anc`. Which shape is legal in a given
# statement is enforced by the surrounding regex, not by this token alone.
_ARG = r"q\[\d+\]|car\[\d+\]|[A-Za-z_]\w*"

_RE_COMMENT = re.compile(r"^//")
_RE_OPENQASM = re.compile(r"^OPENQASM 3\.0;$")
_RE_INCLUDE = re.compile(r'^include "stdgates\.inc";$')
_RE_QUBIT_DECL = re.compile(r"^qubit\[\d+\] (?:q|car);$")
_RE_BIT_DECL = re.compile(r"^bit\[\d+\] disc;$")

_DEF_NAMES = ("temp_and", "unand", "add_phi", "gradient_prep")
_RE_DEF_OPEN = re.compile(r"^def (" + "|".join(_DEF_NAMES) + r")\(.*\)\s*\{$")
_RE_BRACE_CLOSE = re.compile(r"^\}$")

# def-body-only statement forms (control/local-variable statements; these
# never produce an op tuple, they just have to parse without error).
_RE_BIT_LOCAL = re.compile(r"^bit ([A-Za-z_]\w*);$")
_RE_MEASURE_DEF = re.compile(r"^([A-Za-z_]\w*) = measure ([A-Za-z_]\w*);$")
_RE_IF_OPEN = re.compile(r"^if \([A-Za-z_]\w* == 1\) \{$")
_RE_CZ = re.compile(rf"^cz ({_ARG}), ({_ARG});$")

# gate statements: shared by def bodies (bare-name args) and the main body
# (q[N]/car[N] args). The regex accepts either shape; the caller resolves
# each captured token according to which context it is in.
_RE_GATE1 = re.compile(rf"^(x|h|s|sdg) ({_ARG});$")
_RE_GATE2 = re.compile(rf"^(cx|swap) ({_ARG}), ({_ARG});$")
_RE_GATE3 = re.compile(rf"^(ccx) ({_ARG}), ({_ARG}), ({_ARG});$")
_RE_RZ = re.compile(rf"^rz\([^)]*\) ({_ARG});$")
_RE_CALL3 = re.compile(rf"^(temp_and|unand)\(({_ARG}), ({_ARG}), ({_ARG})\);$")

# main-body-only statement forms
_RE_ADD_PHI_CALL = re.compile(r"^add_phi\((.*)\);$")
_RE_GRADIENT_PREP_CALL = re.compile(r"^gradient_prep\((.*)\);$")
_RE_MEASURE_MAIN = re.compile(r"^disc\[\d+\] = measure (q\[\d+\]);$")
_RE_BARRIER = re.compile(r"^barrier q;$")

_RE_Q_ARG = re.compile(r"^q\[(\d+)\]$")
_RE_CAR_ARG = re.compile(r"^car\[(\d+)\]$")


def _q_index(tok, line):
    m = _RE_Q_ARG.match(tok)
    if not m:
        raise VerifyError(f"expected a q[N] argument, got {tok!r}: {line!r}")
    return int(m.group(1))


def _open_lines(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for raw in f:
            line = raw.strip()
            if not line or _RE_COMMENT.match(line):
                continue
            yield line


def _parse_def_body_statement(line):
    """Parse one line from inside a def body. Returns a (gate, args) tuple
    of string argument names for a gate/call statement, or None for a
    control/local-variable statement (bit decl, measure, if-open). Raises
    VerifyError on anything else. Does not consume the closing brace; the
    caller tracks brace depth itself, since `}` closes both `if` blocks
    and the def itself."""
    if _RE_BIT_LOCAL.match(line) or _RE_MEASURE_DEF.match(line):
        return None
    m = _RE_CZ.match(line)
    if m:
        return ("cz", (m.group(1), m.group(2)))
    m = _RE_GATE1.match(line)
    if m:
        return (m.group(1), (m.group(2),))
    if _RE_RZ.match(line):
        return None  # phase gate; not a boolean/adder op, irrelevant to add_phi_body's check
    m = _RE_GATE2.match(line)
    if m:
        return (m.group(1), (m.group(2), m.group(3)))
    m = _RE_GATE3.match(line)
    if m:
        return (m.group(1), (m.group(2), m.group(3), m.group(4)))
    m = _RE_CALL3.match(line)
    if m:
        return (m.group(1), (m.group(2), m.group(3), m.group(4)))
    raise VerifyError(f"unrecognized statement inside def body: {line!r}")


def _parse_def_block(line_iter, def_name):
    """Consume a def body up to (and including) its closing brace (the
    `def NAME(...) {` line has already been matched). Returns the list
    of (gate, args) tuples for gate/call statements in the body, in
    order."""
    depth = 1
    body = []
    for line in line_iter:
        if _RE_DEF_OPEN.match(line):
            raise VerifyError(f"nested def is not part of the accepted subset: {line!r}")
        if _RE_IF_OPEN.match(line):
            depth += 1
            continue
        if _RE_BRACE_CLOSE.match(line):
            depth -= 1
            if depth == 0:
                break
            continue
        stmt = _parse_def_body_statement(line)
        if stmt is not None:
            body.append(stmt)
    else:
        raise VerifyError(f"def {def_name!r} body never closed (ran out of lines)")
    return body


def _parse_add_phi_call(line):
    m = _RE_ADD_PHI_CALL.match(line)
    if not m:
        raise VerifyError(f"malformed add_phi call: {line!r}")
    toks = [t.strip() for t in m.group(1).split(",")]
    if len(toks) != 2 * B + (B - 1):
        raise VerifyError(
            f"add_phi call has {len(toks)} args, expected {2 * B + (B - 1)}: {line!r}"
        )
    operand = tuple(_q_index(t, line) for t in toks[:B])
    phi = tuple(_q_index(t, line) for t in toks[B:2 * B])
    for t in toks[2 * B:]:
        if not _RE_CAR_ARG.match(t):
            raise VerifyError(f"add_phi carry argument must be car[...], got {t!r}: {line!r}")
    return ("add_phi", (operand, phi))


def _parse_gradient_prep_call(line):
    m = _RE_GRADIENT_PREP_CALL.match(line)
    if not m:
        raise VerifyError(f"malformed gradient_prep call: {line!r}")
    toks = [t.strip() for t in m.group(1).split(",")]
    if len(toks) != B:
        raise VerifyError(f"gradient_prep call has {len(toks)} args, expected {B}: {line!r}")
    return ("gradient_prep", tuple(_q_index(t, line) for t in toks))


def _parse_main_body_statement(line):
    """Parse one main-body line into an op tuple. Raises VerifyError on
    anything outside the accepted statement inventory."""
    m = _RE_GATE1.match(line)
    if m:
        return (m.group(1), (_q_index(m.group(2), line),))
    m = _RE_GATE2.match(line)
    if m:
        return (m.group(1), (_q_index(m.group(2), line), _q_index(m.group(3), line)))
    m = _RE_GATE3.match(line)
    if m:
        return (m.group(1), (
            _q_index(m.group(2), line), _q_index(m.group(3), line), _q_index(m.group(4), line),
        ))
    m = _RE_CALL3.match(line)
    if m:
        return (m.group(1), (
            _q_index(m.group(2), line), _q_index(m.group(3), line), _q_index(m.group(4), line),
        ))
    if line.startswith("add_phi("):
        return _parse_add_phi_call(line)
    if line.startswith("gradient_prep("):
        return _parse_gradient_prep_call(line)
    m = _RE_MEASURE_MAIN.match(line)
    if m:
        return ("measure", (_q_index(m.group(1), line),))
    if _RE_BARRIER.match(line):
        return ("barrier", ())
    raise VerifyError(f"unrecognized main-body statement: {line!r}")


def parse_qasm(path):
    """Parse a `.qasm.gz` (or plain `.qasm`) file.

    Returns a dict with:
      - "ops": list of op tuples, one per accepted main-body statement
        (comments dropped, barrier kept, so ops[i] is the i-th
        non-comment statement).
      - "phi": None (the replay engine sets this at gradient_prep).
      - "add_phi_body": list of (gate_name, (argname, ...)) tuples from
        add_phi's def body, for the classical adder check.
      - "b": 38.

    Raises VerifyError on any statement outside the accepted subset.
    """
    lines = _open_lines(path)
    add_phi_body = None

    for line in lines:
        if _RE_OPENQASM.match(line):
            continue
        if _RE_INCLUDE.match(line):
            continue
        if _RE_QUBIT_DECL.match(line):
            continue
        if _RE_BIT_DECL.match(line):
            continue
        m = _RE_DEF_OPEN.match(line)
        if m:
            def_name = m.group(1)
            body = _parse_def_block(lines, def_name)
            if def_name == "temp_and":
                if body != [("ccx", ("a", "b", "anc"))]:
                    raise VerifyError(f"temp_and def body must be exactly one ccx, got {body!r}")
            elif def_name == "add_phi":
                add_phi_body = body
            continue
        # first non-declaration, non-def line: main body begins.
        break
    else:
        line = None

    ops = []
    if line is not None:
        ops.append(_parse_main_body_statement(line))
    for line in lines:
        ops.append(_parse_main_body_statement(line))

    if add_phi_body is None:
        raise VerifyError("file has no add_phi def block")

    return {"ops": ops, "phi": None, "add_phi_body": add_phi_body, "b": B}


# --- sparse replay engine ------------------------------------------------
#
# State is `dict[(key, mask): amp]`: `key` is the wire-i-is-bit-i
# computational-basis key, `mask` is a SEPARATE wire-indexed integer
# carrying the phi-negation toggle (meaningful only on bits that are
# currently known phi wires). `phi` is a frozenset of wire indices, built
# up as `gradient_prep` calls are replayed. It is the accumulated set
# add_phi's mask check reads, distinct from the (operand, phi) tuple
# recorded in an individual add_phi op; the two agree in a
# correctly-built stream, and the accumulated set is authoritative here.


def _check_no_phi(wires, phi, name):
    bad = [w for w in wires if w in phi]
    if bad:
        raise VerifyError(
            "literal op %r touches phi wire(s) %r directly -- only "
            "cx(control, phi_target) may touch phi" % (name, bad))


def _apply_x(aug, w):
    bit = 1 << w
    out = {}
    for (k, m), a in aug.items():
        key = (k ^ bit, m)
        out[key] = out.get(key, 0j) + a
    return out


def _apply_h(aug, w):
    bit = 1 << w
    r = 1 / math.sqrt(2)
    out = {}
    for (k, m), a in aug.items():
        bitval = (k >> w) & 1
        key0 = (k & ~bit, m)
        key1 = (k | bit, m)
        out[key0] = out.get(key0, 0j) + a * r
        out[key1] = out.get(key1, 0j) + (-a * r if bitval else a * r)
    return out


def _apply_phase(aug, w, factor):
    out = {}
    for (k, m), a in aug.items():
        na = a * factor if (k >> w) & 1 else a
        key = (k, m)
        out[key] = out.get(key, 0j) + na
    return out


def _apply_cx(aug, c, t, phi):
    if c in phi:
        raise VerifyError(
            "cx(%d, %d): phi wire %d used as a CONTROL -- only "
            "cx(control, phi_target) is a sanctioned phi touch" % (c, t, c))
    out = {}
    if t in phi:
        for (k, m), a in aug.items():
            bitc = (k >> c) & 1
            key = (k, m ^ (bitc << t))
            out[key] = out.get(key, 0j) + a
    else:
        for (k, m), a in aug.items():
            bitc = (k >> c) & 1
            key = (k ^ (bitc << t), m)
            out[key] = out.get(key, 0j) + a
    return out


def _apply_swap(aug, a_wire, b_wire):
    flip = (1 << a_wire) | (1 << b_wire)
    out = {}
    for (k, m), amp in aug.items():
        bita = (k >> a_wire) & 1
        bitb = (k >> b_wire) & 1
        key = (k ^ flip if bita != bitb else k, m)
        out[key] = out.get(key, 0j) + amp
    return out


def _apply_ccx(aug, c1, c2, t):
    bit = 1 << t
    out = {}
    for (k, m), amp in aug.items():
        nk = k ^ bit if ((k >> c1) & 1) and ((k >> c2) & 1) else k
        key = (nk, m)
        out[key] = out.get(key, 0j) + amp
    return out


def _apply_temp_and(aug, u, v, anc):
    for (k, _m) in aug:
        if (k >> anc) & 1:
            raise VerifyError(
                "temp_and(%d,%d,%d): ancilla wire %d is dirty (bit=1) "
                "before the AND -- dirty-ancilla tripwire" % (u, v, anc, anc))
    return _apply_ccx(aug, u, v, anc)


def _apply_unand(aug, a, b_wire, anc, op_index):
    """Honest-channel unand: H on anc, split the state on the anc bit,
    apply cz(a, b_wire) + clear anc on the anc=1 branch, and require the
    two outcome components to agree elementwise. This mirrors the actual
    def-body gadget (H; measure; classically-controlled cz+reset on the
    1-outcome): the two components agree exactly when anc held a AND b
    on every branch, which is the case the circuit is honestly claiming."""
    aug = _apply_h(aug, anc)
    bit = 1 << anc
    psi0, psi1 = {}, {}
    for (k, m), amp in aug.items():
        if (k >> anc) & 1:
            psi1[(k, m)] = psi1.get((k, m), 0j) + amp
        else:
            psi0[(k, m)] = psi0.get((k, m), 0j) + amp
    merged1 = {}
    for (k, m), amp in psi1.items():
        namp = -amp if ((k >> a) & 1) and ((k >> b_wire) & 1) else amp
        key = (k & ~bit, m)
        merged1[key] = merged1.get(key, 0j) + namp
    for key in set(psi0) | set(merged1):
        if abs(psi0.get(key, 0j) - merged1.get(key, 0j)) >= 1e-12:
            raise VerifyError("unand at op %d: outcomes disagree" % op_index)
    return {key: amp * math.sqrt(2) for key, amp in psi0.items()}


def _apply_measure(aug, w):
    for (k, _m) in aug:
        if (k >> w) & 1:
            raise VerifyError(
                "measure(%d): wire is 1 on some branch -- discarding "
                "measures must be deterministic zero" % w)
    return aug


def _apply_gradient_prep(aug, wires):
    for (k, _m) in aug:
        for w in wires:
            if (k >> w) & 1:
                raise VerifyError(
                    "gradient_prep(%r): wire %d is not |0> before prep"
                    % (wires, w))
    return aug


def _apply_add_phi(aug, operand, phi, b):
    if not phi:
        raise VerifyError("add_phi called before any gradient_prep")
    mod = 1 << b
    out = {}
    for (k, m), amp in aug.items():
        a = 0
        for idx, w in enumerate(operand):
            a |= (((k >> w) & 1) << idx)
        toggled = [(m >> w) & 1 for w in phi]
        if all(t == 0 for t in toggled):
            sign = -1
        elif all(t == 1 for t in toggled):
            sign = +1
        else:
            raise VerifyError(
                "add_phi: partial phi-negation mask on some branch -- "
                "every phi wire must be uniformly toggled or untoggled "
                "at add time")
        phase = cmath.exp(sign * 2j * cmath.pi * a / mod)
        key = (k, m)
        out[key] = out.get(key, 0j) + amp * phase
    return out


def _apply(aug, phi, name, w, op, b, i):
    if name == "x":
        _check_no_phi(w, phi, name)
        return _apply_x(aug, w[0]), phi
    if name == "h":
        _check_no_phi(w, phi, name)
        return _apply_h(aug, w[0]), phi
    if name == "s":
        _check_no_phi(w, phi, name)
        return _apply_phase(aug, w[0], 1j), phi
    if name == "sdg":
        _check_no_phi(w, phi, name)
        return _apply_phase(aug, w[0], -1j), phi
    if name == "cx":
        return _apply_cx(aug, w[0], w[1], phi), phi
    if name == "swap":
        _check_no_phi(w, phi, name)
        return _apply_swap(aug, w[0], w[1]), phi
    if name == "ccx":
        _check_no_phi(w, phi, name)
        return _apply_ccx(aug, w[0], w[1], w[2]), phi
    if name == "temp_and":
        _check_no_phi(w, phi, name)
        return _apply_temp_and(aug, w[0], w[1], w[2]), phi
    if name == "unand":
        _check_no_phi(w, phi, name)
        return _apply_unand(aug, w[0], w[1], w[2], i), phi
    if name == "measure":
        _check_no_phi(w, phi, name)
        return _apply_measure(aug, w[0]), phi
    if name == "gradient_prep":
        return _apply_gradient_prep(aug, w), (phi | frozenset(w))
    if name == "add_phi":
        operand, _declared_phi = w
        return _apply_add_phi(aug, operand, phi, b), phi
    raise VerifyError("unknown op %r" % (name,))


def replay(ops, in_bits, injection, b=38):
    """Replay `ops` (the op-tuple list `parse_qasm` returns under "ops")
    up to and including a `barrier` op (stopping there without executing
    it). `in_bits` is the set of slots the caller wants injected as |1>;
    `injection` maps a parsed op index to the slots to inject immediately
    before executing that op (only slots also in `in_bits` are actually
    injected). Returns the collapsed `{key: amp}` state (phi-negation
    masks discarded). Raises VerifyError on any discipline break."""
    aug = {(0, 0): 1.0 + 0j}
    phi = frozenset()
    for i, op in enumerate(ops):
        for slot in injection.get(i, ()):
            if slot in in_bits:
                bit = 1 << slot
                if any(k & bit for (k, _m) in aug):
                    raise VerifyError(
                        "injection slot %d not clean at op %d" % (slot, i))
                aug = {(k | bit, m): a for (k, m), a in aug.items()}
        name, w = op[0], op[1]
        if name == "barrier":
            break
        aug, phi = _apply(aug, phi, name, w, op, b, i)
        if len(aug) > 1:
            aug = {km: a for km, a in aug.items() if abs(a) >= PRUNE_TOL}
    state = {}
    for (k, _m), a in aug.items():
        state[k] = state.get(k, 0j) + a
    return state


# --- fixture keys + first real-circuit edge ------------------------------
#
# Bridges the b6_compose_sectors_<plane>.json classical fixture (an edge
# s -> t of the block-encoded matrix, with h_st the un-normalized matrix
# element) to the real circuit's wire numbering, via qasm_readout_maps.json
# ("regmap": src_active[12]/src_ctrl[48]/src_mult[4] real wire indices,
# plus the injection schedule). A regmap's three wire lists partition
# exactly the 64 wires the |s> / |t> basis states live on; every other
# wire is implicitly |0> on both the input and (if the circuit is
# correct) the accept branch of the output.


def load_fixture(path):
    """Load a b6_compose_sectors_<plane>.json fixture as-is (a dict with
    "plane", "uhat", "sectors", "meta")."""
    with open(path) as f:
        return json.load(f)


def input_values(edge, sector):
    """The (lam4, ctrl4x4, mult4) input-side values for `edge` drawn from
    `sector`: lam4 = edge's active-link tuple, ctrl4x4 = the sector's
    control pattern (corner-major, then link), mult4[i] = the corner with
    v == i+1's own input multiplicity bit (its "s"[2])."""
    lam4 = edge["lam"]
    ctrl4x4 = sector["ctrl"]
    s2_by_v = {c["v"]: c["s"][2] for c in edge["corners"]}
    mult4 = [s2_by_v[i + 1] for i in range(4)]
    return lam4, ctrl4x4, mult4


# vertex v -> the two active-link slot indices (bi, bj) that corner v's
# "tgt" pair fills.
_VERTEX_LINK_IDX = {1: (0, 3), 2: (0, 1), 3: (2, 1), 4: (2, 3)}


def reconstruct_t_lam(edge):
    """t's full 4-slot active-link tuple, reassembled from the 4 corners'
    "tgt" fields via _VERTEX_LINK_IDX: each of t's 4 links is supplied
    TWICE, by two different corners' overlapping (bi, bj) pairs. The
    overlap agreement doubles as a fixture-consistency check on
    untrusted-until-verified input, so a disagreement raises VerifyError,
    naming the offending corner and slot."""
    t_lam = [None] * 4
    for c in edge["corners"]:
        bi, bj = _VERTEX_LINK_IDX[c["v"]]
        li, lj = c["tgt"]
        for slot, val in ((bi, li), (bj, lj)):
            if t_lam[slot] is None:
                t_lam[slot] = val
            elif t_lam[slot] != val:
                raise VerifyError(
                    "reconstruct_t_lam: corner %d disagrees with an "
                    "earlier corner on slot %d (%r vs %r)"
                    % (c["v"], slot, t_lam[slot], val))
    if any(x is None for x in t_lam):
        raise VerifyError("reconstruct_t_lam: %r left a slot unfilled" % (edge["corners"],))
    return t_lam


def output_values(edge, sector):
    """The (t_lam, ctrl4x4, t_mult) output-side values for `edge`: t_lam
    from reconstruct_t_lam, ctrl4x4 UNCHANGED from the input sector (never
    swapped), t_mult[i] = the corner with v == i+1's own "gp" (target
    multiplicity bit)."""
    t_lam = reconstruct_t_lam(edge)
    ctrl4x4 = sector["ctrl"]
    gp_by_v = {c["v"]: c["gp"] for c in edge["corners"]}
    t_mult = [gp_by_v[i + 1] for i in range(4)]
    return t_lam, ctrl4x4, t_mult


def build_key(values, regmap):
    """Pack (lam4, ctrl4x4, mult4), either `input_values`'s or
    `output_values`'s shape, onto `regmap`'s real wire numbers: each
    lam4[i] is 3 bits, LSB-first, onto src_active[3i:3i+3]; each
    ctrl4x4[c][l] is 3 bits, LSB-first, onto src_ctrl[12c+3l : 12c+3l+3];
    each mult4[i] is 1 bit onto src_mult[i]. Returns the resulting
    computational-basis key as a Python int (bit `w` = wire `w`).

    Raises VerifyError, naming the offending field and value, if any
    lam4/ctrl4x4 entry falls outside 0..7 or any mult4 entry outside 0..1
    -- each is packed into exactly 3 (or 1) bits below, so an out-of-range
    value would otherwise be silently truncated onto the wrong wires
    instead of failing loudly."""
    lam4, ctrl4x4, mult4 = values
    src_active = regmap["src_active"]
    src_ctrl = regmap["src_ctrl"]
    src_mult = regmap["src_mult"]
    for i, v in enumerate(lam4):
        if not (0 <= v <= 7):
            raise VerifyError("build_key: lam4[%d]=%r outside 0..7" % (i, v))
    for c, row in enumerate(ctrl4x4):
        for l, v in enumerate(row):
            if not (0 <= v <= 7):
                raise VerifyError("build_key: ctrl4x4[%d][%d]=%r outside 0..7" % (c, l, v))
    for i, v in enumerate(mult4):
        if not (0 <= v <= 1):
            raise VerifyError("build_key: mult4[%d]=%r outside 0..1" % (i, v))
    key = 0
    for i in range(4):
        v = lam4[i]
        for b in range(3):
            if (v >> b) & 1:
                key |= 1 << src_active[3 * i + b]
    for c in range(4):
        for l in range(4):
            v = ctrl4x4[c][l]
            base = 12 * c + 3 * l
            for b in range(3):
                if (v >> b) & 1:
                    key |= 1 << src_ctrl[base + b]
    for i in range(4):
        if mult4[i] & 1:
            key |= 1 << src_mult[i]
    return key


def _bits_set(n):
    """The set of bit positions (wire indices) set in nonnegative int n."""
    bits = set()
    i = 0
    while n:
        if n & 1:
            bits.add(i)
        n >>= 1
        i += 1
    return bits


def check_edges_grouped(circ, regmap, groups):
    """Replay `circ`["ops"] once per DISTINCT input basis state and read
    off every requested output amplitude from that one replay, instead of
    one replay per (input, output) pair. `groups` is a list of
    `(in_values, [(label, out_values), ...])` pairs: `in_values` and each
    `out_values` are `(lam4, ctrl4x4, mult4)` triples, the shape
    `input_values`/`output_values` return (a hand-built zero-check triple
    of the same shape works too). Labels are caller-chosen and opaque
    here; they must be unique across the WHOLE call, since the result
    pools every group's labels into one flat dict. This is the shared
    engine behind `check_edge`'s dagger swap, `check_zero`'s
    single-replay edge-amp plus zero-amp pairing, and the --full mode's
    per-edge replay loop (replay each edge's input once, read off every
    readout it needs)."""
    injection = {}
    for slot, idx in regmap["injection"]:
        injection.setdefault(idx, []).append(slot)
    out = {}
    for in_values, targets in groups:
        in_key = build_key(in_values, regmap)
        in_bits = _bits_set(in_key)
        state = replay(circ["ops"], in_bits, injection, b=B)
        for label, out_values in targets:
            out[label] = state.get(build_key(out_values, regmap), 0j)
    return out


def check_edge(circ, regmap, edge, sector, uhat, dagger):
    """Replay `circ`["ops"] (from `parse_qasm`) against a basis state
    built from `edge`/`sector`, and read off the amplitude at the paired
    basis state. Compares against the fixture's classical matrix element
    edge["h_st"]/uhat. Returns (amp, delta) with
    delta = abs(amp - edge["h_st"] / uhat).

    dagger=False (forward): prepares `input_values(edge, sector)` and
    reads off at `output_values(edge, sector)`: <t|U_p|s>.

    dagger=True: the caller passes the DAGGER circuit (e.g. `up_dag_12`)
    and its own matching `regmap` block (identical wire numbering to the
    forward block); `check_edge` just swaps which side is prep and which
    is readout. It prepares `output_values(edge, sector)` (edge's
    OUTPUT/t side) and reads off at `input_values(edge, sector)` (edge's
    INPUT/s side): <s|U_p^dagger|t>. h_st is real, so the adjoint matrix
    element equals the forward one; the comparison target is
    edge["h_st"] / uhat either way."""
    if dagger:
        in_values, out_values = output_values(edge, sector), input_values(edge, sector)
    else:
        in_values, out_values = input_values(edge, sector), output_values(edge, sector)
    amp = check_edges_grouped(circ, regmap, [(in_values, [("amp", out_values)])])["amp"]
    delta = abs(amp - edge["h_st"] / uhat)
    return amp, delta


_ZREC_H_ST_TOL = 1e-13  # tolerance for the zrec/edge h_st identity guard


def check_zero(circ, regmap, edge, sector, zrec):
    """Replays `circ` from `edge`'s input basis state (`input_values`)
    and reads off the amplitude at the zero-check target `zrec` names:
    t_lam = zrec["t_lam"], ctrl UNCHANGED from `edge`/`sector`'s own
    input ctrl (never swapped, same convention as `output_values`),
    mult = zrec["t_mult"]. `zrec` is one record naming a zero-check
    target; the caller pairs `zrec` to `edge` BY FIXTURE POSITION. This
    function enforces that pairing itself: it raises VerifyError, naming
    both values, unless `zrec["h_st"]` agrees with `edge["h_st"]` to
    `_ZREC_H_ST_TOL`, before building any key. Returns the amplitude
    alone (this target should read ~0 for an honest circuit).

    Implemented via `check_edges_grouped`: a caller that also wants the
    edge's own forward amp/delta from the SAME replay should call
    `check_edges_grouped` directly with both the edge's `output_values`
    and this zero target in one group, rather than calling `check_edge`
    and `check_zero` separately (which would replay twice)."""
    if abs(zrec["h_st"] - edge["h_st"]) >= _ZREC_H_ST_TOL:
        raise VerifyError(
            "check_zero: zrec[\"h_st\"]=%r does not match edge[\"h_st\"]=%r -- "
            "zrec is not paired to this edge" % (zrec["h_st"], edge["h_st"]))
    in_values = input_values(edge, sector)
    _, ctrl4x4, _ = in_values
    zero_values = (zrec["t_lam"], ctrl4x4, zrec["t_mult"])
    return check_edges_grouped(circ, regmap, [(in_values, [("zero", zero_values)])])["zero"]


# --- classical adder check + CLI -----------------------------------------
#
# check_adder is a purely classical bit-dict simulation of add_phi's def
# body: the boolean ripple-carry adder underlying its phase-kickback
# trick (contrast _apply_add_phi above, the quantum side of the same
# call). parse_qasm's own inventory of add_phi_body, and direct
# inspection of all six shipped circuits, shows it contains only
# x/cx/ccx/unand statements, so that is the whole gate set here.


def check_adder(add_phi_body, b=38, n=200, seed=0):
    """Simulates `add_phi_body` (parse_qasm's "add_phi_body") over a
    dict[argname -> 0/1] for `n` random (a, p) pairs plus four fixed edge
    cases (a=p=0; a=all-ones, p=1, full carry chain; a=1, p=all-ones;
    a=p=all-ones). Argument names follow the exporter's own convention:
    a0..a{b-1} the operand, p0..p{b-1} the register added into, c0..c{b-2}
    the carries. Requires, after the body: p's bits equal (a+p) % 2**b,
    a's bits unchanged, every carry bit back at 0. `unand(a, b_, anc)` is
    interpreted classically as an assertion (anc == a AND b_) followed by
    clearing anc. The honest quantum channel it names in the real
    circuit is a statement about superpositions that a classical
    basis-state trial never exercises. Raises VerifyError naming the
    failing trial and what broke; returns the trial count on success."""
    rng = random.Random(seed)
    mod = 1 << b
    trials = [(rng.getrandbits(b), rng.getrandbits(b)) for _ in range(n)]
    trials += [(0, 0), (mod - 1, 1), (1, mod - 1), (mod - 1, mod - 1)]
    for a, p in trials:
        bits = {}
        for i in range(b):
            bits["a%d" % i] = (a >> i) & 1
            bits["p%d" % i] = (p >> i) & 1
        for i in range(b - 1):
            bits["c%d" % i] = 0
        for gate, args in add_phi_body:
            if gate == "x":
                bits[args[0]] ^= 1
            elif gate == "cx":
                bits[args[1]] ^= bits[args[0]]
            elif gate == "ccx":
                bits[args[2]] ^= bits[args[0]] & bits[args[1]]
            elif gate == "unand":
                u, v, anc = args
                want = bits[u] & bits[v]
                if bits[anc] != want:
                    raise VerifyError(
                        "check_adder: trial a=%d p=%d -- unand(%s,%s,%s) "
                        "ancilla=%d, expected %d" % (a, p, u, v, anc, bits[anc], want))
                bits[anc] = 0
            else:
                raise VerifyError(
                    "check_adder: add_phi body contains unrecognized gate %r "
                    "-- outside the proven x/cx/ccx/unand inventory" % (gate,))
        p_out = sum(bits["p%d" % i] << i for i in range(b))
        a_out = sum(bits["a%d" % i] << i for i in range(b))
        if p_out != (a + p) % mod:
            raise VerifyError(
                "check_adder: trial a=%d p=%d -- got p'=%d, expected %d"
                % (a, p, p_out, (a + p) % mod))
        if a_out != a:
            raise VerifyError("check_adder: trial a=%d p=%d -- a changed to %d" % (a, p, a_out))
        dirty = [i for i in range(b - 1) if bits["c%d" % i]]
        if dirty:
            raise VerifyError(
                "check_adder: trial a=%d p=%d -- carries %r left dirty" % (a, p, dirty))
    return len(trials)


# --- CLI -------------------------------------------------------------------
#
# verify.py sits beside qasm/ (the six .qasm.gz files plus
# qasm_readout_maps.json) and data/ (the compose fixtures and
# zero_checks_b6.json), so the no-flags invocation "python verify.py"
# works out of the box.

_ALL_CIRCUITS = ("up_12", "up_13", "up_23", "up_dag_12", "up_dag_13", "up_dag_23")
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _nonzero_ctrl_count(sector):
    return sum(1 for row in sector["ctrl"] for v in row if v)


def _cost_ordered_edges(fx):
    """Every (fixture_index, sector, edge) in fixture `fx`, sorted by a
    branch-cost proxy (nonzero control entries in the sector's ctrl
    pattern; all-zero, vacuum-like sectors are cheap replays and sort
    first) with ties broken by fixture order (Python's sort is stable,
    and `flat` below is built in fixture order, so equal-cost edges keep
    it). `fixture_index` is the position in the natural sectors-then-edges
    flatten, the order zero_checks_b6.json is aligned to, so a caller can
    pair a selected edge back to its zero record regardless of the cost
    reordering."""
    flat = []
    idx = 0
    for sector in fx["sectors"]:
        for edge in sector["edges"]:
            flat.append((idx, sector, edge))
            idx += 1
    return sorted(flat, key=lambda item: _nonzero_ctrl_count(item[1]))


def _encode_label(label):
    return "%s:%d" % label


def _decode_label(s):
    kind, idx = s.split(":", 1)
    return (kind, int(idx))


def _build_groups(regmap, zc, dagger, selected):
    """Builds replay groups for `selected` (`_cost_ordered_edges`'s
    (fixture_index, sector, edge) triples) against `regmap`. Groups by
    DISTINCT in_values via build_key (injective by regmap construction,
    since the three wire lists partition the 64 basis-state wires), so
    two edges sharing an input basis state ride one replay. `zc` is the
    circuit's zero-check list (None for a dagger circuit: zero_checks_b6
    is keyed to the edge's INPUT/s side, per check_zero's own convention,
    which is not a meaningful target for a circuit whose replay starts
    from the edge's OUTPUT/t side instead) or the per-plane list loaded
    from zero_checks_b6.json. Returns (groups, rows): `groups` maps
    group_key (an int, build_key's output) to (in_values, [(label,
    out_values), ...]); `rows` is one dict per selected edge, naming its
    own edge_label/zero_label (zero_label is None when no zero record
    applies), for the caller to read results back out of the pooled
    replay output. Enforces the zrec/edge h_st identity guard itself,
    since a caller batching through check_edges_grouped directly, as
    this does, must not skip it."""
    groups = {}
    rows = []
    for fidx, sector, edge in selected:
        if dagger:
            in_values, out_values = output_values(edge, sector), input_values(edge, sector)
        else:
            in_values, out_values = input_values(edge, sector), output_values(edge, sector)
        gkey = build_key(in_values, regmap)
        edge_label = ("edge", fidx)
        groups.setdefault(gkey, (in_values, []))[1].append((edge_label, out_values))
        zero_label = None
        if zc is not None and fidx < len(zc) and zc[fidx] is not None:
            zrec = zc[fidx]
            if abs(zrec["h_st"] - edge["h_st"]) >= _ZREC_H_ST_TOL:
                raise VerifyError(
                    "zero record at fixture index %d is not paired to this edge "
                    "(zrec h_st=%r, edge h_st=%r)" % (fidx, zrec["h_st"], edge["h_st"]))
            _, ctrl4x4, _ = in_values
            zero_values = (zrec["t_lam"], ctrl4x4, zrec["t_mult"])
            zero_label = ("zero", fidx)
            groups[gkey][1].append((zero_label, zero_values))
        rows.append({"fidx": fidx, "sector": sector["name"], "edge": edge,
                      "edge_label": edge_label, "zero_label": zero_label})
    return groups, rows


_WORKER_CIRC = None
_WORKER_REGMAP = None


def _worker_init(qasm_path, regmap):
    global _WORKER_CIRC, _WORKER_REGMAP
    _WORKER_CIRC = parse_qasm(qasm_path)
    _WORKER_REGMAP = regmap


def _worker_replay_group(job):
    in_values, targets = job
    return check_edges_grouped(_WORKER_CIRC, _WORKER_REGMAP, [(in_values, targets)])


def _argv_sig(circuits, full, sample):
    """The subset of flags that determine WHICH groups a run computes. A
    checkpoint is only a valid resume point for a run with the same
    signature (order-independent: circuits sorted, sample meaningless
    under --full)."""
    return sorted([
        "circuits=" + ",".join(sorted(circuits)),
        "full=%s" % bool(full),
        "sample=%s" % (None if full else sample),
    ])


def _load_checkpoint(path, argv_sig):
    if not path or not os.path.exists(path):
        return {}
    with open(path) as f:
        payload = json.load(f)
    stored_sig = payload.get("argv_sig")
    if stored_sig != argv_sig:
        raise VerifyError(
            "checkpoint %r was written for a different run (stored argv_sig=%r, "
            "this run's argv_sig=%r) -- resuming would silently overwrite a "
            "possibly multi-hour artifact; delete or move the file, or rerun "
            "with the matching --circuits/--full/--sample flags"
            % (path, stored_sig, argv_sig))
    return payload.get("done", {})


def _save_checkpoint(path, argv_sig, done):
    """Atomic write (temp + os.replace) so a kill mid-write never leaves a
    truncated/corrupt checkpoint."""
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump({"argv_sig": argv_sig, "done": done}, f)
    os.replace(tmp, path)


def _replay_groups(circ, regmap, groups, workers, qasm_path,
                    checkpoint_path, argv_sig, circ_name, checkpoint_done):
    """Replays every group in `groups` not already covered by
    `checkpoint_done[circ_name]`, and returns {label: amp} pooled across
    every group (including ones served from the checkpoint). Runs inline
    when `workers <= 1` (the default single-process path); otherwise
    farms one group per ProcessPoolExecutor job, each worker parsing
    `circ_name`'s own QASM file + regmap once at init (macOS spawns
    fresh processes, and a ~0.4s parse per worker is immaterial next to
    a multi-minute replay). Checkpoints atomically after EVERY group,
    not just at the end, so a kill mid-run loses at most the one group
    in flight."""
    label_amps = {}
    done_for_circuit = checkpoint_done.setdefault(circ_name, {})
    pending = []
    for gkey, (in_values, targets) in groups.items():
        cached = done_for_circuit.get(str(gkey))
        if cached is not None:
            for label_str, re_, im_ in cached:
                label_amps[_decode_label(label_str)] = complex(re_, im_)
            continue
        pending.append((gkey, in_values, targets))

    def _record(gkey, result):
        for label, amp in result.items():
            label_amps[label] = amp
        if checkpoint_path:
            done_for_circuit[str(gkey)] = [
                [_encode_label(label), amp.real, amp.imag] for label, amp in result.items()
            ]
            _save_checkpoint(checkpoint_path, argv_sig, checkpoint_done)

    if workers <= 1:
        for gkey, in_values, targets in pending:
            _record(gkey, check_edges_grouped(circ, regmap, [(in_values, targets)]))
    elif pending:
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=workers, initializer=_worker_init,
                initargs=(qasm_path, regmap)) as ex:
            futures = {ex.submit(_worker_replay_group, (in_values, targets)): gkey
                       for gkey, in_values, targets in pending}
            for fut in concurrent.futures.as_completed(futures):
                _record(futures[fut], fut.result())
    return label_amps


def _build_argparser():
    p = argparse.ArgumentParser(
        description="Verify the shipped OpenQASM 3 U_p streams against classical ground truth.")
    p.add_argument("--circuits", default=",".join(_ALL_CIRCUITS),
                    help="comma-separated circuit names (default: all six)")
    p.add_argument("--sample", type=int, default=None,
                    help="check the N cheapest edges per circuit (default 8)")
    p.add_argument("--full", action="store_true", help="check every edge")
    p.add_argument("--adder-only", action="store_true",
                    help="run only the classical adder check, no edge replays")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--checkpoint", help="JSON path for resumable progress")
    p.add_argument("--workers", type=int, default=os.cpu_count())
    p.add_argument("--qasm-dir", help="default: <script_dir>/qasm")
    p.add_argument("--data-dir", help="default: <script_dir>/data")
    p.add_argument("--maps", help="default: <qasm-dir>/qasm_readout_maps.json")
    return p


def main(argv=None):
    parser = _build_argparser()
    args = parser.parse_args(argv)
    if args.full and args.sample is not None:
        parser.error("--sample and --full are mutually exclusive")
    sample = 8 if args.sample is None else args.sample

    qasm_dir = args.qasm_dir or os.path.join(_SCRIPT_DIR, "qasm")
    data_dir = args.data_dir or os.path.join(_SCRIPT_DIR, "data")
    maps_path = args.maps or os.path.join(qasm_dir, "qasm_readout_maps.json")
    circuits = [c.strip() for c in args.circuits.split(",") if c.strip()]

    maps = None
    if not args.adder_only:
        with open(maps_path) as f:
            maps = json.load(f)
    argv_sig = _argv_sig(circuits, args.full, sample)
    try:
        checkpoint_done = _load_checkpoint(args.checkpoint, argv_sig)
    except VerifyError as exc:
        print("CHECKPOINT ERROR: %s" % exc)
        return 1

    overall_ok = True
    for name in circuits:
        qasm_path = os.path.join(qasm_dir, name + ".qasm.gz")

        # A VerifyError here means this circuit's file doesn't parse (an
        # outside-the-subset statement, a malformed def body, ...), not a
        # bug in the verifier itself, so it must not abort the whole run:
        # catch it, record the failure, and move on to the next circuit.
        try:
            circ = parse_qasm(qasm_path)
        except VerifyError as exc:
            print("[%s] FAILED -- could not parse circuit: %s" % (name, exc))
            overall_ok = False
            continue

        # A VerifyError here means this circuit's adder check failed,
        # not a bug in the verifier itself, so it must not abort the
        # whole run: catch it, record the failure, and keep going, the
        # same discipline every edge check below already has (a bad
        # edge doesn't stop the loop either).
        trials = None
        adder_ok = True
        try:
            trials = check_adder(circ["add_phi_body"], b=circ["b"], seed=args.seed)
            print("[%s] adder check: %d trials OK" % (name, trials))
        except VerifyError as exc:
            adder_ok = False
            print("[%s] adder check: FAILED -- %s" % (name, exc))

        circuit_ok = adder_ok
        n_edges = 0
        max_delta = 0.0
        n_zero = n_zero_pass = 0

        if maps is not None and adder_ok:  # equivalent to "not args.adder_only"
            # Edge replays are skipped when the adder check itself
            # already failed for this circuit: comparing amplitudes
            # against a component already proven broken adds cost
            # without adding information.
            #
            # A VerifyError anywhere in this block (a malformed regmap, a
            # zrec/edge pairing mismatch, a discipline break caught mid-
            # replay, ...) means this circuit's data is broken, not the
            # verifier itself, so it must not abort the whole run: catch
            # it, record the failure, and still print the circuit's
            # summary line below with whatever partial counts were
            # collected before the error.
            try:
                dagger = name.startswith("up_dag_")
                plane = name[len("up_dag_"):] if dagger else name[len("up_"):]
                regmap = maps[name]
                fx = load_fixture(os.path.join(data_dir, "b6_compose_sectors_%s.json" % plane))
                zc = None
                if not dagger:
                    with open(os.path.join(data_dir, "zero_checks_b6.json")) as f:
                        zc = json.load(f).get(plane)

                ordered = _cost_ordered_edges(fx)
                selected = ordered if args.full else ordered[:sample]
                groups, rows = _build_groups(regmap, zc, dagger, selected)
                label_amps = _replay_groups(
                    circ, regmap, groups, args.workers, qasm_path,
                    args.checkpoint, argv_sig, name, checkpoint_done)

                for row in rows:
                    edge = row["edge"]
                    amp = label_amps[row["edge_label"]]
                    delta = abs(amp - edge["h_st"] / fx["uhat"])
                    max_delta = max(max_delta, delta)
                    zero_amp = label_amps.get(row["zero_label"])
                    zero_ok = True
                    zero_str = "-"
                    if zero_amp is not None:
                        n_zero += 1
                        zero_ok = abs(zero_amp) < 1e-9
                        if zero_ok:
                            n_zero_pass += 1
                        zero_str = "%.3e" % abs(zero_amp)
                    edge_ok = delta < 1e-9 and zero_ok
                    circuit_ok = circuit_ok and edge_ok
                    n_edges += 1
                    print("  %s sector=%s edge=%d delta=%.3e zero=%s %s"
                          % (name, row["sector"], row["fidx"], delta, zero_str,
                             "PASS" if edge_ok else "FAIL"))
            except VerifyError as exc:
                circuit_ok = False
                print("[%s] FAILED -- edge replay: %s" % (name, exc))

        if maps is not None:
            adder_field = str(trials) if trials is not None else "FAILED"
            print("[%s] summary: n_edges=%d max_delta=%.3e zero=%d/%d adder_trials=%s -- %s"
                  % (name, n_edges, max_delta, n_zero_pass, n_zero, adder_field,
                     "PASS" if circuit_ok else "FAIL"))
        overall_ok = overall_ok and circuit_ok

    print("RESULT: %s" % ("PASS" if overall_ok else "FAIL"))
    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
