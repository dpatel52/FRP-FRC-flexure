"""
Stage state machine for the fourteen zone solution.

Walks beta and, at each step, picks the zone that matches the current state of the
four subsystems, then updates that state from the strains it just computed.

    T   tension in the matrix,      1 -> 2 -> 3 -> 4 at beta = 1, beta_1, beta_2
    C   compression in the matrix,  1 -> 2 at the compressive yield strain
    R   tension bar,                1 -> 2 at bar yield
    RC  compression bar,            1 -> 2 at yield, and only inside T = 4

Ported from calculateEnvelope_new_2.m and Envelope_Final.m. One deliberate change.
The MATLAB version keeps `hasYieldedBefore4` in a persistent variable that survives
between calls, which meant a second beam solved in the same session inherited the
compression-steel history of the first. Envelope_Final clears it by hand. Here the
flag is explicit state, passed in and returned, so it cannot leak between analyses.
"""

import numpy as np

from . import zones as Z

__all__ = ["StageState", "walk_segment", "assemble_envelope"]

# stage string 'TCRRC' -> the zone that solves it
_ZONE_FOR = {
    "1111": ("zone111", "k111"),
    "2111": ("zone211", "k211"),
    "2121": ("zone212", "k212"),
    "2211": ("zone221", "k221"),
    "2221": ("zone222", "k222"),
    "3111": ("zone311", "k311"),
    "3121": ("zone312", "k312"),
    "3211": ("zone321", "k321"),
    "3221": ("zone322", "k322"),
    "4111": ("zone411", "k411"),
    "4121": ("zone412", "k412"),
    "4211": ("zone421", "k421"),
    "4221": ("zone422", "k422"),
    "4222": ("zone4222", "k4222"),
}


class StageState:
    """The four zone indices plus the compression-steel history flag."""

    def __init__(self, T=1, C=1, R=1, RC=1, has_yielded_before_4=False):
        self.T = T
        self.C = C
        self.R = R
        self.RC = RC
        self.has_yielded_before_4 = has_yielded_before_4

    def stage_string(self):
        return f"{self.T}{self.C}{self.R}{self.RC}"

    def __repr__(self):
        return (f"StageState(T={self.T}, C={self.C}, R={self.R}, RC={self.RC}, "
                f"has_yielded_before_4={self.has_yielded_before_4})")


def walk_segment(beta_all, curves, state, kappa, omega, epsilon_cr,
                 beta_1, beta_2, beta_3, alpha):
    """Walk one beta segment, returning the envelope and the updated state.

    Parameters
    ----------
    beta_all : array
        The beta values of this segment.
    curves : dict
        Maps 'k111', 'M111', 'k211', ... to the arrays returned by the zone solvers,
        each evaluated over this same beta segment.
    state : StageState
        Carried in and mutated, so consecutive segments continue the history.

    Returns
    -------
    (envelope, stage_used) where envelope is an (n, 2) array of [k, M] and
    stage_used records the stage string active at each step.
    """
    beta_all = np.asarray(beta_all, dtype=float)
    n = beta_all.size
    envelope = np.zeros((n, 2))
    stage_used = []
    has_yielded_in_4 = False          # local to this segment, as in the original

    for i in range(n):
        stage = state.stage_string()

        if stage == "4222":
            # 4222 only applies if the compression steel yielded inside T = 4 and
            # had not yielded earlier, otherwise fall back to 4221
            if has_yielded_in_4 and not state.has_yielded_before_4:
                kname, mname = "k4222", "M4222"
            else:
                kname, mname = "k422", "M422"
        elif stage in _ZONE_FOR:
            kname = _ZONE_FOR[stage][1]
            mname = "M" + kname[1:]
        else:
            raise ValueError(f"stage {stage!r} not recognised at beta = {beta_all[i]}")

        k_i = float(curves[kname][i])
        envelope[i, 0] = k_i
        envelope[i, 1] = float(curves[mname][i])
        stage_used.append(stage)

        beta_i = beta_all[i]
        econ = k_i * beta_i * epsilon_cr / (1.0 - k_i)                     # top fibre
        est = ((-alpha + k_i) * beta_i * epsilon_cr) / (k_i - 1.0)          # tension bar
        esc = abs(((k_i - 1.0 + alpha) * beta_i * epsilon_cr) / (k_i - 1.0))  # compression bar

        # tension in the matrix
        if state.T == 1 and beta_i >= 1.0:
            state.T = 2
        elif state.T == 2 and beta_i >= beta_1:
            state.T = 3
        elif state.T == 3 and beta_i >= beta_2:
            state.T = 4

        # compression in the matrix
        if state.C == 1 and econ >= omega * epsilon_cr:
            state.C = 2

        # tension bar
        if state.R == 1 and est >= kappa * epsilon_cr:
            state.R = 2

        # compression bar history
        if state.T < 4 and esc >= kappa * epsilon_cr:
            state.has_yielded_before_4 = True
        if state.T == 4 and esc >= kappa * epsilon_cr:
            has_yielded_in_4 = True

        if state.RC == 1 and esc >= kappa * epsilon_cr:
            if (has_yielded_in_4 and not state.has_yielded_before_4
                    and stage in ("4221", "4222")):
                state.RC = 2

    return envelope, stage_used


# which beta segment each zone is evaluated over, matching the driver. Segment 1
# covers the uncracked range, 2 the first cracked branch, and so on, so each zone
# only ever has to supply values inside the segment where it can be active.
_SEGMENT_OF = {
    "zone111": 0,
    "zone211": 1, "zone212": 1, "zone221": 1, "zone222": 1,
    "zone311": 2, "zone312": 2, "zone321": 2, "zone322": 2,
    "zone411": 3, "zone412": 3, "zone421": 3, "zone422": 3, "zone4222": 3,
}


def assemble_envelope(segments, params, kappa, omega, epsilon_cr,
                      beta_1, beta_2, beta_3, alpha, M_cr):
    """Evaluate the zones and walk the state machine through the beta segments.

    `segments` is the list of four beta arrays, cut at 0, 1, beta_1, beta_2 and
    beta_3. `params` is the argument tuple the zone solvers take after their beta
    argument.

    Each zone is evaluated over its own segment, exactly as the driver wires it,
    and the state machine then indexes position i of whichever zone is active. The
    segmentation is what keeps that indexing consistent, since the tension state
    changes only at a segment boundary.

    Returns the concatenated envelope, the concatenated beta, and the stage used at
    each step.
    """
    segments = [np.asarray(s, dtype=float) for s in segments]
    if len(segments) != 4:
        raise ValueError(f"expected 4 beta segments, got {len(segments)}")

    curves = {}
    for fname, seg_i in _SEGMENT_OF.items():
        kname = "k" + fname[4:]
        beta_seg = segments[seg_i]
        fn = getattr(Z, fname)
        if fname == "zone111":
            k, M = fn(beta_seg, *params)
        else:
            k, M = fn(M_cr, beta_seg, *params)
        curves[kname] = k
        curves["M" + fname[4:]] = M

    state = StageState()
    env_parts, stages = [], []
    for beta_seg in segments:
        if beta_seg.size == 0:
            continue
        env, used = walk_segment(beta_seg, curves, state, kappa, omega, epsilon_cr,
                                 beta_1, beta_2, beta_3, alpha)
        env_parts.append(env)
        stages.extend(used)

    return (np.vstack(env_parts), np.concatenate(segments), stages)
