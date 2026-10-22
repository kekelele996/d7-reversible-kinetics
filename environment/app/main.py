"""Estimate reversible first-order parameters from an absorbance trace."""

from __future__ import annotations

import json
import math
import sys


LAMBDA_MIN = 1.0e-5
LAMBDA_MAX = 2.5e-1
FRACTION_EPS = 1.0e-5


def _response_scale(case):
    path = float(case["path_length_cm"])
    concentration = float(case["total_concentration_mol_l"])
    epsilon_a = float(case["molar_absorptivity_a_l_mol_cm"])
    epsilon_b = float(case["molar_absorptivity_b_l_mol_cm"])
    response_base = case["blank_absorbance"]
    response_span = path * concentration * epsilon_a
    if abs(response_span) < 1.0e-12:
        raise ValueError("the detector response cannot be normalized")
    return float(response_base), float(response_span)


def _read_samples(case):
    rows = []
    for item in case["samples"]:
        rows.append(
            (
                float(item["clock_s"]),
                float(item["absorbance"]),
            )
        )
    rows.sort(key=lambda item: item[0])
    if len(rows) < 6:
        raise ValueError("at least six samples are required")
    if rows[0][0] <= 0.0:
        raise ValueError("sample clocks must be positive")
    return rows


def _linear_fit(rows, base, span, initial, relaxation, dead_time):
    exp_terms = []
    for clock, _ in rows:
        elapsed = clock - dead_time
        if elapsed < -1.0e-10:
            raise ValueError("sample precedes the fitted reaction start")
        exp_terms.append(math.exp(-relaxation * max(0.0, elapsed)))

    numerator = 0.0
    denominator = 0.0
    for (_, absorbance), exponential in zip(rows, exp_terms):
        residual_target = absorbance - base - span * initial * exponential
        weight = 1.0 - exponential
        numerator += weight * residual_target
        denominator += span * weight * weight

    equilibrium = initial
    if abs(denominator) > 1.0e-15:
        equilibrium = numerator / denominator
    equilibrium = min(1.0 - FRACTION_EPS, max(FRACTION_EPS, equilibrium))

    squared = 0.0
    for (_, absorbance), exponential in zip(rows, exp_terms):
        fraction = equilibrium + (initial - equilibrium) * exponential
        prediction = base + span * fraction
        squared += (absorbance - prediction) ** 2
    rms = math.sqrt(squared / len(rows))
    return rms, equilibrium


def _search_grid(rows, case, base, span, relaxation_values, dead_values):
    best = (float("inf"), 0.0, 0.0)
    initial = float(case["feed_fraction_a"])
    for relaxation in relaxation_values:
        for dead_time in dead_values:
            rms, equilibrium = _linear_fit(
                rows, base, span, initial, relaxation, dead_time
            )
            if rms < best[0]:
                best = (rms, relaxation, dead_time)
    return best


def _fit(case):
    rows = _read_samples(case)
    base, span = _response_scale(case)
    dead_upper = 0.0

    coarse_lambda = [
        math.exp(
            math.log(LAMBDA_MIN)
            + (math.log(LAMBDA_MAX) - math.log(LAMBDA_MIN)) * index / 95.0
        )
        for index in range(96)
    ]
    coarse_dead = [dead_upper * index / 72.0 for index in range(73)]
    best = _search_grid(rows, case, base, span, coarse_lambda, coarse_dead)

    lambda_span = 0.42
    dead_span = 0.16
    for iteration in range(24):
        _, relaxation, dead_time = best
        lambda_lo = max(LAMBDA_MIN, relaxation * math.exp(-lambda_span))
        lambda_hi = min(LAMBDA_MAX, relaxation * math.exp(lambda_span))
        lambda_values = [
            math.exp(
                math.log(lambda_lo)
                + (math.log(lambda_hi) - math.log(lambda_lo)) * index / 30.0
            )
            for index in range(31)
        ]
        lambda_candidate = _search_grid(
            rows, case, base, span, lambda_values, [dead_time]
        )
        if lambda_candidate[0] <= best[0]:
            best = lambda_candidate

        _, relaxation, dead_time = best
        dead_lo = max(0.0, dead_time - dead_span)
        dead_hi = min(dead_upper, dead_time + dead_span)
        dead_values = [
            dead_lo + (dead_hi - dead_lo) * index / 30.0 for index in range(31)
        ]
        dead_candidate = _search_grid(
            rows, case, base, span, [relaxation], dead_values
        )
        if dead_candidate[0] <= best[0]:
            best = dead_candidate

        if iteration in (3, 7, 12, 17):
            lambda_span *= 0.5
            dead_span *= 0.5

    rms, equilibrium = _linear_fit(
        rows, base, span, float(case["feed_fraction_a"]), best[1], best[2]
    )
    return rms, equilibrium, best[1], best[2]


def _rates(equilibrium, relaxation):
    reverse_rate = relaxation - relaxation * equilibrium
    forward_rate = relaxation - reverse_rate
    return forward_rate, reverse_rate


def estimate(case):
    _, equilibrium, relaxation, dead_time = _fit(case)
    forward_rate, reverse_rate = _rates(equilibrium, relaxation)
    equilibrium_constant = (1.0 - equilibrium) / equilibrium
    return {
        "equilibrium_fraction_a": round(equilibrium, 6),
        "equilibrium_constant": round(equilibrium_constant, 6),
        "forward_rate_per_s": round(forward_rate, 8),
        "reverse_rate_per_s": round(reverse_rate, 8),
        "dead_time_s": round(dead_time, 3),
    }


def main(argv):
    if len(argv) != 2:
        print("usage: python3 /app/main.py <batch.json>", file=sys.stderr)
        return 2
    with open(argv[1], "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    print(json.dumps(estimate(payload), ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
