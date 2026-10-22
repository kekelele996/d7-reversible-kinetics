"""Fit reversible first-order kinetics from an online absorbance trace."""

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
    response_base = case["blank_absorbance"] + path * concentration * epsilon_b
    response_span = path * concentration * (epsilon_a - epsilon_b)
    if abs(response_span) < 1.0e-12:
        raise ValueError("the two components have indistinguishable response")
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


def _linear_fit(rows, base, span, initial, relaxation, dead_time, deactivation):
    exp_terms = []
    for clock, _ in rows:
        elapsed = clock - dead_time
        if elapsed < -1.0e-10:
            raise ValueError("sample precedes the fitted reaction start")
        elapsed = max(0.0, elapsed)
        if deactivation > 1.0e-12:
            decay_integral = (relaxation / deactivation) * (
                1.0 - math.exp(-deactivation * elapsed)
            )
        else:
            decay_integral = relaxation * elapsed
        exp_terms.append(math.exp(-decay_integral))

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


def _grid_search(
    rows, initial, base, span, relaxation_values, dead_values, deactivation
):
    """Return the lowest-residual point on a two-dimensional search grid."""
    best = (float("inf"), 0.0, 0.0)
    for relaxation in relaxation_values:
        for dead_time in dead_values:
            rms, _ = _linear_fit(
                rows, base, span, initial, relaxation, dead_time, deactivation
            )
            if rms < best[0]:
                best = (rms, relaxation, dead_time)
    return best


def _lambda_window(center, span):
    low = max(LAMBDA_MIN, center * math.exp(-span))
    high = min(LAMBDA_MAX, center * math.exp(span))
    log_low = math.log(low)
    return [
        math.exp(log_low + (math.log(high) - log_low) * index / 30.0)
        for index in range(31)
    ]


def _dead_window(center, span, upper):
    low = max(0.0, center - span)
    high = min(upper, center + span)
    return [low + (high - low) * index / 30.0 for index in range(31)]


def _deactivation_window(center, span):
    low = max(1.0e-6, center / (1.0 + span))
    high = center * (1.0 + span)
    return [
        math.exp(
            math.log(low)
            + (math.log(high) - math.log(low)) * index / 12.0
        )
        for index in range(13)
    ]


def _nelder_mead(case, rows, base, span, allow_decay, start):
    def objective(vector):
        relaxation = math.exp(vector[0])
        dead_time = vector[1]
        deactivation = math.exp(vector[2]) if allow_decay else 0.0
        if not (LAMBDA_MIN <= relaxation <= LAMBDA_MAX):
            return 1.0e6
        upper = max(0.0, min(item[0] for item in rows) - 0.02)
        if dead_time < 0.0 or dead_time > upper:
            return 1.0e6
        rms, _ = _linear_fit(
            rows,
            base,
            span,
            float(case["feed_fraction_a"]),
            relaxation,
            dead_time,
            deactivation,
        )
        return rms * rms

    dimensions = 3 if allow_decay else 2
    initial_decay = max(start[3], 0.001) if allow_decay else 0.0
    initial = [math.log(start[1]), start[2]]
    if allow_decay:
        initial.append(math.log(initial_decay))
    steps = [0.18, max(0.04, 0.08 * max(start[2], 0.1))]
    if allow_decay:
        steps.append(0.55)
    simplex = [list(initial)]
    for index in range(dimensions):
        point = list(initial)
        point[index] += steps[index]
        simplex.append(point)
    values = [objective(point) for point in simplex]

    for _ in range(320):
        order = sorted(range(dimensions + 1), key=lambda index: values[index])
        simplex = [simplex[index] for index in order]
        values = [values[index] for index in order]
        if max(
            abs(simplex[index][axis] - simplex[0][axis])
            for index in range(1, dimensions + 1)
            for axis in range(dimensions)
        ) < 2.0e-9:
            break
        centroid = [
            math.fsum(point[axis] for point in simplex[:-1]) / dimensions
            for axis in range(dimensions)
        ]
        worst = simplex[-1]
        reflected = [centroid[i] + (centroid[i] - worst[i]) for i in range(dimensions)]
        reflected_value = objective(reflected)
        if reflected_value < values[0]:
            expanded = [
                centroid[i] + 2.0 * (reflected[i] - centroid[i])
                for i in range(dimensions)
            ]
            expanded_value = objective(expanded)
            if expanded_value < reflected_value:
                simplex[-1], values[-1] = expanded, expanded_value
            else:
                simplex[-1], values[-1] = reflected, reflected_value
            continue
        if reflected_value < values[-2]:
            simplex[-1], values[-1] = reflected, reflected_value
            continue
        if reflected_value < values[-1]:
            contracted = [
                centroid[i] + 0.5 * (reflected[i] - centroid[i])
                for i in range(dimensions)
            ]
        else:
            contracted = [
                centroid[i] + 0.5 * (worst[i] - centroid[i])
                for i in range(dimensions)
            ]
        contracted_value = objective(contracted)
        if contracted_value < min(values[-1], reflected_value):
            simplex[-1], values[-1] = contracted, contracted_value
        else:
            best = simplex[0]
            for index in range(1, dimensions + 1):
                simplex[index] = [
                    best[axis] + 0.5 * (simplex[index][axis] - best[axis])
                    for axis in range(dimensions)
                ]
                values[index] = objective(simplex[index])

    best_index = min(range(dimensions + 1), key=lambda index: values[index])
    vector = simplex[best_index]
    relaxation = math.exp(vector[0])
    dead_time = vector[1]
    deactivation = math.exp(vector[2]) if allow_decay else 0.0
    rms, equilibrium = _linear_fit(
        rows,
        base,
        span,
        float(case["feed_fraction_a"]),
        relaxation,
        dead_time,
        deactivation,
    )
    return rms, equilibrium, relaxation, dead_time, deactivation


def _fit_mode(case, rows, base, span, allow_decay):
    initial = float(case["feed_fraction_a"])
    dead_upper = max(0.0, min(item[0] for item in rows) - 0.02)

    coarse_lambda = [
        math.exp(
            math.log(LAMBDA_MIN)
            + (math.log(LAMBDA_MAX) - math.log(LAMBDA_MIN)) * index / 95.0
        )
        for index in range(73)
    ]
    coarse_dead = [dead_upper * index / 48.0 for index in range(49)]
    deactivation_values = (
        [0.0005, 0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.12]
        if allow_decay
        else [0.0]
    )
    best = (float("inf"), 0.0, 0.0, 0.0)
    for deactivation in deactivation_values:
        candidate = _grid_search(
            rows,
            initial,
            base,
            span,
            coarse_lambda,
            coarse_dead,
            deactivation,
        )
        candidate = (*candidate, deactivation)
        if candidate[0] < best[0]:
            best = candidate

    return _nelder_mead(case, rows, base, span, allow_decay, best)


def _fit(case):
    rows = _read_samples(case)
    base, span = _response_scale(case)
    stable = _fit_mode(case, rows, base, span, allow_decay=False)
    decay = _fit_mode(case, rows, base, span, allow_decay=True)
    sample_count = len(rows)
    stable_rss = max(stable[0] ** 2 * sample_count, 1.0e-30)
    decay_rss = max(decay[0] ** 2 * sample_count, 1.0e-30)
    stable_bic = sample_count * math.log(stable_rss / sample_count) + 3.0 * math.log(sample_count)
    decay_bic = sample_count * math.log(decay_rss / sample_count) + 4.0 * math.log(sample_count)
    if decay[4] > 1.0e-6 and decay_bic < stable_bic:
        return decay
    return stable


def _rates(equilibrium, relaxation):
    reverse_rate = relaxation * equilibrium
    forward_rate = relaxation - reverse_rate
    return forward_rate, reverse_rate


def estimate(case):
    fit = _fit(case)
    _, equilibrium, relaxation, dead_time, deactivation = fit
    forward_rate, reverse_rate = _rates(equilibrium, relaxation)
    equilibrium_constant = (1.0 - equilibrium) / equilibrium
    return {
        "equilibrium_fraction_a": round(equilibrium, 6),
        "equilibrium_constant": round(equilibrium_constant, 6),
        "forward_rate_per_s": round(forward_rate, 8),
        "reverse_rate_per_s": round(reverse_rate, 8),
        "dead_time_s": round(dead_time, 3),
        "activity_model": "exponential_decay" if deactivation > 1.0e-6 else "stable",
        "catalyst_deactivation_per_s": round(deactivation, 8),
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
