"""Behavioral checks for the reversible isomerization parameter report."""

from __future__ import annotations

import json
import math
import os
import signal
import subprocess
import tempfile
from pathlib import Path


ENTRY = "/app/main.py"
PYTHON = "/usr/local/bin/python3"
CANDIDATE_USER = "runner"
CANDIDATE_ENV = {
    "PATH": "/usr/local/bin:/usr/bin:/bin",
    "HOME": "/home/runner",
    "LANG": "C.UTF-8",
    "TZ": "UTC",
}
CALL_TIMEOUT_SEC = 45


def run_candidate(path=None):
    args = [PYTHON, ENTRY]
    if path is not None:
        args.append(str(path))
    as_root = os.getuid() == 0
    proc = subprocess.Popen(
        args,
        cwd="/app",
        env=CANDIDATE_ENV,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
        user=CANDIDATE_USER if as_root else None,
        group=CANDIDATE_USER if as_root else None,
        extra_groups=[] if as_root else None,
    )
    try:
        out, err = proc.communicate(timeout=CALL_TIMEOUT_SEC)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.communicate()
        raise TimeoutError(f"{ENTRY} exceeded {CALL_TIMEOUT_SEC}s")
    return subprocess.CompletedProcess(proc.args, proc.returncode, out, err)


def make_case(
    *,
    initial,
    concentration,
    epsilon_a,
    epsilon_b,
    blank,
    kf,
    kb,
    dead_time,
    clocks,
    noise=0.0,
    phase=0.0,
):
    equilibrium = kb / (kf + kb)
    samples = []
    for index, clock in enumerate(clocks):
        elapsed = clock - dead_time
        fraction = equilibrium + (initial - equilibrium) * math.exp(
            -(kf + kb) * elapsed
        )
        absorbance = blank + concentration * (
            epsilon_a * fraction + epsilon_b * (1.0 - fraction)
        )
        absorbance += noise * math.sin(1.37 * (index + 1) + phase)
        samples.append(
            {
                "clock_s": clock,
                "absorbance": round(absorbance, 6),
            }
        )
    payload = {
        "feed_fraction_a": initial,
        "total_concentration_mol_l": concentration,
        "molar_absorptivity_a_l_mol_cm": epsilon_a,
        "molar_absorptivity_b_l_mol_cm": epsilon_b,
        "path_length_cm": 1.0,
        "blank_absorbance": blank,
        "samples": samples,
    }
    return payload, {"equilibrium": equilibrium, "kf": kf, "kb": kb, "dead": dead_time}


def call_case(payload):
    folder = Path(tempfile.mkdtemp(prefix="batch-"))
    folder.chmod(0o755)
    path = folder / "case.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    path.chmod(0o644)
    proc = run_candidate(path)
    assert proc.returncode == 0, f"exit={proc.returncode} stderr={proc.stderr[-600:]}"
    return json.loads(proc.stdout)


def relative_error(actual, expected):
    return abs(float(actual) - float(expected)) / abs(float(expected))


# ---------- gate: invocation and output contract ----------

def test_prints_only_the_contract_fields():
    """「程序向标准输出写一个 JSON 对象，只含以下五个字段」"""
    payload, _ = make_case(
        initial=0.8,
        concentration=0.015,
        epsilon_a=1100.0,
        epsilon_b=210.0,
        blank=0.004,
        kf=0.005,
        kb=0.002,
        dead_time=0.0,
        clocks=[0.4, 0.8, 1.5, 2.6, 4.2, 6.5, 9.8, 14.0],
    )
    result = call_case(payload)
    assert set(result) == {
        "equilibrium_fraction_a",
        "equilibrium_constant",
        "forward_rate_per_s",
        "reverse_rate_per_s",
        "dead_time_s",
    }


def test_wrong_argument_count_exits_2():
    """「参数个数不对时退出码为 2，不向标准输出写 JSON」"""
    proc = run_candidate()
    assert proc.returncode == 2
    assert proc.stdout.strip() == ""


# ---------- difficulty_1: mixture response and equilibrium composition ----------

def test_equilibrium_fraction_includes_both_component_absorptivities():
    """「吸收峰同时包含两种异构体的贡献」"""
    payload, truth = make_case(
        initial=0.84,
        concentration=0.0142,
        epsilon_a=1260.0,
        epsilon_b=305.0,
        blank=0.0037,
        kf=0.0052,
        kb=0.0021,
        dead_time=0.0,
        clocks=[0.4, 0.9, 1.7, 2.8, 4.5, 6.9, 9.6, 13.1, 18.4, 24.8],
    )
    result = call_case(payload)
    assert relative_error(result["equilibrium_fraction_a"], truth["equilibrium"]) < 0.008


def test_equilibrium_fraction_handles_reversed_response_direction():
    """「平衡时 A 的摩尔分数，保留六位小数」"""
    payload, truth = make_case(
        initial=0.22,
        concentration=0.019,
        epsilon_a=280.0,
        epsilon_b=980.0,
        blank=0.0061,
        kf=0.0018,
        kb=0.0040,
        dead_time=0.0,
        clocks=[0.3, 0.7, 1.2, 2.0, 3.3, 5.2, 7.9, 11.7, 16.9, 23.5],
    )
    result = call_case(payload)
    assert relative_error(result["equilibrium_fraction_a"], truth["equilibrium"]) < 0.008


def test_equilibrium_fraction_survives_changed_blank_and_loading():
    """「纯溶剂空白由每批任务文件提供」"""
    payload, truth = make_case(
        initial=0.88,
        concentration=0.022,
        epsilon_a=760.0,
        epsilon_b=110.0,
        blank=0.0105,
        kf=0.0027,
        kb=0.00235,
        dead_time=0.0,
        clocks=[0.5, 1.1, 2.3, 4.0, 6.4, 9.3, 13.6, 19.2, 26.8, 35.0],
    )
    result = call_case(payload)
    assert relative_error(result["equilibrium_fraction_a"], truth["equilibrium"]) < 0.008


# ---------- difficulty_2: relaxation rate decomposition ----------

def test_forward_favored_case_splits_relaxation_rate():
    """「A→B 的一级速率常数，单位 s⁻¹」"""
    payload, truth = make_case(
        initial=0.88,
        concentration=0.016,
        epsilon_a=1500.0,
        epsilon_b=0.0,
        blank=0.0,
        kf=0.0070,
        kb=0.0025,
        dead_time=0.0,
        clocks=[0.2, 0.6, 1.3, 2.2, 3.5, 5.1, 7.4, 10.2, 14.1, 18.6],
    )
    result = call_case(payload)
    assert relative_error(result["forward_rate_per_s"], truth["kf"]) < 0.03
    assert relative_error(result["reverse_rate_per_s"], truth["kb"]) < 0.03


def test_reverse_favored_case_splits_relaxation_rate():
    """「B→A 的一级速率常数，单位 s⁻¹」"""
    payload, truth = make_case(
        initial=0.12,
        concentration=0.013,
        epsilon_a=1500.0,
        epsilon_b=0.0,
        blank=0.0,
        kf=0.0015,
        kb=0.00395,
        dead_time=0.0,
        clocks=[0.25, 0.8, 1.6, 2.7, 4.1, 6.0, 8.8, 12.4, 17.5, 24.0],
    )
    result = call_case(payload)
    assert relative_error(result["forward_rate_per_s"], truth["kf"]) < 0.03
    assert relative_error(result["reverse_rate_per_s"], truth["kb"]) < 0.03


def test_rate_ratio_matches_reported_equilibrium_constant():
    """「按 B 摩尔分数除以 A 摩尔分数定义」"""
    payload, truth = make_case(
        initial=0.78,
        concentration=0.017,
        epsilon_a=1500.0,
        epsilon_b=0.0,
        blank=0.0,
        kf=0.0084,
        kb=0.0060,
        dead_time=0.0,
        clocks=[0.18, 0.45, 0.9, 1.5, 2.4, 3.7, 5.5, 8.1, 11.6, 16.5],
    )
    result = call_case(payload)
    expected_ratio = truth["kf"] / truth["kb"]
    reported_ratio = result["forward_rate_per_s"] / result["reverse_rate_per_s"]
    assert relative_error(reported_ratio, expected_ratio) < 0.02
    assert relative_error(result["equilibrium_constant"], expected_ratio) < 0.02


# ---------- difficulty_3: controller clock delay ----------

def test_dead_time_is_estimated_from_early_clock_samples():
    """「控制器时钟与实际反应起点之间存在一个非负且固定的输送延迟」"""
    payload, truth = make_case(
        initial=0.76,
        concentration=0.0155,
        epsilon_a=1380.0,
        epsilon_b=0.0,
        blank=0.0,
        kf=0.0048,
        kb=0.0022,
        dead_time=1.45,
        clocks=[1.70, 1.92, 2.27, 2.74, 3.35, 4.25, 5.60, 7.35, 9.80, 13.20, 18.10],
        noise=0.000018,
        phase=0.4,
    )
    result = call_case(payload)
    assert abs(float(result["dead_time_s"]) - truth["dead"]) < 0.12


def test_dead_time_remains_visible_for_faster_trace():
    """「后段又常把平衡组成判到另一侧」"""
    payload, truth = make_case(
        initial=0.31,
        concentration=0.018,
        epsilon_a=1250.0,
        epsilon_b=0.0,
        blank=0.0,
        kf=0.0065,
        kb=0.0055,
        dead_time=2.20,
        clocks=[2.55, 2.82, 3.18, 3.70, 4.32, 5.10, 6.25, 7.80, 9.90, 12.70, 16.50],
        noise=0.000024,
        phase=1.1,
    )
    result = call_case(payload)
    assert abs(float(result["dead_time_s"]) - truth["dead"]) < 0.12


def test_dead_time_with_irregular_controller_timestamps():
    """「早期采样点与控制器时钟对不上」"""
    payload, truth = make_case(
        initial=0.90,
        concentration=0.011,
        epsilon_a=1710.0,
        epsilon_b=0.0,
        blank=0.0,
        kf=0.0033,
        kb=0.0014,
        dead_time=0.73,
        clocks=[0.95, 1.02, 1.18, 1.41, 1.73, 2.20, 2.81, 3.54, 4.47, 5.76, 7.35, 9.40],
        noise=0.000015,
        phase=2.0,
    )
    result = call_case(payload)
    assert abs(float(result["dead_time_s"]) - truth["dead"]) < 0.12
