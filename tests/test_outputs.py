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
    deactivation=0.0,
):
    equilibrium = kb / (kf + kb)
    samples = []
    for index, clock in enumerate(clocks):
        elapsed = clock - dead_time
        if deactivation > 0.0:
            decay_integral = ((kf + kb) / deactivation) * (
                1.0 - math.exp(-deactivation * elapsed)
            )
        else:
            decay_integral = (kf + kb) * elapsed
        fraction = equilibrium + (initial - equilibrium) * math.exp(-decay_integral)
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
    return payload, {
        "equilibrium": equilibrium,
        "kf": kf,
        "kb": kb,
        "dead": dead_time,
        "deactivation": deactivation,
    }


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
    """「程序向标准输出写一个 JSON 对象，只含以下七个字段」"""
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
        "activity_model",
        "catalyst_deactivation_per_s",
    }


def test_wrong_argument_count_exits_2():
    """「参数个数不对时退出码为 2，不向标准输出写 JSON」"""
    proc = run_candidate()
    assert proc.returncode == 2
    assert proc.stdout.strip() == ""


def test_repeated_run_returns_the_same_json():
    """「同一份批次文件重复执行，输出必须一致」"""
    payload, _ = make_case(
        initial=0.73,
        concentration=0.016,
        epsilon_a=1180.0,
        epsilon_b=320.0,
        blank=0.005,
        kf=0.0046,
        kb=0.0031,
        dead_time=0.35,
        clocks=[0.60, 0.82, 1.15, 1.63, 2.30, 3.25, 4.60, 6.40, 8.90, 12.2],
    )
    folder = Path(tempfile.mkdtemp(prefix="repeat-"))
    folder.chmod(0o755)
    path = folder / "case.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    path.chmod(0o644)
    first = run_candidate(path)
    second = run_candidate(path)
    assert first.returncode == second.returncode == 0
    assert first.stdout == second.stdout


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

def test_strong_activity_decay_is_selected_and_quantified():
    """「活性可能稳定，也可能按一阶指数规律衰减」"""
    payload, truth = make_case(
        initial=0.80,
        concentration=0.0165,
        epsilon_a=1380.0,
        epsilon_b=0.0,
        blank=0.0,
        kf=0.0060,
        kb=0.0030,
        dead_time=1.20,
        clocks=[1.48, 1.72, 2.05, 2.50, 3.10, 3.90, 5.00, 6.50, 8.60, 11.50, 15.50, 21.00, 28.00, 37.00],
        noise=0.000020,
        phase=0.4,
        deactivation=0.018,
    )
    result = call_case(payload)
    assert result["activity_model"] == "exponential_decay"
    assert relative_error(
        result["catalyst_deactivation_per_s"], truth["deactivation"]
    ) < 0.08
    assert abs(float(result["dead_time_s"]) - truth["dead"]) < 0.15


def test_stable_activity_is_not_confused_with_parameter_noise():
    """「取 BIC 较小的模型」"""
    payload, truth = make_case(
        initial=0.43,
        concentration=0.018,
        epsilon_a=1250.0,
        epsilon_b=0.0,
        blank=0.0,
        kf=0.0044,
        kb=0.0021,
        dead_time=1.60,
        clocks=[1.88, 2.12, 2.46, 2.90, 3.52, 4.38, 5.50, 7.05, 9.10, 11.85, 15.50, 20.40, 27.00],
        noise=0.000018,
        phase=1.1,
        deactivation=0.0,
    )
    result = call_case(payload)
    assert result["activity_model"] == "stable"
    assert result["catalyst_deactivation_per_s"] == 0.0
    assert abs(float(result["dead_time_s"]) - truth["dead"]) < 0.12


def test_weak_decay_is_retained_when_bic_supports_the_extra_parameter():
    """「程序需要根据整段曲线在“稳定活性”和“指数失活”之间作出选择」"""
    payload, truth = make_case(
        initial=0.91,
        concentration=0.0105,
        epsilon_a=1680.0,
        epsilon_b=0.0,
        blank=0.0,
        kf=0.0090,
        kb=0.0040,
        dead_time=0.85,
        clocks=[1.02, 1.15, 1.34, 1.62, 2.00, 2.55, 3.25, 4.20, 5.45, 7.05, 9.20, 12.10, 16.00, 21.20, 28.00, 37.00, 49.00],
        noise=0.000012,
        phase=2.0,
        deactivation=0.0045,
    )
    result = call_case(payload)
    assert result["activity_model"] == "exponential_decay"
    assert relative_error(
        result["catalyst_deactivation_per_s"], truth["deactivation"]
    ) < 0.12
    assert abs(float(result["dead_time_s"]) - truth["dead"]) < 0.15


def test_exact_stable_trace_does_not_report_decay():
    """「活性可能稳定，也可能按一阶指数规律衰减」"""
    payload, truth = make_case(
        initial=0.03,
        concentration=0.025448327989052225,
        epsilon_a=1559.3853889859956,
        epsilon_b=0.0,
        blank=0.01154116322640787,
        kf=0.034539705511752575,
        kb=0.009708832597949815,
        dead_time=1.2601320078482245,
        clocks=[
            1.3401320078482246,
            3.798998644218323,
            7.198028607151327,
            8.514031129490501,
            8.934921766260189,
            10.72800605235682,
            11.759563942691573,
            12.446254138463106,
            13.631094662592233,
            14.393012298700214,
            14.541260388402733,
            15.612403120855346,
            21.37670627721609,
        ],
        deactivation=0.0,
    )
    result = call_case(payload)
    assert result["activity_model"] == "stable"
    assert result["catalyst_deactivation_per_s"] == 0.0
    assert abs(float(result["dead_time_s"]) - truth["dead"]) < 0.12


def test_weak_decay_near_pure_feed_keeps_equilibrium_identifiable():
    """「程序需要根据整段曲线在“稳定活性”和“指数失活”之间作出选择」"""
    payload, truth = make_case(
        initial=0.995,
        concentration=0.012238292330501251,
        epsilon_a=1248.8271190103612,
        epsilon_b=0.0,
        blank=0.03752451609579396,
        kf=0.0015969734398884605,
        kb=0.002792729637222259,
        dead_time=8.329148887343214,
        clocks=[
            8.439281705579264,
            8.52070630147967,
            9.26904917308594,
            10.047743784943096,
            10.201916282441882,
            10.976685700031314,
            11.085502649354325,
            11.334497030349207,
            12.965806866051938,
            14.073216748481933,
            14.189380924201721,
            16.626506707102067,
            19.589607740611147,
            20.588171151273478,
            20.702070285924627,
            26.57506546873707,
            29.147327687611618,
            29.18078192523167,
            30.905330691704833,
            34.34965054946141,
            34.41601561892689,
            35.17560805013328,
            35.275523211488846,
            35.478617813173344,
            36.722594448537656,
        ],
        noise=0.000020,
        phase=1.085338759863179,
        deactivation=0.0015992714349186835,
    )
    result = call_case(payload)
    assert result["activity_model"] == "exponential_decay"
    assert relative_error(result["equilibrium_fraction_a"], truth["equilibrium"]) < 0.008
    assert relative_error(
        result["catalyst_deactivation_per_s"], truth["deactivation"]
    ) < 0.12
    assert abs(float(result["dead_time_s"]) - truth["dead"]) < 0.12


def test_noisy_decay_trace_keeps_core_parameters_coupled():
    """「两条候选模型的参数个数不同」"""
    payload, truth = make_case(
        initial=0.43470210194374836,
        concentration=0.025481694252694403,
        epsilon_a=975.1470410028222,
        epsilon_b=0.0,
        blank=0.0032105504860038938,
        kf=0.0062393012042914345,
        kb=0.0015965501129257171,
        dead_time=2.399668074463084,
        clocks=[
            4.227127565934367,
            4.317169821234588,
            4.5132577798565805,
            5.11676844167763,
            5.261997955093118,
            5.606929457791726,
            7.064390539533694,
            7.265006844154274,
            7.376197976035774,
            8.402607130623265,
            9.12629635846363,
            9.298262785300414,
            9.337036491470132,
            10.655630704936968,
            10.821907136796856,
            10.823567365483997,
            11.267958283562498,
            12.079719407766106,
            12.410051158848972,
            13.00155970927245,
            14.897810524527634,
            15.013059671760049,
            15.64770836454195,
            16.092551648649092,
            16.19432696245362,
            17.073217839327754,
            17.49083172976544,
            18.760967932543743,
            19.81673775143293,
            22.616568365596294,
            23.260918831187873,
            23.408933458989573,
            23.951800671035542,
            24.739828369346455,
            26.33035958180031,
            27.20307428229671,
            27.221387157151433,
            27.238854440694105,
            27.460679988052526,
            27.580551604045716,
            28.521783949363517,
            29.090624142044938,
            29.251918701914825,
            30.032583855871408,
            31.498820223678898,
            31.752772124252374,
            32.30690870067497,
            32.69394002072747,
            33.001273074633836,
            33.21911190077168,
            33.755467905454495,
            34.142637952765625,
            34.585363801971006,
            34.88773134815577,
        ],
        noise=0.000080,
        phase=1.3490371604869287,
        deactivation=0.0050929405942706125,
    )
    result = call_case(payload)
    assert relative_error(result["equilibrium_fraction_a"], truth["equilibrium"]) < 0.008
    assert abs(float(result["dead_time_s"]) - truth["dead"]) < 0.12
