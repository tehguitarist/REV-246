"""End-to-end test of harness stages A–C on a synthetic library built from the reference IR.

    .venv/bin/python -m pytest tests -q
"""
import os
import subprocess
import sys

import pandas as pd
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "irs", "eugtone_EMT246_CLA_shortRT.wav")
PY = sys.executable
BASE = "EMT246_S_RT1.0_LF0.9_HF0.9_D30_A10"


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    if not os.path.exists(SRC):
        pytest.skip("reference IR not present")
    lib = tmp_path_factory.mktemp("synthlib")
    subprocess.run([PY, os.path.join(ROOT, "tools/harness/synth.py"), SRC, str(lib)], check=True, capture_output=True)
    name = "pytest_synthetic"
    subprocess.run([PY, os.path.join(ROOT, "tools/harness/run.py"), str(lib), "--name", name], check=True, capture_output=True)
    out = os.path.join(ROOT, "analysis", "library", name)
    files = pd.read_csv(os.path.join(out, "files.csv"))
    files["flags"] = files["flags"].fillna("").str.split(";")
    yield files.set_index("relpath"), out


def row(files, rel):
    return files.loc[rel]


def has(files, rel, flag):
    return flag in row(files, rel)["flags"]


def test_base_file(run):
    f, _ = run
    r = row(f, BASE + ".wav")
    assert r["cluster_status"] == "ok"
    assert r["program_measured"] == "S"
    assert abs(r["rt_mid"] - 0.99) < 0.03
    assert abs(r["lr_offset_us"] - 18.6) < 0.5


@pytest.mark.parametrize("suffix,status", [("_swapped", "swapped"), ("_invR", "inverted_R")])
def test_channel_fingerprint(run, suffix, status):
    f, _ = run
    assert row(f, f"qc/{BASE}{suffix}.wav")["cluster_status"] == status


@pytest.mark.parametrize("suffix,flag", [
    ("_faded", "fade_out"), ("_preecho", "pre_onset_energy"), ("_mono", "mono"), ("_clipped", "clipped"),
    ("_trimmed", "no_pre_roll"), ("_skew150ppm", "clock_outlier"), ("_mislabelled", "label_mismatch_rt"),
    ("_swapped", "label_group_inconsistent")])
def test_qc_flags(run, suffix, flag):
    f, _ = run
    assert has(f, f"qc/{BASE}{suffix}.wav", flag)


def test_duplicates(run):
    f, _ = run
    assert has(f, "qc/EMT246_S_RT2.0_LF0.9_HF0.9_D30_A10_dup_other_label.wav", "duplicate_with_different_label")
    assert has(f, "qc/random_unparseable_name.wav", "unlabelled_duplicate_of_labelled_file")
    assert has(f, f"repeat/{BASE}.wav", "has_repeat_capture")
    assert has(f, f"48k/{BASE}.wav", "other_rate_version")
    assert not has(f, BASE + ".wav", "duplicate_with_different_label")  # delay-only label differences are fine


def test_settings_inference(run):
    f, _ = run
    assert abs(row(f, "EMT246_S_RT0.5_LF0.9_HF0.9_D30_A10.wav")["rt_mid"] - 0.5) < 0.05
    assert row(f, "EMT246_S_RT1.0_LF0.9_HF0.9_D30_A0.wav")["cluster_status"] == "no_cluster"
    assert row(f, "EMT246_S_RT1.0_LF0.9_HF0.9_D30_A0.wav")["program_measured"] == "S"
    a5 = row(f, "EMT246_S_RT1.0_LF0.9_HF0.9_D30_A5.wav")["cluster_gain_db"]
    a10 = row(f, BASE + ".wav")["cluster_gain_db"]
    assert abs((a10 - a5) - 6.02) < 0.2
    assert row(f, "EMT246_S_RT1.0_LF0.9_HF0.5_D30_A10.wav")["hf_ratio"] < 0.7


def test_clock_and_latency(run):
    f, out = run
    skew = row(f, f"qc/{BASE}_skew150ppm.wav")["clock_ppm"] - row(f, BASE + ".wav")["clock_ppm"]
    assert abs(skew + 150) < 5
    import json
    lib = json.load(open(os.path.join(out, "library.json")))
    assert lib["latency"]["status"] == "ok"
    assert abs(lib["latency"]["slope"] - 1) < 0.001
    assert abs(lib["latency"]["intercept_ms"] + 20.0) < 0.1  # synthetic labels say D30, first tap at ~10 ms
    assert os.path.exists(os.path.join(out, "report.md"))


def test_kienast_naming():
    sys.path.insert(0, os.path.join(ROOT, "tools", "harness"))
    from core import load_config
    from naming import parse_name
    cfg = load_config()["naming"]
    s = parse_name("96kHz/Long RT/12_0_80_08_100_01_dc.wav", cfg)
    assert (s["label_program"], s["label_amp"], s["label_lf"], s["label_rt"], s["label_hf"]) == ("L", 80, 0.8, 10.0, 0.1)
    n = parse_name("48kHz/Non-Lin/12_0_05_-30_20_dc.wav", cfg)
    assert (n["label_program"], n["label_span"], n["label_early"], n["label_late"]) == ("N", 0.5, -3.0, 2.0)
