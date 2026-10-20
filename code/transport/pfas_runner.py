#!/usr/bin/env python3
"""Live, resumable runner for the four PFAS transport and inverse-SSL cases."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


CASES = {
    1: {"name": "PFOS_Langmuir", "cycles": 35, "seed": 0.004,
        "initial_low": 0.00129427701300129, "initial_high": 0.00134710464618502},
    2: {"name": "PFOS_Kralchevsky", "cycles": 10, "seed": 0.014,
        "initial_low": 0.09, "initial_high": 0.095},
    3: {"name": "PFOA_Langmuir", "cycles": 5, "seed": 0.13},
    4: {"name": "PFOA_Kralchevsky", "cycles": 10, "seed": 0.40},
}

SOLVER_VERSION = "PFAS_transport_SSL_2026-08_SIparams_concentration_v5_bounded_retry"


def say(message: str) -> None:
    print(message, flush=True)


def matlab_executable() -> str:
    configured = os.environ.get("PFAS_MATLAB_CMD")
    if configured:
        return configured
    found = shutil.which("matlab")
    if found:
        return found
    windows_fallback = Path(r"D:\Matlab\bin\matlab.exe")
    if windows_fallback.exists():
        return str(windows_fallback)
    raise RuntimeError(
        "MATLAB was not found. Load the HPC MATLAB module or set PFAS_MATLAB_CMD."
    )


def read_one_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise RuntimeError(f"Expected one data row in {path}, found {len(rows)}")
    return rows[0]


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def bool_field(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def validate_response(
    row: dict[str, str], case_id: int, response_path: Path,
    maximum_inventory_error: float | None = None,
) -> None:
    # PFOS-Kralchevsky is best conditioned with bounded Gamma-logit. At the
    # much lower PFOA SSL concentrations, aqueous C is the equivalent but
    # markedly better-conditioned primary variable.
    expected_primary = "Gamma_logit" if case_id == 2 else "aqueous_concentration"
    checks = {
        "Case_ID": str(case_id),
        "Interface_Primary": expected_primary,
        "PFAS_Inventory_Assembly": "single_total_inventory_BE",
        "Solver_Formulation_Version": SOLVER_VERSION,
    }
    for field, expected in checks.items():
        actual = row.get(field)
        if actual != expected:
            raise RuntimeError(
                f"Invalid solver provenance in {response_path}: "
                f"{field}={actual!r}, expected {expected!r}"
            )
    if maximum_inventory_error is None:
        maximum_inventory_error = float(os.environ.get("PFAS_MAX_INVENTORY_ERROR", "1e-3"))
    pfas_error = float(row["PFAS_Mass_Error_Inventory_Normalized"])
    water_error = float(row["Water_Mass_Error_Flux_Normalized"])
    if not math.isfinite(pfas_error) or pfas_error > maximum_inventory_error:
        raise RuntimeError(
            f"PFAS inventory error {pfas_error:.3e} exceeds "
            f"{maximum_inventory_error:.3e}: {response_path}"
        )
    if not math.isfinite(water_error) or water_error > 1e-3:
        raise RuntimeError(f"Water balance error {water_error:.3e} is unacceptable: {response_path}")
    daf = float(row["Groundwater_DAF"])
    af = float(row["Groundwater_AF"])
    leachate = float(row["Peak_Leachate_ug_L"])
    receptor = float(row["Peak_Receptor_ug_L"])
    if not (daf > 1 and 0 < af < 1 and math.isclose(daf * af, 1.0, rel_tol=1e-10)):
        raise RuntimeError(f"Groundwater DAF/AF are not reciprocal in {response_path}")
    if not math.isclose(receptor, leachate * af, rel_tol=1e-10, abs_tol=1e-15):
        raise RuntimeError(f"Peak receptor concentration is inconsistent with AF in {response_path}")


def scan_log(log_path: Path) -> dict[str, int]:
    text = log_path.read_text(encoding="utf-8", errors="replace").lower()
    return {
        "singular_or_badly_scaled": (
            text.count("singular") + text.count("badly scaled") + text.count("nearly singular")
        ),
        "newton_failures": text.count("newton solves did not converge"),
        "matlab_errors": text.count("error using"),
    }


def run_matlab_live(package: Path, env: dict[str, str], log_path: Path) -> float:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command = [matlab_executable()]
    if env.get("PFAS_MATLAB_SINGLE_THREAD", "false").strip().lower() in {"1", "true", "yes"}:
        command.append("-singleCompThread")
    command.extend(["-batch", "run('PFAS_transport_SSL_forward.m')"])
    say(f"MATLAB log: {log_path}")
    started = time.time()
    with log_path.open("w", encoding="utf-8", buffering=1) as log:
        process = subprocess.Popen(
            command,
            cwd=package,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
        return_code = process.wait()
    elapsed = time.time() - started
    if return_code != 0:
        raise RuntimeError(f"MATLAB exited with code {return_code}; see {log_path}")
    return elapsed


def request_matches(
    row: dict[str, str], case_id: int, scale: float, cycles: int,
    short_days: int | None, target: float,
) -> bool:
    expected_years = short_days / 365.25 if short_days else cycles * 3652 / 365.25
    try:
        return (
            int(row["Case_ID"]) == case_id
            and math.isclose(float(row["Initial_Profile_Scale"]), scale, rel_tol=1e-10)
            and math.isclose(float(row["Simulated_Years"]), expected_years, rel_tol=1e-8)
            and math.isclose(float(row["Target_MCL_ug_L"]), target, rel_tol=1e-12)
            and row.get("Solver_Formulation_Version") == SOLVER_VERSION
        )
    except (KeyError, TypeError, ValueError):
        return False


def run_forward(
    package: Path,
    case_id: int,
    scale: float,
    cycles: int,
    output_dir: Path,
    response_only: bool,
    target: float,
    output_interval_days: int,
    progress_interval_days: float,
    short_days: int | None = None,
    allow_cache: bool = True,
) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    response_path = output_dir / "ssl_forward_response.csv"
    if allow_cache and response_path.exists() and response_path.stat().st_size > 0:
        row = read_one_row(response_path)
        if request_matches(row, case_id, scale, cycles, short_days, target):
            validate_response(row, case_id, response_path)
            row["Runner_Wall_Seconds"] = "0"
            row["Cache_Source"] = "trial_directory"
            say(f"RESUME cached case {case_id}, scale={scale:.8g}: {response_path}")
            return row
        raise RuntimeError(
            f"Existing output does not match this request: {response_path}. "
            "Choose a new --output-root or remove only that incomplete trial directory."
        )

    env = os.environ.copy()
    env.update({
        "PFAS_CASE_ID": str(case_id),
        "PFAS_INITIAL_SCALE": f"{scale:.17g}",
        "PFAS_SIMULATION_CYCLES": str(cycles),
        "PFAS_OUTPUT_INTERVAL_DAYS": str(output_interval_days),
        "PFAS_PROGRESS_INTERVAL_DAYS": f"{progress_interval_days:.17g}",
        "PFAS_OUTPUT_ROOT": str(output_dir),
        "PFAS_SSL_RESPONSE_ONLY": "true" if response_only else "false",
        "PFAS_RECEPTOR_MCL_UG_L": f"{target:.17g}",
        "PFAS_AFFF_ACTIVE_YEARS": env.get("PFAS_AFFF_ACTIVE_YEARS", "30"),
        "PFAS_CONSERVATIVE_INVENTORY": "true",
        "PFAS_CONSERVATIVE_STEP_MASS_TOL": env.get(
            "PFAS_CONSERVATIVE_STEP_MASS_TOL", "1e-8"
        ),
        "PFAS_MINIMUM_DT_SECONDS": env.get("PFAS_MINIMUM_DT_SECONDS", "1e-6"),
        "PFAS_VERBOSE_TIMESTEPS": env.get("PFAS_VERBOSE_TIMESTEPS", "false"),
    })
    if case_id == 2:
        env["PFAS_INTERFACE_PRIMARY"] = "gamma"
    else:
        env.pop("PFAS_INTERFACE_PRIMARY", None)
    if short_days:
        env["PFAS_SHORT_TEST_DAYS"] = str(short_days)
    else:
        env.pop("PFAS_SHORT_TEST_DAYS", None)

    say(
        f"START case={case_id} ({CASES[case_id]['name']}) scale={scale:.8g} "
        f"horizon={'%d days' % short_days if short_days else '%d cycles' % cycles}"
    )
    elapsed = run_matlab_live(package, env, output_dir / "matlab.log")
    if not response_path.exists():
        raise RuntimeError(f"MATLAB finished without producing {response_path}")
    row = read_one_row(response_path)
    validate_response(row, case_id, response_path)
    warnings = scan_log(output_dir / "matlab.log")
    if any(warnings.values()):
        raise RuntimeError(f"Numerical warning/error text detected in MATLAB log: {warnings}")
    row["Runner_Wall_Seconds"] = f"{elapsed:.6f}"
    row["Cache_Source"] = "computed"
    row.update({f"Log_{key}": str(value) for key, value in warnings.items()})
    say(
        f"PASS case={case_id} wall={elapsed:.1f}s "
        f"PFAS_error={float(row['PFAS_Mass_Error_Inventory_Normalized']):.3e} "
        f"water_error={float(row['Water_Mass_Error_Flux_Normalized']):.3e}"
    )
    return row


def run_transport(args: argparse.Namespace) -> int:
    package = Path(__file__).resolve().parent
    root = args.output_root.resolve()
    status_root = root / "status"
    status_root.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, str]] = []
    for case_id in args.cases:
        case_name = CASES[case_id]["name"]
        done_path = status_root / f"case_{case_id}_{case_name}.done"
        case_dir = root / f"case_{case_id}_{case_name}"
        if done_path.exists() and not args.force:
            say(f"SKIP completed transport case {case_id}: {done_path}")
            response = case_dir / "ssl_forward_response.csv"
            if response.exists():
                summaries.append(read_one_row(response))
            continue
        row = run_forward(
            package, case_id, args.scale, args.cycles, case_dir, False,
            args.target, args.output_interval_days, args.progress_interval_days,
            args.short_days, allow_cache=not args.force,
        )
        row["Run_Type"] = "transport_smoke" if args.short_days else "transport"
        summaries.append(row)
        done_path.write_text(
            f"case_id={case_id}\ncase_name={case_name}\ncompleted={time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n",
            encoding="utf-8",
        )
        write_rows(root / "transport_summary.csv", summaries)
    say(f"Transport workflow complete: {root}")
    return 0


def run_ssl(args: argparse.Namespace) -> int:
    package = Path(__file__).resolve().parent
    case_id = args.case
    config = CASES[case_id]
    cycles = args.cycles or config["cycles"]
    case_root = args.output_root.resolve() / f"case_{case_id}_{config['name']}"
    case_root.mkdir(parents=True, exist_ok=True)
    done_path = case_root / "SSL.done"
    if done_path.exists() and not args.force:
        say(f"SKIP completed SSL case {case_id}: {done_path}")
        return 0

    if args.smoke_days:
        row = run_forward(
            package, case_id, config["seed"], 1, case_root / "smoke", False,
            args.target, 1, 1, args.smoke_days, allow_cache=not args.force,
        )
        row["Run_Type"] = "ssl_forward_smoke"
        write_rows(case_root / "smoke_summary.csv", [row])
        say(json.dumps(row, indent=2))
        return 0

    history_path = case_root / "inverse_history.csv"
    history: list[dict[str, str]] = []
    trial_number = 0

    def trial(scale: float, label: str) -> dict[str, str]:
        nonlocal trial_number
        trial_number += 1
        row = run_forward(
            package, case_id, scale, cycles,
            case_root / "trials" / f"{trial_number:02d}_{label}_s{scale:.8g}",
            True, args.target, 365, args.progress_interval_days,
            allow_cache=not args.force,
        )
        row["Trial"] = str(trial_number)
        row["Trial_Label"] = label
        history.append(row)
        write_rows(history_path, history)
        if not bool_field(row["Horizon_Peak_Covered"]):
            raise RuntimeError(
                f"Case {case_id} horizon did not demonstrate a passed receptor peak. "
                "Increase --cycles and use a new --output-root."
            )
        say(
            f"INVERSE case={case_id} trial={trial_number} scale={scale:.8g} "
            f"peak={float(row['Peak_Receptor_ug_L']):.8g} ug/L target={args.target:.8g}"
        )
        return row

    def peak(row: dict[str, str]) -> float:
        return float(row["Peak_Receptor_ug_L"])

    low_scale = config.get("initial_low", config["seed"] / 4.0)
    high_scale = config.get("initial_high", config["seed"] * 4.0)
    low_row = trial(low_scale, "bracket_low")
    high_row = trial(high_scale, "bracket_high")

    for _ in range(12):
        if peak(low_row) <= args.target:
            break
        high_scale, high_row = low_scale, low_row
        low_scale /= 4.0
        low_row = trial(low_scale, "expand_low")
    else:
        raise RuntimeError("Could not find a scale below the receptor target")

    for _ in range(12):
        if peak(high_row) > args.target:
            break
        low_scale, low_row = high_scale, high_row
        high_scale *= 4.0
        high_row = trial(high_scale, "expand_high")
    else:
        raise RuntimeError("Could not find a scale above the receptor target")

    if not (low_scale < high_scale and peak(low_row) <= args.target < peak(high_row)):
        raise RuntimeError("Invalid inverse bracket or non-monotone endpoint responses")

    if high_scale / low_scale > 1.0 + args.tolerance:
        for iteration in range(1, args.max_iterations + 1):
            mid_scale = math.sqrt(low_scale * high_scale)
            mid_row = trial(mid_scale, f"bisect_{iteration}")
            if peak(mid_row) <= args.target:
                low_scale, low_row = mid_scale, mid_row
            else:
                high_scale, high_row = mid_scale, mid_row
            if high_scale / low_scale <= 1.0 + args.tolerance:
                break

    final_row = run_forward(
        package, case_id, low_scale, cycles, case_root / "final", False,
        args.target, 30, args.progress_interval_days, allow_cache=not args.force,
    )
    final_row.update({
        "Run_Type": "accepted_final",
        "Safe_Lower_Scale": f"{low_scale:.17g}",
        "Unsafe_Upper_Scale": f"{high_scale:.17g}",
        "Scale_Bracket_Ratio": f"{high_scale / low_scale:.17g}",
        "Inverse_Trials": str(trial_number),
        "Requested_Relative_Tolerance": f"{args.tolerance:.17g}",
    })
    if not bool_field(final_row["Horizon_Peak_Covered"]):
        raise RuntimeError("Final run did not demonstrate a passed receptor peak")
    write_rows(case_root / "SSL_summary.csv", [final_row])
    done_path.write_text(
        f"case_id={case_id}\ncompleted={time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n",
        encoding="utf-8",
    )
    say(json.dumps(final_row, indent=2))
    return 0


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest="workflow", required=True)

    transport = sub.add_parser("transport", help="Run fixed-initial-profile transport cases")
    transport.add_argument("--cases", type=int, nargs="+", choices=CASES, default=list(CASES))
    transport.add_argument("--output-root", type=Path, default=Path("Results_transport_SI"))
    transport.add_argument("--cycles", type=int, default=5)
    transport.add_argument("--scale", type=float, default=1.0)
    transport.add_argument("--target", type=float, default=0.004)
    transport.add_argument("--output-interval-days", type=int, default=30)
    transport.add_argument("--progress-interval-days", type=float, default=30)
    transport.add_argument("--short-days", type=int)
    transport.add_argument("--force", action="store_true")
    transport.set_defaults(func=run_transport)

    ssl = sub.add_parser("ssl", help="Run one inverse soil-screening-level case")
    ssl.add_argument("--case", type=int, choices=CASES, required=True)
    ssl.add_argument("--output-root", type=Path, default=Path("Results_SSL_SI"))
    ssl.add_argument("--cycles", type=int)
    ssl.add_argument("--target", type=float, default=0.004)
    ssl.add_argument("--tolerance", type=float, default=0.05)
    ssl.add_argument("--max-iterations", type=int, default=6)
    ssl.add_argument("--progress-interval-days", type=float, default=365)
    ssl.add_argument("--smoke-days", type=int)
    ssl.add_argument("--force", action="store_true")
    ssl.set_defaults(func=run_ssl)
    return root


def main() -> int:
    args = parser().parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
