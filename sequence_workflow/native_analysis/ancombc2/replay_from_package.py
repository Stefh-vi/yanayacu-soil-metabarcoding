"""Portable entry point for the frozen ANCOM-BC2 fit; no implicit execution.

--check-only verifies archive inputs without R, writing, network access or fits.
An explicit --rscript and fresh --output are required to fit the same models.
"""
from pathlib import Path
import argparse
import concurrent.futures
import datetime
import hashlib
import json
import os
import subprocess

ROOT = Path(__file__).resolve().parent
EXPECTED = {
    "frozen_inputs/FROZEN_CONFIG.json": "c26e0f6d7ad01cd923da725ef7958cc936de60ecd9c1c8f10d19be7a72c8362a",
    "numerical_extension/NUMERICAL_EXTENSION_CONFIG.json": "ec9cf8ba918854e306b27672b61b3c190e85cf56c6599c70b85f40637992d57a",
    "numerical_extension/run_extended_ancombc2.R": "c54c6ee27fe9e837c9a28065cc41f6b19e8ca4ce044b2bad5ba6a57a14a204ae",
    "executed_run/CONDA_INSTALLED_HASH_LOCK.json": "7b19cca8866efec62c570aa2829b3adeeaa3e7ab07c852e08db2cc5ad6889f36",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def preflight():
    checked = []
    for name, expected in EXPECTED.items():
        path = ROOT / name
        require(path.is_file() and sha(path) == expected, "Hash mismatch: " + name)
        checked.append({"path": name, "sha256": expected})
    cfg = json.loads((ROOT / "frozen_inputs/FROZEN_CONFIG.json").read_text())
    for row in cfg["inputs"]:
        path = ROOT / "frozen_inputs" / row["path"]
        require(path.is_file() and sha(path) == row["sha256"], "Input hash mismatch: " + row["path"])
        checked.append({"path": "frozen_inputs/" + row["path"], "sha256": row["sha256"]})
    require(cfg["formula"] == "locality + zone", "Unexpected formula")
    require([cfg["marker_families"][m]["n_hypotheses"] for m in ["16S", "ITS"]] == [356, 179], "Unexpected families")
    ext = json.loads((ROOT / "numerical_extension/NUMERICAL_EXTENSION_CONFIG.json").read_text())
    require(ext["base_config_sha256"] == EXPECTED["frozen_inputs/FROZEN_CONFIG.json"], "Extension/config mismatch")
    require(ext["change"] == {"em_max_iter_from": 1000, "em_max_iter_to": 10000, "tolerance_unchanged": 1e-5}, "Unexpected extension")
    lock = json.loads((ROOT / "executed_run/CONDA_INSTALLED_HASH_LOCK.json").read_text())
    require(len(lock["packages"]) == 228, "Unexpected package count")
    explicit = [s.strip() for s in (ROOT / "executed_run/CONDA_EXPLICIT_WITH_MD5.txt").read_text().splitlines() if s.strip() and not s.startswith(("#", "@"))]
    require(set(explicit) == {p["url"] + "#" + p["md5"] for p in lock["packages"]}, "Explicit lock differs from hashed package lock")
    return {"status": "PASS_INPUT_AND_REPLAY_CODE_PREFLIGHT", "checked_files": checked,
            "explicit_packages": len(explicit), "fit_executed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--rscript", help="Absolute Rscript path in the recreated pinned Linux environment")
    parser.add_argument("--output", help="Fresh directory outside code/ancombc2")
    parser.add_argument("--marker", choices=["16S", "ITS", "both"], default="both")
    args = parser.parse_args()
    verified = preflight()
    if args.check_only:
        print(json.dumps(verified, indent=2))
        return
    require(args.rscript and args.output, "Fitting requires explicit --rscript and --output")
    # Keep the supplied prefix/bin path even if Rscript is a symlink.
    rscript = Path(args.rscript).expanduser().absolute()
    require(rscript.is_file(), "Rscript does not exist")
    out = Path(args.output).resolve()
    require(not out.exists(), "Choose a fresh output directory; existing results are never overwritten")
    require(ROOT != out and ROOT not in out.parents, "Output must be outside the archived code/ancombc2 directory")
    env = os.environ.copy()
    r_library = rscript.parent.parent / "lib/R/library"
    require(r_library.is_dir(), "Expected the pinned Conda R library layout")
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LC_ALL="C.UTF-8",
               R_LIBS_USER=str(r_library), R_LIBS_SITE=str(r_library))
    version_command = [str(rscript), "--vanilla", "-e",
                       'stopifnot(as.character(getRversion()) == "4.5.3", as.character(packageVersion("ANCOMBC")) == "2.14.0"); sessionInfo()']
    version = subprocess.run(version_command, env=env, capture_output=True, text=True, check=True)
    out.mkdir(parents=True, exist_ok=False)
    (out / "PREFLIGHT_INPUT_HASHES.json").write_text(json.dumps(verified, indent=2) + "\n")
    (out / "REPLAY_ENVIRONMENT_PREFLIGHT.txt").write_text(version.stdout + version.stderr)
    markers = ["16S", "ITS"] if args.marker == "both" else [args.marker]

    def run(marker):
        dest = out / marker
        dest.mkdir()
        cmd = [str(rscript), "--vanilla", str(ROOT / "numerical_extension/run_extended_ancombc2.R"),
               str(ROOT / "frozen_inputs"), str(dest / "results"), marker]
        start = utc()
        with (dest / "stdout.log").open("wb") as stdout, (dest / "stderr.log").open("wb") as stderr:
            result = subprocess.run(cmd, env=env, stdout=stdout, stderr=stderr)
        status_path = dest / "results" / marker / "STATUS.json"
        status = json.loads(status_path.read_text()) if status_path.is_file() else None
        return {"marker": marker, "argv": cmd, "returncode": result.returncode,
                "started_utc": start, "finished_utc": utc(), "native_status": status}

    with concurrent.futures.ThreadPoolExecutor(max_workers=len(markers)) as pool:
        results = list(pool.map(run, markers))
    receipt = {"status": "EXECUTION_FINISHED_REQUIRES_RESULT_REVIEW", "fit_executed": True,
               "input_config_sha256": EXPECTED["frozen_inputs/FROZEN_CONFIG.json"],
               "native_R_script_sha256": EXPECTED["numerical_extension/run_extended_ancombc2.R"],
               "results": results}
    (out / "REPRODUCTION_EXECUTION_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))
    if any(r["returncode"] or r["native_status"] is None or r["native_status"]["status"] == "FIT_FAILED" for r in results):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
