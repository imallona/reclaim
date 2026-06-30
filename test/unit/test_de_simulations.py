"""
Pytest wrapper for the R-side unit tests of workflow/scripts/de_simulations.R.

Skips when Rscript or any required Bioconductor package is missing so the
suite still runs in CI environments without R installed. Locally, with the
edger conda env active (or system-wide R + edgeR + RUVSeq + EDASeq available),
this calls the simulator helpers end to end.
"""
import os
import shutil
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "de_simulations.R"
R_TESTS = REPO_ROOT / "test" / "unit" / "test_de_simulations.R"
REQUIRED_R_PACKAGES = ("edgeR", "RUVSeq", "EDASeq", "Biobase", "matrixStats",
                       "ggplot2", "dplyr", "tidyr", "viridis")


def _rscript_available():
    return shutil.which("Rscript") is not None


def _r_packages_available():
    if not _rscript_available():
        return False
    expr = (
        "pkgs <- c(" + ", ".join(repr(p) for p in REQUIRED_R_PACKAGES) + ");"
        "missing <- pkgs[!pkgs %in% rownames(installed.packages())];"
        "if (length(missing) > 0) { cat(missing, sep='\\n'); quit(status=1) }"
    )
    try:
        result = subprocess.run(
            ["Rscript", "-e", expr],
            capture_output=True, text=True, timeout=60
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


@pytest.mark.skipif(not _rscript_available(),
                    reason="Rscript not in PATH")
@pytest.mark.skipif(not _r_packages_available(),
                    reason="required R Bioconductor packages not installed")
def test_de_simulations_r_unit_tests():
    env = os.environ.copy()
    env["DE_SIM_SCRIPT"] = str(SCRIPT)
    result = subprocess.run(
        ["Rscript", str(R_TESTS)],
        cwd=str(REPO_ROOT), env=env,
        capture_output=True, text=True, timeout=600
    )
    assert result.returncode == 0, (
        "R unit tests failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )
    assert "OK: all de_simulations unit tests passed" in result.stdout


def test_de_simulations_script_exists():
    assert SCRIPT.exists(), f"Missing R script at {SCRIPT}"


def test_de_simulations_module_exists():
    module = REPO_ROOT / "workflow" / "modules" / "de_simulations.snmk"
    assert module.exists(), f"Missing snakemake module at {module}"
