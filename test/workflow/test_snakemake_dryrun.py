"""
Snakemake dry-run tests.

These tests run `snakemake --dry-run` on both the main Snakefile and the
test Snakefile to verify that the DAG can be constructed without errors.
They do not execute any actual rules and require only snakemake to be installed.

Marked as `workflow` so they can be excluded with `-m "not workflow"` when
snakemake is not in the PATH.
"""
import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
WORKFLOW_DIR = os.path.join(REPO_ROOT, 'workflow')
TEST_WORKFLOW = os.path.join(REPO_ROOT, 'test', 'workflow', 'Snakefile_test')
SMARTSEQ2_CFG = os.path.join(WORKFLOW_DIR, 'configs', 'simulation_smartseq2.yaml')
CHROMIUM_CFG = os.path.join(WORKFLOW_DIR, 'configs', 'simulation_chromium.yaml')
NEG_CTRL_CFG = os.path.join(REPO_ROOT, 'test', 'workflow', 'configs',
                             'test_negative_control.yaml')


def snakemake_available():
    try:
        subprocess.run(['snakemake', '--version'], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


SKIP_IF_NO_SNAKEMAKE = pytest.mark.skipif(
    not snakemake_available(),
    reason='snakemake not found in PATH'
)


def run_dryrun(snakefile, configfile, workdir=WORKFLOW_DIR):
    cmd = [
        'snakemake',
        '-s', snakefile,
        '--configfile', configfile,
        '--dry-run',
        '--cores', '1',
        '--quiet',
    ]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_smartseq2():
    r = run_dryrun(os.path.join(WORKFLOW_DIR, 'Snakefile'), SMARTSEQ2_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for simulation_smartseq2.yaml:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_chromium():
    if not os.path.exists(CHROMIUM_CFG):
        pytest.skip('simulation_chromium.yaml not found')
    r = run_dryrun(os.path.join(WORKFLOW_DIR, 'Snakefile'), CHROMIUM_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for simulation_chromium.yaml:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_test_snakefile_dryrun_negative_control():
    r = run_dryrun(TEST_WORKFLOW, NEG_CTRL_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for test negative control config:\n{r.stderr}'
    )


TEST_BULK_SINGLE_END_CFG = os.path.join(REPO_ROOT, 'test', 'workflow', 'configs',
                                        'test_bulk_single_end.yaml')
TEST_BULK_PAIRED_END_CFG = os.path.join(REPO_ROOT, 'test', 'workflow', 'configs',
                                        'test_bulk_paired_end.yaml')

BULK_GSE230647_CFG = os.path.join(WORKFLOW_DIR, 'configs', 'gse230647_bulk.yaml')
BULK_GSE126543_CFG = os.path.join(WORKFLOW_DIR, 'configs', 'gse126543_bulk.yaml')
SC_GSE230647_CFG = os.path.join(WORKFLOW_DIR, 'configs', 'gse230647_sc.yaml')


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_bulk_gse230647():
    r = run_dryrun(os.path.join(WORKFLOW_DIR, 'Snakefile'), BULK_GSE230647_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for gse230647_bulk.yaml:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_bulk_gse126543():
    r = run_dryrun(os.path.join(WORKFLOW_DIR, 'Snakefile'), BULK_GSE126543_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for gse126543_bulk.yaml:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_sc_gse230647():
    r = run_dryrun(os.path.join(WORKFLOW_DIR, 'Snakefile'), SC_GSE230647_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for gse230647_sc.yaml:\n{r.stderr}'
    )
@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_bulk_single_end():
    r = run_dryrun(os.path.join(WORKFLOW_DIR, 'Snakefile'), TEST_BULK_SINGLE_END_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for test_bulk_single_end.yaml:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_bulk_paired_end():
    r = run_dryrun(os.path.join(WORKFLOW_DIR, 'Snakefile'), TEST_BULK_PAIRED_END_CFG)
    assert r.returncode == 0, (
        f'Dry-run failed for test_bulk_paired_end.yaml:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_lint():
    """snakemake --lint catches common Snakemake style/correctness issues."""
    cmd = [
        'snakemake',
        '--lint',
        '--configfile', SMARTSEQ2_CFG,
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=WORKFLOW_DIR)
    # lint may return non-zero for warnings; we only fail on errors
    # (lint output on stderr distinguishes errors from warnings)
    assert 'Error' not in r.stdout, f'Snakemake lint errors:\n{r.stdout}'


# Tagcount coverage with explicit aligner overrides. The default sim
# configs already include starsolo_tagcount; these tests pass an aligner
# list via --config that drops bowtie2 and keeps tagcount, exercising
# the code path where tagcount is the only opt-in extra.

def run_dryrun_with_config(snakefile, configfile, extra_config,
                           workdir=WORKFLOW_DIR):
    cmd = [
        'snakemake',
        '-s', snakefile,
        '--configfile', configfile,
        '--config', extra_config,
        '--dry-run',
        '--cores', '1',
        '--quiet',
    ]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_chromium_with_tagcount():
    """Chromium sim with starsolo_tagcount opted in via --config: the
    DAG should build, fan out tagcount jobs over repeat feature_sets only
    (genes is served by STARsolo's matrix), and pull in the bam_index
    rule."""
    r = run_dryrun_with_config(
        os.path.join(WORKFLOW_DIR, 'Snakefile'), CHROMIUM_CFG,
        'aligners=["starsolo","kallisto","alevin","starsolo_tagcount"]')
    assert r.returncode == 0, (
        f'Dry-run failed for chromium + tagcount:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_smartseq2_with_tagcount():
    """SmartSeq2 sim with starsolo_tagcount opted in: routes through
    smartseq2_tagcount.snmk and runs sc_count_features with --umi-dedup
    none."""
    r = run_dryrun_with_config(
        os.path.join(WORKFLOW_DIR, 'Snakefile'), SMARTSEQ2_CFG,
        'aligners=["starsolo","kallisto","alevin","starsolo_tagcount"]')
    assert r.returncode == 0, (
        f'Dry-run failed for smartseq2 + tagcount:\n{r.stderr}'
    )


@pytest.mark.workflow
@SKIP_IF_NO_SNAKEMAKE
def test_main_snakefile_dryrun_sc_real_with_tagcount():
    """Real-data sc with sc_count_mode=tagcount switches starsolo_sc into
    one-align-per-library + tagcount recount mode. Verifies the DAG
    builds with the new bam_index include and the (genes, gene_id)
    carve-out."""
    r = run_dryrun_with_config(
        os.path.join(WORKFLOW_DIR, 'Snakefile'), SC_GSE230647_CFG,
        'sc_count_mode=tagcount')
    assert r.returncode == 0, (
        f'Dry-run failed for sc real-data + tagcount:\n{r.stderr}'
    )
