import argparse
import os
import unittest
from pathlib import Path
from sys import path

import vcf

from juno_mapping import JunoMapping


@unittest.skipIf(
    (not Path("tests/gordonia_s_mutated_R1.fastq.gz").exists())
    or (not Path("tests/gordonia_s_mutated_R2.fastq.gz").exists()),
    "Skipped because test data is missing",
)
class TestJunoMappingPipelineSingularity(unittest.TestCase):
    """Testing the JunoMapping class"""

    reference_sample_dict = {
        "gordonia_s_mutated": {
            "R1": str(Path("tests/gordonia_s_mutated_R1.fastq.gz").resolve()),
            "R2": str(Path("tests/gordonia_s_mutated_R2.fastq.gz").resolve()),
        }
    }

    expected_files = [
        "multiqc/multiqc.html",
        "audit_trail/log_git.yaml",
        "audit_trail/log_pipeline.yaml",
        "audit_trail/log_conda.txt",
        "audit_trail/snakemake_report.html",
        "audit_trail/sample_sheet.yaml",
        "audit_trail/user_parameters.yaml",
    ]

    vcf_dict = {
        820: {"REF": "T", "ALT": "C"},
        5659: {"REF": "A", "ALT": "G"},
        15162: {"REF": "T", "ALT": "A"},
        46679: {"REF": "C", "ALT": "G"},
    }

    output_dir = Path("pipeline_test_output_singularity")
    input_dir = "tests"

    @classmethod
    def setUpClass(cls, output_dir=output_dir) -> None:
        os.system(f"rm -rf {output_dir}")

    @classmethod
    def tearDownClass(cls, output_dir=output_dir) -> None:
        os.system(f"rm -rf {output_dir}")

    def test_010_junomapping_run_in_singularity(self) -> None:
        """Testing the pipeline runs properly with real samples"""
        # Check if running in Github actions
        if not Path("/home/runner").exists():
            kraken_db = Path.home().joinpath("kraken-database")
            assert kraken_db.exists(), "Kraken database not found"
        else:
            kraken_db = Path("/home/runner/kraken-database")

        pipeline = JunoMapping(
            argv=[
                "-i",
                self.input_dir,
                "-o",
                str(self.output_dir),
                "--species",
                "mycobacterium_tuberculosis",
                "--reference",
                "tests/gordonia.fasta",
                "--disable-mask",
                "--local",
                "--snakemake-args",
                "cores=1",
                "nodes=2",
                "--db-dir",
                str(kraken_db),
                "--prefix",
                "/home/runner/sing_containers",
            ]
        )
        pipeline.run()

        global pipeline_sample_dict
        pipeline_sample_dict = pipeline.sample_dict

    def test_020_sample_dict(self):
        self.assertDictEqual(
            pipeline_sample_dict,
            self.reference_sample_dict,
        )

    def test_030_expected_files(self):
        for file_ in self.expected_files:
            self.assertTrue(self.output_dir.joinpath(file_).exists())

    def test_040_mutations(self):
        reader = vcf.Reader(open(f"{self.output_dir}/variants/gordonia_s_mutated.vcf"))

        for var in reader:
            self.assertEqual(self.vcf_dict[var.POS]["REF"], var.REF)
            self.assertEqual(self.vcf_dict[var.POS]["ALT"], var.ALT[0])
