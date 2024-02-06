# type: ignore

import argparse
import os
import unittest
from pathlib import Path
from sys import path

from juno_mapping import JunoMapping


def make_non_empty_file(file_path: Path, num_lines: int = 1000) -> None:
    content = "a\n" * num_lines
    with open(file_path, "w") as file_:
        file_.write(content)


class TestJunoMappingDryRun(unittest.TestCase):
    """Testing the junomapping class (code specific for this pipeline)"""

    @classmethod
    def setUpClass(cls) -> None:
        fake_dirs = ["fake_dir_wsamples", "fake_empty_dir"]

        fake_files = [
            Path("fake_dir_wsamples/sample1_R1.fastq"),
            Path("fake_dir_wsamples/sample1_R2.fastq.gz"),
            Path("fake_dir_wsamples/sample2_R1_filt.fq"),
            Path("fake_dir_wsamples/sample2_R2_filt.fq.gz"),
            Path("fake_dir_wsamples/1234_R1.fastq.gz"),
            Path("fake_dir_wsamples/1234_R2.fastq.gz"),
            Path("reference.fasta"),
        ]

        for folder in fake_dirs:
            Path(folder).mkdir(exist_ok=True)
        for file_ in fake_files:
            make_non_empty_file(file_)

    @classmethod
    def tearDownClass(cls) -> None:
        fake_dirs = [
            "fake_dir_wsamples",
            "fake_empty_dir",
            "test_output",
            "reference.fasta",
        ]
        for folder in fake_dirs:
            os.system("rm -rf {}".format(str(folder)))

    def test_fails_with_empty_dir(self) -> None:
        """Testing the pipeline fails if input_dir is empty"""
        with self.assertRaisesRegex(
            ValueError,
            "does not contain any files with the expected format/naming",
        ):
            argv = [
                "-i",
                "fake_empty_dir",
                "-o",
                "test_output",
                "-n",
                "--reference",
                "reference.fasta",
                "-s",
                "mycobacterium_tuberculosis",
            ]
            pipeline = JunoMapping(argv=argv)
            pipeline.setup()

    def test_dryrun(self) -> None:
        """Testing the pipeline runs properly as a dry run"""
        full_input_dir = Path("fake_dir_wsamples").resolve()
        pipeline = JunoMapping(
            argv=[
                "-i",
                "fake_dir_wsamples",
                "-o",
                "test_output",
                "-n",
                "--reference",
                "reference.fasta",
                "-s",
                "mycobacterium_tuberculosis",
            ]
        )
        pipeline.run()
        expected_sample_sheet = {
            "sample1": {
                "R1": str(full_input_dir.joinpath("sample1_R1.fastq")),
                "R2": str(full_input_dir.joinpath("sample1_R2.fastq.gz")),
            },
            "sample2": {
                "R1": str(full_input_dir.joinpath("sample2_R1_filt.fq")),
                "R2": str(full_input_dir.joinpath("sample2_R2_filt.fq.gz")),
            },
            "1234": {
                "R1": str(full_input_dir.joinpath("1234_R1.fastq.gz")),
                "R2": str(full_input_dir.joinpath("1234_R2.fastq.gz")),
            },
        }
        self.assertEqual(pipeline.sample_dict, expected_sample_sheet)

    def test_junomapping_dryrun_user_params(self) -> None:
        """Check whether user parameters are correct"""
        pipeline = JunoMapping(
            argv=[
                "-i",
                "fake_dir_wsamples",
                "-o",
                "test_output",
                "-n",
                "-s",
                "mycobacterium_tuberculosis",
                "--reference",
                "reference.fasta",
            ]
        )
        full_input_dir = Path("fake_dir_wsamples").resolve()
        pipeline.run()
        # reference genome cannot be left at default because snakemake dry run will look for it as an input file and fail if absent
        expected_user_param_values = {
            "species": "mycobacterium_tuberculosis",
            "reference": "reference.fasta",
            "minimum_allele_frequency": 0.8,
            "minimum_depth": 10,
            "disable_mask": "False",
            "mask_bed": "files/RLC_Farhat.bed",
        }
        for key, value in expected_user_param_values.items():
            self.assertTrue(
                pipeline.user_parameters[key] == value,
                msg=f"{key} is {pipeline.user_parameters[key]} in user_parameters, expected {value}",
            )

    def test_junomapping_dryrun_changed_user_params(self) -> None:
        """Check whether user parameters can be changed correctly"""
        pipeline = JunoMapping(
            argv=[
                "-i",
                "fake_dir_wsamples",
                "-o",
                "test_output",
                "-n",
                "-s",
                "mycobacterium_tuberculosis",
                "--reference",
                "reference.fasta",
                "--disable-mask",
                "-smaf",
                "0.5",
                "-md",
                "20",
            ]
        )
        full_input_dir = Path("fake_dir_wsamples").resolve()
        pipeline.run()
        # reference genome cannot be left at default because snakemake dry run will look for it as an input file and fail if absent
        expected_user_param_values = {
            "species": "mycobacterium_tuberculosis",
            "reference": "reference.fasta",
            "soft_filter_minimum_allele_frequency": 0.5,
            "minimum_depth": 20,
            "disable_mask": "True",
            "mask_bed": "None",
        }
        for key, value in expected_user_param_values.items():
            self.assertTrue(
                pipeline.user_parameters[key] == value,
                msg=f"{key} is {pipeline.user_parameters[key]} in user_parameters, expected {value}",
            )

    def test_junomapping_fails_with_unsupported_species(self) -> None:
        """Testing the pipeline fails if a fake species is provided"""
        with self.assertRaises(argparse.ArgumentError):
            pipeline = JunoMapping(
                argv=[
                    "-i",
                    "fake_dir_wsamples",
                    "-o",
                    "test_output",
                    "-s",
                    "fake_species",
                    "--reference",
                    "reference.fasta",
                ]
            )
            pipeline.parser.exit_on_error = False  # type: ignore
            pipeline.setup()


if __name__ == "__main__":
    unittest.main()
