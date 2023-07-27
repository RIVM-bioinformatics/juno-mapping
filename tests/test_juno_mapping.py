import argparse

import csv
import os
from pathlib import Path
from sys import path
import unittest
import vcf

main_script_path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
path.insert(0, main_script_path)
from juno_mapping import JunoMapping


def make_non_empty_file(file_path: Path, num_lines: int = 1000) -> None:
    content = "a\n" * num_lines
    with open(file_path, "w") as file_:
        file_.write(content)

def fh(fname, mode='rt'):
    return open(os.path.join(os.path.dirname(__file__), fname), mode)



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

        # with open("fake_dir_wsamples/fake_metadata.csv", mode="w") as metadata_file:
        #     metadata_writer = csv.writer(
        #         metadata_file, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
        #     )
        #     metadata_writer.writerow(["sample", "genus", "species"])
        #     metadata_writer.writerow(["sample1", "Mycobacterium", "tuberculosis"])
        #     metadata_writer.writerow(["sample2", "Streptococcus", "pyogenes"])
        #     metadata_writer.writerow(["1234", "Mycobacterium", "tuberculosis"])

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
                "-maf",
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
            "minimum_allele_frequency": 0.5,
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


@unittest.skipIf(
    not Path("/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly").exists(),
    "Skipped in non-RIVM environments (because test data is needed)",
)
class TestJunoAssemblyPipeline(unittest.TestCase):
    """Testing the junoassembly class (code specific for this pipeline)"""

    @classmethod
    def setUpClass(cls) -> None:
        os.system("rm -rf test_output")

    @classmethod
    def tearDownClass(cls) -> None:
        os.system("rm -rf test_output")
        os.system("rm -rf test_output_sing")
        os.system("rm -rf test_output_sing_prefix")
        os.system("rm -rf sing_containers")

    def test_junoassembly_run_wMetadata_in_conda(self) -> None:
        """
        Testing the pipeline runs properly with real samples when providing
        a metadata file
        """
        output_dir = Path("test_output")
        input_dir = "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly"
        metadata_file = (
            "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/metadata.csv"
        )
        pipeline = JunoMapping(
            argv=[
                "-i",
                input_dir,
                "-o",
                str(output_dir),
                "-m",
                metadata_file,
                "--no-containers",
            ]
        )

        pipeline.run()
        expected_sample_sheet = {
            "sample1": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R2_001.fastq.gz",
                "genus": "salmonella",
            },
            "sample2": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R2.fastq.gz",
                "genus": "escherichia",
            },
            "sample3": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R2_001.fastq.gz",
                "genus": "streptococcus",
            },
            "sample4": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R2.fastq.gz",
                "genus": "shigella",
            },
        }

        self.assertDictEqual(
            pipeline.sample_dict,
            expected_sample_sheet,
            pipeline.sample_dict,
        )
        self.assertTrue(output_dir.joinpath("multiqc", "multiqc.html").exists())
        self.assertTrue(
            output_dir.joinpath(
                "identify_species", "top1_species_multireport.csv"
            ).exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_git.yaml").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "log_pipeline.yaml").exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_conda.txt").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "snakemake_report.html").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "sample_sheet.yaml").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "user_parameters.yaml").exists()
        )

    def test_junoassembly_run_in_singularity(self) -> None:
        """Testing the pipeline runs properly with real samples when providing
        a metadata file
        """
        output_dir = Path("test_output_sing")
        input_dir = "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly"

        pipeline = JunoMapping(
            argv=[
                "-i",
                input_dir,
                "-o",
                str(output_dir),
            ]
        )
        pipeline.run()

        expected_sample_sheet = {
            "sample1": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R2_001.fastq.gz",
                "genus": None,
            },
            "sample2": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R2.fastq.gz",
                "genus": None,
            },
            "sample3": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R2_001.fastq.gz",
                "genus": None,
            },
            "sample4": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R2.fastq.gz",
                "genus": None,
            },
        }

        self.assertEqual(
            pipeline.sample_dict,
            expected_sample_sheet,
            pipeline.sample_dict,
        )
        self.assertTrue(output_dir.joinpath("multiqc", "multiqc.html").exists())
        self.assertTrue(
            output_dir.joinpath(
                "identify_species", "top1_species_multireport.csv"
            ).exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_git.yaml").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "log_pipeline.yaml").exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_conda.txt").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "snakemake_report.html").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "sample_sheet.yaml").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "user_parameters.yaml").exists()
        )

    def test_junoassembly_wsingularity_prefix(self) -> None:
        """Testing the pipeline runs properly with real samples when providing
        a metadata file
        """
        output_dir = Path("test_output_sing_prefix")
        input_dir = "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly"

        pipeline = JunoMapping(
            argv=[
                "-i",
                input_dir,
                "-o",
                str(output_dir),
                "--prefix",
                "sing_containers",
            ]
        )
        pipeline.run()

        expected_sample_sheet = {
            "sample1": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample1_S14_R2_001.fastq.gz",
                "genus": None,
            },
            "sample2": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample2_R2.fastq.gz",
                "genus": None,
            },
            "sample3": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R1_001.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample3_R2_001.fastq.gz",
                "genus": None,
            },
            "sample4": {
                "R1": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R1.fastq.gz",
                "R2": "/data/BioGrid/hernanda/test_data_per_pipeline/Juno_assembly/sample4_R2.fastq.gz",
                "genus": None,
            },
        }

        self.assertEqual(
            pipeline.sample_dict,
            expected_sample_sheet,
            pipeline.sample_dict,
        )
        self.assertTrue(output_dir.joinpath("multiqc", "multiqc.html").exists())
        self.assertTrue(
            output_dir.joinpath(
                "identify_species", "top1_species_multireport.csv"
            ).exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_git.yaml").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "log_pipeline.yaml").exists()
        )
        self.assertTrue(output_dir.joinpath("audit_trail", "log_conda.txt").exists())
        self.assertTrue(
            output_dir.joinpath("audit_trail", "snakemake_report.html").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "sample_sheet.yaml").exists()
        )
        self.assertTrue(
            output_dir.joinpath("audit_trail", "user_parameters.yaml").exists()
        )

class test_mutation_calls(unittest.TestCase):

    def test_mutations(self):
        reader = vcf.Reader(fh('AL123456.3_mutated_1_depth30_AF0.5.vcf'))

        for var in reader:
            is_snp = var.is_snp
            if var.POS == 820:
                self.assertEqual(True, is_snp)
            if var.POS == 5659:
                self.assertEqual(False, is_snp)
            if var.POS == 15162:
                self.assertEqual(True, is_snp)
            if var.POS == 46679:
                self.assertEqual(True, is_snp)

if __name__ == "__main__":
    unittest.main()
