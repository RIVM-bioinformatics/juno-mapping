"""
Juno mapping
Authors: Boas van der Putten
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 25-05-2023   
"""

import argparse
import pathlib
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional, Union

import yaml
from juno_library import Pipeline

from version import __description__, __package_name__, __version__


def main() -> None:
    juno_mapping = JunoMapping()
    juno_mapping.run()


def check_number_within_range(
    minimum: float = 0, maximum: float = 1
) -> Union[Callable[[str], str], argparse.FileType]:
    """
    Creates a function to check whether a numeric value is within a range, inclusive.

    The generated function can be used by the `type` parameter in argparse.ArgumentParser.
    See https://stackoverflow.com/a/53446178.

    Args:
        value: the numeric value to check.
        minimum: minimum of allowed range, inclusive.
        maximum: maximum of allowed range, inclusive.

    Returns:
        A function which takes a single argument and checks this against the range.

    Raises:
        argparse.ArgumentTypeError: if the value is outside the range.
        ValueError: if the value cannot be converted to float.
    """

    def generated_func_check_range(value: str) -> str:
        value_f = float(value)
        if (value_f < minimum) or (value_f > maximum):
            raise argparse.ArgumentTypeError(
                f"Supplied value {value} is not within expected range {minimum} to {maximum}."
            )
        return str(value)

    return generated_func_check_range


@dataclass
class JunoMapping(Pipeline):
    pipeline_name: str = __package_name__
    pipeline_version: str = __version__
    input_type: str = "fastq"
    species_options = ["mycobacterium_tuberculosis"]

    def _add_args_to_parser(self) -> None:
        super()._add_args_to_parser()

        self.parser.description = (
            "Juno-mapping pipeline for reference mapping of bacterial genomes"
        )

        self.add_argument(
            "-s",
            "--species",
            dest="species",
            type=str.lower,
            required=True,
            metavar="STR",
            help=f"Species to use, choose from: {self.species_options}",
            choices=self.species_options,
        )
        self.add_argument(
            "--reference",
            type=Path,
            metavar="FILE",
            dest="custom_reference",
            help="Reference genome to use. Default is chosen based on species argument, defaults per species can be found in: /mnt/db/juno/mapping/[species]",
            required=False,
        )
        self.add_argument(
            "--mask",
            type=Path,
            metavar="FILE",
            dest="custom_mask",
            help="Mask file to use, defaults per species can be found in: /mnt/db/juno/mapping/[species]",
            required=False,
        )
        self.add_argument(
            "--disable-mask",
            action="store_true",
            dest="disable_mask",
            help="Disable masking, use at your own risk: this might cause the appearance of low quality variants in the final VCF file.",
            required=False,
        )
        self.add_argument(
            "--db-dir",
            type=Path,
            default="/mnt/db/juno/kraken2_db",
            metavar="DIR",
            help="Kraken2 database directory.",
        )
        self.add_argument(
            "-mpt",
            "--mean-quality-threshold",
            type=check_number_within_range(minimum=1, maximum=36),
            metavar="INT",
            default=28,
            help="Phred score to be used as threshold for cleaning (filtering) fastq files.",
        )
        self.add_argument(
            "-ws",
            "--window-size",
            type=check_number_within_range(minimum=1, maximum=1000),
            metavar="INT",
            default=5,
            help="Window size to use for cleaning (filtering) fastq files.",
        )
        self.add_argument(
            "-ml",
            "--minimum-length",
            type=check_number_within_range(minimum=0, maximum=500),
            metavar="INT",
            default=50,
            help="Minimum length for fastq reads to be kept after trimming.",
        )
        self.add_argument(
            "-md",
            "--minimum-depth-variant",
            type=check_number_within_range(minimum=0, maximum=9999),
            metavar="INT",
            default=10,
            help="Minimum length for fastq reads to be kept after trimming.",
        )
        self.add_argument(
            "-maf",
            "--minimum-allele-frequency",
            type=check_number_within_range(minimum=0, maximum=1),
            metavar="FLOAT",
            default=0.8,
            help="Minimum allele frequency to filter variants on.",
        )

    def _parse_args(self) -> argparse.Namespace:
        args = super()._parse_args()

        # Optional arguments are loaded into self here
        self.db_dir: Path = args.db_dir
        self.mean_quality_threshold: int = args.mean_quality_threshold
        self.window_size: int = args.window_size
        self.min_read_length: int = args.minimum_length
        self.reference: Optional[Path] = None
        self.custom_reference: Path = args.custom_reference
        self.mask: Optional[Path] = None
        self.custom_mask: Path = args.custom_mask
        self.disable_mask: bool = args.disable_mask
        self.time_limit: int = args.time_limit
        self.species: str = args.species
        self.minimum_depth_variant: int = args.minimum_depth_variant
        self.minimum_allele_frequency: float = args.minimum_allele_frequency

        return args

    # # Extra class methods for this pipeline can be defined here
    # def example_class_method(self):
    #     print(f"example option is set to {self.example}")

    def setup(self) -> None:
        super().setup()

        if self.snakemake_args["use_singularity"]:
            self.snakemake_args["singularity_args"] = " ".join(
                [
                    self.snakemake_args["singularity_args"],
                    f"--bind {self.db_dir}:{self.db_dir}",
                ]  # paths that singularity should be able to read from can be bound by adding to the above list
            )

        # # Extra class methods for this pipeline can be invoked here
        # if self.example:
        #     self.example_class_method()

        # select a reference based on species:
        # self.ref_dir = "/mnt/db/apollo/mapping/"
        # replace by dictionary if expansion is needed, could also support acronyms (mtb, tb)
        if self.species == "mycobacterium_tuberculosis":
            self.reference = Path(
                "/mnt/db/juno/mapping/mycobacterium_tuberculosis/AL123456.3.fasta"
            )
            self.mask = Path("files/RLC_Farhat.bed")
        else:
            self.reference = None
            self.mask = None
        # elif self.species == "speciesX":
        #     self.reference = "/mnt/db/juno/mapping/speciesX/awesome_assembly.fasta"

        if self.custom_reference is not None:
            print(
                "A reference genome was specified by the user, which may not be the default reference genome for this species."
            )
            self.reference = self.custom_reference

        if self.custom_mask is not None:
            print(
                "A mask file was specified by the user, which may not be the default mask file for this species."
            )
            self.mask = self.custom_mask

        if self.disable_mask:
            print(
                "Masking was disabled by the user, which may not be the default for this species."
            )
            self.mask = None

        print(f"Running pipeline for {self.species} with reference: {self.reference}.")
        if self.mask is None:
            print(f"Masking will not be performed.")
        else:
            print(f"Masking will be performed using {self.mask}.")

        with open(
            Path(__file__).parent.joinpath("config/pipeline_parameters.yaml")
        ) as f:
            parameters_dict = yaml.safe_load(f)
        self.snakemake_config.update(parameters_dict)

        self.user_parameters = {
            "input_dir": str(self.input_dir),
            "output_dir": str(self.output_dir),
            "exclusion_file": str(self.exclusion_file),
            "db_dir": str(self.db_dir),
            "mean_quality_threshold": int(self.mean_quality_threshold),
            "window_size": int(self.window_size),
            "min_read_length": int(self.min_read_length),
            "reference": str(self.reference),
            "mask_bed": str(self.mask),
            "disable_mask": str(self.disable_mask),
            "use_singularity": str(self.snakemake_args["use_singularity"]),
            "species": str(self.species),
            "minimum_depth": int(self.minimum_depth_variant),
            "minimum_allele_frequency": float(self.minimum_allele_frequency),
        }


if __name__ == "__main__":
    main()
