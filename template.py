"""
Juno template
Authors: Karim Hajji, Roxanne Wolthuis
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 05-04-2023   
"""

from pathlib import Path
import pathlib
import yaml
import argparse
import sys
from dataclasses import dataclass, field
from juno_library import Pipeline
from typing import Optional
from version import __package_name__, __version__, __description__

def main() -> None:
    juno_template = JunoTemplate()
    juno_template.run()

@dataclass
class JunoTemplate(Pipeline):
    pipeline_name: str = __package_name__
    pipeline_version: str = __version__
    input_type: str = "fastq"

    def _add_args_to_parser(self) -> None:
        super()._add_args_to_parser()

        self.parser.description = "Template juno pipeline. If you see this message please change it to something appropriate"
        
        self.add_argument(
            "--example-option",
            dest="example",
            type=str,
            required=False,
            metavar="STR",
            help="This is an optional argument, specific for this pipeline. General options are included in juno-library.",
        )
        
    def _parse_args(self) -> argparse.Namespace:
        args = super()._parse_args()

        # Optional arguments are loaded into self here
        self.example: bool = args.example

        return args
    
    # Extra class methods for this pipeline can be defined here
    def example_class_method(self):
        print(f"example option is set to {self.example}")

    def setup(self) -> None:
        super().setup()

        if self.snakemake_args["use_singularity"]:
            self.snakemake_args["singularity_args"] = " ".join(
                [
                    self.snakemake_args["singularity_args"]
                ] # paths that singularity should be able to read from can be bound by adding to the above list
            )

        # Extra class methods for this pipeline can be invoked here
        if self.example:
            self.example_class_method()

        with open(
            Path(__file__).parent.joinpath("config/pipeline_parameters.yaml")
        ) as f:
            parameters_dict = yaml.safe_load(f)
        self.snakemake_config.update(parameters_dict)

        self.user_parameters = {
            "input_dir": str(self.input_dir),
            "output_dir": str(self.output_dir),
            "exclusion_file": str(self.exclusion_file),
            "example": str(self.example), # other user parameters can be included in user_parameters.yaml here
        }


if __name__ == "__main__":
    main()
