import pandas as pd

def parse_wgs_metrics(wgs_metrics_file, output_mean_coverage):
    """
    Parse the WgsMetrics file and return the mean coverage.
    """
    # Read the WgsMetrics file
    df = pd.read_csv(wgs_metrics_file, sep='\t', comment='#', skiprows=6, nrows=1)

    # Extract the mean coverage
    mean_coverage = df['MEAN_COVERAGE'].values[0]

    # Write the mean coverage to the output file
    with open(output_mean_coverage, "w") as f:
        f.write(str(mean_coverage))


# Extract file paths from snakemake input and output
wgs_metrics_file = str(snakemake.input[0])  # type: ignore[name-defined]
output_mean_coverage = str(snakemake.output[0])  # type: ignore[name-defined]

# Call the function with the extracted file paths
parse_wgs_metrics(wgs_metrics_file, output_mean_coverage)

# import pandas as pd
# import sys

# input_wgs_metrics = sys.argv[1]
# output_mean_coverage = sys.argv[2]

# # Read the WgsMetrics file
# df = pd.read_csv(input_wgs_metrics, sep='\t', comment='#', skiprows=6, nrows=1)

# # Extract the mean coverage
# mean_coverage = df['MEAN_COVERAGE'].values[0]

# # Write the mean coverage to the output file
# with open(output_mean_coverage, "w") as f:
#     f.write(str(mean_coverage))