"""Console script for validate_meta."""
import argparse
import sys
from pathlib import Path
import pandas as pd
import yaml
from validate_meta.classes import DataFrameValidator


def main():
    """Console script for validate_meta."""
    parser = argparse.ArgumentParser(description="Validates a GMS-mikro metadata file")
    parser.add_argument("-i", dest="filename", required=True,
                        help="input file with two matrices",
                        metavar="CSV DATA FILE")

    parser.add_argument("-s", dest="definition", required=True,
                        help="yaml definition file",
                        metavar="DEFINITION YAML FILE")

    args = parser.parse_args()

    datafile = Path(args.filename)
    definitionfile = Path(args.definition)

    if datafile.exists() and definitionfile.exists():
        df = pd.read_csv(datafile)
        with definitionfile.open(encoding='utf8') as fp:
            definition = yaml.safe_load(fp)

    v = DataFrameValidator(df, definition)
    errors = v.validate()
    for e in errors:
        print(e)

    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
