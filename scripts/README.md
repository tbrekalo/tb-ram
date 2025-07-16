# Directory containing scripts for evaluating ram performance

## Evaluate reference mapping

```bash
uv run eval.py -h

usage: eval [-h] [-n NAME] [-m MIN_RATIO] [-o OUT] overlaps origins

Evaluate mapper precision for mapping simulated reads onto the reference

positional arguments:
  overlaps
  origins

options:
  -h, --help            show this help message and exit
  -n, --name NAME       Evaluation name.
  -m, --min_ratio MIN_RATIO
                        How much read mapping has to overlap with the origin
                        on the reference to be considered a true positive.
  -o, --out OUT         Csv output file.

```
