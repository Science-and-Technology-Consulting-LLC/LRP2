"""paths.py

Where a run's results are, set once. The same `paths.json` as paths.R reads:
the run folder, the comparison name, the dataset name, and the folders inside the run
as short relative names. To point the notebooks at another run, edit `paths.json` only.
The environment variable LRP2_PATHS can name a different config file.
"""

import json
import os


def run_locations(config=None):
    """Return a dict: comparison, dataset_name, and the full path of every folder in `folders`."""
    config = config or os.environ.get("LRP2_PATHS", "paths.json")
    with open(config) as fh:
        cfg = json.load(fh)

    def fill(rel):
        return os.path.join(cfg["run_dir"], rel.replace("{comparison}", cfg["comparison"]))

    loc = {"run_dir": cfg["run_dir"], "comparison": cfg["comparison"], "dataset_name": cfg["dataset_name"]}
    loc.update({name: fill(rel) for name, rel in cfg["folders"].items()})
    return loc
