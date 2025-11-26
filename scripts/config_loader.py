import yaml
from pathlib import Path

def load_config(path="config/config.yml"):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def load_species_params(path="config/species_params.yml"):
    with open(path, "r") as f:
        return yaml.safe_load(f)
