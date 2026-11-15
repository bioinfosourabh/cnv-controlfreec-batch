import importlib.util
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="session")
def mfc():
    return _load("make_freec_configs")


@pytest.fixture(scope="session")
def summary():
    return _load("freec_cohort_summary")


@pytest.fixture
def template_path():
    return ROOT / "assets" / "config_template.txt"


@pytest.fixture
def sheet(tmp_path):
    def _write(text):
        p = tmp_path / "samples.csv"
        p.write_text(text)
        return p
    return _write
