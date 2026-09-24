"""Regression tests for the refactored public package surface."""

from app.cli import main
from biology import PrimerDesigner, RestrictionMapper
from core.sequence import clean_sequence, gc_content, reverse_complement, translate, validate_sequence


def test_core_sequence_helpers():
    seq = "ATGCGTACGTA"
    assert clean_sequence(seq) == seq
    assert reverse_complement(seq) == "TACGTACGCAT"
    assert gc_content(seq) == 50.0
    assert translate(seq) == "MRT"
    assert validate_sequence(seq) == seq


def test_public_package_imports():
    assert callable(main)
    assert PrimerDesigner is not None
    assert RestrictionMapper("GAATTCAAA").find_sites(["EcoRI"])["EcoRI"]["positions"] == [1]


def test_primer_designer_rejects_invalid_bounds():
    try:
        PrimerDesigner("A" * 40, min_len=25, max_len=18)
    except ValueError:
        return
    raise AssertionError("invalid primer bounds should raise ValueError")
