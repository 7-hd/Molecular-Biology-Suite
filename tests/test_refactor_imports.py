"""Regression tests for the refactored package surface."""

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
