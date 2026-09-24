from core.sequence import clean_sequence, gc_content, reverse_complement, translate, validate_sequence
from app.cli import main


def test_core_sequence_helpers():
    seq = "ATGCGTACGTA"
    assert clean_sequence(seq) == seq
    assert reverse_complement(seq) == "TACGTACGCAT"
    assert gc_content(seq) == 50.0
    assert translate(seq) == "MRT"
    assert validate_sequence(seq) == seq


def test_cli_entrypoint_is_defined():
    assert callable(main)
