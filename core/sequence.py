"""Shared DNA, RNA, and protein sequence helpers."""

from __future__ import annotations

import re

GENETIC_CODE = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M", "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T", "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*", "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K", "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W", "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R", "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}

AMBIGUOUS_MAP = {"N":"A", "R":"A", "Y":"C", "S":"G", "W":"A", "K":"G", "M":"A", "B":"C", "D":"A", "H":"A", "V":"A"}
DNA_BASES = frozenset("ACGT")


def clean_sequence(seq: str, strip_ambiguous: bool = True) -> str:
    """Normalize a DNA sequence, removing FASTA whitespace and numbering."""
    normalized = re.sub(r"[\s\d]", "", seq.upper())
    normalized = re.sub(r"^>[^\n]*", "", normalized)
    if strip_ambiguous:
        return "".join(base if base in DNA_BASES else AMBIGUOUS_MAP.get(base, "") for base in normalized)
    invalid = set(normalized) - DNA_BASES
    if invalid:
        raise ValueError(f"Invalid DNA characters: {sorted(invalid)}")
    return normalized


def validate_sequence(seq: str, minimum_length: int = 10) -> str:
    """Normalize and validate a DNA sequence for downstream analysis."""
    cleaned = clean_sequence(seq)
    if len(cleaned) < minimum_length:
        raise ValueError(f"Sequence too short (min {minimum_length} bp).")
    return cleaned


def complement(base: str) -> str:
    return {"A":"T", "T":"A", "G":"C", "C":"G"}.get(base.upper(), "N")


def reverse_complement(seq: str) -> str:
    return "".join(complement(base) for base in reversed(seq.upper()))


def gc_content(seq: str) -> float:
    seq = seq.upper()
    return 0.0 if not seq else (seq.count("G") + seq.count("C")) / len(seq) * 100


def translate(seq: str) -> str:
    """Translate a DNA sequence in reading frame zero."""
    protein = []
    seq = seq.upper()
    for index in range(0, len(seq) - 2, 3):
        amino_acid = GENETIC_CODE.get(seq[index:index + 3], "?")
        protein.append(amino_acid)
        if amino_acid == "*":
            break
    return "".join(protein)


def dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def rna_complement(base: str) -> str:
    return {"A":"U", "U":"A", "G":"C", "C":"G"}.get(base.upper(), "N")


__all__ = ["AMBIGUOUS_MAP", "DNA_BASES", "GENETIC_CODE", "clean_sequence", "complement", "dna_to_rna", "gc_content", "reverse_complement", "rna_complement", "translate", "validate_sequence"]
