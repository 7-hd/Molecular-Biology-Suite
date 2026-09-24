"""Shared sequence helpers and validation for the bioinformatics suite."""

from __future__ import annotations

import re

GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AMBIGUOUS_MAP = {
    "N": "A", "R": "A", "Y": "C", "S": "G", "W": "A",
    "K": "G", "M": "A", "B": "C", "D": "A", "H": "A", "V": "A",
}


def clean_sequence(seq: str, strip_ambiguous: bool = True) -> str:
    """Strip whitespace, FASTA headers, and non-DNA characters where needed."""
    seq = re.sub(r"[\s\d>]", "", seq).upper()
    seq = re.sub(r">.*", "", seq)

    if strip_ambiguous:
        cleaned = []
        for base in seq:
            if base in "ACGT":
                cleaned.append(base)
            elif base in AMBIGUOUS_MAP:
                cleaned.append(AMBIGUOUS_MAP[base])
        return "".join(cleaned)

    invalid = set(seq) - set("ACGT")
    if invalid:
        raise ValueError(f"Invalid characters: {sorted(invalid)}")
    return seq


def validate_sequence(seq: str, minimum_length: int = 10) -> str:
    """Ensure a DNA sequence is usable for downstream analysis."""
    cleaned = clean_sequence(seq)
    if len(cleaned) < minimum_length:
        raise ValueError(f"Sequence too short (min {minimum_length} bp).")
    return cleaned


def complement(base: str) -> str:
    """Return the complementary DNA base."""
    return {"A": "T", "T": "A", "G": "C", "C": "G"}.get(base.upper(), "N")


def reverse_complement(seq: str) -> str:
    """Compute the reverse complement of a DNA sequence."""
    return "".join(complement(base) for base in reversed(seq.upper()))


def gc_content(seq: str) -> float:
    """Return the percent GC content."""
    seq = seq.upper()
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def translate(seq: str) -> str:
    """Translate a DNA sequence to a protein sequence."""
    seq = seq.upper()
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = GENETIC_CODE.get(codon, "?")
        protein.append(aa)
        if aa == "*":
            break
    return "".join(protein)


def dna_to_rna(seq: str) -> str:
    """Convert DNA to RNA."""
    return seq.upper().replace("T", "U")


def rna_complement(base: str) -> str:
    """Return the complementary RNA base."""
    return {"A": "U", "U": "A", "G": "C", "C": "G"}.get(base.upper(), "N")
