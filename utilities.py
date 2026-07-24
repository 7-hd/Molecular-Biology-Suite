# ================================================================
#  DBS  –  Digital BioScience Suite
#  utilities.py  |  Shared engines, parsers, thermodynamics
#
#  Supports all modules across Phase 1, 2 and 3
# ================================================================

import os, re, math, gzip

# ──────────────────────────────────────────────────────────────────
#  SantaLucia (1998) nearest-neighbor DNA parameters
#  (ΔH kcal/mol , ΔS cal/mol·K)
# ──────────────────────────────────────────────────────────────────
NN_PARAMS = {
    "AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "GT": (-8.4, -22.4), "CT": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "CG": (-10.6,-27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "AC": (-7.8, -21.0), "TC": (-7.8, -21.0),
    "TG": (-8.5, -22.7), "AG": (-7.8, -21.0), "TT": (-7.9, -22.2),
    "CC": (-8.0, -19.9),
}
INIT_PARAMS = {"dH": 0.2, "dS": -5.7}

# ──────────────────────────────────────────────────────────────────
#  RNA Turner energy parameters (simplified, 2004 rules)
#  Used for RNA secondary structure prediction
#  Stacking pairs: (5'X-Y3') energy in kcal/mol
# ──────────────────────────────────────────────────────────────────
RNA_STACK = {
    ("A","U"): -0.9, ("U","A"): -0.9, ("G","C"): -2.9, ("C","G"): -2.9,
    ("G","U"): -0.5, ("U","G"): -0.5, ("A","A"): 0.0,  ("U","U"): 0.0,
    ("G","A"): -1.1, ("A","G"): -1.1, ("C","U"): -1.8, ("U","C"): -1.8,
    ("G","G"): -1.6, ("C","C"): -1.6, ("A","C"): -0.8, ("C","A"): -0.8,
}
RNA_HAIRPIN_INIT = {3: 5.4, 4: 4.1, 5: 4.3, 6: 4.1, 7: 4.1}  # kcal/mol

# ──────────────────────────────────────────────────────────────────
#  IUPAC ambiguous base mapping
# ──────────────────────────────────────────────────────────────────
AMBIGUOUS_MAP = {
    "N":"A","R":"A","Y":"C","S":"G","W":"A",
    "K":"G","M":"A","B":"C","D":"A","H":"A","V":"A",
}

# ──────────────────────────────────────────────────────────────────
#  Standard genetic code (DNA triplet → amino acid, 1-letter)
# ──────────────────────────────────────────────────────────────────
GENETIC_CODE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CAT":"H","CAC":"H","AAT":"N","AAC":"N",
    "AAA":"K","AAG":"K","GAT":"D","GAC":"D",
    "GAA":"E","GAG":"E","TGT":"C","TGC":"C",
    "TGA":"*","TGG":"W","CGT":"R","CGC":"R",
    "CGA":"R","CGG":"R","AGT":"S","AGC":"S",
    "AGA":"R","AGG":"R","GGT":"G","GGC":"G",
    "GGA":"G","GGG":"G",
}

AMINO_ACIDS = {
    "A":"Ala","C":"Cys","D":"Asp","E":"Glu","F":"Phe",
    "G":"Gly","H":"His","I":"Ile","K":"Lys","L":"Leu",
    "M":"Met","N":"Asn","P":"Pro","Q":"Gln","R":"Arg",
    "S":"Ser","T":"Thr","V":"Val","W":"Trp","Y":"Tyr",
    "*":"Stop",
}

# ──────────────────────────────────────────────────────────────────
#  Codon usage tables (fraction of synonymous codons used)
#  Source: Kazusa Codon Usage Database
#  Format: {organism: {codon: frequency}}
# ──────────────────────────────────────────────────────────────────
CODON_TABLES = {
    "E. coli K-12": {
        "TTT":0.58,"TTC":0.42,"TTA":0.14,"TTG":0.13,"CTT":0.12,
        "CTC":0.10,"CTA":0.04,"CTG":0.47,"ATT":0.51,"ATC":0.39,
        "ATA":0.11,"ATG":1.00,"GTT":0.28,"GTC":0.20,"GTA":0.17,
        "GTG":0.35,"TCT":0.17,"TCC":0.15,"TCA":0.14,"TCG":0.14,
        "CCT":0.18,"CCC":0.13,"CCA":0.20,"CCG":0.49,"ACT":0.19,
        "ACC":0.40,"ACA":0.17,"ACG":0.25,"GCT":0.18,"GCC":0.26,
        "GCA":0.23,"GCG":0.33,"TAT":0.59,"TAC":0.41,"TAA":0.61,
        "TAG":0.09,"TGA":0.30,"CAT":0.57,"CAC":0.43,"CAA":0.34,
        "CAG":0.66,"AAT":0.49,"AAC":0.51,"AAA":0.74,"AAG":0.26,
        "GAT":0.63,"GAC":0.37,"GAA":0.68,"GAG":0.32,"TGT":0.46,
        "TGC":0.54,"TGG":1.00,"CGT":0.36,"CGC":0.36,"CGA":0.07,
        "CGG":0.11,"AGT":0.16,"AGC":0.25,"AGA":0.07,"AGG":0.04,
        "GGT":0.35,"GGC":0.37,"GGA":0.13,"GGG":0.15,
    },
    "S. cerevisiae": {
        "TTT":0.59,"TTC":0.41,"TTA":0.28,"TTG":0.29,"CTT":0.13,
        "CTC":0.06,"CTA":0.14,"CTG":0.11,"ATT":0.46,"ATC":0.26,
        "ATA":0.27,"ATG":1.00,"GTT":0.39,"GTC":0.21,"GTA":0.21,
        "GTG":0.19,"TCT":0.26,"TCC":0.16,"TCA":0.21,"TCG":0.10,
        "CCT":0.31,"CCC":0.15,"CCA":0.42,"CCG":0.12,"ACT":0.35,
        "ACC":0.22,"ACA":0.30,"ACG":0.14,"GCT":0.38,"GCC":0.22,
        "GCA":0.29,"GCG":0.11,"TAT":0.56,"TAC":0.44,"TAA":0.47,
        "TAG":0.22,"TGA":0.31,"CAT":0.64,"CAC":0.36,"CAA":0.69,
        "CAG":0.31,"AAT":0.59,"AAC":0.41,"AAA":0.58,"AAG":0.42,
        "GAT":0.65,"GAC":0.35,"GAA":0.70,"GAG":0.30,"TGT":0.63,
        "TGC":0.37,"TGG":1.00,"CGT":0.14,"CGC":0.06,"CGA":0.07,
        "CGG":0.04,"AGT":0.16,"AGC":0.11,"AGA":0.48,"AGG":0.21,
        "GGT":0.47,"GGC":0.19,"GGA":0.22,"GGG":0.12,
    },
    "H. sapiens": {
        "TTT":0.45,"TTC":0.55,"TTA":0.07,"TTG":0.13,"CTT":0.13,
        "CTC":0.20,"CTA":0.07,"CTG":0.41,"ATT":0.36,"ATC":0.47,
        "ATA":0.17,"ATG":1.00,"GTT":0.18,"GTC":0.24,"GTA":0.11,
        "GTG":0.47,"TCT":0.15,"TCC":0.22,"TCA":0.15,"TCG":0.06,
        "CCT":0.28,"CCC":0.33,"CCA":0.27,"CCG":0.11,"ACT":0.25,
        "ACC":0.36,"ACA":0.28,"ACG":0.11,"GCT":0.26,"GCC":0.40,
        "GCA":0.23,"GCG":0.11,"TAT":0.43,"TAC":0.57,"TAA":0.28,
        "TAG":0.20,"TGA":0.52,"CAT":0.41,"CAC":0.59,"CAA":0.27,
        "CAG":0.73,"AAT":0.46,"AAC":0.54,"AAA":0.42,"AAG":0.58,
        "GAT":0.46,"GAC":0.54,"GAA":0.42,"GAG":0.58,"TGT":0.45,
        "TGC":0.55,"TGG":1.00,"CGT":0.08,"CGC":0.19,"CGA":0.11,
        "CGG":0.21,"AGT":0.15,"AGC":0.24,"AGA":0.20,"AGG":0.20,
        "GGT":0.16,"GGC":0.34,"GGA":0.25,"GGG":0.25,
    },
    "B. subtilis": {
        "TTT":0.57,"TTC":0.43,"TTA":0.17,"TTG":0.19,"CTT":0.14,
        "CTC":0.07,"CTA":0.06,"CTG":0.15,"ATT":0.52,"ATC":0.32,
        "ATA":0.16,"ATG":1.00,"GTT":0.37,"GTC":0.17,"GTA":0.20,
        "GTG":0.27,"TCT":0.15,"TCC":0.14,"TCA":0.17,"TCG":0.11,
        "CCT":0.15,"CCC":0.08,"CCA":0.24,"CCG":0.19,"ACT":0.24,
        "ACC":0.29,"ACA":0.25,"ACG":0.21,"GCT":0.25,"GCC":0.21,
        "GCA":0.27,"GCG":0.26,"TAT":0.57,"TAC":0.43,"TAA":0.56,
        "TAG":0.08,"TGA":0.37,"CAT":0.60,"CAC":0.40,"CAA":0.50,
        "CAG":0.50,"AAT":0.57,"AAC":0.43,"AAA":0.64,"AAG":0.36,
        "GAT":0.64,"GAC":0.36,"GAA":0.65,"GAG":0.35,"TGT":0.46,
        "TGC":0.54,"TGG":1.00,"CGT":0.25,"CGC":0.18,"CGA":0.10,
        "CGG":0.07,"AGT":0.19,"AGC":0.23,"AGA":0.14,"AGG":0.08,
        "GGT":0.33,"GGC":0.28,"GGA":0.24,"GGG":0.16,
    },
    "P. aeruginosa": {
        "TTT":0.36,"TTC":0.64,"TTA":0.04,"TTG":0.12,"CTT":0.10,
        "CTC":0.22,"CTA":0.04,"CTG":0.48,"ATT":0.32,"ATC":0.55,
        "ATA":0.13,"ATG":1.00,"GTT":0.14,"GTC":0.22,"GTA":0.14,
        "GTG":0.50,"TCT":0.09,"TCC":0.18,"TCA":0.08,"TCG":0.16,
        "CCT":0.09,"CCC":0.20,"CCA":0.14,"CCG":0.57,"ACT":0.10,
        "ACC":0.44,"ACA":0.10,"ACG":0.36,"GCT":0.11,"GCC":0.37,
        "GCA":0.12,"GCG":0.40,"TAT":0.34,"TAC":0.66,"TAA":0.22,
        "TAG":0.14,"TGA":0.65,"CAT":0.37,"CAC":0.63,"CAA":0.23,
        "CAG":0.77,"AAT":0.31,"AAC":0.69,"AAA":0.37,"AAG":0.63,
        "GAT":0.37,"GAC":0.63,"GAA":0.37,"GAG":0.63,"TGT":0.27,
        "TGC":0.73,"TGG":1.00,"CGT":0.19,"CGC":0.35,"CGA":0.06,
        "CGG":0.21,"AGT":0.10,"AGC":0.29,"AGA":0.06,"AGG":0.13,
        "GGT":0.15,"GGC":0.44,"GGA":0.16,"GGG":0.25,
    },
}

# ──────────────────────────────────────────────────────────────────
#  Restriction enzyme database (name → recognition sequence)
#  Sourced from REBASE (Roberts et al.)
# ──────────────────────────────────────────────────────────────────
RESTRICTION_ENZYMES = {
    # Class II — most common in cloning
    "EcoRI":  ("GAATTC",  1),   # cut after pos 1 on top strand
    "BamHI":  ("GGATCC",  1),
    "HindIII":("AAGCTT",  1),
    "NcoI":   ("CCATGG",  1),
    "NdeI":   ("CATATG",  2),
    "XhoI":   ("CTCGAG",  1),
    "XbaI":   ("TCTAGA",  1),
    "SalI":   ("GTCGAC",  1),
    "SphI":   ("GCATGC",  1),
    "KpnI":   ("GGTACC",  1),
    "SacI":   ("GAGCTC",  1),
    "AvaI":   ("CYCGRG",  1),
    "BglII":  ("AGATCT",  1),
    "ClaI":   ("ATCGAT",  2),
    "EcoRV":  ("GATATC",  3),   # blunt cutter
    "MluI":   ("ACGCGT",  1),
    "NheI":   ("GCTAGC",  1),
    "NotI":   ("GCGGCCGC",2),
    "PstI":   ("CTGCAG",  4),
    "PvuII":  ("CAGCTG",  3),   # blunt
    "ScaI":   ("AGTACT",  3),   # blunt
    "SmaI":   ("CCCGGG",  3),   # blunt
    "SpeI":   ("ACTAGT",  1),
    "StuI":   ("AGGCCT",  3),   # blunt
    "SwaI":   ("ATTTAAAT",4),
    "AscI":   ("GGCGCGCC",2),
    "FseI":   ("GGCCGGCC",6),
    "PacI":   ("TTAATTAA",4),
    "SfiI":   ("GGCCNNNNNGGCC",8),
    # Type IIS (Golden Gate)
    "BsaI":   ("GGTCTCN",7),
    "BsmBI":  ("CGTCTCN",7),
    "BbsI":   ("GAAGACNN",8),
    "SapI":   ("GCTCTTCN",8),
}

# ──────────────────────────────────────────────────────────────────
#  Antibiotic Resistance Gene (ARG) database
#  Format: {gene_name: (drug_class, mechanism, conserved_motif)}
# ──────────────────────────────────────────────────────────────────
ARG_DATABASE = {
    # Beta-lactams
    "blaTEM":    ("Beta-lactam",      "Enzymatic hydrolysis (beta-lactamase)",
                  "MSIQHFRVALIPFFAAFCLPVFA"),
    "blaSHV":    ("Beta-lactam",      "Enzymatic hydrolysis (ESBL)",
                  "MSLLTLVVVTIVCGSGLA"),
    "blaCTX-M":  ("Beta-lactam",      "Extended-spectrum beta-lactamase",
                  "MMTAFKRLALLPLLLPLA"),
    "blaOXA":    ("Beta-lactam",      "Oxacillinase",
                  "MNNFIPLAIFCLSTLSG"),
    "blaKPC":    ("Carbapenem",       "Carbapenemase (class A)",
                  "MSLYRLLPLLSLLPQPLA"),
    "blaNDM":    ("Carbapenem",       "Metallo-beta-lactamase",
                  "MQKKALLIAALITLSGMSAHA"),
    "blaVIM":    ("Carbapenem",       "Metallo-beta-lactamase",
                  "MKHSIALLAAGLLLSAPA"),
    "blaIMP":    ("Carbapenem",       "Metallo-beta-lactamase",
                  "MLNLKTLLSLAFFCS"),
    # Aminoglycosides
    "aac6":      ("Aminoglycoside",   "Acetyltransferase (AAC)",
                  "MSGIHASELLK"),
    "aph3":      ("Aminoglycoside",   "Phosphotransferase (APH)",
                  "MSHIQRETSCS"),
    "ant2":      ("Aminoglycoside",   "Nucleotidyltransferase (ANT)",
                  "MAVQHILNLSK"),
    "rmtB":      ("Aminoglycoside",   "16S rRNA methyltransferase",
                  "MRIVTSEGPN"),
    # Tetracyclines
    "tetA":      ("Tetracycline",     "Efflux pump (MFS family)",
                  "MKSNNALIVILGTVTLDAVGIGLVMPVLPDFY"),
    "tetB":      ("Tetracycline",     "Efflux pump",
                  "MKSNNALVLWLGTVTLDTIGI"),
    "tetM":      ("Tetracycline",     "Ribosomal protection protein",
                  "MSKNKIFKEFSENR"),
    "tetX":      ("Tetracycline",     "FAD-dependent monooxygenase",
                  "MDRQFELTPNYLQ"),
    # Chloramphenicol
    "cat":       ("Chloramphenicol",  "Acetyltransferase (CAT)",
                  "MEKKITGYTTVDISQWHRKEHF"),
    "cfr":       ("Chloramphenicol",  "23S rRNA methyltransferase",
                  "MLNELLERIKKKLRE"),
    # Fluoroquinolones
    "qnrA":      ("Fluoroquinolone",  "Quinolone resistance protein",
                  "MTEQQNLLNKGKPLDGLM"),
    "qnrB":      ("Fluoroquinolone",  "Quinolone resistance protein",
                  "MSEHTQSIFLDQFLN"),
    "qnrS":      ("Fluoroquinolone",  "Quinolone resistance protein",
                  "MTEKEQQILRGKKPL"),
    "aac6-aph2": ("Aminoglycoside/Fluoroquinolone","Bifunctional enzyme",
                  "MSHIQRETSA"),
    # Macrolides
    "ermA":      ("Macrolide",        "23S rRNA methyltransferase",
                  "MNNNKNNKPTNKLK"),
    "ermB":      ("Macrolide",        "23S rRNA methyltransferase",
                  "MENIAKNINNN"),
    "ermC":      ("Macrolide",        "23S rRNA methyltransferase",
                  "MPIFTKKINDQT"),
    "mefA":      ("Macrolide",        "Efflux pump",
                  "MKKNNTLYNLGIFG"),
    # Glycopeptides
    "vanA":      ("Glycopeptide",     "D-Ala-D-Lac ligase (VanA)",
                  "MSEVQIFNENQKEQ"),
    "vanB":      ("Glycopeptide",     "D-Ala-D-Lac ligase (VanB)",
                  "MSEIQIFNED"),
    # Trimethoprim
    "dhfr":      ("Trimethoprim",     "Dihydrofolate reductase",
                  "MISLIAALAVDRVIGMENAMPWN"),
    "dfrA":      ("Trimethoprim",     "Dihydrofolate reductase variant",
                  "MNLAQLIAA"),
    # Sulfonamides
    "sul1":      ("Sulfonamide",      "Dihydropteroate synthase",
                  "MNKSLIIFAVNRTPLAK"),
    "sul2":      ("Sulfonamide",      "Dihydropteroate synthase",
                  "MNIFSLIILSNR"),
    # Colistin (last resort)
    "mcr1":      ("Colistin",         "Phosphoethanolamine transferase",
                  "MKLKSSKYLTLGLALL"),
    "mcr2":      ("Colistin",         "Phosphoethanolamine transferase",
                  "MKLKASKYLTLGLALL"),
}

# ──────────────────────────────────────────────────────────────────
#  MHC-I binding motifs (simplified anchor residues per allele)
#  Format: {allele: {pos: [preferred residues]}}
#  Based on SYFPEITHI and NetMHCpan anchor data
# ──────────────────────────────────────────────────────────────────
MHC_I_MOTIFS = {
    "HLA-A*02:01": {2: list("LMV"),  9: list("LV")},
    "HLA-A*01:01": {3: list("DE"),   9: list("Y")},
    "HLA-A*03:01": {2: list("LMV"),  9: list("KR")},
    "HLA-B*07:02": {2: list("P"),    9: list("L")},
    "HLA-B*44:02": {2: list("E"),    9: list("LY")},
}

# ──────────────────────────────────────────────────────────────────
#  Eisenberg hydrophobicity scale (amino acids)
#  Reference: Eisenberg et al., 1982
# ──────────────────────────────────────────────────────────────────
HYDROPHOBICITY = {
    "A": 0.25, "R":-1.76, "N":-0.64, "D":-0.72, "C": 0.04,
    "Q":-0.69, "E":-0.62, "G": 0.16, "H":-0.40, "I": 0.73,
    "L": 0.53, "K":-1.10, "M": 0.26, "F": 0.61, "P":-0.07,
    "S":-0.26, "T":-0.18, "W": 0.37, "Y": 0.02, "V": 0.54,
}

# ──────────────────────────────────────────────────────────────────
#  CRISPR systems (carried over from original)
# ──────────────────────────────────────────────────────────────────
CRISPR_SYSTEMS = {
    "1": {"name":"SpCas9 (S. pyogenes)","pam":"NGG","guide_len":20,
          "pam_pos":"3prime","cut_offset":-3,"enzyme":"Cas9",
          "organisms":"Bacteria, Fungi, most microbes"},
    "2": {"name":"SaCas9 (S. aureus)","pam":"NNGRRT","guide_len":21,
          "pam_pos":"3prime","cut_offset":-3,"enzyme":"Cas9",
          "organisms":"Gram-positive bacteria"},
    "3": {"name":"Cas12a/Cpf1 (Acidaminococcus)","pam":"TTTV","guide_len":23,
          "pam_pos":"5prime","cut_offset":18,"enzyme":"Cas12a",
          "organisms":"AT-rich genomes, Mycoplasma, Clostridium"},
    "4": {"name":"CjCas9 (C. jejuni)","pam":"NNNNRYAC","guide_len":22,
          "pam_pos":"3prime","cut_offset":-3,"enzyme":"Cas9",
          "organisms":"Campylobacter, gut microbiome"},
    "5": {"name":"NmCas9 (N. meningitidis)","pam":"NNNNGATT","guide_len":24,
          "pam_pos":"3prime","cut_offset":-3,"enzyme":"Cas9",
          "organisms":"Neisseria, AAV-compatible"},
}

DEGEN = {"N":"[ACGT]","R":"[AG]","Y":"[CT]","V":"[ACG]",
         "W":"[AT]","S":"[GC]","M":"[AC]","K":"[GT]"}


# ════════════════════════════════════════════════════════════════
#  ──  SEQUENCE UTILITIES  ──
# ════════════════════════════════════════════════════════════════

def clean_sequence(seq: str, strip_ambiguous: bool = True) -> str:
    seq = re.sub(r"[\s\d>]", "", seq).upper()
    seq = re.sub(r">.*", "", seq)
    if strip_ambiguous:
        cleaned = []
        for b in seq:
            if b in "ACGT":
                cleaned.append(b)
            elif b in AMBIGUOUS_MAP:
                cleaned.append(AMBIGUOUS_MAP[b])
        seq = "".join(cleaned)
    else:
        bad = set(seq) - set("ACGT")
        if bad:
            raise ValueError(f"Invalid characters: {bad}")
    return seq


def validate_sequence(seq: str) -> str:
    seq = clean_sequence(seq)
    if len(seq) < 10:
        raise ValueError("Sequence too short (min 10 bp).")
    return seq


def complement(base: str) -> str:
    return {"A":"T","T":"A","G":"C","C":"G"}.get(base.upper(),"N")


def reverse_complement(seq: str) -> str:
    return "".join(complement(b) for b in reversed(seq.upper()))


def gc_content(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def translate(seq: str) -> str:
    """Translate a DNA sequence to amino acids."""
    seq  = seq.upper()
    prot = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa    = GENETIC_CODE.get(codon, "?")
        prot.append(aa)
        if aa == "*":
            break
    return "".join(prot)


def dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def rna_complement(base: str) -> str:
    return {"A":"U","U":"A","G":"C","C":"G"}.get(base.upper(),"N")


# ════════════════════════════════════════════════════════════════
#  ──  FILE I/O  ──
# ════════════════════════════════════════════════════════════════

def _open_file(filepath: str):
    filepath = filepath.strip().strip('"').strip("'")
    if not os.path.isfile(filepath):
        raise FileNotFoundError(
            f"File not found: {filepath}\nCurrent dir: {os.getcwd()}"
        )
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt", encoding="utf-8", errors="ignore")
    return open(filepath, "r", encoding="utf-8", errors="ignore")


def parse_fasta(filepath: str) -> dict:
    records, current_header, current_seq = {}, None, []
    with _open_file(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records[current_header] = "".join(current_seq).upper()
                current_header = line[1:].split()[0]
                current_seq    = []
            else:
                current_seq.append(line.replace(" ", ""))
    if current_header is not None:
        records[current_header] = "".join(current_seq).upper()
    if not records:
        raise ValueError("No FASTA records found.")
    return records


def parse_genbank(filepath: str) -> dict:
    records, accession, in_origin, seq_parts = {}, None, False, []
    with _open_file(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("ACCESSION"):
                accession = line.split()[1]; in_origin = False; seq_parts = []
            elif line.startswith("ORIGIN"):
                in_origin = True
            elif line.startswith("//"):
                if accession and seq_parts:
                    seq = re.sub(r"[^ACGTacgt]", "", "".join(seq_parts)).upper()
                    records[accession] = seq
                in_origin = False; seq_parts = []; accession = None
            elif in_origin:
                seq_parts.append(re.sub(r"[\d\s]", "", line))
    if not records:
        raise ValueError("No sequences found in GenBank file.")
    return records


def detect_and_parse(filepath: str) -> dict:
    fp = filepath.strip().strip('"').strip("'").lower()
    fasta_exts   = (".fasta",".fa",".fna",".ffn",
                    ".fasta.gz",".fa.gz",".fna.gz",".ffn.gz")
    genbank_exts = (".gbff",".gb",".gbk",".gbff.gz",".gb.gz",".gbk.gz")
    if any(fp.endswith(e) for e in fasta_exts):
        return parse_fasta(filepath)
    elif any(fp.endswith(e) for e in genbank_exts):
        return parse_genbank(filepath)
    else:
        try:
            r = parse_fasta(filepath)
            if r: return r
        except Exception:
            pass
        return parse_genbank(filepath)


def handle_large_sequence(seq: str, record_id: str) -> str:
    from rich.console import Console
    from rich.prompt  import Prompt
    console = Console()
    if len(seq) <= 100_000:
        return seq
    MB = len(seq) / 1_000_000
    console.print(
        f"\n[yellow]⚠ Large sequence:[/yellow] [bold]{record_id}[/bold] "
        f"[dim]({MB:.2f} Mb)[/dim]"
    )
    console.print("  [1] Full sequence  [2] Extract coordinates  [3] First N bp")
    choice = Prompt.ask("How to proceed?", choices=["1","2","3"], default="2")
    if choice == "1":
        return seq
    elif choice == "2":
        try:
            s = int(Prompt.ask("Start (1-based)")) - 1
            e = int(Prompt.ask("End"))
            if s < 0 or e > len(seq) or s >= e:
                raise ValueError
            return seq[s:e]
        except ValueError:
            return seq
    else:
        try:
            n = int(Prompt.ask("Number of bp", default="5000"))
        except ValueError:
            n = 5000
        return seq[:min(n, len(seq))]


def get_sequence_from_user(prompt_label: str = "target") -> tuple:
    from rich.console import Console
    from rich.prompt  import Prompt
    console = Console()
    console.print(f"\n[bold cyan]── {prompt_label.upper()} SEQUENCE INPUT ──[/bold cyan]")
    console.print("  [1] Paste sequence directly")
    console.print("  [2] Load from file  (.fasta .fna .fa .gbff .gz)")
    choice = Prompt.ask("Choice", choices=["1","2"], default="1")
    if choice == "1":
        console.print("[dim]Paste sequence then press Enter twice:[/dim]")
        lines = []
        while True:
            line = input()
            if line == "" and lines:
                break
            lines.append(line)
        seq = validate_sequence(" ".join(lines))
        return seq, "pasted_sequence"
    else:
        console.print("\n[dim]Tip: drag-and-drop file into terminal to fill path.[/dim]")
        path    = Prompt.ask("File path")
        with console.status("[bold green]Reading file...[/bold green]", spinner="dots"):
            records = detect_and_parse(path)
        console.print(f"\n[green]✔ {len(records)} record(s) found[/green]")
        if len(records) == 1:
            rid, raw = next(iter(records.items()))
        else:
            ids = list(records.keys())
            for i, rid in enumerate(ids, 1):
                console.print(f"  [bold yellow][{i}][/bold yellow] {rid} "
                              f"[dim]({len(records[rid]):,} bp)[/dim]")
            try:
                idx = int(Prompt.ask("Select record", default="1")) - 1
            except ValueError:
                idx = 0
            rid = ids[max(0, min(idx, len(ids)-1))]
            raw = records[rid]
        raw = clean_sequence(raw)
        seq = handle_large_sequence(raw, rid)
        seq = validate_sequence(seq)
        console.print(f"[green]✔ Ready:[/green] [bold]{rid}[/bold] "
                      f"[dim]({len(seq):,} bp)[/dim]")
        return seq, rid


# ════════════════════════════════════════════════════════════════
#  ──  THERMODYNAMICS  ──
# ════════════════════════════════════════════════════════════════

def tm_nearest_neighbor(seq: str, conc_primer: float = 250e-9,
                         conc_salt: float = 0.05) -> float:
    seq = seq.upper()
    dH  = INIT_PARAMS["dH"]
    dS  = INIT_PARAMS["dS"]
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if dinuc in NN_PARAMS:
            h, s  = NN_PARAMS[dinuc]
            dH   += h
            dS   += s
    R          = 1.987
    dS_total   = dS + R * math.log(conc_primer / 4.0)
    tm_celsius = (dH * 1000.0) / dS_total - 273.15
    if conc_salt > 0:
        tm_celsius /= (1.0 + 16.6 * math.log10(conc_salt) * (-1 / tm_celsius))
    return round(tm_celsius, 1)


def check_gc_clamp(primer: str, window: int = 3) -> bool:
    return any(b in "GC" for b in primer[-window:].upper())


def self_complementarity_score(primer: str) -> int:
    rc = reverse_complement(primer)
    max_run = run = 0
    for a, b in zip(primer, rc):
        if complement(a) == b:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0
    return max_run


def primer_dimer_check(fwd: str, rev: str, window: int = 4) -> bool:
    return reverse_complement(fwd[-window:]) == rev[-window:]


def rna_stem_energy(seq_i: str, seq_j: str) -> float:
    pair = (seq_i, seq_j)
    return RNA_STACK.get(pair, 0.0)


# ════════════════════════════════════════════════════════════════
#  ──  DISPLAY HELPERS  ──
# ════════════════════════════════════════════════════════════════

def truncate(seq: str, maxlen: int = 40) -> str:
    return seq if len(seq) <= maxlen else seq[:maxlen] + "…"


def format_sequence_block(seq: str, width: int = 60) -> str:
    return "\n".join(f"{i+1:>6} {seq[i:i+width]}"
                     for i in range(0, len(seq), width))


def expand_pam(pam: str) -> str:
    return "".join(DEGEN.get(c, c) for c in pam.upper())