# ================================================================
#  DBS  –  Digital BioScience Suite
#  structure.py  |  All biological engine classes
#
#  Phase 1: PrimerDesigner, CRISPRDesigner, CRISPRAnimator
#           RestrictionMapper, CodonOptimizer, PlasmidDesigner
#
#  Phase 2: siRNADesigner, RNAStructurePredictor,
#           PhylogeneticsBuilder
#
#  Phase 3: VaccineEpitopeDesigner, AMPDesigner, ARGScanner
# ================================================================

import re, time, math, os, random
from collections import defaultdict

from rich.console import Console
from rich.table   import Table
from rich.panel   import Panel
from rich.text    import Text
from rich.align   import Align
from rich        import box

from utilities import (
    reverse_complement, gc_content, tm_nearest_neighbor,
    check_gc_clamp, self_complementarity_score, primer_dimer_check,
    truncate, expand_pam, translate, dna_to_rna, rna_complement,
    rna_stem_energy, complement,
    GENETIC_CODE, CODON_TABLES, RESTRICTION_ENZYMES,
    ARG_DATABASE, MHC_I_MOTIFS, HYDROPHOBICITY,
    CRISPR_SYSTEMS, AMINO_ACIDS,
)

console = Console()


# ════════════════════════════════════════════════════════════════
#  ── PHASE 1 · A  ──  PRIMER DESIGNER
# ════════════════════════════════════════════════════════════════

class PrimerDesigner:
    DEFAULT = dict(min_len=18, max_len=25, gc_min=40.0, gc_max=60.0,
                   tm_min=52.0, tm_max=68.0, tm_delta=5.0, top_n=5)

    def __init__(self, sequence: str, **kwargs):
        self.seq = sequence.upper()
        self.cfg = {**self.DEFAULT, **kwargs}

    def _score_candidate(self, primer: str) -> dict | None:
        gc    = gc_content(primer)
        tm    = tm_nearest_neighbor(primer)
        clamp = check_gc_clamp(primer)
        sc    = self_complementarity_score(primer)
        fails = []
        if not (self.cfg["gc_min"] <= gc <= self.cfg["gc_max"]): fails.append("GC%")
        if not (self.cfg["tm_min"] <= tm <= self.cfg["tm_max"]): fails.append("Tm")
        if not clamp: fails.append("GC-clamp")
        if sc > 6:    fails.append("self-complement")
        score = (100 - abs(gc-50)*0.8 - abs(tm-60)*0.5
                 + (5 if clamp else 0) - sc * 2)
        return {"primer":primer,"len":len(primer),"gc":round(gc,1),
                "tm":tm,"clamp":clamp,"sc":sc,"score":round(score,1),
                "warnings":fails}

    def _scan_candidates(self, strand_seq: str, label: str) -> list:
        candidates, cfg = [], self.cfg
        for start in range(0, len(strand_seq) - cfg["min_len"] + 1):
            for length in range(cfg["min_len"], cfg["max_len"] + 1):
                end = start + length
                if end > len(strand_seq): break
                primer = strand_seq[start:end]
                result = self._score_candidate(primer)
                if result and not result["warnings"]:
                    result.update({"start":start+1,"end":end,"label":label})
                    candidates.append(result)
        seen, unique = set(), []
        for c in sorted(candidates, key=lambda x: -x["score"]):
            if c["primer"] not in seen:
                seen.add(c["primer"]); unique.append(c)
        return unique

    def design(self) -> dict:
        console.print("\n[bold green]⚙  Scanning for primer candidates...[/bold green]")
        fwd = self._scan_candidates(self.seq, "Forward")
        rev = self._scan_candidates(reverse_complement(self.seq), "Reverse")
        pairs = []
        for f in fwd[:20]:
            for r in rev[:20]:
                dt = abs(f["tm"] - r["tm"])
                if dt <= self.cfg["tm_delta"] and not primer_dimer_check(f["primer"],r["primer"]):
                    pairs.append({"fwd":f,"rev":r,"delta_tm":round(dt,1),
                                  "pair_score":f["score"]+r["score"]-dt*2})
        pairs.sort(key=lambda x: -x["pair_score"])
        return {"forward":fwd[:self.cfg["top_n"]],
                "reverse":rev[:self.cfg["top_n"]], "pairs":pairs[:5]}

    @staticmethod
    def display_results(results: dict, seq: str):
        for direction, color in [("forward","cyan"),("reverse","magenta")]:
            cands = results[direction]
            if not cands:
                console.print(f"[red]No valid {direction} primers found.[/red]"); continue
            tbl = Table(title=f"[bold {color}]{direction.upper()} PRIMERS[/bold {color}]",
                        box=box.ROUNDED, show_lines=True, header_style=f"bold {color}")
            for col in ["#","Sequence (5'→3')","Len","Tm","GC%","Clamp","SC","Score"]:
                tbl.add_column(col, justify="center")
            for i, c in enumerate(cands, 1):
                cc = "green" if c["clamp"] else "red"
                sc = "green" if c["sc"]<=4 else "yellow" if c["sc"]<=6 else "red"
                tbl.add_row(str(i),f"[{color}]{c['primer']}[/{color}]",str(c["len"]),
                            str(c["tm"]),str(c["gc"]),
                            f"[{cc}]{'✔' if c['clamp'] else '✘'}[/{cc}]",
                            f"[{sc}]{c['sc']}[/{sc}]",f"[bold]{c['score']}[/bold]")
            console.print(tbl)
        if not results["pairs"]:
            console.print("[red]⚠ No valid primer pairs found.[/red]"); return
        tbl2 = Table(title="[bold yellow]★ RECOMMENDED PAIRS[/bold yellow]",
                     box=box.DOUBLE, show_lines=True, header_style="bold yellow")
        for col in ["#","Forward (5'→3')","Fwd Tm","Reverse (5'→3')","Rev Tm","ΔTm","Score"]:
            tbl2.add_column(col, justify="center")
        for i, p in enumerate(results["pairs"], 1):
            tbl2.add_row(str(i),f"[cyan]{p['fwd']['primer']}[/cyan]",
                         str(p["fwd"]["tm"]),f"[magenta]{p['rev']['primer']}[/magenta]",
                         str(p["rev"]["tm"]),f"[green]{p['delta_tm']}[/green]",
                         f"[bold]{round(p['pair_score'],1)}[/bold]")
        console.print(tbl2)
        best = results["pairs"][0]
        console.print(Panel(
            f"Forward : [cyan]{best['fwd']['primer']}[/cyan]\n"
            f"Reverse : [magenta]{best['rev']['primer']}[/magenta]\n"
            f"Fwd Tm  : [yellow]{best['fwd']['tm']} °C[/yellow]\n"
            f"Rev Tm  : [yellow]{best['rev']['tm']} °C[/yellow]\n"
            f"Anneal  : [bold green]{round(min(best['fwd']['tm'],best['rev']['tm'])-5,1)} °C[/bold green]\n"
            f"ΔTm     : {best['delta_tm']} °C",
            title="[bold yellow]Primer Summary[/bold yellow]", border_style="yellow"))


# ════════════════════════════════════════════════════════════════
#  ── PHASE 1 · B  ──  CRISPR DESIGNER + ANIMATOR
# ════════════════════════════════════════════════════════════════

class CRISPRDesigner:
    def __init__(self, sequence: str, system_key: str):
        self.seq = sequence.upper()
        self.sys = CRISPR_SYSTEMS[system_key]

    def _build_regex(self) -> str:
        gl  = self.sys["guide_len"]
        pam = expand_pam(self.sys["pam"])
        if self.sys["pam_pos"] == "3prime":
            return rf"(?=(.{{{gl}}})({pam}))"
        return rf"(?=({pam})(.{{{gl}}}))"

    def _score_guide(self, guide: str) -> dict:
        gc = gc_content(guide)
        gc_score = 40 if 40<=gc<=65 else max(0,40-abs(gc-52.5)*1.5)
        homo_pen = 20 if any(b*4 in guide for b in "ACGT") else 0
        sc       = self_complementarity_score(guide)
        sc_score = max(0, 20 - sc*3)
        seed     = guide[-12:] if self.sys["pam_pos"]=="3prime" else guide[:12]
        sgc      = gc_content(seed)
        seed_sc  = 20 if 40<=sgc<=65 else max(0,20-abs(sgc-52.5))
        return {"gc":round(gc,1),"homo":homo_pen>0,"sc":sc,
                "seed_gc":round(sgc,1),
                "score":round(gc_score+(20-homo_pen)+sc_score+seed_sc,1)}

    def find_guides(self) -> list:
        pattern = self._build_regex()
        results = []
        for label, strand in [("sense (+)", self.seq),
                               ("antisense (−)", reverse_complement(self.seq))]:
            try:
                for m in re.finditer(pattern, strand):
                    if self.sys["pam_pos"] == "3prime":
                        proto, pam = m.group(1), m.group(2)
                        pos = m.start() + 1
                    else:
                        pam, proto = m.group(1), m.group(2)
                        pos = m.start() + len(pam) + 1
                    if len(proto) != self.sys["guide_len"]: continue
                    if set(proto) - set("ACGT"):           continue
                    sc = self._score_guide(proto)
                    results.append({"guide":proto,"pam":pam,"strand":label,
                                    "position":pos,
                                    "cut_site":pos+self.sys["cut_offset"],
                                    **sc})
            except re.error:
                pass
        results.sort(key=lambda x: -x["score"])
        return results

    @staticmethod
    def display_results(guides: list, system: dict, top_n: int = 10):
        if not guides:
            console.print("[red]⚠ No guide RNAs found.[/red]"); return
        show = guides[:top_n]
        tbl = Table(title=f"[bold green]CRISPR Guides – {system['name']}[/bold green]",
                    box=box.ROUNDED, show_lines=True, header_style="bold green")
        for col in ["#","Guide (5'→3')","PAM","Strand","Pos","Cut","GC%","Seed","Homo","Score"]:
            tbl.add_column(col, justify="center")
        for i, g in enumerate(show, 1):
            sc = "[bold green]" if g["score"]>=70 else "[bold yellow]" if g["score"]>=45 else "[bold red]"
            ec = sc.replace("[","").replace("]","")
            tbl.add_row(str(i),f"[bold cyan]{g['guide']}[/bold cyan]",
                        f"[yellow]{g['pam']}[/yellow]",g["strand"],
                        str(g["position"]),str(g["cut_site"]),
                        str(g["gc"]),str(g["seed_gc"]),
                        "[red]✘[/red]" if g["homo"] else "[green]✔[/green]",
                        f"{sc}{g['score']}[/{ec}]")
        console.print(tbl)
        best = show[0]
        console.print(Panel(
            f"Protospacer : [cyan]{best['guide']}[/cyan]\n"
            f"PAM         : [yellow]{best['pam']}[/yellow]\n"
            f"Strand      : {best['strand']}\n"
            f"Position    : {best['position']}\n"
            f"Cut Site    : ~{best['cut_site']}\n"
            f"GC%%         : {best['gc']}\n"
            f"Score       : [bold green]{best['score']}[/bold green]",
            title="[bold green]★ Best Guide[/bold green]", border_style="green"))


class CRISPRAnimator:
    CAS9_ART = [
        "      ╔═══════════════╗",
        "   ┌──╢  ✦  Cas9  ✦  ╟──┐",
        "   │  ╚═══════════════╝  │",
        "   │   [HNH] [RuvC]     │",
        "   └─────────────────────┘",
    ]
    CAS12A_ART = [
        "      ╔════════════════╗",
        "   ┌──╢  ✦  Cas12a ✦  ╟──┐",
        "   │  ╚════════════════╝  │",
        "   │      [RuvC]         │",
        "   └──────────────────────┘",
    ]

    def __init__(self, sequence: str, guide: dict, system: dict):
        self.seq    = sequence.upper()
        self.guide  = guide
        self.system = system
        self.pos    = guide["position"] - 1
        self.enzyme = system["enzyme"]

    def _clear(self): os.system("cls" if os.name=="nt" else "clear")

    def _print_header(self, phase: str, color: str = "green"):
        console.print(Panel(f"[bold {color}]{phase}[/bold {color}]",
                            title=f"[bold]DBS CRISPR Animation – {self.system['name']}[/bold]",
                            border_style=color))

    def _dna_line(self, seq, cut=-1, hs=-1, he=-1) -> Text:
        t = Text()
        for i, b in enumerate(seq):
            if hs <= i < he:    t.append(b, style="bold cyan on blue")
            elif i in (cut, cut+1): t.append(b, style="bold red")
            else:               t.append(b, style="white")
        return t

    def _show_duplex(self, title, hl_s=-1, hl_e=-1, cut=-1, scanner=-1):
        disp  = min(70, len(self.seq))
        off   = max(0, min(self.pos-10, len(self.seq)-disp))
        win   = self.seq[off:off+disp]
        rc    = reverse_complement(win)
        ahs   = hl_s - off; ahe = hl_e - off; ac = cut - off if cut>=0 else -1
        console.print(f"[dim]{title}[/dim]")
        console.print("5'  ", end=""); console.print(self._dna_line(win,ac,ahs,ahe), end=""); console.print("  3'")
        console.print(f"    {'|'*len(win)}")
        console.print("3'  ", end=""); console.print(self._dna_line(rc,ac,ahs,ahe), end=""); console.print("  5'")
        if 0 <= scanner-off < disp:
            console.print("    "+" "*(scanner-off)+"[bold yellow]▲[/bold yellow]")
            console.print("    "+" "*(scanner-off)+"[bold yellow]gRNA[/bold yellow]")
        console.print()

    def _print_enzyme(self, bound=False, cutting=False):
        art   = self.CAS12A_ART if self.enzyme=="Cas12a" else self.CAS9_ART
        color = "red" if cutting else ("green" if bound else "blue")
        for l in art: console.print(f"[{color}]{l}[/{color}]")

    def run(self, speed: float = 0.07):
        guide_seq = self.guide["guide"]
        guide_len = len(guide_seq)
        cut_site  = self.guide["cut_site"] - 1
        # Phase 1 – Introduction
        self._clear(); self._print_header("Phase 1/6 – DNA Duplex","blue")
        self._show_duplex("Genomic DNA:"); self._print_enzyme(); time.sleep(1.8)
        # Phase 2 – Scanning
        step = max(1, guide_len//8)
        for sp in list(range(0, self.pos, step)) + [self.pos]:
            self._clear(); self._print_header("Phase 2/6 – Guide RNA Scanning","yellow")
            console.print(f"[yellow]Scanning for PAM: [bold]{self.system['pam']}[/bold]...[/yellow]\n")
            self._show_duplex("Scanning:", scanner=sp); self._print_enzyme()
            console.print(f"\n[dim]Guide: 5'-[cyan]{guide_seq}[/cyan]-3'[/dim]"); time.sleep(speed)
        # Phase 3 – PAM recognition
        for _ in range(4):
            self._clear(); self._print_header("Phase 3/6 – PAM Recognition","cyan")
            console.print(f"[bold cyan]✦ PAM [yellow]{self.guide['pam']}[/yellow] found at position {self.guide['position']}[/bold cyan]\n")
            self._show_duplex("PAM recognised:", hl_s=self.pos, hl_e=self.pos+guide_len, scanner=self.pos)
            self._print_enzyme(); time.sleep(speed*4)
        # Phase 4 – R-loop
        self._clear(); self._print_header("Phase 4/6 – R-loop Formation","green")
        target = self.seq[self.pos:self.pos+guide_len]
        pairing = Text()
        for g, d in zip(guide_seq, target):
            pairing.append("|" if g==d else "x", style="green" if g==d else "red")
        console.print(f"  Guide 5'-[cyan]{guide_seq}[/cyan]-3'")
        console.print("        3' ", end=""); console.print(pairing, end=""); console.print(" 5'")
        console.print(f"  DNA   3'-[white]{target}[/white]-5'\n")
        self._show_duplex("Enzyme bound:", hl_s=self.pos, hl_e=self.pos+guide_len)
        self._print_enzyme(bound=True); time.sleep(1.8)
        # Phase 5 – Cleavage
        for flash in range(5):
            self._clear()
            color = "red" if flash%2==0 else "yellow"
            self._print_header("Phase 5/6 – DOUBLE-STRAND BREAK", color)
            console.print(f"[bold {color}]⚡ {self.enzyme} nuclease domains activate![/bold {color}]\n")
            self._show_duplex(f"DSB at ~{self.guide['cut_site']}:",
                              hl_s=self.pos, hl_e=self.pos+guide_len, cut=cut_site)
            self._print_enzyme(bound=True, cutting=(flash%2==0))
            if flash%2==0:
                lf = self.seq[max(0,cut_site-12):cut_site]
                rf = self.seq[cut_site:cut_site+12]
                console.print(f"\n  5'─[white]{lf}[/white][bold red]↓[/bold red][white]{rf}[/white]─3'")
                console.print(f"  3'─[white]{reverse_complement(lf)[::-1]}[/white][bold red]↑[/bold red][white]{reverse_complement(rf)[::-1]}[/white]─5'")
            time.sleep(0.45)
        # Phase 6 – Summary
        self._clear(); self._print_header("Phase 6/6 – Edit Complete ✔","bright_green")
        lf = self.seq[:cut_site]; rf = self.seq[cut_site:]
        console.print(Panel(
            f"[bold green]Double-strand break introduced![/bold green]\n\n"
            f"  Fragment A : 5'-[cyan]...{truncate(lf,28)}[/cyan]-3'  [dim]({len(lf)} bp)[/dim]\n"
            f"  Fragment B : 5'-[magenta]{truncate(rf,28)}...[/magenta]-3'  [dim]({len(rf)} bp)[/dim]\n\n"
            f"  Guide  : [cyan]{guide_seq}[/cyan]\n"
            f"  PAM    : [yellow]{self.guide['pam']}[/yellow]\n"
            f"  System : [green]{self.system['name']}[/green]\n\n"
            f"[dim]Cell may repair via NHEJ (indels) or HDR (precise edit).[/dim]",
            title="[bold bright_green]✦ CRISPR Edit Summary[/bold bright_green]",
            border_style="bright_green"))
        console.print(Align.center(
            "\n[bold yellow]✦ ─────────────────────────────── ✦[/bold yellow]\n"
            "[bold bright_green]    CRISPR Editing Event Complete!   [/bold bright_green]\n"
            "[bold yellow]✦ ─────────────────────────────── ✦[/bold yellow]"))


# ════════════════════════════════════════════════════════════════
#  ── PHASE 1 · C  ──  RESTRICTION ENZYME MAPPER
# ════════════════════════════════════════════════════════════════

class RestrictionMapper:
    """
    Scans a DNA sequence for restriction enzyme recognition sites.
    Supports REBASE-derived database, degenerate bases, virtual digest.
    """

    def __init__(self, sequence: str):
        self.seq = sequence.upper()

    def _match_site(self, seq: str, pattern: str) -> list:
        """Find all positions where pattern matches seq (IUPAC degenerate)."""
        regex = "".join({
            "N":"[ACGT]","R":"[AG]","Y":"[CT]","S":"[GC]",
            "W":"[AT]","K":"[GT]","M":"[AC]","B":"[CGT]",
            "D":"[AGT]","H":"[ACT]","V":"[ACG]",
        }.get(c, c) for c in pattern.upper())
        return [m.start() for m in re.finditer(f"(?={regex})", seq)]

    def find_sites(self, enzymes: list = None) -> dict:
        """
        Returns {enzyme_name: [cut_positions]} for all enzymes with ≥1 hit.
        enzymes = None → scan all known enzymes.
        """
        db  = {k:v for k,v in RESTRICTION_ENZYMES.items()
               if enzymes is None or k in enzymes}
        hits = {}
        for name, (site, cut_offset) in db.items():
            positions = self._match_site(self.seq, site)
            if positions:
                hits[name] = {
                    "site"       : site,
                    "positions"  : [p+1 for p in positions],   # 1-based
                    "cut_offset" : cut_offset,
                    "cut_type"   : self._overhang_type(site, cut_offset),
                    "count"      : len(positions),
                }
        return hits

    def _overhang_type(self, site: str, cut_offset: int) -> str:
        site_len = len(re.sub(r"[^ACGT]","",site))
        mid      = len(site) // 2
        if cut_offset == mid: return "Blunt"
        if cut_offset < mid:  return "5' overhang"
        return "3' overhang"

    def virtual_digest(self, enzyme1: str, enzyme2: str = None) -> list:
        """Simulate single or double digest, return fragment sizes."""
        cuts = []
        for name in ([enzyme1] if enzyme2 is None else [enzyme1, enzyme2]):
            if name not in RESTRICTION_ENZYMES: continue
            site, off = RESTRICTION_ENZYMES[name]
            for p in self._match_site(self.seq, site):
                cuts.append(p + off)
        cuts = sorted(set(cuts))
        if not cuts: return []
        frags = [cuts[0]] + [cuts[i]-cuts[i-1] for i in range(1,len(cuts))]
        frags.append(len(self.seq) - cuts[-1])
        return frags

    def draw_map(self, hits: dict, width: int = 60) -> str:
        """ASCII linear restriction map."""
        lines = []
        seq_line = "─" * width
        lines.append(f"5'──{seq_line}──3'  ({len(self.seq)} bp)")
        for name, data in sorted(hits.items(), key=lambda x: x[1]["positions"][0]):
            for pos in data["positions"]:
                frac   = (pos-1) / len(self.seq)
                tick   = int(frac * width)
                arrow  = " "*tick + "│"
                label  = f" {name} ({pos})"
                lines.append(f"    {arrow}{label}")
        return "\n".join(lines)

    @staticmethod
    def display_results(hits: dict, digest_frags: list = None, map_str: str = ""):
        if not hits:
            console.print("[red]No restriction sites found.[/red]"); return
        tbl = Table(title="[bold cyan]RESTRICTION ENZYME SITES[/bold cyan]",
                    box=box.ROUNDED, show_lines=True, header_style="bold cyan")
        for col in ["Enzyme","Recognition Site","Positions","Cut Type","#Cuts"]:
            tbl.add_column(col, justify="center")
        for name, d in sorted(hits.items(), key=lambda x: -x[1]["count"]):
            pos_str = ", ".join(str(p) for p in d["positions"][:6])
            if len(d["positions"]) > 6: pos_str += "…"
            oc = {"Blunt":"dim","5' overhang":"cyan","3' overhang":"magenta"}.get(d["cut_type"],"white")
            tbl.add_row(f"[bold]{name}[/bold]",
                        f"[green]{d['site']}[/green]",
                        pos_str,
                        f"[{oc}]{d['cut_type']}[/{oc}]",
                        f"[bold]{d['count']}[/bold]")
        console.print(tbl)
        if map_str:
            console.print(Panel(f"[dim]{map_str}[/dim]",
                                title="[bold cyan]Linear Restriction Map[/bold cyan]",
                                border_style="cyan"))
        if digest_frags:
            gel = "\n".join(f"  {f:>6} bp  {'█'*min(40,f//20)}" for f in sorted(digest_frags, reverse=True))
            console.print(Panel(f"[bold]Virtual Gel (largest→smallest):[/bold]\n{gel}",
                                title="[bold cyan]Virtual Digest[/bold cyan]",
                                border_style="cyan"))


# ════════════════════════════════════════════════════════════════
#  ── PHASE 1 · D  ──  CODON OPTIMIZER
# ════════════════════════════════════════════════════════════════

class CodonOptimizer:
    """
    Back-translates or re-codes a DNA/protein sequence for maximum
    expression in a target host using real Kazusa codon usage tables.
    Calculates Codon Adaptation Index (CAI) before and after.
    """

    def __init__(self, sequence: str, host: str = "E. coli K-12",
                 input_type: str = "dna"):
        self.input_type = input_type
        self.host       = host
        self.table      = CODON_TABLES.get(host, CODON_TABLES["E. coli K-12"])
        if input_type == "dna":
            self.dna_seq  = sequence.upper()
            self.prot_seq = translate(sequence)
        else:
            self.prot_seq = sequence.upper()
            self.dna_seq  = None

    # ── build reverse codon map: aa → {codon: freq} ──
    def _rev_codon_map(self) -> dict:
        rc = defaultdict(dict)
        for codon, aa in GENETIC_CODE.items():
            if codon in self.table:
                rc[aa][codon] = self.table[codon]
        return rc

    def _best_codon(self, aa: str, rcm: dict) -> str:
        codons = rcm.get(aa, {})
        if not codons: return "NNN"
        return max(codons, key=codons.get)

    def optimize(self) -> dict:
        rcm = self._rev_codon_map()
        opt_codons, orig_codons = [], []

        if self.dna_seq:
            for i in range(0, len(self.dna_seq)-2, 3):
                codon = self.dna_seq[i:i+3]
                aa    = GENETIC_CODE.get(codon, "?")
                if aa == "*": break
                orig_codons.append(codon)
                opt_codons.append(self._best_codon(aa, rcm))
        else:
            for aa in self.prot_seq:
                if aa == "*": break
                opt_codons.append(self._best_codon(aa, rcm))

        opt_seq    = "".join(opt_codons)
        orig_cai   = self._calc_cai(orig_codons) if orig_codons else 0.0
        opt_cai    = self._calc_cai(opt_codons)
        orig_gc    = gc_content(self.dna_seq) if self.dna_seq else 0.0
        opt_gc     = gc_content(opt_seq) if opt_seq else 0.0

        # Flag rare codons (freq < 0.15) in original
        rare = [f"{c}({GENETIC_CODE.get(c,'?')}) @{i*3+1}"
                for i, c in enumerate(orig_codons)
                if self.table.get(c, 1.0) < 0.15]

        return {
            "original"      : self.dna_seq or "",
            "optimized"     : opt_seq,
            "protein"       : translate(opt_seq),
            "orig_cai"      : round(orig_cai, 3),
            "opt_cai"       : round(opt_cai, 3),
            "orig_gc"       : round(orig_gc, 1),
            "opt_gc"        : round(opt_gc, 1),
            "host"          : self.host,
            "length"        : len(opt_seq),
            "rare_codons"   : rare[:10],
            "num_rare"      : len(rare),
        }

    def _calc_cai(self, codons: list) -> float:
        """Codon Adaptation Index — geometric mean of relative synonymous codon usage."""
        if not codons: return 0.0
        rcm    = self._rev_codon_map()
        scores = []
        for codon in codons:
            aa = GENETIC_CODE.get(codon, "?")
            if aa in ("?","*"): continue
            synonyms = rcm.get(aa, {})
            if not synonyms: continue
            max_freq = max(synonyms.values())
            freq     = self.table.get(codon, 0.001)
            scores.append(math.log(freq / max_freq + 1e-9))
        return math.exp(sum(scores)/len(scores)) if scores else 0.0

    @staticmethod
    def display_results(result: dict):
        console.print(Panel(
            f"[bold]Host organism :[/bold] [green]{result['host']}[/green]\n"
            f"[bold]Sequence length:[/bold] {result['length']} bp\n\n"
            f"[bold]CAI (original) :[/bold] [{'red' if result['orig_cai']<0.6 else 'yellow'}]{result['orig_cai']}[/{'red' if result['orig_cai']<0.6 else 'yellow'}]\n"
            f"[bold]CAI (optimized):[/bold] [bold green]{result['opt_cai']}[/bold green]  "
            f"[dim](≥0.8 = high expression expected)[/dim]\n\n"
            f"[bold]GC% (original) :[/bold] {result['orig_gc']}%\n"
            f"[bold]GC% (optimized):[/bold] {result['opt_gc']}%\n\n"
            f"[bold]Rare codons    :[/bold] [{'red' if result['num_rare']>3 else 'green'}]{result['num_rare']}[/{'red' if result['num_rare']>3 else 'green'}]\n"
            + (f"  [dim]{', '.join(result['rare_codons'])}[/dim]\n" if result["rare_codons"] else "") +
            f"\n[bold]Optimized sequence (5'→3'):[/bold]\n[cyan]{result['optimized'][:120]}{'…' if len(result['optimized'])>120 else ''}[/cyan]\n\n"
            f"[bold]Translated protein:[/bold]\n[magenta]{result['protein'][:80]}{'…' if len(result['protein'])>80 else ''}[/magenta]",
            title="[bold yellow]★ Codon Optimization Results[/bold yellow]",
            border_style="yellow"))


# ════════════════════════════════════════════════════════════════
#  ── PHASE 1 · E  ──  PLASMID DESIGNER
# ════════════════════════════════════════════════════════════════

class PlasmidDesigner:
    """
    Assembles a virtual expression plasmid from modular parts.
    Checks reading frames, generates circular ASCII map, predicts expression.
    """

    PROMOTERS = {
        "1": ("T7",        "Very strong (IPTG-inducible)", "Prokaryote",  "TAATACGACTCACTATA"),
        "2": ("Ptrc",      "Strong (IPTG-inducible)",      "Prokaryote",  "TTGACAATTAATCATCGGCTCGTATAATGTGTGGA"),
        "3": ("Plac",      "Medium (IPTG-inducible)",      "Prokaryote",  "TTACACTTTATGCTTCCGGCTCGTATGTTGTGTGG"),
        "4": ("Para",      "Arabinose-inducible",          "Prokaryote",  "TTTTCTCCTTTAGTGAGGGTTAATTCCGAGCTCAG"),
        "5": ("Ptac",      "Strong hybrid promoter",       "Prokaryote",  "TTGACAATTAATCATCGGCTCGTATAATGTGTGGA"),
        "6": ("CMV",       "Very strong constitutive",     "Eukaryote",   "GACATTGATTATTGACTAGTTATTAATAGTAATC"),
        "7": ("EF1alpha",  "Strong constitutive",          "Eukaryote",   "GCGAGGGCCTATTTCCCATGATTCCTTCATATTTG"),
        "8": ("PGK",       "Medium constitutive",          "Eukaryote",   "GCGGTTCACTCGCCGTATTCGGAGCCGGCTCGCTG"),
    }

    RBS_OPTIONS = {
        "1": ("Shine-Dalgarno (strong)",  "AAGGAGG", "prokaryote", 0.95),
        "2": ("Shine-Dalgarno (medium)", "AAGGAG",  "prokaryote", 0.65),
        "3": ("Kozak (optimal)",         "GCCACCATG","eukaryote", 0.90),
        "4": ("Kozak (minimal)",         "ACCATG",  "eukaryote",  0.60),
        "5": ("IRES (cap-independent)",  "GGATCC",  "eukaryote",  0.40),
    }

    TAGS = {
        "1": ("His6",       "CACCACCACCACCACCAC", "N-term/C-term", "Ni-NTA purification"),
        "2": ("His10",      "CACCACCACCACCACCACCACCACCAC","N-term/C-term","High-affinity Ni-NTA"),
        "3": ("FLAG",       "GACTACAAAGACCATGACGGTGATTATAAAGATCATGATATCGATTACAAGGATGACGATGACAAG","N-term","Anti-FLAG antibody"),
        "4": ("Strep-II",   "TGGAGCCACCCGCAGTTCGAAAAATAA","C-term","Strep-Tactin purification"),
        "5": ("GST",        "ATGTCCCCTATACTAGGTTAT","N-term","Glutathione purification"),
        "6": ("GFP",        "ATGGTGAGCAAGGGCGAGGAG","C-term","Fluorescence"),
        "7": ("None",       "","","No tag"),
    }

    TERMINATORS = {
        "1": ("T7 Te",      "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG","Strong"),
        "2": ("rrnB T1T2",  "GCTTTTTTTATACTAAGTTGGCATTATAAAAAAGCATTGCTTATCAATTTGTTGCAA","Strong dual"),
        "3": ("CYC1",       "CATGTAATTAGTTATGTCACGCTTACATTCACGCCCTCCCCCACATCC","Yeast"),
        "4": ("BGH polyA",  "CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCC","Mammalian"),
    }

    ORIGINS = {
        "1": ("ColE1/pUC",  "High copy (~500-700)",  "E. coli"),
        "2": ("p15A",       "Medium copy (~10-12)",  "E. coli compatible"),
        "3": ("pSC101",     "Low copy (~5)",         "E. coli stable"),
        "4": ("f1/M13",     "Phagemid",              "E. coli ssDNA production"),
        "5": ("2micron",    "High copy (yeast)",     "S. cerevisiae"),
        "6": ("CEN/ARS",    "Low copy (yeast)",      "S. cerevisiae stable"),
    }

    RESISTANCES = {
        "1": ("AmpR/bla",   "Ampicillin",    100),
        "2": ("KanR/aph",   "Kanamycin",     50),
        "3": ("CmR/cat",    "Chloramphenicol",25),
        "4": ("TetR/tetA",  "Tetracycline",   15),
        "5": ("SpecR/aadA", "Spectinomycin", 100),
        "6": ("ZeoR/Sh-ble","Zeocin",        25),
    }

    def __init__(self, insert_seq: str = ""):
        self.insert = insert_seq.upper() if insert_seq else ""

    def interactive_build(self) -> dict:
        from rich.prompt import Prompt
        console.print("\n[bold cyan]── PLASMID BUILDER ──[/bold cyan]")

        # Promoter
        self._show_choice_table("PROMOTER", self.PROMOTERS,
                                ["#","Name","Strength","Host","Sequence (partial)"])
        pk = Prompt.ask("Select promoter", choices=list(self.PROMOTERS), default="1")

        # RBS
        self._show_choice_table("RIBOSOME BINDING SITE / KOZAK", self.RBS_OPTIONS,
                                ["#","Name","Sequence","Host","Efficiency"])
        rk = Prompt.ask("Select RBS", choices=list(self.RBS_OPTIONS), default="1")

        # Tag
        self._show_choice_table("PURIFICATION TAG", self.TAGS,
                                ["#","Tag","Sequence","Position","Use"])
        tk = Prompt.ask("Select tag", choices=list(self.TAGS), default="1")

        # Terminator
        self._show_choice_table("TERMINATOR", self.TERMINATORS,
                                ["#","Name","Sequence (partial)","Strength"])
        ttk = Prompt.ask("Select terminator", choices=list(self.TERMINATORS), default="1")

        # Origin
        self._show_choice_table("ORIGIN OF REPLICATION", self.ORIGINS,
                                ["#","Origin","Copy Number","Host"])
        ok = Prompt.ask("Select origin", choices=list(self.ORIGINS), default="1")

        # Resistance
        self._show_choice_table("ANTIBIOTIC RESISTANCE", self.RESISTANCES,
                                ["#","Gene","Antibiotic","Concentration (μg/mL)"])
        rsk = Prompt.ask("Select resistance", choices=list(self.RESISTANCES), default="1")

        # Assemble
        prom  = self.PROMOTERS[pk]
        rbs   = self.RBS_OPTIONS[rk]
        tag   = self.TAGS[tk]
        term  = self.TERMINATORS[ttk]
        ori   = self.ORIGINS[ok]
        res   = self.RESISTANCES[rsk]

        insert_w_tag = (tag[1] + self.insert) if tag[1] else self.insert
        expression_score = round(
            ({"Very strong":0.95,"Strong":0.80,"Medium":0.55,"Arabinose-inducible":0.60,
              "Very strong constitutive":0.90,"Strong constitutive":0.80,
              "Medium constitutive":0.60,"Strong hybrid promoter":0.85,
              "Strong (IPTG-inducible)":0.80,"Weak":0.30}.get(prom[1],0.50))
            * rbs[3], 2)

        plasmid_seq = (prom[3] + rbs[1] + insert_w_tag + term[1])
        total_size  = len(plasmid_seq) + 2000   # approx backbone

        return {
            "promoter"         : prom,
            "rbs"              : rbs,
            "tag"              : tag,
            "terminator"       : term,
            "origin"           : ori,
            "resistance"       : res,
            "insert"           : self.insert,
            "plasmid_seq"      : plasmid_seq,
            "total_size"       : total_size,
            "expression_score" : expression_score,
        }

    @staticmethod
    def _show_choice_table(title: str, data: dict, cols: list):
        tbl = Table(title=f"[bold cyan]{title}[/bold cyan]",
                    box=box.SIMPLE, header_style="bold cyan")
        for c in cols: tbl.add_column(c)
        for key, vals in data.items():
            tbl.add_row(key, *[str(v)[:50] for v in vals])
        console.print(tbl)

    @staticmethod
    def display_results(plan: dict):
        size = plan["total_size"]
        circ = PlasmidDesigner._draw_circle(size, plan)
        exp  = plan["expression_score"]
        exp_color = "green" if exp>=0.7 else "yellow" if exp>=0.4 else "red"
        console.print(Panel(
            f"{circ}\n\n"
            f"[bold]Promoter    :[/bold] [cyan]{plan['promoter'][0]}[/cyan]  "
            f"({plan['promoter'][1]})\n"
            f"[bold]RBS/Kozak   :[/bold] [green]{plan['rbs'][0]}[/green]\n"
            f"[bold]Insert      :[/bold] [white]{len(plan['insert'])} bp[/white]\n"
            f"[bold]Tag         :[/bold] [magenta]{plan['tag'][0]}[/magenta]  "
            f"– {plan['tag'][3]}\n"
            f"[bold]Terminator  :[/bold] {plan['terminator'][0]}\n"
            f"[bold]Origin      :[/bold] {plan['origin'][0]}  "
            f"[dim]({plan['origin'][1]})[/dim]\n"
            f"[bold]Resistance  :[/bold] [red]{plan['resistance'][0]}[/red]  "
            f"({plan['resistance'][1]} {plan['resistance'][2]} μg/mL)\n"
            f"[bold]Total size  :[/bold] ~{size:,} bp\n"
            f"[bold]Expression  :[/bold] [{exp_color}]{exp} / 1.0[/{exp_color}]",
            title="[bold yellow]★ Plasmid Design Summary[/bold yellow]",
            border_style="yellow"))

    @staticmethod
    def _draw_circle(size: int, plan: dict) -> str:
        lines = [
            "          ╭──── AmpR ──╮",
            "      ╭───╯            ╰───╮",
            f"      │    pDBS Custom     │",
            f"      │     ~{size//1000}.{(size%1000)//100}kb        │",
            f"  ori─┤                   ├─ MCS/Insert",
            f"      │                   │",
            "      ╰───╮            ╭───╯",
            f"          ╰── {plan['promoter'][0]} / {plan['terminator'][0][:4]} ──╯",
        ]
        return "\n".join(f"[dim]{l}[/dim]" for l in lines)


# ════════════════════════════════════════════════════════════════
#  ── PHASE 2 · A  ──  siRNA / shRNA DESIGNER
# ════════════════════════════════════════════════════════════════

class siRNADesigner:
    """
    Designs 21-nt siRNA duplexes from an mRNA/DNA target.
    Applies Reynolds (2004) 8-criterion rational design rules.
    Also generates shRNA hairpin sequences for vector cloning.
    """

    REYNOLDS_WEIGHTS = {
        "A_at_pos19"      : 1,   # A at position 19 of sense
        "A_at_pos3"       : 1,   # A at position 3
        "U_at_pos10"      : 1,   # U at position 10
        "no_G_at_pos19"   : 1,   # not G at 19
        "no_C_at_pos19"   : 1,   # not C at 19
        "no_G_stretch"    : 1,   # no run ≥ 4G
        "gc_between_4_14" : 1,   # low GC in 4-14 window
        "low_5prime_dg"   : 1,   # 5' end less stable
    }

    HAIRPIN_LOOP = "TTCAAGAGA"   # common shRNA loop sequence

    def __init__(self, sequence: str, guide_len: int = 21):
        self.seq       = sequence.upper().replace("U","T")
        self.guide_len = guide_len

    def _score_siRNA(self, sense: str) -> int:
        score = 0
        rna   = dna_to_rna(sense)
        if rna[18]  == "A": score += 1
        if rna[2]   == "A": score += 1
        if rna[9]   == "U": score += 1
        if rna[18] != "G":  score += 1
        if rna[18] != "C":  score += 1
        if "GGGG"  not in rna: score += 1
        mid_gc = gc_content(sense[4:14])
        if mid_gc < 52:     score += 1
        t5_gc  = gc_content(sense[:7])
        t3_gc  = gc_content(sense[14:])
        if t5_gc < t3_gc:   score += 1
        return score

    def design(self, top_n: int = 8) -> dict:
        results = []
        for i in range(0, len(self.seq) - self.guide_len):
            sense     = self.seq[i:i+self.guide_len]
            if len(sense) < self.guide_len: break
            if set(sense) - set("ACGT"):    continue
            antisense = reverse_complement(sense)
            score     = self._score_siRNA(sense)
            gc        = gc_content(sense)
            if not (30 <= gc <= 65): continue     # hard GC filter
            # Build two-base 3' overhang (TT standard)
            sense_overhang     = sense     + "TT"
            antisense_overhang = antisense + "TT"
            # shRNA
            shrna = sense + self.HAIRPIN_LOOP + antisense
            results.append({
                "sense"      : sense,
                "antisense"  : antisense,
                "sense_rna"  : dna_to_rna(sense_overhang),
                "anti_rna"   : dna_to_rna(antisense_overhang),
                "shrna"      : shrna,
                "position"   : i+1,
                "gc"         : round(gc,1),
                "score"      : score,
                "max_score"  : 8,
            })
        results.sort(key=lambda x: -x["score"])
        return {"candidates": results[:top_n], "total_found": len(results)}

    @staticmethod
    def display_results(data: dict):
        cands = data["candidates"]
        if not cands:
            console.print("[red]No valid siRNAs found.[/red]"); return
        tbl = Table(title="[bold magenta]siRNA DESIGN RESULTS[/bold magenta]",
                    box=box.ROUNDED, show_lines=True, header_style="bold magenta")
        for col in ["#","Sense strand (5'→3'  +TT)","Pos","GC%","Reynolds Score"]:
            tbl.add_column(col, justify="center")
        for i, c in enumerate(cands, 1):
            sc = "[bold green]" if c["score"]>=6 else "[bold yellow]" if c["score"]>=4 else "[bold red]"
            ec = sc.replace("[","").replace("]","")
            tbl.add_row(str(i),f"[magenta]{c['sense_rna']}[/magenta]",
                        str(c["position"]),str(c["gc"]),
                        f"{sc}{c['score']}/8[/{ec}]")
        console.print(tbl)
        best = cands[0]
        console.print(Panel(
            f"[bold]Best siRNA (Reynolds score {best['score']}/8)[/bold]\n\n"
            f"  Sense    5'-[magenta]{best['sense_rna']}[/magenta]-3'\n"
            f"  Antisense5'-[cyan]{best['anti_rna']}[/cyan]-3'\n\n"
            f"  [bold]shRNA hairpin (for U6 vector cloning):[/bold]\n"
            f"  5'-[green]{best['shrna']}[/green]-3'\n\n"
            f"  GC%%      : {best['gc']}%%\n"
            f"  Position  : {best['position']}",
            title="[bold magenta]★ Best siRNA[/bold magenta]",
            border_style="magenta"))


# ════════════════════════════════════════════════════════════════
#  ── PHASE 2 · B  ──  RNA SECONDARY STRUCTURE PREDICTOR
# ════════════════════════════════════════════════════════════════

class RNAStructurePredictor:
    """
    Predicts minimum free energy (MFE) RNA secondary structure.
    Uses Nussinov base-pair maximization + simplified Turner energetics.
    Outputs ASCII stem-loop representation.
    """

    MIN_LOOP_SIZE = 3

    def __init__(self, sequence: str):
        self.rna = dna_to_rna(sequence.upper().replace("T","U"))
        self.n   = len(self.rna)

    def _can_pair(self, a: str, b: str) -> bool:
        return (a,b) in {("A","U"),("U","A"),("G","C"),("C","G"),
                         ("G","U"),("U","G")}

    def _fill_dp(self) -> list:
        """Nussinov dynamic programming — maximize base pairs."""
        n   = self.n
        dp  = [[0]*n for _ in range(n)]
        bt  = [[None]*n for _ in range(n)]
        for span in range(self.MIN_LOOP_SIZE+1, n):
            for i in range(n-span):
                j = i + span
                best = dp[i][j-1]
                move = ("del_j", None)
                if self._can_pair(self.rna[i], self.rna[j]):
                    v = dp[i+1][j-1] + 1
                    if v > best:
                        best = v; move = ("pair", (i,j))
                for k in range(i, j):
                    v = dp[i][k] + dp[k+1][j]
                    if v > best:
                        best = v; move = ("split", k)
                dp[i][j] = best
                bt[i][j] = move
        return dp, bt

    def _traceback(self, bt, i, j, pairs):
        if i >= j: return
        move = bt[i][j]
        if move is None: return
        kind = move[0]
        if kind == "del_j":
            self._traceback(bt, i, j-1, pairs)
        elif kind == "pair":
            pairs.append((i, j))
            self._traceback(bt, i+1, j-1, pairs)
        elif kind == "split":
            k = move[1]
            self._traceback(bt, i, k, pairs)
            self._traceback(bt, k+1, j, pairs)

    def predict(self) -> dict:
        if self.n < 6:
            return {"rna":self.rna,"pairs":[],"mfe":0.0,
                    "structure":"."*self.n,"dot_bracket":"."*self.n}
        _, bt = self._fill_dp()
        pairs = []
        self._traceback(bt, 0, self.n-1, pairs)
        pairs.sort()
        # Build dot-bracket
        db = ["."] * self.n
        for i, j in pairs:
            db[i] = "("; db[j] = ")"
        # Estimate MFE (simplified Turner)
        mfe = 0.0
        for i, j in pairs:
            mfe += rna_stem_energy(self.rna[i], self.rna[j])
        return {"rna"        : self.rna,
                "pairs"      : pairs,
                "mfe"        : round(mfe, 2),
                "dot_bracket": "".join(db),
                "num_pairs"  : len(pairs),
                "structure"  : self._ascii_structure(pairs)}

    def _ascii_structure(self, pairs: list) -> str:
        """Generate a simple ASCII stem-loop representation."""
        paired = {i for p in pairs for i in p}
        lines  = []
        lines.append("5'-" + self.rna + "-3'")
        lines.append("   " + "".join("|" if i in paired else " "
                                     for i in range(self.n)))
        # show pair numbers underneath for stems
        lines.append("   " + "".join(
            str(self.n - 1 - i)[-1] if i in paired else " "
            for i in range(self.n)))
        return "\n".join(lines)

    @staticmethod
    def display_results(result: dict):
        mfe_color = "green" if result["mfe"] < -5 else "yellow" if result["mfe"] < 0 else "red"
        console.print(Panel(
            f"[bold]RNA sequence ({len(result['rna'])} nt):[/bold]\n"
            f"[white]{result['rna']}[/white]\n\n"
            f"[bold]Dot-bracket notation:[/bold]\n"
            f"[cyan]{result['dot_bracket']}[/cyan]\n\n"
            f"[bold]MFE:[/bold] [{mfe_color}]{result['mfe']} kcal/mol[/{mfe_color}]\n"
            f"[bold]Base pairs formed:[/bold] {result['num_pairs']}\n\n"
            f"[bold]Structure:[/bold]\n[dim]{result['structure']}[/dim]",
            title="[bold cyan]RNA Secondary Structure[/bold cyan]",
            border_style="cyan"))
        if result["pairs"]:
            tbl = Table(title="[bold cyan]Base Pairs[/bold cyan]",
                        box=box.SIMPLE, header_style="bold cyan")
            tbl.add_column("Pair #"); tbl.add_column("Position i")
            tbl.add_column("Base i"); tbl.add_column("Position j")
            tbl.add_column("Base j"); tbl.add_column("Type")
            for n, (i,j) in enumerate(result["pairs"][:20], 1):
                bi = result["rna"][i]; bj = result["rna"][j]
                pt = "Watson-Crick" if (bi,bj) in {("A","U"),("U","A"),("G","C"),("C","G")} else "Wobble"
                tbl.add_row(str(n),str(i+1),bi,str(j+1),bj,pt)
            console.print(tbl)


# ════════════════════════════════════════════════════════════════
#  ── PHASE 2 · C  ──  PHYLOGENETICS BUILDER
# ════════════════════════════════════════════════════════════════

class PhylogeneticsBuilder:
    """
    Multiple sequence alignment (progressive) + Neighbor-Joining tree.
    Calculates Jukes-Cantor evolutionary distance.
    Draws ASCII phylogenetic tree with bootstrap labels.
    """

    def __init__(self, sequences: dict):
        """sequences: {name: dna_string}"""
        self.seqs = {k: v.upper() for k,v in sequences.items()}
        self.names = list(self.seqs.keys())

    # ── Pairwise alignment (Needleman-Wunsch) ──
    def _align_pair(self, s1: str, s2: str,
                    match=2, mismatch=-1, gap=-2) -> tuple:
        m, n   = len(s1), len(s2)
        dp     = [[0]*(n+1) for _ in range(m+1)]
        tb     = [[""]*(n+1) for _ in range(m+1)]
        for i in range(1, m+1): dp[i][0] = dp[i-1][0] + gap; tb[i][0]="U"
        for j in range(1, n+1): dp[0][j] = dp[0][j-1] + gap; tb[0][j]="L"
        for i in range(1, m+1):
            for j in range(1, n+1):
                sc = match if s1[i-1]==s2[j-1] else mismatch
                opts = [(dp[i-1][j-1]+sc,"D"),(dp[i-1][j]+gap,"U"),(dp[i][j-1]+gap,"L")]
                best = max(opts, key=lambda x: x[0])
                dp[i][j] = best[0]; tb[i][j] = best[1]
        a1,a2 = [],[]
        i,j   = m, n
        while i>0 or j>0:
            d = tb[i][j]
            if d=="D": a1.append(s1[i-1]); a2.append(s2[j-1]); i-=1; j-=1
            elif d=="U": a1.append(s1[i-1]); a2.append("-"); i-=1
            else:        a1.append("-"); a2.append(s2[j-1]); j-=1
        return "".join(reversed(a1)), "".join(reversed(a2))

    def _jc_distance(self, s1: str, s2: str) -> float:
        """Jukes-Cantor distance between two aligned sequences."""
        pairs = [(a,b) for a,b in zip(s1,s2) if a!="-" and b!="-"]
        if not pairs: return 1.0
        p = sum(1 for a,b in pairs if a!=b) / len(pairs)
        p = min(p, 0.74)    # avoid log(0)
        try:
            return -0.75 * math.log(1 - 4*p/3)
        except ValueError:
            return 1.0

    def build_distance_matrix(self) -> dict:
        names = self.names
        dist  = {n:{m:0.0 for m in names} for n in names}
        for i, n1 in enumerate(names):
            for n2 in names[i+1:]:
                a1, a2 = self._align_pair(self.seqs[n1][:300], self.seqs[n2][:300])
                d      = self._jc_distance(a1, a2)
                dist[n1][n2] = dist[n2][n1] = d
        return dist

    def neighbor_joining(self, dist: dict) -> dict:
        """Classic Neighbor-Joining algorithm (Saitou & Nei, 1987)."""
        names = list(dist.keys())
        if len(names) < 2:
            return {}

        # Work on a copy
        d = {n:{m:dist[n][m] for m in dist[n]} for n in dist}
        nodes  = {n:{"name":n,"children":[],"length":0.0} for n in names}
        n_taxa = len(names)
        node_counter = 0

        while len(names) > 2:
            n = len(names)
            # Q-matrix
            totals = {ni: sum(d[ni][nj] for nj in names if nj!=ni) for ni in names}
            min_q  = float("inf"); best_pair = (None,None)
            for i in range(len(names)):
                for j in range(i+1, len(names)):
                    ni,nj = names[i], names[j]
                    q = (n-2)*d[ni][nj] - totals[ni] - totals[nj]
                    if q < min_q: min_q=q; best_pair=(ni,nj)

            ni, nj = best_pair
            # Branch lengths
            dij = d[ni][nj]
            li  = dij/2 + (totals[ni]-totals[nj]) / (2*(n-2))
            lj  = dij - li
            # New node
            node_counter += 1
            new_name = f"node{node_counter}"
            new_node = {"name":new_name,"children":[
                {**nodes[ni],"length":round(li,4)},
                {**nodes[nj],"length":round(lj,4)}
            ],"length":0.0}
            nodes[new_name] = new_node

            # Update distances
            remaining = [x for x in names if x not in (ni,nj)]
            d[new_name] = {}
            for nk in remaining:
                d[new_name][nk] = (d[ni][nk]+d[nj][nk]-dij)/2
                d[nk][new_name] = d[new_name][nk]
            for x in (ni,nj): del d[x]; [d[k].pop(x,None) for k in d]
            names = remaining + [new_name]

        # Connect final two nodes
        if len(names)==2:
            node_counter += 1
            root = {"name":f"root{node_counter}","children":[
                {**nodes[names[0]],"length":round(d[names[0]][names[1]]/2,4)},
                {**nodes[names[1]],"length":round(d[names[0]][names[1]]/2,4)},
            ],"length":0.0}
        else:
            root = nodes[names[0]]
        return root

    def _draw_tree(self, node: dict, prefix: str = "",
                   is_last: bool = True) -> list:
        lines    = []
        connector = "└── " if is_last else "├── "
        name      = node.get("name","")
        length    = node.get("length", 0.0)
        label     = f"{name}  [dim][:{length:.4f}][/dim]"
        lines.append(prefix + connector + label)
        children = node.get("children",[])
        for i, child in enumerate(children):
            ext = "    " if is_last else "│   "
            lines.extend(self._draw_tree(
                child, prefix+ext, i==len(children)-1))
        return lines

    @staticmethod
    def display_results(tree: dict, dist: dict, names: list):
        # Distance matrix table
        tbl = Table(title="[bold green]DISTANCE MATRIX (Jukes-Cantor)[/bold green]",
                    box=box.SIMPLE, header_style="bold green")
        tbl.add_column("Taxon", style="bold")
        for n in names: tbl.add_column(n[:12], justify="center")
        for n1 in names:
            row = [n1[:12]] + [f"{dist[n1][n2]:.4f}" for n2 in names]
            tbl.add_row(*row)
        console.print(tbl)
        # Tree
        builder = PhylogeneticsBuilder({})
        tree_lines = builder._draw_tree(tree, "", True)
        tree_str   = "\n".join(tree_lines)
        console.print(Panel(
            tree_str,
            title="[bold green]Neighbor-Joining Phylogenetic Tree[/bold green]",
            border_style="green"))
        console.print(Panel(
            "[dim]Branch lengths are Jukes-Cantor substitutions per site.\n"
            "Algorithm: Neighbor-Joining (Saitou & Nei, 1987)\n"
            "Distance model: Jukes-Cantor (1969)[/dim]",
            title="[bold green]Tree Notes[/bold green]", border_style="dim"))


# ════════════════════════════════════════════════════════════════
#  ── PHASE 3 · A  ──  VACCINE EPITOPE DESIGNER
# ════════════════════════════════════════════════════════════════

class VaccineEpitopeDesigner:
    """
    Predicts B-cell and T-cell epitopes from a protein/DNA sequence.
    Designs a multi-epitope vaccine construct with linkers.
    Scores antigenicity using Vaxijen-inspired propensity model.
    """

    CTL_LEN   = 9    # MHC-I (CD8+) epitope length
    HTL_LEN   = 15   # MHC-II (CD4+) epitope length
    BCELL_LEN = 12   # B-cell linear epitope length

    CTL_LINKER   = "AAY"
    HTL_LINKER   = "GPGPG"
    BCELL_LINKER = "KK"

    def __init__(self, sequence: str, input_type: str = "protein"):
        if input_type == "dna":
            self.prot = translate(sequence)
        else:
            self.prot = sequence.upper()
        self.prot = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]","",self.prot)

    def _antigenicity(self, peptide: str) -> float:
        """
        Vaxijen-inspired score: sum of auto-cross-covariance
        descriptors (simplified via amino acid propensity table).
        """
        phy = {
            "A":0.62,"R":-2.53,"N":-0.78,"D":-0.90,"C":0.29,
            "Q":-0.85,"E":-0.74,"G":0.48,"H":-0.40,"I":1.38,
            "L":1.06,"K":-1.50,"M":0.64,"F":1.19,"P":0.12,
            "S":-0.18,"T":-0.05,"V":1.08,"W":0.81,"Y":0.26,
        }
        scores = [phy.get(aa, 0.0) for aa in peptide]
        return round(sum(scores)/len(scores) + 0.5, 3) if scores else 0.0

    def _predict_ctl(self, allele: str = "HLA-A*02:01") -> list:
        motif   = MHC_I_MOTIFS.get(allele, {})
        results = []
        for i in range(len(self.prot) - self.CTL_LEN + 1):
            pep  = self.prot[i:i+self.CTL_LEN]
            ok   = all(pep[p-1] in res
                       for p,res in motif.items()
                       if p <= len(pep))
            score = self._antigenicity(pep) + (0.3 if ok else 0.0)
            results.append({"peptide":pep,"position":i+1,
                            "score":round(score,3),"allele":allele,"anchor":ok})
        return sorted(results, key=lambda x:-x["score"])

    def _predict_htl(self) -> list:
        """MHC-II (HTL): sliding 15-mer, score by antigenicity + hydrophilicity."""
        results = []
        for i in range(len(self.prot) - self.HTL_LEN + 1):
            pep   = self.prot[i:i+self.HTL_LEN]
            score = self._antigenicity(pep)
            hyd   = sum(HYDROPHOBICITY.get(aa,0) for aa in pep) / self.HTL_LEN
            score += (0.2 if -0.5 < hyd < 0.5 else -0.1)
            results.append({"peptide":pep,"position":i+1,"score":round(score,3)})
        return sorted(results, key=lambda x:-x["score"])

    def _predict_bcell(self) -> list:
        """B-cell linear epitopes: Parker hydrophilicity-based."""
        parker = {"A":-0.5,"R":3.0,"N":0.2,"D":3.0,"C":-1.0,
                  "Q":0.2,"E":3.0,"G":0.0,"H":-0.5,"I":-1.8,
                  "L":-1.8,"K":3.0,"M":-1.3,"F":-2.5,"P":0.0,
                  "S":0.3,"T":-0.4,"V":-1.5,"W":-3.4,"Y":-2.3}
        results = []
        for i in range(len(self.prot) - self.BCELL_LEN + 1):
            pep   = self.prot[i:i+self.BCELL_LEN]
            score = sum(parker.get(aa,0) for aa in pep) / self.BCELL_LEN
            results.append({"peptide":pep,"position":i+1,"score":round(score,3)})
        return sorted(results, key=lambda x:-x["score"])

    def _allergenicity(self, pep: str) -> str:
        """
        Simplified allergenicity check based on known IgE-binding motifs.
        Returns 'Non-allergen' or 'Putative allergen'.
        """
        allergen_motifs = ["VHFHV","GDSVHH","VRIVI","QLQSF","ILLPK"]
        return ("Putative allergen" if any(m in pep for m in allergen_motifs)
                else "Non-allergen")

    def design(self, allele: str = "HLA-A*02:01") -> dict:
        ctl   = self._predict_ctl(allele)[:3]
        htl   = self._predict_htl()[:2]
        bcell = self._predict_bcell()[:2]

        construct_parts = (
            [c["peptide"] for c in ctl]  +
            [self.CTL_LINKER] +
            [h["peptide"] for h in htl]  +
            [self.HTL_LINKER] +
            [b["peptide"] for b in bcell] +
            [self.BCELL_LINKER]
        )
        construct = "".join(construct_parts)
        antigen   = self._antigenicity(construct)
        allergy   = self._allergenicity(construct)

        return {
            "ctl"        : ctl,
            "htl"        : htl,
            "bcell"      : bcell,
            "construct"  : construct,
            "antigenicity": antigen,
            "allergenicity": allergy,
            "allele"     : allele,
            "length"     : len(construct),
        }

    @staticmethod
    def display_results(data: dict):
        console.print(Panel(
            f"[bold]MHC-I allele:[/bold] [cyan]{data['allele']}[/cyan]\n"
            f"[bold]Construct length:[/bold] {data['length']} aa\n"
            f"[bold]Antigenicity:[/bold] [{'green' if data['antigenicity']>0.5 else 'red'}]{data['antigenicity']}[/{'green' if data['antigenicity']>0.5 else 'red'}]  [dim](threshold >0.5)[/dim]\n"
            f"[bold]Allergenicity:[/bold] [{'green' if 'Non' in data['allergenicity'] else 'red'}]{data['allergenicity']}[/{'green' if 'Non' in data['allergenicity'] else 'red'}]",
            title="[bold yellow]Vaccine Construct Overview[/bold yellow]",
            border_style="yellow"))

        for section, color, label in [
            ("ctl","cyan","CTL (MHC-I / CD8+) Epitopes"),
            ("htl","green","HTL (MHC-II / CD4+) Epitopes"),
            ("bcell","magenta","B-cell Linear Epitopes"),
        ]:
            tbl = Table(title=f"[bold {color}]{label}[/bold {color}]",
                        box=box.SIMPLE, header_style=f"bold {color}")
            tbl.add_column("#"); tbl.add_column("Peptide"); tbl.add_column("Pos"); tbl.add_column("Score")
            for i, ep in enumerate(data[section], 1):
                tbl.add_row(str(i),f"[{color}]{ep['peptide']}[/{color}]",
                            str(ep["position"]),str(ep["score"]))
            console.print(tbl)

        console.print(Panel(
            f"[bold]Multi-epitope vaccine construct:[/bold]\n\n"
            f"[cyan]{'─'*3}CTL{'─'*3}[/cyan] [yellow]AAY[/yellow] "
            f"[green]{'─'*3}HTL{'─'*3}[/green] [yellow]GPGPG[/yellow] "
            f"[magenta]{'─'*2}B-cell{'─'*2}[/magenta] [yellow]KK[/yellow]\n\n"
            f"[white]{data['construct']}[/white]",
            title="[bold yellow]★ Vaccine Construct[/bold yellow]",
            border_style="yellow"))


# ════════════════════════════════════════════════════════════════
#  ── PHASE 3 · B  ──  AMP DESIGNER
# ════════════════════════════════════════════════════════════════

class AMPDesigner:
    """
    Evaluates or designs Antimicrobial Peptides (AMPs).
    Calculates amphipathicity, Boman index, net charge.
    Draws ASCII helical wheel.
    """

    # Boman index (protein-binding potential, kcal/mol)
    BOMAN = {
        "A":0.17,"R":0.81,"N":0.42,"D":1.23,"C":-0.24,"Q":0.58,
        "E":2.02,"G":0.01,"H":0.96,"I":-1.56,"L":-1.27,"K":0.99,
        "M":-0.23,"F":-1.71,"P":0.45,"S":0.13,"T":0.14,"V":-0.97,
        "W":-2.09,"Y":-0.71,
    }

    CHARGE_TABLE = {
        "R":1,"K":1,"H":0.5,"D":-1,"E":-1,
    }

    def __init__(self, peptide: str):
        self.pep = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]","",peptide.upper())

    def analyze(self) -> dict:
        pep    = self.pep
        gc_net = sum(self.CHARGE_TABLE.get(aa, 0) for aa in pep)
        hydro  = sum(HYDROPHOBICITY.get(aa, 0) for aa in pep) / len(pep)
        boman  = sum(self.BOMAN.get(aa, 0) for aa in pep) / len(pep)

        # Amphipathicity — helical wheel method (3.6 residues per turn)
        angle_step = 100  # degrees per residue
        x_sum = sum(HYDROPHOBICITY.get(aa,0) * math.cos(math.radians(i*angle_step))
                    for i,aa in enumerate(pep))
        y_sum = sum(HYDROPHOBICITY.get(aa,0) * math.sin(math.radians(i*angle_step))
                    for i,aa in enumerate(pep))
        amphip = math.sqrt(x_sum**2 + y_sum**2) / len(pep)

        # AMP potential scoring (APD3-inspired)
        amp_score = 0
        if len(pep) >= 10:      amp_score += 20
        if gc_net >= 2:         amp_score += 20
        if 0.3 <= hydro <= 0.6: amp_score += 20
        if amphip >= 0.4:       amp_score += 20
        if boman <= 0.0:        amp_score += 20
        if len(pep) <= 50:      amp_score += 10

        return {
            "peptide"       : pep,
            "length"        : len(pep),
            "net_charge"    : round(gc_net, 1),
            "hydrophobicity": round(hydro, 3),
            "amphipathicity": round(amphip, 3),
            "boman_index"   : round(boman, 3),
            "amp_score"     : amp_score,
            "wheel"         : self._helical_wheel(pep),
        }

    def _helical_wheel(self, pep: str, n: int = 18) -> str:
        """Generate a simplified ASCII helical wheel projection."""
        wheel = pep[:min(n, len(pep))]
        rows  = ["   Helical Wheel (N→C, top→right, 100° rotation)"]
        rows.append("   " + "─"*30)
        # Place residues on a 3-row ASCII grid
        positions = []
        for i, aa in enumerate(wheel):
            angle = math.radians(i * 100 % 360)
            x     = int(14 + 10*math.cos(angle))
            y     = int(4  +  3*math.sin(angle))
            positions.append((x, y, aa))
        grid = [[" "]*30 for _ in range(9)]
        for x, y, aa in positions:
            if 0<=y<9 and 0<=x<30:
                grid[y][x] = aa
        for row in grid:
            rows.append("   " + "".join(row))
        return "\n".join(rows)

    @staticmethod
    def display_results(data: dict):
        sc = data["amp_score"]
        sc_color = "green" if sc>=80 else "yellow" if sc>=50 else "red"
        console.print(Panel(
            f"[bold]Peptide  :[/bold] [cyan]{data['peptide']}[/cyan]\n"
            f"[bold]Length   :[/bold] {data['length']} aa\n"
            f"[bold]Net charge:[/bold] [{'green' if data['net_charge']>=2 else 'red'}]+{data['net_charge']}[/{'green' if data['net_charge']>=2 else 'red'}]  [dim](≥+2 ideal for AMPs)[/dim]\n"
            f"[bold]Hydrophobicity:[/bold] {data['hydrophobicity']}  [dim](0.3–0.6 ideal)[/dim]\n"
            f"[bold]Amphipathicity:[/bold] [{'green' if data['amphipathicity']>=0.4 else 'red'}]{data['amphipathicity']}[/{'green' if data['amphipathicity']>=0.4 else 'red'}]  [dim](≥0.4 ideal)[/dim]\n"
            f"[bold]Boman index:[/bold]   {data['boman_index']}  [dim](≤0 = membrane-active)[/dim]\n\n"
            f"[bold]AMP Score  :[/bold] [{sc_color}]{sc}/110[/{sc_color}]\n\n"
            f"{data['wheel']}",
            title="[bold magenta]★ AMP Analysis[/bold magenta]",
            border_style="magenta"))


# ════════════════════════════════════════════════════════════════
#  ── PHASE 3 · C  ──  ARG SCANNER
# ════════════════════════════════════════════════════════════════

class ARGScanner:
    """
    Scans a DNA/protein sequence for known Antibiotic Resistance Genes.
    Uses conserved motif database aligned to WHO priority drug classes.
    """

    def __init__(self, sequence: str, input_type: str = "dna"):
        if input_type == "dna":
            self.prot = translate(sequence)
            self.dna  = sequence.upper()
        else:
            self.prot = sequence.upper()
            self.dna  = ""
        self.input_type = input_type

    def _motif_similarity(self, query: str, motif: str) -> float:
        """
        Simple sliding-window identity score.
        Returns best fraction of motif found in query.
        """
        if not motif or not query: return 0.0
        best = 0.0
        mlen = min(len(motif), len(query))
        for i in range(len(query) - mlen + 1):
            window = query[i:i+mlen]
            match  = sum(a==b for a,b in zip(window, motif[:mlen]))
            best   = max(best, match/mlen)
        return round(best, 3)

    def scan(self, threshold: float = 0.60) -> list:
        """
        Scan both DNA translation and protein directly.
        Returns list of hits above threshold.
        """
        hits = []
        for gene, (drug_class, mechanism, motif) in ARG_DATABASE.items():
            sim = self._motif_similarity(self.prot, motif)
            if sim >= threshold:
                hits.append({
                    "gene"       : gene,
                    "drug_class" : drug_class,
                    "mechanism"  : mechanism,
                    "similarity" : sim,
                    "risk"       : ("CRITICAL" if sim>=0.90 else
                                    "HIGH"     if sim>=0.75 else "MODERATE"),
                })
        hits.sort(key=lambda x: -x["similarity"])
        return hits

    @staticmethod
    def display_results(hits: list):
        if not hits:
            console.print(Panel(
                "[bold green]✔ No known resistance genes detected.[/bold green]\n"
                "[dim]Sequence does not match ARG motifs at ≥60% identity.[/dim]",
                title="[bold green]ARG Scan – Clean[/bold green]",
                border_style="green"))
            return
        risk_color = {"CRITICAL":"bold red","HIGH":"bold orange3","MODERATE":"yellow"}
        tbl = Table(
            title=f"[bold red]⚠ ANTIBIOTIC RESISTANCE GENES DETECTED ({len(hits)})[/bold red]",
            box=box.ROUNDED, show_lines=True, header_style="bold red")
        for col in ["Gene","Drug Class","Mechanism","Identity","Risk"]:
            tbl.add_column(col, justify="center")
        for h in hits:
            rc = risk_color.get(h["risk"],"white")
            tbl.add_row(
                f"[bold]{h['gene']}[/bold]",
                f"[cyan]{h['drug_class']}[/cyan]",
                h["mechanism"],
                f"[yellow]{h['similarity']*100:.1f}%[/yellow]",
                f"[{rc}]{h['risk']}[/{rc}]")
        console.print(tbl)

        # WHO priority alert
        critical = [h for h in hits if h["risk"]=="CRITICAL"]
        if critical:
            console.print(Panel(
                "[bold red]⚠ CRITICAL RESISTANCE GENES FOUND[/bold red]\n\n"
                + "\n".join(f"  • [red]{h['gene']}[/red] – {h['drug_class']} ({h['mechanism']})"
                            for h in critical) +
                "\n\n[dim]These genes are on the WHO critical priority pathogen resistance list.\n"
                "Notify infection control if this is a clinical isolate.[/dim]",
                title="[bold red]WHO Priority Alert[/bold red]",
                border_style="red"))