# ================================================================
#  DBS  –  Digital BioScience Suite
#  main.py  |  Terminal GUI  –  Main menu and module router
#
#  Phase 1: Primer Design, CRISPR, Restriction Mapping,
#            Codon Optimization, Plasmid Design
#
#  Phase 2: siRNA/shRNA Design, RNA Secondary Structure,
#            Phylogenetics
#
#  Phase 3: Vaccine Epitope Design, AMP Analysis, ARG Scanner
#
#  Install: pip install rich
#  Run    : python main.py
# ================================================================

import sys, time

from rich.console import Console
from rich.panel   import Panel
from rich.prompt  import Prompt, Confirm
from rich.table   import Table
from rich.align   import Align
from rich        import box

from utilities import get_sequence_from_user, truncate, CRISPR_SYSTEMS, CODON_TABLES
from structure  import (
    # Phase 1
    PrimerDesigner, CRISPRDesigner, CRISPRAnimator,
    RestrictionMapper, CodonOptimizer, PlasmidDesigner,
    # Phase 2
    siRNADesigner, RNAStructurePredictor, PhylogeneticsBuilder,
    # Phase 3
    VaccineEpitopeDesigner, AMPDesigner, ARGScanner,
)

console = Console()

# ──────────────────────────────────────────────────────────────
DBS_LOGO = """
[bold bright_cyan]
  ██████╗ ██████╗ ███████╗
  ██╔══██╗██╔══██╗██╔════╝
  ██║  ██║██████╔╝███████╗
  ██║  ██║██╔══██╗╚════██║
  ██████╔╝██████╔╝███████║
  ╚═════╝ ╚═════╝ ╚══════╝
[/bold bright_cyan]
[bold white]  Digital BioScience Suite[/bold white]
[dim]  Three-phase in silico biology platform for microbiology research[/dim]
"""
# ──────────────────────────────────────────────────────────────


def splash():
    console.clear()
    console.print(DBS_LOGO)
    console.print(Align.center(
        "[dim]Phase 1: Primer · CRISPR · Restriction · Codon Opt · Plasmid[/dim]\n"
        "[dim]Phase 2: siRNA  · RNA Structure · Phylogenetics            [/dim]\n"
        "[dim]Phase 3: Vaccine · AMP · ARG Scanner                       [/dim]"
    ))
    console.print()
    time.sleep(0.4)


# ════════════════════════════════════════════════════════════════
#  ── PHASE 1 MODULE RUNNERS ──
# ════════════════════════════════════════════════════════════════

def run_primer_design():
    console.print(Panel(
        "[bold cyan]PRIMER DESIGN[/bold cyan]\n"
        "[dim]SantaLucia NN thermodynamics · GC-clamp · ΔTm pairing[/dim]",
        border_style="cyan"))
    try:
        seq, rid = get_sequence_from_user("primer target")
    except Exception as e:
        console.print(f"[red]{e}[/red]"); return
    console.print(f"\n[green]✔[/green] [bold]{rid}[/bold] [dim]({len(seq)} bp)[/dim]")
    try:
        min_len = int(Prompt.ask("Min primer length", default="18"))
        max_len = int(Prompt.ask("Max primer length", default="25"))
        tm_min  = float(Prompt.ask("Min Tm °C",       default="52"))
        tm_max  = float(Prompt.ask("Max Tm °C",       default="68"))
        top_n   = int(Prompt.ask("Candidates shown",  default="5"))
    except ValueError:
        min_len, max_len, tm_min, tm_max, top_n = 18, 25, 52, 68, 5
    with console.status("[bold green]Designing...[/bold green]", spinner="dots"):
        results = PrimerDesigner(seq, min_len=min_len, max_len=max_len,
                                  tm_min=tm_min, tm_max=tm_max,
                                  top_n=top_n).design()
    PrimerDesigner.display_results(results, seq)
    if results["pairs"] and Confirm.ask("\nSave results?", default=False):
        _save_text(f"primers_{rid[:15]}.txt",
                   "\n".join(f"Fwd: {p['fwd']['primer']}\nRev: {p['rev']['primer']}\n"
                             for p in results["pairs"]))


def run_crispr():
    console.print(Panel(
        "[bold green]CRISPR GUIDE DESIGN + ANIMATION[/bold green]\n"
        "[dim]SpCas9 · SaCas9 · Cas12a · CjCas9 · NmCas9[/dim]",
        border_style="green"))
    # System selection
    tbl = Table(title="[bold green]CRISPR Systems[/bold green]",
                box=box.ROUNDED, header_style="bold green", show_lines=True)
    for col in ["#","System","PAM","Guide (nt)","Best for"]:
        tbl.add_column(col, justify="center")
    for k, s in CRISPR_SYSTEMS.items():
        tbl.add_row(k, s["name"], s["pam"], str(s["guide_len"]), s["organisms"])
    console.print(tbl)
    sk = Prompt.ask("Select system", choices=list(CRISPR_SYSTEMS.keys()), default="1")
    system = CRISPR_SYSTEMS[sk]
    try:
        seq, rid = get_sequence_from_user("CRISPR target")
    except Exception as e:
        console.print(f"[red]{e}[/red]"); return
    with console.status("[bold green]Scanning...[/bold green]", spinner="dots"):
        guides = CRISPRDesigner(seq, sk).find_guides()
    console.print(f"[green]✔ {len(guides)} guide(s) found[/green]")
    top_n = int(Prompt.ask("Top guides to show", default="10"))
    CRISPRDesigner.display_results(guides, system, top_n)
    if not guides: return
    if Confirm.ask("Run CRISPR animation?", default=True):
        try:
            idx = int(Prompt.ask("Animate guide #", default="1")) - 1
        except ValueError:
            idx = 0
        speed_map = {"1":0.04,"2":0.08,"3":0.15}
        console.print("Speed: [1] Fast  [2] Normal  [3] Slow")
        spd = speed_map.get(Prompt.ask("Speed",choices=["1","2","3"],default="2"),0.08)
        time.sleep(1.5)
        try:
            CRISPRAnimator(seq, guides[min(idx,len(guides)-1)], system).run(spd)
        except KeyboardInterrupt:
            console.print("\n[yellow]Animation cancelled.[/yellow]")
    if Confirm.ask("Save guide results?", default=False):
        lines = [f"{g['guide']}\t{g['pam']}\t{g['strand']}\t{g['score']}"
                 for g in guides[:top_n]]
        _save_text(f"crispr_{rid[:15]}.txt", "\n".join(lines))


def run_restriction():
    console.print(Panel(
        "[bold cyan]RESTRICTION ENZYME MAPPING[/bold cyan]\n"
        "[dim]REBASE-derived database · virtual digest · ASCII map[/dim]",
        border_style="cyan"))
    try:
        seq, rid = get_sequence_from_user("restriction target")
    except Exception as e:
        console.print(f"[red]{e}[/red]"); return
    mapper = RestrictionMapper(seq)
    console.print("\n[1] Scan all enzymes")
    console.print("[2] Scan specific enzymes")
    choice = Prompt.ask("Choice", choices=["1","2"], default="1")
    if choice == "2":
        names_input = Prompt.ask("Enter enzyme names (comma-separated)")
        enzymes = [e.strip() for e in names_input.split(",")]
    else:
        enzymes = None
    with console.status("[bold cyan]Scanning...[/bold cyan]", spinner="dots"):
        hits = mapper.find_sites(enzymes)
    map_str = mapper.draw_map(hits)
    # Optional virtual digest
    frags = None
    if hits and Confirm.ask("Run virtual digest?", default=True):
        e1 = Prompt.ask("Enzyme 1")
        if Confirm.ask("Double digest?", default=False):
            e2 = Prompt.ask("Enzyme 2")
        else:
            e2 = None
        frags = mapper.virtual_digest(e1, e2)
    RestrictionMapper.display_results(hits, frags, map_str)


def run_codon_opt():
    console.print(Panel(
        "[bold yellow]CODON OPTIMIZATION[/bold yellow]\n"
        "[dim]Kazusa codon tables · CAI · rare codon detection[/dim]",
        border_style="yellow"))
    console.print("\n[1] Optimize DNA sequence")
    console.print("[2] Back-translate protein sequence")
    choice = Prompt.ask("Input type", choices=["1","2"], default="1")
    input_type = "dna" if choice=="1" else "protein"
    if input_type == "dna":
        try:
            seq, rid = get_sequence_from_user("coding sequence (CDS)")
        except Exception as e:
            console.print(f"[red]{e}[/red]"); return
    else:
        console.print("[dim]Paste protein sequence (1-letter AA, Enter twice):[/dim]")
        lines = []
        while True:
            line = input()
            if line=="" and lines: break
            lines.append(line)
        seq = "".join(lines).upper()
        rid = "protein_input"
    # Host selection
    hosts = list(CODON_TABLES.keys())
    tbl = Table(title="[bold yellow]Select Host Organism[/bold yellow]",
                box=box.SIMPLE, header_style="bold yellow")
    tbl.add_column("#"); tbl.add_column("Organism")
    for i, h in enumerate(hosts, 1): tbl.add_row(str(i), h)
    console.print(tbl)
    try:
        hidx = int(Prompt.ask("Host #", default="1")) - 1
        host = hosts[max(0,min(hidx,len(hosts)-1))]
    except ValueError:
        host = hosts[0]
    with console.status("[bold yellow]Optimizing...[/bold yellow]", spinner="dots"):
        result = CodonOptimizer(seq, host, input_type).optimize()
    CodonOptimizer.display_results(result)
    if Confirm.ask("Save optimized sequence?", default=False):
        _save_text(f"optimized_{rid[:15]}.txt",
                   f">Optimized for {host}\n{result['optimized']}")


def run_plasmid():
    console.print(Panel(
        "[bold yellow]PLASMID DESIGNER[/bold yellow]\n"
        "[dim]Modular assembly · promoter/RBS/tag/terminator/ori/resistance[/dim]",
        border_style="yellow"))
    if Confirm.ask("Do you have an insert sequence to clone?", default=True):
        try:
            seq, _ = get_sequence_from_user("insert sequence")
        except Exception as e:
            console.print(f"[red]{e}[/red]"); seq = ""
    else:
        seq = ""
    designer = PlasmidDesigner(seq)
    plan     = designer.interactive_build()
    PlasmidDesigner.display_results(plan)
    if Confirm.ask("Save plasmid design?", default=False):
        _save_text("plasmid_design.txt",
                   f"Promoter: {plan['promoter'][0]}\n"
                   f"RBS: {plan['rbs'][0]}\nTag: {plan['tag'][0]}\n"
                   f"Terminator: {plan['terminator'][0]}\n"
                   f"Origin: {plan['origin'][0]}\nResistance: {plan['resistance'][0]}\n"
                   f"Size: ~{plan['total_size']} bp\n"
                   f"Expression score: {plan['expression_score']}")


# ════════════════════════════════════════════════════════════════
#  ── PHASE 2 MODULE RUNNERS ──
# ════════════════════════════════════════════════════════════════

def run_sirna():
    console.print(Panel(
        "[bold magenta]siRNA / shRNA DESIGN[/bold magenta]\n"
        "[dim]Reynolds (2004) 8-criterion scoring · hairpin shRNA output[/dim]",
        border_style="magenta"))
    try:
        seq, rid = get_sequence_from_user("mRNA target")
    except Exception as e:
        console.print(f"[red]{e}[/red]"); return
    top_n = int(Prompt.ask("Top candidates to show", default="8"))
    with console.status("[bold magenta]Designing siRNAs...[/bold magenta]", spinner="dots"):
        data = siRNADesigner(seq).design(top_n=top_n)
    console.print(f"[green]✔ {data['total_found']} valid siRNA(s) found[/green]")
    siRNADesigner.display_results(data)
    if data["candidates"] and Confirm.ask("Save results?", default=False):
        lines = [f">{i+1} pos={c['position']} score={c['score']}\n{c['sense_rna']}"
                 for i,c in enumerate(data["candidates"])]
        _save_text(f"sirna_{rid[:15]}.txt","\n".join(lines))


def run_rna_structure():
    console.print(Panel(
        "[bold cyan]RNA SECONDARY STRUCTURE PREDICTOR[/bold cyan]\n"
        "[dim]Nussinov algorithm · Turner energetics · dot-bracket output[/dim]",
        border_style="cyan"))
    try:
        seq, rid = get_sequence_from_user("RNA/DNA sequence")
    except Exception as e:
        console.print(f"[red]{e}[/red]"); return
    if len(seq) > 200:
        console.print("[yellow]⚠ Sequence trimmed to 200 nt for structure prediction.[/yellow]")
        seq = seq[:200]
    with console.status("[bold cyan]Predicting structure...[/bold cyan]", spinner="dots"):
        result = RNAStructurePredictor(seq).predict()
    RNAStructurePredictor.display_results(result)


def run_phylogenetics():
    console.print(Panel(
        "[bold green]PHYLOGENETICS BUILDER[/bold green]\n"
        "[dim]Needleman-Wunsch alignment · Jukes-Cantor distance · "
        "Neighbor-Joining tree[/dim]",
        border_style="green"))
    console.print("\nLoad sequences from:")
    console.print("  [1] Multiple FASTA file (recommended)")
    console.print("  [2] Paste individual sequences")
    choice = Prompt.ask("Choice", choices=["1","2"], default="1")
    sequences = {}
    if choice == "1":
        try:
            _, _ = get_sequence_from_user  # just to confirm import
            path = Prompt.ask("FASTA file path (multi-record)")
            from utilities import detect_and_parse, clean_sequence
            raw_recs = detect_and_parse(path)
            for rid, raw in raw_recs.items():
                try:
                    sequences[rid] = clean_sequence(raw)
                except Exception:
                    pass
        except Exception as e:
            console.print(f"[red]{e}[/red]"); return
    else:
        try:
            n = int(Prompt.ask("Number of sequences to enter", default="4"))
        except ValueError:
            n = 4
        for i in range(n):
            name = Prompt.ask(f"Name for sequence {i+1}")
            console.print(f"[dim]Paste sequence {i+1} (Enter twice):[/dim]")
            lines = []
            while True:
                line = input()
                if line=="" and lines: break
                lines.append(line)
            try:
                from utilities import clean_sequence
                sequences[name] = clean_sequence(" ".join(lines))
            except Exception as ex:
                console.print(f"[red]{ex}[/red]")
    if len(sequences) < 3:
        console.print("[red]Need at least 3 sequences for phylogenetics.[/red]"); return
    console.print(f"\n[green]✔ {len(sequences)} sequences loaded[/green]")
    with console.status("[bold green]Computing distances & tree...[/bold green]", spinner="dots"):
        builder = PhylogeneticsBuilder(sequences)
        dist    = builder.build_distance_matrix()
        tree    = builder.neighbor_joining(dist)
    PhylogeneticsBuilder.display_results(tree, dist, list(sequences.keys()))


# ════════════════════════════════════════════════════════════════
#  ── PHASE 3 MODULE RUNNERS ──
# ════════════════════════════════════════════════════════════════

def run_vaccine():
    console.print(Panel(
        "[bold yellow]VACCINE EPITOPE DESIGNER[/bold yellow]\n"
        "[dim]B-cell · MHC-I CTL · MHC-II HTL · multi-epitope construct[/dim]",
        border_style="yellow"))
    console.print("\nInput type:")
    console.print("  [1] DNA / CDS sequence")
    console.print("  [2] Protein sequence (1-letter)")
    choice = Prompt.ask("Choice", choices=["1","2"], default="1")
    if choice == "1":
        try:
            seq, rid = get_sequence_from_user("antigen CDS")
        except Exception as e:
            console.print(f"[red]{e}[/red]"); return
        input_type = "dna"
    else:
        console.print("[dim]Paste protein sequence (Enter twice):[/dim]")
        lines = []
        while True:
            line = input()
            if line=="" and lines: break
            lines.append(line)
        seq = " ".join(lines); rid = "protein"; input_type = "protein"
    # Allele selection
    alleles = ["HLA-A*02:01","HLA-A*01:01","HLA-A*03:01",
               "HLA-B*07:02","HLA-B*44:02"]
    console.print("\nCommon MHC-I alleles:")
    for i,a in enumerate(alleles,1): console.print(f"  [{i}] {a}")
    try:
        ai = int(Prompt.ask("Select allele", default="1")) - 1
        allele = alleles[max(0,min(ai,len(alleles)-1))]
    except ValueError:
        allele = alleles[0]
    with console.status("[bold yellow]Predicting epitopes...[/bold yellow]", spinner="dots"):
        data = VaccineEpitopeDesigner(seq, input_type).design(allele)
    VaccineEpitopeDesigner.display_results(data)
    if Confirm.ask("Save construct?", default=False):
        _save_text(f"vaccine_{rid[:15]}.txt",
                   f">Multi-epitope construct\n{data['construct']}")


def run_amp():
    console.print(Panel(
        "[bold magenta]ANTIMICROBIAL PEPTIDE DESIGNER[/bold magenta]\n"
        "[dim]Amphipathicity · Boman index · net charge · helical wheel[/dim]",
        border_style="magenta"))
    console.print("[dim]Paste your peptide sequence (1-letter AA):[/dim]")
    pep = input().strip().upper()
    if not pep:
        console.print("[red]Empty input.[/red]"); return
    with console.status("[bold magenta]Analysing...[/bold magenta]", spinner="dots"):
        data = AMPDesigner(pep).analyze()
    AMPDesigner.display_results(data)


def run_arg():
    console.print(Panel(
        "[bold red]ANTIBIOTIC RESISTANCE GENE SCANNER[/bold red]\n"
        "[dim]WHO priority ARGs · 33-gene database · beta-lactam/AMG/tet/mcr[/dim]",
        border_style="red"))
    console.print("\nInput type:")
    console.print("  [1] DNA sequence / genome")
    console.print("  [2] Protein sequence")
    choice = Prompt.ask("Choice", choices=["1","2"], default="1")
    if choice == "1":
        try:
            seq, rid = get_sequence_from_user("genome / gene")
        except Exception as e:
            console.print(f"[red]{e}[/red]"); return
        input_type = "dna"
    else:
        console.print("[dim]Paste protein sequence (Enter twice):[/dim]")
        lines = []
        while True:
            line = input()
            if line=="" and lines: break
            lines.append(line)
        seq = " ".join(lines); rid = "protein"; input_type = "protein"
    try:
        threshold = float(Prompt.ask("Similarity threshold (0–1)", default="0.60"))
    except ValueError:
        threshold = 0.60
    with console.status("[bold red]Scanning for ARGs...[/bold red]", spinner="dots"):
        hits = ARGScanner(seq, input_type).scan(threshold)
    console.print(f"[green]✔ Scan complete – {len(hits)} hit(s)[/green]")
    ARGScanner.display_results(hits)
    if hits and Confirm.ask("Save report?", default=False):
        lines = [f"{h['gene']}\t{h['drug_class']}\t{h['mechanism']}\t{h['similarity']}\t{h['risk']}"
                 for h in hits]
        _save_text(f"arg_report_{rid[:15]}.txt", "\n".join(lines))


# ════════════════════════════════════════════════════════════════
#  ── HELPERS ──
# ════════════════════════════════════════════════════════════════

def _save_text(filename: str, content: str):
    with open(filename, "w") as fh:
        fh.write(content)
    console.print(f"[green]✔ Saved to [bold]{filename}[/bold][/green]")


# ════════════════════════════════════════════════════════════════
#  ── MAIN MENU ──
# ════════════════════════════════════════════════════════════════

PHASE_MENUS = {
    "1": {
        "title" : "Phase 1  –  Core Molecular Biology",
        "color" : "cyan",
        "items"  : {
            "1": ("Primer Design",                          run_primer_design),
            "2": ("CRISPR Guide Design + Animation",        run_crispr),
            "3": ("Restriction Enzyme Mapping",             run_restriction),
            "4": ("Codon Optimization",                     run_codon_opt),
            "5": ("Plasmid / Expression Vector Designer",   run_plasmid),
        }
    },
    "2": {
        "title" : "Phase 2  –  RNA Biology & Genomics",
        "color" : "green",
        "items"  : {
            "1": ("siRNA / shRNA Design",                   run_sirna),
            "2": ("RNA Secondary Structure Prediction",     run_rna_structure),
            "3": ("Phylogenetics & Tree Builder",           run_phylogenetics),
        }
    },
    "3": {
        "title" : "Phase 3  –  Translational & Clinical",
        "color" : "yellow",
        "items"  : {
            "1": ("Vaccine Epitope Designer",               run_vaccine),
            "2": ("Antimicrobial Peptide (AMP) Analyser",  run_amp),
            "3": ("Antibiotic Resistance Gene (ARG) Scanner", run_arg),
        }
    },
}


def phase_menu(phase_key: str):
    phase = PHASE_MENUS[phase_key]
    color = phase["color"]
    while True:
        console.print(Panel(
            "\n".join(f"  [bold {color}][{k}][/bold {color}] {v[0]}"
                      for k,v in phase["items"].items()) +
            f"\n  [bold red][B][/bold red] Back to main menu",
            title=f"[bold {color}]{phase['title']}[/bold {color}]",
            border_style=color, width=60))
        choices = list(phase["items"].keys()) + ["b","B"]
        ch = Prompt.ask("Select module", choices=choices,
                        default=list(phase["items"].keys())[0])
        if ch.upper() == "B": break
        console.print()
        try:
            phase["items"][ch][1]()
        except KeyboardInterrupt:
            console.print("\n[yellow]Module interrupted.[/yellow]")
        except Exception as exc:
            console.print(f"\n[bold red]Error:[/bold red] {exc}")
            import traceback; console.print_exception()
        console.print()
        if not Confirm.ask("Continue in this phase?", default=True):
            break


def main_menu():
    splash()
    while True:
        console.print(Panel(
            "  [bold cyan][1][/bold cyan]  Phase 1  –  Core Molecular Biology\n"
            "             Primer · CRISPR · Restriction · Codon Opt · Plasmid\n\n"
            "  [bold green][2][/bold green]  Phase 2  –  RNA Biology & Genomics\n"
            "             siRNA · RNA Structure · Phylogenetics\n\n"
            "  [bold yellow][3][/bold yellow]  Phase 3  –  Translational & Clinical\n"
            "             Vaccine · AMP · ARG Scanner\n\n"
            "  [bold red][Q][/bold red]  Quit",
            title="[bold white]DBS  –  Digital BioScience Suite[/bold white]",
            border_style="bright_white", width=62))
        ch = Prompt.ask("Select phase", choices=["1","2","3","q","Q"],
                        default="1").upper()
        if ch == "Q":
            console.print("\n[bold cyan]Goodbye from DBS!  Happy researching 🧬[/bold cyan]\n")
            sys.exit(0)
        console.print()
        try:
            phase_menu(ch)
        except KeyboardInterrupt:
            console.print("\n[yellow]Returned to main menu.[/yellow]")
        console.clear(); splash()


if __name__ == "__main__":
    main_menu()