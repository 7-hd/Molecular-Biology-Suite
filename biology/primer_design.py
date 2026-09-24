"""Primer design domain module."""

from __future__ import annotations

from rich.console import Console
from rich.table import Table
from rich import box

from utilities import (
    check_gc_clamp,
    gc_content,
    primer_dimer_check,
    reverse_complement,
    self_complementarity_score,
    tm_nearest_neighbor,
)

console = Console()


class PrimerDesigner:
    """Design forward/reverse primer pairs using configurable constraints."""

    DEFAULT = {
        "min_len": 18,
        "max_len": 25,
        "gc_min": 40.0,
        "gc_max": 60.0,
        "tm_min": 52.0,
        "tm_max": 68.0,
        "tm_delta": 5.0,
        "top_n": 5,
    }

    def __init__(self, sequence: str, **kwargs):
        self.seq = sequence.upper()
        self.cfg = {**self.DEFAULT, **kwargs}
        if self.cfg["min_len"] <= 0 or self.cfg["max_len"] < self.cfg["min_len"]:
            raise ValueError("Primer lengths must be positive and min_len <= max_len.")

    def _score_candidate(self, primer: str) -> dict:
        gc = gc_content(primer)
        tm = tm_nearest_neighbor(primer)
        clamp = check_gc_clamp(primer)
        self_complementarity = self_complementarity_score(primer)
        warnings = []
        if not self.cfg["gc_min"] <= gc <= self.cfg["gc_max"]:
            warnings.append("GC%")
        if not self.cfg["tm_min"] <= tm <= self.cfg["tm_max"]:
            warnings.append("Tm")
        if not clamp:
            warnings.append("GC-clamp")
        if self_complementarity > 6:
            warnings.append("self-complement")
        score = (
            100 - abs(gc - 50) * 0.8 - abs(tm - 60) * 0.5
            + (5 if clamp else 0) - self_complementarity * 2
        )
        return {
            "primer": primer,
            "len": len(primer),
            "gc": round(gc, 1),
            "tm": tm,
            "clamp": clamp,
            "sc": self_complementarity,
            "score": round(score, 1),
            "warnings": warnings,
        }

    def _scan_candidates(self, strand_seq: str, label: str) -> list[dict]:
        candidates = []
        for start in range(0, len(strand_seq) - self.cfg["min_len"] + 1):
            for length in range(self.cfg["min_len"], self.cfg["max_len"] + 1):
                end = start + length
                if end > len(strand_seq):
                    break
                result = self._score_candidate(strand_seq[start:end])
                if not result["warnings"]:
                    result.update({"start": start + 1, "end": end, "label": label})
                    candidates.append(result)

        unique = []
        seen = set()
        for candidate in sorted(candidates, key=lambda item: -item["score"]):
            if candidate["primer"] not in seen:
                seen.add(candidate["primer"])
                unique.append(candidate)
        return unique

    def design(self) -> dict:
        forward = self._scan_candidates(self.seq, "Forward")
        reverse = self._scan_candidates(reverse_complement(self.seq), "Reverse")
        pairs = []
        for fwd in forward[:20]:
            for rev in reverse[:20]:
                delta_tm = abs(fwd["tm"] - rev["tm"])
                if delta_tm <= self.cfg["tm_delta"] and not primer_dimer_check(fwd["primer"], rev["primer"]):
                    pairs.append({
                        "fwd": fwd,
                        "rev": rev,
                        "delta_tm": round(delta_tm, 1),
                        "pair_score": fwd["score"] + rev["score"] - delta_tm * 2,
                    })
        pairs.sort(key=lambda item: -item["pair_score"])
        return {
            "forward": forward[: self.cfg["top_n"]],
            "reverse": reverse[: self.cfg["top_n"]],
            "pairs": pairs[:5],
        }

    @staticmethod
    def display_results(results: dict, seq: str) -> None:
        """Render candidates in the existing Rich terminal format."""
        for direction, color in (("forward", "cyan"), ("reverse", "magenta")):
            candidates = results[direction]
            if not candidates:
                console.print(f"[red]No valid {direction} primers found.[/red]")
                continue
            table = Table(
                title=f"[bold {color}]{direction.upper()} PRIMERS[/bold {color}]",
                box=box.ROUNDED,
                show_lines=True,
                header_style=f"bold {color}",
            )
            for column in ("#", "Sequence (5'→3')", "Len", "Tm", "GC%", "Clamp", "SC", "Score"):
                table.add_column(column, justify="center")
            for index, candidate in enumerate(candidates, 1):
                clamp_color = "green" if candidate["clamp"] else "red"
                table.add_row(
                    str(index),
                    f"[{color}]{candidate['primer']}[/{color}]",
                    str(candidate["len"]),
                    str(candidate["tm"]),
                    str(candidate["gc"]),
                    f"[{clamp_color}]{'✔' if candidate['clamp'] else '✘'}[/{clamp_color}]",
                    str(candidate["sc"]),
                    f"[bold]{candidate['score']}[/bold]",
                )
            console.print(table)

        if not results["pairs"]:
            console.print("[red]⚠ No valid primer pairs found.[/red]")
            return
        table = Table(
            title="[bold yellow]★ RECOMMENDED PAIRS[/bold yellow]",
            box=box.DOUBLE,
            show_lines=True,
            header_style="bold yellow",
        )
        for column in ("#", "Forward (5'→3')", "Fwd Tm", "Reverse (5'→3')", "Rev Tm", "ΔTm", "Score"):
            table.add_column(column, justify="center")
        for index, pair in enumerate(results["pairs"], 1):
            table.add_row(
                str(index), pair["fwd"]["primer"], str(pair["fwd"]["tm"]),
                pair["rev"]["primer"], str(pair["rev"]["tm"]),
                str(pair["delta_tm"]), str(round(pair["pair_score"], 1)),
            )
        console.print(table)
