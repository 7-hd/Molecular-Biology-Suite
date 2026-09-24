# Molecular Biology Suite

A terminal-based Python suite for sequence analysis and in-silico molecular biology workflows.

## Install

```bash
python -m pip install -e .
```

## Run

```bash
python main.py
# or
molecular-biology-suite
```

## Architecture

- `app/` — CLI entry points
- `biology/` — domain-facing public modules
- `core/` — shared sequence utilities
- `structure.py` — legacy engine retained during the incremental migration
- `utilities.py` — legacy shared data and file/thermodynamic helpers
- `tests/` — regression coverage

The current refactor keeps the existing interactive behavior while providing
stable package imports. The legacy engine will be migrated module-by-module in
subsequent changes.

> Results are computational estimates and should be independently validated
> before research or clinical use.
