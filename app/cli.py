"""CLI entry point for the refactored application."""

from __future__ import annotations


def main() -> None:
    """Run the existing interactive suite through the new app package."""
    from main import main_menu

    main_menu()


if __name__ == "__main__":
    main()
