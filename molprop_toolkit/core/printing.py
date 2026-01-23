"""Printing helpers for CLI analyzers."""

from __future__ import annotations


def print_banner(title: str, width: int = 80, char: str = "=") -> None:
    print(char * width)
    print(title)
    print(char * width)


def print_section(title: str, width: int = 80, char: str = "-") -> None:
    print(char * width)
    print(title)
    print(char * width)

