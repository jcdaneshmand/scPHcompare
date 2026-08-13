#!/usr/bin/env python3
"""Normalize native dependency-tool output for stable candidate manifests."""

from __future__ import annotations

import argparse
import os
import re


def linux_dependencies(lines: list[str]) -> list[str]:
    address = re.compile(r"\s+\(0x[0-9a-fA-F]+\)$")
    return sorted({address.sub("", line.strip()) for line in lines if line.strip()})


def macos_dependencies(lines: list[str]) -> list[str]:
    values = [line.strip() for line in lines if line.strip()]
    if values and values[0].endswith(":"):
        values = values[1:]
    return sorted(set(values))


def windows_dependencies(lines: list[str]) -> list[str]:
    pattern = re.compile(r"^\s*([A-Za-z0-9_.+-]+\.dll)\s*$", re.IGNORECASE)
    values = {match.group(1).lower() for line in lines if (match := pattern.match(line))}
    if not values:
        raise ValueError("dumpbin output contained no DLL dependencies")
    return sorted(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--platform", required=True, choices=("Linux", "macOS", "Windows"))
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    with open(args.input, encoding="utf-8", errors="replace") as handle:
        lines = handle.read().splitlines()
    normalizers = {
        "Linux": linux_dependencies,
        "macOS": macos_dependencies,
        "Windows": windows_dependencies,
    }
    dependencies = normalizers[args.platform](lines)
    if not dependencies:
        raise ValueError("native dependency inventory is empty")
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    temporary = args.output + ".tmp"
    with open(temporary, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(dependencies) + "\n")
    os.replace(temporary, args.output)


if __name__ == "__main__":
    main()
