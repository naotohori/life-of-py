#!/usr/bin/env python3

'''
This script find all jupyter notebook files (.ipynb) under the current direcotry,
including any sub directories, and convert them to a html file under ./html_output directory.
Each file name of the html file will include the directory name in which the notebook was found.

To use this, you need jupytext and pretty-jupyter. If not,
 $ mamba install jupytext
 $ pip install pretty-jupyter
'''

from pathlib import Path
import subprocess
import sys
import os

OUTDIR = Path("./html_output").resolve()
TEMPLATE = "pj"

EXCLUDE_DIRS = {
    ".ipynb_checkpoints",
    OUTDIR.name,
    ".git",
    "__pycache__",
}

def is_hidden_or_metadata_path(path: Path) -> bool:
    # Exclude files staring with "."
    for part in path.parts:
        if part.startswith("._"):
            return True
        if part.startswith(".") and part not in (".", ".."):
            return True
    return False


def find_notebooks():
    for root, dirs, files in os.walk(".", topdown=True):
        root_path = Path(root)

        # Exclude specific directories.
        dirs[:] = [
            d for d in dirs
            if d not in EXCLUDE_DIRS
            and not d.startswith(".")
            and not d.startswith("._")
        ]

        for f in files:
            if not f.endswith(".ipynb"):
                continue
            if f.startswith(".") or f.startswith("._"):
                continue

            nb = root_path / f

            # Exclude any hidden elements.
            if is_hidden_or_metadata_path(nb):
                continue

            yield nb.resolve()


def make_flat_stem(nb: Path) -> str:
    rel = nb.relative_to(Path(".").resolve())
    return "__".join(rel.with_suffix("").parts)


def needs_update(nb: Path, html: Path) -> bool:
    if not html.exists():
        return True
    return nb.stat().st_mtime > html.stat().st_mtime


def convert_notebook(nb: Path, outdir: Path, template: str) -> tuple[bool, str]:
    flat_stem = make_flat_stem(nb)
    outfile = outdir / f"{flat_stem}.html"

    if not needs_update(nb, outfile):
        return True, f"Skipped (up-to-date): {nb} -> {outfile}"

    cmd = [
        "jupyter",
        "nbconvert",
        "--to",
        "html",
        "--template",
        template,
        "--output-dir",
        str(outdir),
        "--output",
        flat_stem,
        nb.name,
    ]

    try:
        print(f"Converting: {nb} -> {outfile}")
        subprocess.run(
            cmd,
            check=True,
            cwd=nb.parent,
            text=True,
        )
        return True, f"Converted: {nb} -> {outfile}"

    except subprocess.CalledProcessError as e:
        return False, f"Failed: {nb} (return code: {e.returncode})"


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    notebooks = sorted(find_notebooks())

    if not notebooks:
        print("No .ipynb files found.")
        return 0

    failures: list[str] = []
    converted = 0
    skipped = 0

    print("Notebooks to process:")
    for nb in notebooks:
        print(f"  {nb}")

    for nb in notebooks:
        ok, msg = convert_notebook(nb, OUTDIR, TEMPLATE)
        print(msg)

        if ok:
            if msg.startswith("Converted:"):
                converted += 1
            elif msg.startswith("Skipped"):
                skipped += 1
        else:
            failures.append(msg)

    print("\nSummary")
    print(f"  total notebooks: {len(notebooks)}")
    print(f"  converted:       {converted}")
    print(f"  skipped:         {skipped}")
    print(f"  failed:          {len(failures)}")

    if failures:
        print("\nFailures")
        for i, failure in enumerate(failures, start=1):
            print(f"  [{i}] {failure}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
