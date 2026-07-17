# Changelog

All notable changes to **MTS-Finder** are documented in this file.

## [3.7] - 2026-07-16

Abundance-channel detection rewritten so that the underscore header style
(`Abundance_126` / `Abundances_Normalized_126`) is supported and the normalized
and non-normalized groups can no longer be confused with one another. A
normalized run now says so in the output workbook. Adds `requirements.txt`.

### Added (`functions.py`) especially important for MassCore output
- **Underscore abundance headers supported.** Files written with
  `Abundance_126` / `Abundances_Normalized_126` style columns matched none of
  the known patterns. On a normalized run no channel was found at all; on a
  non-normalized run the bare `'Abundance'` fallback matched the raw *and* the
  normalized columns (28 instead of 14) and the run died with an `IndexError`.
  Channel detection now lives in one `findChannels()` helper covering all three
  header styles, and the two selections can no longer overlap.
- **Normalization recorded in the output.** When *Normalized values* is
  selected the condition columns are written as `WT (Normalized)`, so the
  Excel file states which abundances it was built from. Ratio names
  (`Log2(KO/WT)`) are unchanged; non-normalized output is unchanged.

### Fixed (`functions.py`, `main.py`)
- **Derived columns mistaken for channels.** The last-resort `'Abundance'`
  match also picked up Proteome Discoverer's `Abundance Ratio`,
  `Abundances Count`, `Grouped`, `Scaled` and `CV` columns. These are now
  excluded, so bare `Abundance F1 126 Sample` headers resolve to real channels
  only.
- **Silent/cryptic failures.** A missing channel group, or a condition count
  that does not match the channel count, now raises a `ValueError` naming the
  columns it did find, and `main.py` shows it in a dialog and re-enables the
  RUN button instead of letting the worker thread die unnoticed.

### Added (project)
- **`requirements.txt`.** The four third-party dependencies (`pandas`, `scipy`,
  `requests`, `openpyxl`) are now declared in one place and installable with
  `pip install -r requirements.txt`. `tkinter` is deliberately excluded — it is
  part of the standard library and is not available from PyPI.

### Changed (project)
- Bumped the window title to `MTS Finder App v3.7 by S. Bozkurt @2026`.
- README: install steps now use `requirements.txt` (with virtual-environment
  and verification commands), the three supported abundance header styles are
  documented in a table, and the `(Normalized)` output suffix is described.

## [3.6] - 2026-06-02

This is the first release that runs end-to-end correctly. It combines the
six bug fixes made to `functions.py` with follow-up refinements to the
mitochondrial-targeting-sequence (MTS) detection logic and the file dialogs in
`main.py`. The v3.4 release from 2023 did not work; v3.6 supersedes it.

### Fixed (`functions.py`)
- **Off-by-one cut-site lookup.** The cleavage-site index was read as
  `entryNumberList[j-1]` instead of `entryNumberList[j]`, so every matched
  accession was compared against the *previous* protein's presequence
  location. Now uses the correct index `j`.
- **Accession-not-found guard.** When an accession was not present in the
  Uniprot / TargetP MTS list, the loop counter `j` ran past the end of the
  list and produced a wrong (or out-of-range) match. Added an explicit
  `if j >= len(...): return None` bounds check in both `uniprot()` and
  `targetP()`.
- **`None` concatenation.** `out_df` was concatenated even when the
  `uniprot()` / `targetP()` lookup returned `None`, polluting the results
  with empty/garbage rows. The `pd.concat(...)` call is now inside the
  `if cache != None:` block.
- **Multi-accession duplicate rows.** In the multi-accession branch
  (`;`-separated IDs) execution fell through into the single-accession code
  path, so those peptides were processed twice. Added `num += 1; continue`
  to terminate the branch cleanly.
- **Normalization flag never matched.** The Tkinter radio button passes the
  selection as the *string* `'True'`/`'False'`, so `if normalization:` was
  always truthy and non-normalized runs silently used normalized columns.
  Changed to an explicit `if normalization == 'True':` comparison.
- **Unhandled calculation errors.** `Log2`, `pvalue`, and `-Log10 pvalue`
  computations crashed the whole run on zero/negative ratios or degenerate
  t-tests. Each is now wrapped in `try/except (ValueError, ZeroDivisionError)`
  and writes `NaN` for the offending row instead of aborting.

### Changed (`functions.py`)
- **MTS membership criterion corrected.** A peptide is now flagged as part of
  the presequence only when its **last** residue ends at or before the
  cleavage site (`locationLast <= entryNumber`), replacing the earlier check
  on the peptide's first residue. The position parse is wrapped in
  `try/except` and returns `None` on malformed `[start-end]` strings.
- Migrated the deprecated `DataFrame.append(...)` pattern to
  `pd.concat([...], ignore_index=True)` for pandas 2.x compatibility.

### Changed (`main.py`)
- Rewrote the **Browse** dialogs (peptides, conditions, pairs) to use
  `filedialog.askopenfilename()` with `os.path` handling instead of
  `askopenfile()` plus fragile string parsing of the file-handle `repr`.
  This fixes broken output paths on Windows and with paths containing
  quotes/special characters, and lets the dialogs reopen in the last-used
  directory.
- Bumped the window title to `MTS Finder App v3.6 by S. Bozkurt @2026`.

## [3.4] - 2023

- Initial public release. Known to be non-functional due to the bugs listed
  under 3.6.
