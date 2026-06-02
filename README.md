# MTS-Finder

**MTS-Finder** is a desktop application for discovering **mitochondrial targeting sequences (MTS / presequences)** in peptide-level shotgun (bottom-up) proteomics data. It cross-references identified peptides against curated Uniprot presequence annotations and TargetP-2.0 cleavage-site predictions, flags peptides that fall within a mitochondrial presequence, fetches gene symbols from Uniprot, and computes quantitative comparisons between experimental conditions.

**Current version: v3.6 (2026)** — see [CHANGELOG.md](CHANGELOG.md) for the full release history. v3.6 is the first fully working release; the original v3.4 (2023) was non-functional.

---

## Features

- **MTS / presequence detection** — matches each peptide's master-protein accession against two reference sets:
  - **Uniprot** curated presequence (transit-peptide) annotations (`files/Uniprot_MTS.xlsx`)
  - **TargetP-2.0** mitochondrial cleavage-site predictions for the MitoCarta3 proteome (`files/TargetP2_0_prediction_mitocarta3.xlsx`)
- **Cleavage-site aware filtering** — a peptide is reported only when it lies *within* the presequence, i.e. its last residue ends at or before the predicted/annotated cleavage site.
- **Automatic gene-symbol lookup** — resolves accessions to gene symbols via the [EBI Proteins API](https://www.ebi.ac.uk/proteins/api/doc/), with a retry and a graceful fallback to the accession ID when offline.
- **Flexible input** — reads Proteome Discoverer–style peptide exports as either Excel (`.xlsx`) or tab-delimited text (`.txt`).
- **Normalized or raw abundances** — pick "Abundances (Normalized)" or raw "Abundance" channels with a single toggle.
- **Channel-to-condition mapping** — a simple comma-separated list assigns each abundance channel to an experimental condition; `skip` and `Boost` channels are dropped automatically.
- **Built-in quantification** — for each requested comparison pair it computes:
  - `Log2(ratio)` of the group means
  - `pvalue` from an unpaired (Welch-style, `equal_var=True`) two-sample *t*-test
  - `-Log10 pvalue` (volcano-plot ready)
  - Degenerate rows (zero/negative ratios, empty groups) yield `NaN` instead of crashing the run.
- **Excel output** — results are written to an `.xlsx` workbook with one row per matched peptide, and an **Open** button launches the file when done.
- **Simple GUI** — a Tkinter interface with file browsers and editable Conditions/Pairs boxes; long-running analysis executes on a background thread so the window stays responsive.

---

## Requirements

- **Python 3.8+**
- Python packages:
  - `pandas`
  - `scipy`
  - `requests`
  - `openpyxl` (Excel read/write)
  - `tkinter` (ships with most Python installations; on Linux install `python3-tk`)
- Internet access (optional but recommended) for gene-symbol lookups via the EBI Proteins API.

Install the dependencies:

```bash
pip install pandas scipy requests openpyxl
```

---

## Installation

```bash
git clone https://github.com/science64/MTS-Finder.git
cd MTS-Finder
pip install pandas scipy requests openpyxl
```

Make sure the `files/` folder (shipped with the repository) stays next to `main.py` — it contains the reference databases and the application icon that the tool loads at runtime.

---

## Usage

Launch the application from the project root:

```bash
python main.py
```

Then, in the window:

1. **Browse** and select your peptide file (`.xlsx` or tab-delimited `.txt`).
2. Choose **Normalized values** or **Non-Normalized values**.
3. Enter an **Output Name** for the result workbook.
4. Review/edit the **Conditions** and **Pairs** boxes (these are pre-filled from `condtions.txt` and `pairs.txt`; you can also load them from a file with their **Browse** buttons).
5. Click **RUN**. Progress is shown in the status box.
6. When finished, click **Open** to view the generated `<Output Name>.xlsx`, saved next to your input file.

> Note: the chosen Conditions and Pairs are also written back as `condtions.txt` and `pairs.txt` in the output folder, so your last configuration is preserved.

### Conditions

A single comma-separated line with **one label per abundance channel**, in the same order as the channels appear in the file. Special labels are dropped before analysis:

- `skip` — ignore this channel
- `Boost` — ignore the carrier/boost channel

Example:

```
Light,WT,WT,WT,WT,KO,KO,KO,KO,skip,Boost
```

### Pairs

Semicolon-separated comparisons, each written as `numerator/denominator`. Each pair produces a `Log2`, `pvalue`, and `-Log10 pvalue` column. For a `KO/WT` comparison, `KO` is the numerator (group 2) and `WT` is the denominator (group 1):

```
KO/WT
```

Multiple comparisons:

```
KO/WT; Rotenone/0DMSO; Antimycin/0DMSO
```

### Expected input columns

The peptide file should contain the typical Proteome Discoverer peptide-export columns, including:

- `Positions in Master Proteins` (e.g. `P49189 [275-293]`, or multiple `; `-separated accessions)
- `Modifications`
- Abundance channels — either `Abundances (Normalized)...` / `Abundances Normalized...` or `Abundance: ...` / `Abundance...`

### Output columns

| Column | Description |
| --- | --- |
| `Accession` | Uniprot accession of the master protein |
| `Gene Symbol` | Gene symbol resolved from Uniprot (falls back to accession) |
| `Positions in Master Proteins` | Gene symbol with peptide position range |
| `Modifications` | Peptide modifications |
| `Uniprot` / `Uniprot Location` | Whether matched via Uniprot annotation, and the cleavage-site position |
| `TargetP` / `TargetP Location` | Whether matched via TargetP-2.0 prediction, and the cleavage-site position |
| *(condition channels)* | Abundance values, renamed to your condition labels |
| `Log2(num/den)` | Log2 ratio of group means, per pair |
| `pvalue(num/den)` | Unpaired two-sample *t*-test p-value, per pair |
| `-Log10 pvalue(num/den)` | −log10 of the p-value, per pair |

---

## Reference data

The `files/` directory contains the bundled reference databases used by the matcher:

- `Uniprot_MTS.xlsx` — Uniprot presequence (transit-peptide) locations
- `TargetP2_0_prediction_mitocarta3.xlsx` — TargetP-2.0 mitochondrial cleavage-site predictions
- `Total_accession_precurser.xlsx` — combined accession list used for fast pre-filtering
- `icon.ico` — application icon

---

## Project structure

```
MTS-Finder/
├── main.py          # Tkinter GUI and analysis orchestration
├── functions.py     # MTS-finding engine, Uniprot/TargetP matching, statistics
├── condtions.txt    # Default conditions line
├── pairs.txt        # Default comparison pairs
├── files/           # Reference databases and icon
├── CHANGELOG.md
├── LICENSE
└── README.md
```

---

## License

Released under the **MIT License**. Copyright © 2023 Süleyman Bozkurt. See [LICENSE](LICENSE) for the full text.

---

## Citation / contact

Developed by **Süleyman Bozkurt**. For questions or issues, please open an issue on the project repository.
