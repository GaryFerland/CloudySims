# Mao et al. (2022) → Cloudy STOUT conversion

**Project documentation — 16 October 2026**  
**Scope:** prepare H-like and He-like atomic datasets supplied in `Mao22.zip` as three-file STOUT datasets for **subsequent review by the Cloudy team**. **This workflow does not install data into Cloudy or choose which atomic model Cloudy should use.**

## 1. Quick start

Keep these files together in a working directory:

```text
Mao22/
├── Mao22.zip                         # Original, unmodified source archive
├── ConvertToStoutAll_Mao22.py        # CONVERSION: needed
├── validate_Mao22_STOUT.py           # VALIDATION: needed
└── Mao22_STOUT/                      # Created by converter
```

From inside the directory:

```bash
# Convert all ions from the original folder:
python3 ConvertToStoutAll_Mao22.py . -o Mao22_STOUT

# Validate every generated ion against the source archive:
python3 validate_Mao22_STOUT.py Mao22.zip -o Mao22_STOUT

# To convert or validate just one ion (example: Fe XXVI):
python3 ConvertToStoutAll_Mao22.py Mao22.zip --ion fe_26 -o Mao22_STOUT
python3 validate_Mao22_STOUT.py Mao22.zip --ion fe_26 -o Mao22_STOUT
```

If the original ZIP has already been extracted, use the **directory that contains `adf04_ecs`**, or its parent directory:

```bash
python3 ConvertToStoutAll_Mao22.py . -o Mao22_STOUT
python3 validate_Mao22_STOUT.py . -o Mao22_STOUT
```

Here `.` means **the current directory**. The positional `source` argument is mandatory. Omitting it causes the error `the following arguments are required: source`.

Both new scripts use the **Python standard library only** (no pandas or NumPy required).

## 2. What are we converting, and why?

`Mao22.zip` is the original archive. For this particular workflow, the converter selects files ending in `.dat` under **`Mao22/adf04_ecs/`**; it does not convert the alternative `adf04_ocs/`, `adf04_str/`, or `data/` tables. Those other sources may still be useful for independent scientific comparisons.

Each input ADF04 file combines several kinds of atomic information:

- **Energy levels:** level index, electron configuration, term, angular momentum `J`, and level energy.
- **Radiative transitions:** transition indices and Einstein coefficient `A` (s⁻¹).
- **Electron-impact excitation:** effective collision strengths **Υ(T)** listed on an input temperature grid.

For an ion, the output is three STOUT-format files:

| Output file | Atomic information | Primary use |
|---|---|---|
| `<ion>.nrg` | Energies in cm⁻¹, statistical weights `g`, level labels | Atomic levels / transition wavelengths |
| `<ion>.tp` | Einstein A-values in s⁻¹ | Radiative rates |
| `<ion>.coll` | Electron effective collision strengths vs. temperature | Collisional excitation rates |

**Important:** effective collision strength `Υ(T)` is a dimensionless thermally averaged collision quantity, *not* a rate coefficient in cm³ s⁻¹. The conversion copies Υ(T); it does **not** calculate excitation rate coefficients. The converter treats the temperature header as Kelvin, as intended for its STOUT `TEMP` row; verify the original ADF04 temperature-unit convention independently before adopting the data scientifically.

## 3. Folder naming and ion identification

Input examples (the number in these ADF04 filenames denotes **ion charge**, not the spectroscopic stage):

```text
hlike_jm22#fe25.dat   -> Fe+25, one remaining electron -> fe_26/
helike_jm22#fe24.dat  -> Fe+24, two remaining electrons -> fe_25/
hlike_jm22#o7.dat     -> O+7  -> o_8/
helike_jm22#o6.dat    -> O+6  -> o_7/
```

The Cloudy/STOUT file convention uses **spectroscopic ionization stage = charge + 1**. `fe_26` therefore means Fe XXVI, not charge +26.

A completed directory looks like:

```text
Mao22_STOUT/
├── conversion_report.json
├── fe_26/
│   ├── fe_26.nrg
│   ├── fe_26.tp
│   └── fe_26.coll
├── fe_25/
│   ├── fe_25.nrg
│   ├── fe_25.tp
│   └── fe_25.coll
└── ...
```

The converter checks that a source filename is consistent with its purported H-like or He-like sequence by comparing atomic number, charge, and remaining electron count. Unknown and inconsistent names trigger errors instead of producing misleading filenames.

## 4. Main script: `ConvertToStoutAll_Mao22.py`

**Purpose:** transform Mao+22 ADF04 ECS records into **STOUT text files**, preserving original level indices and the tabulated quantities. The converter **does not** replace level energies or A-values with NIST values, adjust wavelengths, interpolate/extrapolate collision strengths, or modify Cloudy's C++ code.

### 4.1 Entry point and command-line arguments

The `main()` function uses `argparse`:

- `source` (**required**): path to `Mao22.zip` or an extracted archive directory.
- `-o`, `--outdir`: output directory; defaults to `Mao22_STOUT`.
- `--ion`: optional filter for one STOUT ion name, e.g. `fe_26`.

It loops through inputs, converts each ion, prints `OK` or `ERROR`, and writes `conversion_report.json`. It exits with a nonzero status if any ion fails or nothing matches the input/filter.

### 4.2 `input_files(path)` — find and read original files

- If `source` is a ZIP, reads `Mao22/adf04_ecs/*.dat` **directly from the ZIP** via `zipfile`.
- If `source` is a directory, recursively locates `.dat` files with `adf04_ecs` in their path.
- Ignores macOS `__MACOSX` metadata.
- Yields `(source_filename, file_text)` for each candidate file.

This is why a simple `Mao22.zip` command works without manually extracting the archive.

### 4.3 `parse_ion_name(filename)` — determine the STOUT ion name

The pattern `hlike_jm22#fe25.dat` is read as iron, charge +25, H-like. The code uses its built-in periodic-table symbols to identify atomic number and verifies that **atomic number − charge** equals 1 (H-like) or 2 (He-like), then produces `fe_26`.

### 4.4 `parse_adf04(content, filename)` — read the atomic data

The parser separates the original ADF04 file into **energy-level section**, **temperature header**, and **transition records**.

**A. Energy levels:** after the ADF04 header, it reads each level up to the `-1` sentinel. It obtains the fixed-position level index and configuration, parses the spectroscopic term and `J`, and calculates statistical weight:

\[
g=2J+1.
\]

For example, `J = 1.5` gives `g = 4`. Level labels are assembled from configuration, spin multiplicity, angular-momentum letter (`S`, `P`, `D`, …), and `J`. Energy strings are preserved as supplied, rather than converted through a rounded intermediate number.

**B. Temperature grid:** immediately after the level section, the first two fields of the temperature record are treated as metadata and the remaining entries are read as temperature values. The code requires positive, finite, increasing/non-decreasing temperatures.

**C. Transition records:** each record supplies **upper-level index, lower-level index, A-value**, then the temperature-dependent effective collision strengths. The parser requires the indices to correspond to known levels (`lower < upper`) and rejects negative/nonfinite A-values or Υ values.

**D. ADAS number formatting:** ADAS data can use `1.41-03` instead of `1.41E-03`, and numeric fields can touch. `number()` inserts the missing `E` when required; the `NUM` regular expression extracts scientific-notation quantities even when no blank separates adjacent values.

**E. Extra trailing field:** the code accepts either **N** or **N+1** numerical fields after the A-value, where **N** is the number of temperatures. If an extra value is present, it uses the **first N** as Υ(T), **ignores the final field**, and records how often this occurs in the report. This is a deliberate behavior of the code, **not yet an independently confirmed interpretation of the extra ADF04 field**.

**F. Error checks:** missing sections, duplicate level numbers, invalid indices, duplicate transition pairs, inconsistent temperature counts, or nonfinite values cause a conversion failure.

### 4.5 `write_ion(...)` — produce the three STOUT files

The shared file structure begins with `17 09 05` (the format's magic/version header), contains the appropriate data lines, and ends in a `**************` reference footer.

**Energy file (`fe_26.nrg`, illustrative structure):**

```text
17 09 05
1    0.0000          2    "1s1 2S0.5"
2    56097682.1693   2    "2p1 2P0.5"
...
**************
#Reference:
# Source: Mao et al. (2022), original Mao22 ADF04 ECS data; input=...
```

Columns: **level index | energy (cm⁻¹) | statistical weight | quoted configuration/term**. Input indices are **not renumbered**.

**Radiative file (`fe_26.tp`, illustrative structure):**

```text
17 09 05
A
1    2    2.830000000e+14
...
**************
#Reference:
...
```

After the line `A`, data are **lower index | upper index | A-value (s⁻¹)**. Only entries with `A > 1.0001e-30` are written. Rows are sorted by **upper**, then **lower** index.

**Collision file (`fe_26.coll`, illustrative structure):**

```text
17 09 05
CS ELECTRON
TEMP    1.350000000e+05    ...
1    2    1.410000000e-03    ...
...
**************
#Reference:
...
```

After `TEMP` (the temperature grid), rows are **lower index | upper index | Υ(T₁) | Υ(T₂) | … | Υ(Tₙ)**. Rows are sorted by **upper**, then **lower**. Each collision row has exactly one Υ value for each temperature. The output uses 9 digits after the decimal point in scientific notation (`.9e`) for A-values, temperatures, and effective collision strengths.

The converter labels the source **Mao et al. (2022)**, with the original filename, rather than incorrectly carrying forward references to previous APAP datasets. **A fully formatted journal citation/DOI should be added only after verification against the source publication.**

### 4.6 `conversion_report.json` — what it contains

The output JSON has `converted` and `failed` lists. Successful ion records include:

- `ion`, `source`, `output`;
- number of `levels`, `transitions`, `a_values`, and `temperatures`;
- `extra_adf04_fields_ignored`.

Failures include the source path and an error message. This report is useful for inventory and troubleshooting, **but a successful report alone does not establish scientific accuracy**.

## 5. Validator: `validate_Mao22_STOUT.py`

**Purpose:** check that the exported files reproduce what our ADF04 parser read, and that each ion's output is internally consistent. It is an **export/consistency audit**, not an independent atomic-physics validation.

### 5.1 What the validator imports

```python
from ConvertToStoutAll_Mao22 import input_files, parse_ion_name, parse_adf04
```

This matters: the validator **shares the converter's interpretation of the input format**. If `parse_adf04()` misinterprets a source field, both conversion and validation could agree on the same wrong value. Therefore its `PASS` result is not independent confirmation that every ADF04 field was understood correctly.

### 5.2 Its main functions

- **`load_data(path, kind)`**: reads generated `.nrg`, `.tp`, `.coll`; checks the `17 09 05` magic string, `A` or `CS ELECTRON` header as appropriate, collision `TEMP` row, and reference footer; extracts numeric data rows.
- **`close(a, b)`**: compares floating-point values using `math.isclose` with relative tolerance `5e-9` and absolute tolerance `1e-12`. This allows for the output's fixed numerical formatting.
- **`check(source, folder, wanted=None)`**: reparses the source ADF04 files, reads the generated files, compares all rows, prints per-ion `PASS`/`FAIL` and a final summary.

### 5.3 What is compared for every ion

| File | Checks made |
|---|---|
| `.nrg` | Same number of levels, level indices in output sequence, energy values, statistical weights |
| `.tp` | Same number of retained A-values, unique lower/upper pairs, A-values match those parsed from ADF04 |
| `.coll` | Same number of transition rows, same temperature grid, unique lower/upper pairs, correct number and values of Υ(T) entries |
| All | Required header/footer forms, presence of files and parseable numerical data |

The validator **does not independently check level-label strings against a spectroscopic reference**, decide whether the temperature-unit interpretation is correct, or test the Cloudy C++ reader.

### 5.4 Actual saved validation result

The provided `Mao22_validation_summary.txt` ends with:

```text
Summary: 50/50 ions passed; 0 failed
NOTE: source interpretation is not independent: validator imports the converter parser.
NOTE: this does not verify Cloudy C++ reader compatibility or collision physics.
```

The set contains **25 H-like ions and 25 He-like ions**. In this dataset's conversion, H-like files each have 36 levels and 630 collision transitions, while He-like files each have 71 levels and 2,485 collision transitions. The number of retained A-value rows differs by ion due to the A-value threshold.

Those counts demonstrate consistent conversion **of the selected `adf04_ecs/` files**, not coverage of every different dataset within `Mao22.zip`.

## 6. What about the five older scripts?

These were supplied from previous STOUT/APAP work. They are **not called by the current Mao+22 converter or validator**:

| Older script | What it does | Needed for this Mao+22 workflow? |
|---|---|---|
| `ConvertToStoutAll.py` | Original ADF04 parser that inspired the new converter. Processes `.dat` in current directory, using more rigid formatting assumptions and generic input-based names. | **No**; superseded by the Mao+22-specific converter for this dataset. |
| `runSubFolders.py` | Finds `ConvertToStoutAll.py` and executes it in immediate subfolders via `subprocess`. | **No**; new converter walks the input archive/directory automatically. |
| `sortAll.py` | Post-processes pre-existing `.coll` blocks: sorts CS rows, splits by `TEMP`, adds headers/refs. Also has `.tp` sorting support. | **No** for initial Mao+22 conversion. Use only if a separate, tested STOUT reader requires splitting collision files. |
| `fix_v2_order.py` | Reorders existing `.nrg` levels by energy and remaps `.tp` and `.coll` indices. | **No** unless an independently confirmed level-order issue is found. Any reindexing must remain consistent across all three files. |
| `Full-best.py` | Previous project's NIST level-energy matching and A-value processing; uses hard-coded references/assumptions. | **No**. Must not apply directly to Mao+22 without adapting and auditing the mapping. |
| `Full-TP.py` | Earlier version of NIST energy/A-value adjustment. | **No**; should not be mixed into raw Mao+22 conversion. |

**Minimal reproducible Mao+22 workflow = two scripts:** `ConvertToStoutAll_Mao22.py` + `validate_Mao22_STOUT.py`.

## 7. Validation status and limitations (important for the Cloudy team)

### Confirmed by the provided scripts and report

1. The conversion selects the original **Mao+22 ADF04 ECS** files and creates STOUT-style `.nrg`, `.tp`, `.coll` files for each recognized ion.
2. It preserves source level indices, level-energy numbers, statistical weights, non-thresholded A-values, and the first N effective collision-strength entries for N tabulated temperatures (subject to numeric output formatting).
3. The validator confirms the generated numerical data agree with those read by the **same** parser: **50/50 pass, no failures** in the saved audit.
4. The scripts require no changes to Cloudy or selection among the existing H-like/He-like atomic models.

### Not yet established — future scientific/integration checks

1. **Interpret the extra ADF04 field** following the N temperature-dependent collision strengths, against authoritative ADF04 format documentation. The code currently discards it.
2. **Independently check** representative levels, A-values, temperature grids, and Υ(T) using separately formatted Mao+22 `data/` tables or a trusted independent parser (for example, Fe XXV, Fe XXVI, O VII, and O VIII). Ideally also compare transition counts and the original source metadata.
3. **Confirm the original ADF04 temperature units and conventions** for these files; the STOUT output labels the grid with `TEMP` and assumes Kelvin.
4. **Test acceptance by Cloudy's current STOUT reader**, including the single `CS ELECTRON` block, headers, footer, and any database indexing requirements. This has **not** been demonstrated by the Python audit.
5. **Resolve precise provenance and publication metadata** (reference/DOI, source version, theoretical versus experimental level provenance) before distributing as production-grade atomic data.
6. **Consider H-like/He-like treatment in Cloudy separately**. The Cloudy team will decide whether and how these STOUT datasets can or should be used; this repository does not enable or override those atomic models.

## 8. Suggested handoff note

> We converted Mao et al. (2022) H-like and He-like ADF04 effective-collision-strength files into STOUT `.nrg`, `.tp`, and `.coll` datasets. The dedicated Python converter produces 50 ion folders and a JSON inventory. A second script validates exported numerical values against the same ADF04 parser and reports 50/50 ions passing. This is a successful structural/numerical export audit, **not yet** an independent physics or Cloudy-reader validation. The remaining items are independent ADF04/atomic-data comparisons, interpretation of an extra trailing ADF04 field, confirmation of temperature conventions, and testing with the STOUT reader. We have not changed Cloudy's treatment of H-like or He-like ions.

## 9. Reproducibility checklist

- [ ] Keep original `Mao22.zip` unchanged.
- [ ] Save exact versions of both Python scripts alongside the ZIP (ideally commit to version control).
- [ ] Convert from a clean output directory and save `conversion_report.json`.
- [ ] Run `validate_Mao22_STOUT.py`; archive the console report.
- [ ] Compare independently to original data for selected ions.
- [ ] Verify the discarded field and the source temperature convention.
- [ ] Have the Cloudy team test one or more ions against the STOUT reader.
- [ ] Only after independent checks, decide whether any scientific data substitutions, NIST anchoring, or Cloudy integration are appropriate.

---

**The guiding principle:** converting atomic data is a *representation* task. We preserve the published numbers first, validate the conversion, and leave model adoption and any scientifically motivated replacements for a separate review.
