The Stout Data Format is described in the [StoutData documentation](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/StoutData).


Every `*.stout` file begins with a **magic number** that identifies:

- Which *version* of the format it follows (pre-C25 or C25+)

This allows tools to reliably detect incompatible or outdated stout files before attempting to parse them.

---

## Pre-C25 Stout Files (Legacy Format)

Older stout files (generated before the **C25** release):

- Used the **original magic number**: 17 09 05

---

## C25+ Stout Files (New Format)

Starting with **C25**, stout files use a **new header format** and **new magic number** (25 10 15).
---

# Notes on `*.coll` Files in C25+ versions

`coll` files describe collision data.  
A given dataset may include **multiple files**, typically named with priority-coded suffixes such as:

...blk-10.coll
...blk-20.coll
...blk-30.coll


These files **do not represent separate datasets**—instead, they form a **layered override system**.

---

## 1. Multiple Files, Ordered by Priority

When loading collision data:

1. The loader reads the **lowest-priority file first** (e.g., `blk_10.coll`).
2. It then loads the next file (e.g., `blk_20.coll`), **replacing/overriding** any entries that appear again.
3. This process continues upward through all available `blk_*` files.

Thus, later (higher-numbered) files **replace** earlier (lower-numbered) values.

---

## 2. Meaning of Priority Levels

### **`blk_10` — Lowest priority**
- Usually produced by **R-matrix calculations**.
- Provides a complete baseline dataset.
- May include theoretical or approximate values.

### **`blk_20`, `blk_30`, … — Higher priority**
- Contain **targeted calculations**, **refined models**, or **observational data**.
- Whenever a dataset in these files overlaps with one from a lower-priority file:
  - The higher-priority entry **overrides** the earlier one.
- This ensures the final merged collision dataset always uses:
  1. Observational/targeted data (if available)  
  2. Else refined calculations  
  3. Else the R-matrix baseline

---

## 3. Final Outcome

The final in-memory collision dataset is a **merged and priority-resolved** combination of all `blk_*` files:

- Start with **blk_10**  
- Override with **blk_20**  
- Override with **blk_30**  
- … and so on

This layering mechanism guarantees that:

- We always have *complete* data (from `blk_10`)
- But we also get the *best available* values (from higher-priority files)

---

This behavior is automatic and built into the loader—users only need to supply the available `blk_*` files, and the system will combine them according to the rules described above.

