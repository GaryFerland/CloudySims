readme for molecules
# How to Add Molecules, Chemical Reactions, and Molecular Lines in CLOUDY

## Gargi Shaw

### 1. Adding Molecular Reactions
- Add molecular reactions, document the source and date in the `mole_base_CO.dat` file located in the `data` directory.
- Use `hmrate`, `crnurate`, and `th85rate` for two-body reactions, `CRPHOT`, and `PHOTON` reactions, respectively.
- Note that these reactions are not (V,J) resolved.

**Example:**
```plaintext
O,C2=>CO,C:hmrate:2.e-10:-.12:0 # UMIST2012, updated on 2023may08,GS .
OCS,PHOTON=>S,CO:th85rate:3.7e-9:0.:0. # UMIST
OCS,CRPHOT=>S,CO:crnurate:5360 :0.:0. # UMIST

## 1. Updating Molecular Species
- Ensure the molecular species is included in the `chem_species.dat` file in the `data` directory. 
- Add the formation enthalpy at 0K in kJ/mol. 

---

## 2. Testing for Errors
- Run `test.in` to check for duplicate reactions, charge imbalances, or elemental imbalances. 

---

## 3. Handling Special Elements
- Additional reaction files in the `data` directory include: 
  - `mole_ti.dat` (Ti-related reactions) 
  - `mole_lithium.dat` (Li-related reactions) 
  - `mole_deuterium.dat` (D-related reactions) 
  - `mole_misc.dat` (miscellaneous reactions, including Ar) 

---

## 4. Adding Molecular Lines
- Provide energy levels and collisional rates for a species in **LAMDA** format. 
- Store the file in the `lamda` subdirectory of the `data` directory. 
- List the filename in `Lamda.ini`, located in `data/lamda/masterlist`. 

---

## 5. Adding Grain Reactions
- Include grain reactions in the `mole_co_base.dat` file in the `data` directory. 
- Ensure the grain species are listed in `chem_species_grn.dat` in the `data` directory. 
