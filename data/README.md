# Data

## What ships in this repository

`derived/paper1_data/` — the **derived** data products: per-epoch aggregated Moran's *I*
values (`I_gnd`, `I_120`, `I_300`, `Δ = I_120 − I_300`, auroral fraction `f`, station
counts, ...) for every event in the paper, plus quiet-time records and a cross-event
summary. These are sufficient to regenerate the figures and summary statistics
(e.g. `code/make_fig1b_clean.py` reads `derived/paper1_data/halloween_2003_altitude_data.json`).

These products are released under CC-BY-4.0.

## What you must fetch (raw inputs are NOT redistributed here)

The analysis loaders (`code/loader_supermag.py`, `code/loader_tec.py`) expect raw inputs
under this `data/` directory. Fetch them from the original providers:

| Source | Used for | Where the loader expects it | How to obtain |
|--------|----------|-----------------------------|---------------|
| **SuperMAG** ground magnetometer data | `I_gnd`, `I_120` inputs | `data/supermag/` (`*_YYYYMM_full.npz`, `station_coords.json`, `station_ranges.json`) | Register at https://supermag.jhuapl.edu and download per its data policy. **SuperMAG data may not be redistributed**, which is why it is not bundled here. |
| **CODE** Global Ionosphere Maps (GPS-TEC) | `I_300`, `I_120` (Chapman split) | `data/code_gims/` | Freely available from AIUB/CODE (ftp://ftp.aiub.unibe.ch/CODE/ or the IGS ionosphere products). Cite Schaer (1999). |
| **Swarm** (ESA) | Greenland cross-validation (FAC, AOB, Langmuir) | `data/satellite/swarm/` | ESA VirES (https://vires.services), open license. |
| **Van Allen Probes** (NASA) | Greenland cross-validation (HOPE composition) | `data/satellite/rbsp/` | NASA CDAWeb (https://cdaweb.gsfc.nasa.gov), public domain. |

Once the raw data are in place, `code/export_altitude_data.py` regenerates the derived
products, and `code/test_coordinate.py` / `test_persector.py` / `test_robustness.py`
reproduce the robustness results.

## Acknowledgment requirements

Please honor each provider's acknowledgment text. In particular, SuperMAG requests a
specific acknowledgment of the collaborators and funding for the ground magnetometer
arrays it aggregates; see the SuperMAG website for the current wording.
