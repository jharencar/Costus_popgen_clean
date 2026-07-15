# Chemistry QTL — trait prioritization

Source: [`../data/Chem_QTL_mass_WATERS.csv`](../data/Chem_QTL_mass_WATERS.csv) (146 Waters mass features + NPC ontology + collaborator Notes).

**QTL mapping is on hold** until you confirm which tier/shortlist to use. Do not start phenotype merge / Farm scantwo until then.

## Priority tables

| File | Contents |
|------|----------|
| [`../data/Chem_QTL_trait_priority.csv`](../data/Chem_QTL_trait_priority.csv) | All 146 features with `tier` A/B/C, pathway annotations, LOD |
| [`../data/Chem_QTL_TierA_alkaloids_flavonoids.csv`](../data/Chem_QTL_TierA_alkaloids_flavonoids.csv) | Tier A only (54): Alkaloids + Flavonoids |
| [`../data/Chem_QTL_recommended_shortlist.csv`](../data/Chem_QTL_recommended_shortlist.csv) | Suggested ~17: top-10 LOD imidazole alkaloids + 7 noted flavonoids |

Trust hierarchy (collaborator): **Pathway > Superclass > Class**. Do not treat Class as a firm compound ID.

## Tiers (summary)

- **A (54):** Alkaloids (41) + Flavonoids (13; 7 with structural Notes)
- **B (5):** Amino acids and Peptides (alkaloid-adjacent optional)
- **C (87):** Unclassified + other pathways (fatty acids, carbs, terpenoids, etc.)

## Sample size (expected at merge)

~130 chemistry samples → ~104 Kay Lab_ID → **~83** in F3phe / cross (subset of 385). Reconfirm when prep runs.

## After confirmation

Update prep + QTL scripts under this folder to keep only confirmed `row_ID`s; Farm scantwo for significant traits only.
