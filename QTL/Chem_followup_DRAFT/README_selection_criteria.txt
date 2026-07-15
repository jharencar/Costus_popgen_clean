Chem QTL follow-up trait selection (DRAFT for review)
=====================================================

Criteria
--------
1. Seed: Alkaloid OR Flavonoid with scanone max LOD > approximate 15% perm threshold
   (more lenient than prior 10% list).
2. Flavonoids: ONLY if they meet that 15% threshold (no longer keep all flavonoids;
   no coloc rescue for flavonoids).
3. Extra alkaloids: peak within +/- 5 cM of any seed peak (was +/- 15 cM).

Approximate 15% threshold
-------------------------
Chem_260703 did not save alpha=0.15 cutoffs (only 10% for traits that passed, and 50%
used internally for shared peaks).

Draft approximation:
  cutoff_15 ~= 0.9274 * cutoff_10
0.9274 = median (Gumbel-interpolated c15 / c10) from prelim UV chemistry cutoffs
(10% and 50%, N_PERM=100).

- Traits in significant_10pct.csv: use saved cutoff_10.
- Other alk/flav: provisional cutoff_10 = median WatersMass 10% cutoff (3.4372)
  -> cutoff_15 ~= 3.1877.

Flagged as approximate; recompute summary(perm, alpha=0.15) before locking for mapping.

Counts
------
Total: 32  (prior draft was 38)
  Alkaloids: 24
  Flavonoids: 8
  Seeds (~15%): 32
  Coloc+/-5cM alkaloids: 7
