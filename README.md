# Loggerhead Sea Turtle Population Dynamics

A stage-structured matrix population model analyzing the long-term dynamics of loggerhead sea turtles. The model represents the population as a 7×7 projection matrix and uses eigenvalue analysis to determine whether the population grows or declines over time, then quantifies which life stages most influence that outcome.

**Language:** Python · **Methods:** Linear algebra (eigenvalue analysis, diagonalization), population modeling

---

## Summary

Loggerhead sea turtles pass through several distinct life stages — hatchlings, juveniles, subadults, and multiple breeding stages — each with its own survival and reproduction rates. A single-stage growth model cannot capture this structure, so the population is modeled as a **stage-structured (Lefkovitch) matrix**, where a state vector tracks the number of turtles in each stage and a projection matrix advances the population one year at a time.

The central question: **how do stage-specific survival and reproduction rates shape long-term population growth?**

**Key finding.** The baseline projection matrix has a dominant eigenvalue of **λ ≈ 0.945**. Because this is below 1, the population declines over the long run regardless of its starting composition. A sensitivity analysis shows that raising survival in the **late-juvenile and subadult stages** lifts the dominant eigenvalue above 1 (to ≈ 1.02) — turning decline into growth. This identifies those stages, rather than hatchling abundance, as the highest-leverage targets for conservation.

## Results

| Scenario | Dominant eigenvalue (λ) | Long-term outcome |
|---|---|---|
| Baseline | 0.945 | Decline |
| Improved late-juvenile / subadult survival | 1.02 | Growth |

## Methods

The analysis applies core linear-algebra techniques to a real ecological dataset:

- **Projection matrix.** A 7×7 matrix encodes fecundity (first row), survival within a stage (diagonal), and survival with advancement to the next stage (subdiagonal).
- **Population projection.** The recurrence `x_{t+1} = A·x_t` advances the population one year; `x_{t+k} = Aᵏ·x_t` projects k years ahead.
- **Eigenvalue analysis.** The dominant eigenvalue of A is the asymptotic growth factor — λ < 1 implies decline, λ > 1 implies growth.
- **Stable stage distribution.** The dominant eigenvector gives the fixed proportion of turtles across stages that the population converges to.
- **Diagonalization.** Writing `A = SΛS⁻¹` shows why long-term behavior is governed by the dominant eigenvalue and is independent of the initial population.
- **Sensitivity analysis.** Perturbing stage-specific survival rates quantifies each stage's leverage on λ.

Matrix entries are drawn from Crouse, Crowder & Caswell (1987); all values reflect the published demographic data.

## Repository Contents

```
├── model.py           # Projection matrix, eigenvalue and sensitivity analysis
├── report.pdf         # Full written report
├── requirements.txt   # Dependencies
└── README.md
```

## Running the Model

```bash
pip install -r requirements.txt
python model.py
```

Output includes the baseline dominant eigenvalue, the stable stage distribution, short-term population projections, and the improved-survival sensitivity comparison — reproducing the figures reported in `report.pdf`.

## Report

The full write-up — mathematical formulation, worked projections, diagonalization, and ecological discussion — is available in **[report.pdf](report.pdf)**.

## Contributors

Completed for APPM 3310 (Matrix Methods), Section 002.

| Contributor | Primary focus |
|---|---|
| **Sanvitha Vallem** | Eigenvalue analysis and long-term behavior |
| **Sanjitha Chakka** | Background research and source data |
| **Mahi Chalicham** | Numerical analysis and matrix computations |

All contributors participated in writing, editing, verifying the mathematics, and developing the code.

## Reference

Crouse, D. T., L. B. Crowder, and H. Caswell (1987). "A Stage-Based Population Model for Loggerhead Sea Turtles and Implications for Conservation." *Ecology* 68(5): 1412–1423.
