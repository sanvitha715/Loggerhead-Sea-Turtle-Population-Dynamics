"""
Loggerhead Sea Turtle Population Dynamics
=========================================

A stage-structured (Lefkovitch) matrix model of a loggerhead sea turtle
population. The population is divided into seven life stages and advanced one
year at a time by a 7x7 projection matrix A:

    x_{t+1} = A @ x_t

The dominant eigenvalue of A is the asymptotic growth factor: values below 1
imply long-term decline, values above 1 imply growth. This script computes the
baseline dominant eigenvalue, the stable stage distribution, short-term
projections, and a sensitivity analysis on stage-specific survival.

Matrix entries are from:
    Crouse, D. T., L. B. Crowder, and H. Caswell (1987).
    "A Stage-Based Population Model for Loggerhead Sea Turtles and
    Implications for Conservation." Ecology 68(5): 1412-1423.

Baseline result: dominant eigenvalue lambda ~= 0.945 (long-term decline).
"""

from __future__ import annotations

import numpy as np


STAGE_NAMES: list[str] = [
    "Hatchlings",
    "Small juveniles",
    "Large juveniles",
    "Subadults",
    "Novice breeders",
    "First-year remigrants",
    "Mature breeders",
]


def build_projection_matrix() -> np.ndarray:
    """Construct the 7x7 baseline projection matrix.

    Matrix structure:
        - Row 0 (first row): fecundities; only breeding stages 5-7 reproduce.
        - Diagonal:          probability of surviving and staying in a stage.
        - Subdiagonal:       probability of surviving and advancing a stage.
        - All other entries: zero (turtles cannot skip stages).

    Returns
    -------
    np.ndarray
        The 7x7 projection matrix A.
    """
    return np.array([
        [0,      0,      0,      0,      127,    4,      80    ],
        [0.6747, 0.7370, 0,      0,      0,      0,      0     ],
        [0,      0.0486, 0.6610, 0,      0,      0,      0     ],
        [0,      0,      0.0147, 0.6907, 0,      0,      0     ],
        [0,      0,      0,      0.0518, 0,      0,      0     ],
        [0,      0,      0,      0,      0.8091, 0,      0     ],
        [0,      0,      0,      0,      0,      0.8091, 0.8089],
    ])


def dominant_eigen(A: np.ndarray) -> tuple[float, np.ndarray]:
    """Return the dominant eigenvalue and stable stage distribution of A.

    The stable stage distribution is the dominant eigenvector expressed as
    non-negative proportions summing to 1.

    Parameters
    ----------
    A : np.ndarray
        A square projection matrix.

    Returns
    -------
    tuple[float, np.ndarray]
        The dominant (largest-magnitude) eigenvalue and the corresponding
        normalized stage distribution.
    """
    eigenvalues, eigenvectors = np.linalg.eig(A)
    dominant_index = int(np.argmax(np.abs(eigenvalues)))
    eigenvalue = float(eigenvalues[dominant_index].real)

    distribution = np.abs(eigenvectors[:, dominant_index].real)
    distribution /= distribution.sum()
    return eigenvalue, distribution


def project(A: np.ndarray, x0: np.ndarray, years: int) -> np.ndarray:
    """Project the population forward, returning each yearly state as a row.

    Parameters
    ----------
    A : np.ndarray
        The projection matrix.
    x0 : np.ndarray
        The initial population vector.
    years : int
        Number of years (matrix applications) to project.

    Returns
    -------
    np.ndarray
        An array of shape (years + 1, len(x0)); row k is the population at year k.
    """
    states = [np.asarray(x0, dtype=float)]
    for _ in range(years):
        states.append(A @ states[-1])
    return np.array(states)


def improved_survival_scenario(A: np.ndarray, boost: float = 0.20) -> np.ndarray:
    """Return a copy of A with late-juvenile and subadult survival increased.

    Scales the large-juvenile and subadult survival and advancement terms by
    ``(1 + boost)``. A 20% boost raises the dominant eigenvalue above 1,
    demonstrating that these stages have the greatest leverage on long-term
    growth.

    Parameters
    ----------
    A : np.ndarray
        The baseline projection matrix.
    boost : float, optional
        Fractional increase applied to the targeted survival terms.

    Returns
    -------
    np.ndarray
        The modified projection matrix.
    """
    modified = A.copy()
    modified[2, 2] *= (1 + boost)  # large juvenile: stay
    modified[3, 2] *= (1 + boost)  # large juvenile -> subadult
    modified[3, 3] *= (1 + boost)  # subadult: stay
    modified[4, 3] *= (1 + boost)  # subadult -> novice breeder
    return modified


def main() -> None:
    A = build_projection_matrix()

    print("Loggerhead Sea Turtle Population Model")
    print("=" * 55)

    # --- Baseline -----------------------------------------------------------
    eigenvalue, distribution = dominant_eigen(A)
    outcome = "decline" if eigenvalue < 1 else "growth"
    print("\nBaseline scenario")
    print("-" * 55)
    print(f"Dominant eigenvalue (lambda): {eigenvalue:.4f}  ->  long-term {outcome}")

    print("\nStable stage distribution:")
    for name, proportion in zip(STAGE_NAMES, distribution):
        print(f"  {name:<24} {proportion:7.2%}")

    # --- Short-term projection ---------------------------------------------
    print("\nShort-term projection from x0 = [1, 1, 1, 1, 1, 1, 1]:")
    states = project(A, np.ones(7), years=2)
    print(f"  A  x0 = {np.round(states[1], 4)}")
    print(f"  A^2 x0 = {np.round(states[2], 4)}")

    # --- Sensitivity analysis ----------------------------------------------
    modified = improved_survival_scenario(A, boost=0.20)
    improved_eigenvalue, _ = dominant_eigen(modified)
    improved_outcome = "growth" if improved_eigenvalue > 1 else "decline"
    print("\nSensitivity: improved late-juvenile / subadult survival")
    print("-" * 55)
    print(f"  Baseline lambda: {eigenvalue:.4f}  ->  decline")
    print(f"  Improved lambda: {improved_eigenvalue:.4f}  ->  {improved_outcome}")
    print(
        "\nConservation implication: survival gains in the late-juvenile and\n"
        "subadult stages have the largest effect on long-term growth,\n"
        "consistent with Crouse, Crowder & Caswell (1987)."
    )


if __name__ == "__main__":
    main()
