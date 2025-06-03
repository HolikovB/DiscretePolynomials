# Discrete polynomials

This project defines and manipulates functions from the finite abelian group **G** (mainly **ℤ/4ℤ × ℤ/4ℤ**) to the circle group **ℝ/ℤ** (represented as `Fraction % 1` in Python). In order to investigate discrete polynomials.

It includes tools to:
- Define arbitrary functions (not necessarily homomorphisms)
- Take **directional derivatives**
- Check if a function is a **discrete polynomial** of degree ≤ *k*

---

## 📚 Requirements

- Python 3.x
- No external libraries (just `fractions` and `itertools`)

---

## 🚀 Getting Started

### 1. Clone the repo
```bash
git clone https://github.com/HolikovB/DiscretePolynomials.git
```

## Example of usage 
example_z4x4.ipynb contains Jypiter notebook with some examples of usage of package.


## Current experiments
find_all_poly.ipynb contains list of all 24 polynomials (up to constant) of degree $2$ such that $4f \neq 0$.
