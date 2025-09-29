from sympy import Matrix, latex
from fractions import Fraction
import re


def latex_scalar(scalar: Fraction) -> str:
    if scalar.denominator != 1:
        if scalar.numerator < 0:
            return f"-\\tfrac{{{-scalar.numerator}}}{{{scalar.denominator}}}"
        return f"\\tfrac{{{scalar.numerator}}}{{{scalar.denominator}}}"
    return f"{scalar}"


def latex_ero3(B: Matrix, row: int, pivot_row: int, scalar: Fraction) -> None:
    if scalar < 0:
        if scalar == -1:
            scalar = "+ "
        else:
            scalar = f"+ {latex_scalar(-scalar)}"
    else:
        if scalar == 1:
            scalar = "- "
        else:
            scalar = f"- {latex_scalar(scalar)}"
    return f"& \\sim {latex(B)} && r_{row + 1} \\rightarrow r_{row + 1} {scalar}r_{pivot_row + 1}\\\\\n"


def row_echelon_form(
    A: Matrix, begin_matrix: str = r"\\begin{bmatrix}", end_matrix=r"\\end{bmatrix}"
) -> tuple[Matrix, str]:
    B = A.copy()
    r, c = B.shape
    pivot_col = 0
    result = latex(B)
    for pivot_row in range(r):
        found_pivot = False
        while not found_pivot and pivot_col < c and B[pivot_row, pivot_col] == 0:
            for i in range(pivot_row + 1, r):
                if B[i, pivot_col] != 0:
                    B.row_swap(pivot_row, i)
                    result += f"& \\sim {latex(B)} && r_{pivot_row + 1} \\leftrightarrow r_{i + 1} \\\\\n"
                    found_pivot = True
                    break
            else:
                pivot_col += 1
        if pivot_col == c:
            return B
        pivot = B[pivot_row, pivot_col]
        if pivot == 0:
            return B
        for i in range(pivot_row + 1, r):
            if B[i, pivot_col] != 0:
                scalar = Fraction(B[i, pivot_col], pivot)
                for j in range(pivot_col, c):
                    B[i, j] -= scalar * B[pivot_row, j]
                result += latex_ero3(B, i, pivot_row, scalar)
    result = re.sub(r"\\left\[\\begin{matrix}", begin_matrix, result)
    result = re.sub(r"\\end{matrix}\\right]", end_matrix, result)
    return (B, result)


def reduced_row_echelon_form(
    A: Matrix, begin_matrix: str = r"\\begin{bmatrix}", end_matrix=r"\\end{bmatrix}"
) -> tuple[Matrix, str]:
    B, result = row_echelon_form(A, begin_matrix, end_matrix)
    r, c = B.shape
    for pivot_row in range(r - 1, -1, -1):
        for pivot_col in range(0, c):
            pivot = B[pivot_row, pivot_col]
            if pivot == 0:
                continue
            for i in range(pivot_col, c):
                B[pivot_row, i] *= Fraction(1, pivot)
            if pivot != 1:
                if pivot == -1:
                    scalar = "-"
                else:
                    scalar = latex_scalar(Fraction(1, pivot))
                result += f"& \\sim {latex(B)} && r_{pivot_row + 1} \\rightarrow {scalar}r_{pivot_row + 1}\\\\\n"

            for i in range(pivot_row - 1, -1, -1):
                scalar = B[i, pivot_col]
                if scalar != 0:
                    for j in range(pivot_col, c):
                        B[i, j] -= B[pivot_row, j] * scalar
                    result += latex_ero3(B, i, pivot_row, scalar)

            break
    result = re.sub(r"\\left\[\\begin{matrix}", begin_matrix, result)
    result = re.sub(r"\\end{matrix}\\right]", end_matrix, result)
    return (B, result)
