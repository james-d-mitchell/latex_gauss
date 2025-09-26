# -*- coding: utf-8 -*-

# Copyright (c) 2025, J. D. Mitchell
#
# Distributed under the terms of the GPL license version 3.
#
# The full license is in the file LICENSE, distributed with this software.

"""
This module contains some tests for the latex_gauss package.
"""

from sympy import Matrix, latex
from fractions import Fraction

from latex_gauss import row_echelon_form, reduced_row_echelon_form


def test_latex_gauss():
    A = Matrix([[1, 2, 1], [2, 3, 4], [3, 1, 0]])
    assert row_echelon_form(A) == (
        Matrix([[1, 2, 1], [0, -1, 2], [0, 0, -13]]),
        "\\left[\\begin{matrix}1 & 2 & 1\\\\2 & 3 & 4\\\\3 & 1 & 0\\end{matrix}\\right]& \\sim \\left[\\begin{matrix}1 & 2 & 1\\\\0 & -1 & 2\\\\3 & 1 & 0\\end{matrix}\\right] && r_2 \\rightarrow r_2 - 2r_1\\\\& \\sim \\left[\\begin{matrix}1 & 2 & 1\\\\0 & -1 & 2\\\\0 & -5 & -3\\end{matrix}\\right] && r_3 \\rightarrow r_3 - 3r_1\\\\& \\sim \\left[\\begin{matrix}1 & 2 & 1\\\\0 & -1 & 2\\\\0 & 0 & -13\\end{matrix}\\right] && r_3 \\rightarrow r_3 - 5r_2\\\\",
    )
