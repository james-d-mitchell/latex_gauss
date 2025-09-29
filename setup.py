from setuptools import setup

setup(
    name="latex_gauss",
    version="0.1.1",
    description="Package for producing the latex of Gaussian elimination",
    url="http://github.com/james-d-mitchell/latex_gauss",
    author="James D. Mitchell",
    author_email="jdm3@st-andrews.ac.uk",
    license="MIT",
    packages=["latex_gauss"],
    install_requires=[
        "sympy",
    ],
    zip_safe=False,
)
