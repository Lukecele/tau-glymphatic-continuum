from setuptools import setup, find_packages

setup(
    name="tau-glymphatic-continuum",
    version="1.1.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    entry_points={
        "console_scripts": [
            "tau-calc=tau_glymphatic.cli:main",
        ],
    },
)
