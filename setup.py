# setup.py  –  RadSmina meta-package

from pathlib import Path
import sys
import tempfile

import numpy as np
import pybind11
from setuptools import Extension, setup, find_packages
from setuptools.command.build_ext import build_ext

################################################################################
# Build helpers (same logic the original rad project used)
################################################################################
def _has_flag(compiler, flag):
    """Return True if the given compiler supports `flag`."""
    with tempfile.NamedTemporaryFile("w", suffix=".cpp", delete=False) as f:
        f.write("int main(){return 0;}")
    try:
        compiler.compile([f.name], extra_postargs=[flag])
    except Exception:
        return False
    finally:
        Path(f.name).unlink(missing_ok=True)
    return True


def _cpp_std_flag(compiler):
    for flag in ("-std=c++17", "-std=c++14", "-std=c++11"):
        if _has_flag(compiler, flag):
            return flag
    raise RuntimeError("A compiler with C++11 support or better is required")


class _BuildExt(build_ext):
    """Add compiler-specific flags automatically."""

    c_opts = {
        "msvc": ["/EHsc", "/O2", "/openmp"],
        "unix": ["-O3", "-march=native"],
    }
    link_opts = {"msvc": [], "unix": []}

    if sys.platform == "darwin":
        c_opts["unix"] += ["-stdlib=libc++", "-mmacosx-version-min=10.13"]
        # (drop OpenMP on macOS unless you know libomp is present)
    else:
        c_opts["unix"].append("-fopenmp")
        link_opts["unix"] += ["-fopenmp", "-pthread"]

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = list(self.c_opts.get(ct, []))

        if ct == "unix":
            opts += [
                _cpp_std_flag(self.compiler),
                f'-DVERSION_INFO="{self.distribution.get_version()}"',
            ]
            if _has_flag(self.compiler, "-fvisibility=hidden"):
                opts.append("-fvisibility=hidden")
        elif ct == "msvc":
            opts += [
                f'/DVERSION_INFO=\\"{self.distribution.get_version()}\\"',
            ]

        for ext in self.extensions:
            ext.extra_compile_args += opts
            ext.extra_link_args += self.link_opts.get(ct, [])

        super().build_extensions()


################################################################################
# Extension module definition
################################################################################
ext_modules = [
    Extension(
        "hnswlib_rad",
        sources=["rad/hnswlib/python_bindings/bindings.cpp"],
        include_dirs=[
            pybind11.get_include(),
            np.get_include(),
            "rad/hnswlib/hnswlib",
        ],
        language="c++",
    )
]

################################################################################
# Setup metadata
################################################################################
ROOT = Path(__file__).parent

setup(
    name="radsmina",                # distribution name on PyPI / pip
    version="0.1.0",
    description="Retrieval-Augmented Docking plus Smina utilities",
    long_description=(ROOT / "README.md").read_text(encoding="utf8")
    if (ROOT / "README.md").exists()
    else "",
    long_description_content_type="text/markdown",
    author="Betty Li",
    url="https://github.com/bl521/RadSmina",
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.19",
        "pybind11>=2.10",

    ],
    packages=find_packages(
        include=[
            "rad",
            "rad.*",
            "radsmina",
            "radsmina.*",
        ]
    ),
    ext_modules=ext_modules,
    cmdclass={"build_ext": _BuildExt},
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
    ],
)
