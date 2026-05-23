from __future__ import annotations

import shlex
import subprocess

from Cython.Build import cythonize
import numpy as np
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as build_ext_orig


def _pkg_config_flags(package: str) -> list[str]:
    try:
        output = subprocess.check_output(
            ["pkg-config", "--cflags", "--libs", package],
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return []
    return shlex.split(output)


class build_ext(build_ext_orig):
    def build_extensions(self) -> None:
        for ext in self.extensions:
            for flag in _pkg_config_flags("fftw3"):
                if flag.startswith("-I"):
                    ext.include_dirs.append(flag[2:])
                elif flag.startswith("-L"):
                    ext.library_dirs.append(flag[2:])
                elif flag.startswith("-l"):
                    ext.libraries.append(flag[2:])
                else:
                    ext.extra_link_args.append(flag)
        super().build_extensions()


extensions = [
    Extension(
        "pygft._core",
        ["pygft/_core.pyx", "gft.c"],
        include_dirs=[np.get_include(), "."],
        libraries=["fftw3"],
    )
]


setup(
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
    cmdclass={"build_ext": build_ext},
)
