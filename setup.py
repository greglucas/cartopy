# Copyright Cartopy Contributors
#
# This file is part of Cartopy and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import os
from pathlib import Path

import numpy as np
from setuptools import Extension, setup


# The existence of a PKG-INFO directory is enough to tell us whether this is a
# source installation or not (sdist).
HERE = Path(__file__).parent
IS_SDIST = (HERE / 'PKG-INFO').exists()
FORCE_CYTHON = os.environ.get('FORCE_CYTHON', False)

if not IS_SDIST or FORCE_CYTHON:
    import Cython
    if Cython.__version__ < '0.29':
        raise ImportError(
            "Cython 0.29+ is required to install cartopy from source.")

    from Cython.Distutils import build_ext as cy_build_ext
    ext = '.pyx'
    cmdclass = {'build_ext': cy_build_ext}
else:
    ext = '.cpp'
    cmdclass = {}

# General extension paths
compiler_directives = {}
define_macros = []

compiler_directives['profile'] = True
compiler_directives['linetrace'] = True


cython_coverage_enabled = os.environ.get('CYTHON_COVERAGE', False)
if cython_coverage_enabled:
    define_macros.append(('CYTHON_TRACE_NOGIL', '1'))

extensions = [
    Extension(
        'cartopy.trace',
        [f'lib/cartopy/trace{ext}'],
        include_dirs=[np.get_include()],
        language='c++',
        define_macros=define_macros),
]


if cython_coverage_enabled:
    # We need to explicitly cythonize the extension in order
    # to control the Cython compiler_directives.
    from Cython.Build import cythonize
    extensions = cythonize(extensions, compiler_directives=compiler_directives)


# Main setup
# ==========
setup(
    ext_modules=extensions,
    cmdclass=cmdclass,

)
