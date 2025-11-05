from skbuild import setup

import sys
_cmake_args = []

# if 'win32' in sys.platform:
#    _cmake_args += ['-G', 'MinGW Makefiles']
# _cmake_args += ['-G', 'Visual Studio 17 2022']
# _cmake_args += ['-G', 'Visual Studio 16 2019']
    
# setup(
#     name="ASCsoft",
#     version="0.0.2",
#     author="J. Schoeberl",
#     license="MIT",
#     packages=["ASCsoft"],
#     cmake_args=_cmake_args
# )

setup(
    name="ASCsoft",
    version="0.0.2",
    description="A demo package",
    author="Joachim Schoeberl",
    license="LGPL-2.1",
    packages=['ASCsoft'],
    package_dir={'ASCsoft': 'ASCsoft'},
    cmake_install_dir="ASCsoft",
    python_requires=">=3.6",
)