import platform
import numpy as np
from setuptools import Extension, setup


is_arm = "arm" in platform.machine().lower() or "aarch64" in platform.machine().lower()
is_macos = platform.system() == "Darwin"

extra_cflags = []

# macOS on arm64 uses fused multiply adds which changes the results of evaluating 
# a benchmark function enough to make our tests fail. Disable FMA for now.
if is_macos and is_arm:
    extra_cflags.append("-ffp-contract=off")


extensions = []
extensions.append(Extension(name="cocoex.interface",
                            sources=["src/cocoex/coco.c", "src/cocoex/interface.pyx"],
                            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
                            include_dirs=[np.get_include()],
                            extra_compile_args=extra_cflags
                            ))
extensions.append(Extension(name="cocoex.bare_problem",
                            sources=["src/cocoex/coco.c", "src/cocoex/bare_problem.pyx"],
                            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
                            include_dirs=[np.get_include()],
                            extra_compile_args=extra_cflags
                            ))

setup(ext_modules=extensions)
