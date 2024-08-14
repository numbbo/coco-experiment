import sys
import numpy as np

from pathlib import Path
from setuptools import Extension, setup


## Ensure source files are present:
BASE_DIR = Path(__file__).parent.absolute()
COCOEX_DIR = Path("src") / "cocoex"
MISSING = []
for fname in ["coco.c", "coco.h"]:
    full_fname = BASE_DIR / COCOEX_DIR / fname
    if not full_fname.exists():
        MISSING.append(str(COCOEX_DIR / fname))

if len(MISSING) > 0:
    print(f"""BUILD ERROR

Missing source files: {", ".join(MISSING)}
    
Could not find some generated source files. These are required to build the 
Python package. If you are building in-tree from a git checkout of the 
'numbbo/coco-experiment' repository, did you forget to run the 
`scripts/fabricate` script to generate the COCO source bundles?

If you are not building in-tree or ran `scripts/fabricate`, please file an 
issue on Github: https://github.com/numbbo/coco-experiment/issues
""")
    sys.exit(-1)


## Configure Cython build
extensions = []
extensions.append(Extension(name="cocoex.interface",
                            sources=["src/cocoex/coco.c", "src/cocoex/interface.pyx"],
                            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
                            include_dirs=[np.get_include()]))
extensions.append(Extension(name="cocoex.function",
                            sources=["src/cocoex/coco.c", "src/cocoex/function.pyx"],
                            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
                            include_dirs=[np.get_include()]))

setup(ext_modules=extensions)
