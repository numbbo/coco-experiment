project(
  'coco-experiment',
  'c',
  'cpp',
  'cython',
  version: run_command(['python', 'src/cocoex/_version.py'], check: true).stdout().strip(),
  license: 'BSD-3',
  meson_version: '>= 1.5.0',
  default_options: [
    'buildtype=debugoptimized',
    'b_ndebug=if-release',
    'c_std=c17',
    'cpp_std=c++17',
  ],
)

################################################################################
## Find Python
py3 = import('python').find_installation(pure: false)
py3_dep = py3.dependency()

################################################################################
## Configure Cython
cy = meson.get_compiler('cython')
cython = find_program(cy.cmd_array()[0])

################################################################################
# Configure C compiler
cc = meson.get_compiler('c')
c_args = cc.get_supported_arguments('-Wno-unused-but-set-variable', '-Wno-unused-function', '-ffp-contract=off')
add_project_arguments(c_args, language: 'c')

m_dep = cc.find_library('m', required: false)
if m_dep.found()
  add_project_link_arguments('-lm', language: 'c')
endif

################################################################################
## Configure C++ Compiler
##
## NOTE: Currently unused, but will be once we merge the cpp branch
cpp = meson.get_compiler('cpp')
cpp_args = cpp.get_supported_arguments('-Wno-unused-but-set-variable', '-Wno-unused-function', '-ffp-contract=off')
add_project_arguments(cpp_args, language: 'cpp')

################################################################################
## Find numpy headers
np_incdir = meson.get_external_property('numpy-include-dir', 'not-given')
if np_incdir == 'not-given'
  np_incdir = run_command(
    py3,
    [
      '-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())',
    ],
    check: true,
  ).stdout().strip()
endif
np_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_23_API_VERSION'
np_dep = declare_dependency(include_directories: np_incdir, compile_args: np_nodepr_api)

################################################################################
## Build cocoex module

python_sources = [
  'src/cocoex/__init__.py',
  'src/cocoex/_version.py',
  'src/cocoex/exceptions.py',
  'src/cocoex/solvers.py',
  'src/cocoex/utilities.py',
]

py3.install_sources(python_sources, subdir: 'cocoex')

inc_coco = include_directories('src/cocoex')

py3.extension_module(
  'interface',
  ['src/cocoex/interface.pyx', 'src/cocoex/coco.c'],
  include_directories: [inc_coco],
  dependencies: np_dep,
  install: true,
  subdir: 'cocoex',
)

py3.extension_module(
  'bare_problem',
  ['src/cocoex/bare_problem.pyx', 'src/cocoex/coco.c'],
  include_directories: [inc_coco],
  dependencies: np_dep,
  install: true,
  subdir: 'cocoex',
)
