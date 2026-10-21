project('coco', 'c', 'cpp',
  default_options: ['warning_level=3', 'buildtype=release', 'cpp_std=c++17']
  )

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)
# Disable floating-point contraction (FMA), see src/coco_fp_contract.h
add_project_arguments(cc.get_supported_arguments('-ffp-contract=off'), language: 'c')
cpp = meson.get_compiler('cpp')
add_project_arguments(cpp.get_supported_arguments('-ffp-contract=off'), language: 'cpp')

coco_lib = static_library('coco', 
  sources: 'coco.cpp',
  dependencies: m_dep
  )

executable('example_experiment', 
  sources: 'example_experiment.c',
  link_with: coco_lib,
  dependencies: m_dep
  )

executable('test_coco', 
  sources: 'test_coco.c',
  link_with: coco_lib,
  dependencies: m_dep
  )

