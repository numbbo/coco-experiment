project('coco', 'c', 'cpp',
  default_options: ['warning_level=3', 'buildtype=release', 'cpp_std=c++20']
  )

cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')

m_dep = cc.find_library('m', required : false)

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

