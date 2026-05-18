# =============================================================================
# meson.build — Meson build description for the gwswex package
# =============================================================================
#
# Build pipeline overview
# -----------------------
#  1. constants.f08 → types.f08 → physics.f08 → ... → kernel.f08
#                    compiled into a *static* library (gwswex_core) so that
#                    all generated .mod files exist before f2py runs.
#  2. wrapper.f90   → processed by f2py to generate C glue code, then
#                    compiled together with the static lib into the Python
#                    extension f_gwswex (.pyd on Windows, .so on Linux/macOS).
#  3. Pure Python sources are installed alongside the extension inside
#     the gwswex/ package directory.
# =============================================================================

project('gwswex', 'c', 'fortran',
  version : run_command(
    find_program('python3', 'python'), ['-c', 'import tomllib; data = tomllib.loads(open("' + meson.project_source_root() / 'pyproject.toml' + '").read()); print(data["project"]["version"], end="")'],
    check : true
  ).stdout(),
  default_options : ['warning_level=1']
)

py_mod = import('python')
py = py_mod.find_installation(required : true, pure : false)

incdir_numpy = run_command(py, ['-c', 'import numpy; print(numpy.get_include())'], check : true).stdout().strip()
incdir_f2py  = run_command(py, ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'], check : true).stdout().strip()
fortranobject_c = incdir_f2py + '/fortranobject.c'

# Pass debug flags when building with --buildtype=debug.
debug_fortran_args = []
debug_c_args = []
if get_option('buildtype') == 'debug'
  debug_fortran_args += ['-g', '-O0', '-fbacktrace']
  debug_c_args += ['-g', '-O0']
endif

# OpenMP is optional.
openmp_dep = dependency('openmp', required : false)

# ---- Step 1: Build the pure-Fortran static library ---------------------------
lib_core = static_library('gwswex_core',
  [
    'gwswex/src/shared/constants.f08',
    'gwswex/src/shared/types.f08',
    'gwswex/src/shared/physics.f08',
    'gwswex/src/shared/geometry.f08',
    'gwswex/src/shared/lateral.f08',
    'gwswex/src/explicit/processes.f08',
    'gwswex/src/explicit/timestep.f08',
    'gwswex/src/explicit/solver.f08',
    'gwswex/src/implicit/solver.f08',
    'gwswex/src/shared/mass_balance.f08',
    'gwswex/src/kernel.f08',
  ],
  dependencies : openmp_dep,
  fortran_args : debug_fortran_args,
  install : false
)

# ---- Step 2: Run f2py on the wrapper to generate C glue code -----------------
f2py_c_source = custom_target('f2py_gwswex_c',
  input   : ['gwswex/wrapper.f90'],
  output  : ['f_gwswexmodule.c', 'f_gwswex-f2pywrappers2.f90'],
  depends : lib_core,
  command : [py, '-m', 'numpy.f2py',
             '@INPUT@',
             '-m', 'f_gwswex',
             '--lower',
             '--build-dir', meson.current_build_dir()
            ]
)

# ---- Step 3: Compile everything into the Python extension --------------------
# macOS-specific linker flags
macos_link_args = host_machine.system() == 'darwin' ? ['-Wl,-no_compact_unwind'] : []

py.extension_module('f_gwswex',
  [f2py_c_source, 'gwswex/wrapper.f90', fortranobject_c],
  link_with           : lib_core,
  dependencies        : openmp_dep,
  include_directories : [
    include_directories('gwswex'),
  ],
  c_args      : debug_c_args + ['-I' + incdir_numpy, '-I' + incdir_f2py, '-Wno-misleading-indentation'],
  fortran_args: debug_fortran_args + ['-I' + meson.current_build_dir()],
  link_args   : macos_link_args,
  subdir      : 'gwswex',
  install     : true
)

# ---- Step 4: Install pure-Python sources alongside the extension -------------
py.install_sources([
  'gwswex/__init__.py',
  'gwswex/version.py',
  'gwswex/config.py',
  'gwswex/model.py',
  'gwswex/io.py',
], subdir : 'gwswex')
