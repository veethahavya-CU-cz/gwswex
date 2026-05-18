! kernel_bmi.f08 -- BMI 2.0 Fortran shim for the GWSWEX kernel.
!
!! Forward-looking stub. None of the procedures below are wired up to the
!! kernel yet; each `error stop`s with a message identifying it as a stub.
!! The contract follows the CSDMS Basic Model Interface 2.0 Fortran
!! specification (https://bmi.csdms.io/en/stable/bmi.spec.html). When
!! implemented, this module will provide a parallel entry surface to
!! `gwswex_wrapper` that BMI-based coupling frameworks (NextGen, PyMT,
!! BMI-Tester, WRF-Hydro coupler, etc.) can drive without going through
!! the bespoke f2py API.
!!
!! Naming: the Fortran BMI spec uses `bmi_<verb>` for control functions and
!! `bmi_get_<thing>` / `bmi_set_<thing>` for variable / grid accessors. We
!! follow that convention here so the eventual production version drops
!! into a `bmi-fortran` derived-type extension cleanly.
!!
!! Implementation plan:
!!   1. Wire `bmi_initialize` to read a YAML/TOML config and call the
!!      existing `kernel_init` chain (will require a small Fortran YAML
!!      reader or, more likely, restrict input to a flat namelist file).
!!   2. Wire `bmi_update` to a single `kernel_step` call using forcing
!!      arrays held in module-private state.
!!   3. Map BMI variable names (CSDMS Standard Names) to existing
!!      `kernel_get_*` / `kernel_set_*` getters/setters.
!!   4. Expose two grids: id 0 = element-only `(ne)`; id 1 = layered
!!      `(nl, ne)`.
!!
!! See `gwswex/coupler.py::BmiGwswex` for the parallel Python adapter and
!! `todo.md` for the prioritised work breakdown.
module gwswex_kernel_bmi
  use gwswex_constants, only: dp
  use gwswex_kernel,    only: kernel_init, kernel_step, kernel_deinit, &
                              kernel_get_gw, kernel_get_sw, kernel_get_uz
  implicit none
  private

  ! Public BMI 2.0 control surface
  public :: bmi_initialize, bmi_update, bmi_update_until, bmi_finalize

  ! Public BMI 2.0 info surface
  public :: bmi_get_component_name
  public :: bmi_get_input_item_count, bmi_get_output_item_count
  public :: bmi_get_input_var_names,  bmi_get_output_var_names

  ! Public BMI 2.0 variable-info surface
  public :: bmi_get_var_grid, bmi_get_var_type, bmi_get_var_units
  public :: bmi_get_var_itemsize, bmi_get_var_nbytes, bmi_get_var_location

  ! Public BMI 2.0 time surface
  public :: bmi_get_current_time, bmi_get_start_time, bmi_get_end_time
  public :: bmi_get_time_units, bmi_get_time_step

  ! Public BMI 2.0 getter/setter surface
  public :: bmi_get_value, bmi_get_value_ptr, bmi_get_value_at_indices
  public :: bmi_set_value, bmi_set_value_at_indices

  ! Public BMI 2.0 grid surface
  public :: bmi_get_grid_rank, bmi_get_grid_size, bmi_get_grid_type
  public :: bmi_get_grid_shape, bmi_get_grid_spacing, bmi_get_grid_origin
  public :: bmi_get_grid_x, bmi_get_grid_y, bmi_get_grid_z
  public :: bmi_get_grid_node_count, bmi_get_grid_edge_count, bmi_get_grid_face_count
  public :: bmi_get_grid_edge_nodes, bmi_get_grid_face_edges
  public :: bmi_get_grid_face_nodes, bmi_get_grid_nodes_per_face

  character(len=*), parameter :: COMPONENT_NAME = "GWSWEX"

contains

  ! ───────────────────────────────────────────────────────────────────────────
  ! Control functions
  ! ───────────────────────────────────────────────────────────────────────────

  subroutine bmi_initialize(config_file, ierr)
    !! Read `config_file`, allocate kernel state, set IC, attach forcings.
    character(len=*), intent(in)  :: config_file
    integer,          intent(out) :: ierr
    error stop "bmi_initialize: stub. See gwswex/src/kernel_bmi.f08 and todo.md."
  end subroutine bmi_initialize

  subroutine bmi_update(ierr)
    !! Advance one macro time-step using forcing held in module-private state.
    integer, intent(out) :: ierr
    error stop "bmi_update: stub."
  end subroutine bmi_update

  subroutine bmi_update_until(time_target, ierr)
    !! Advance until model time reaches `time_target` (seconds).
    real(dp), intent(in)  :: time_target
    integer,  intent(out) :: ierr
    error stop "bmi_update_until: stub."
  end subroutine bmi_update_until

  subroutine bmi_finalize(ierr)
    !! Release all kernel resources (delegates to `kernel_deinit`).
    integer, intent(out) :: ierr
    error stop "bmi_finalize: stub."
  end subroutine bmi_finalize

  ! ───────────────────────────────────────────────────────────────────────────
  ! Info functions
  ! ───────────────────────────────────────────────────────────────────────────

  subroutine bmi_get_component_name(name)
    character(len=*), intent(out) :: name
    name = COMPONENT_NAME
  end subroutine bmi_get_component_name

  subroutine bmi_get_input_item_count(count)
    integer, intent(out) :: count
    count = 5  ! precip, pet, ptt, lat_gw, lat_sw
  end subroutine bmi_get_input_item_count

  subroutine bmi_get_output_item_count(count)
    integer, intent(out) :: count
    count = 6  ! theta, SW, GWH, UZ_storage, recharge_rate, runoff_rate
  end subroutine bmi_get_output_item_count

  subroutine bmi_get_input_var_names(names, ierr)
    !! Populate `names` with CSDMS Standard Names (see Python `BmiGwswex`).
    character(len=*), intent(out) :: names(:)
    integer,          intent(out) :: ierr
    error stop "bmi_get_input_var_names: stub."
  end subroutine bmi_get_input_var_names

  subroutine bmi_get_output_var_names(names, ierr)
    character(len=*), intent(out) :: names(:)
    integer,          intent(out) :: ierr
    error stop "bmi_get_output_var_names: stub."
  end subroutine bmi_get_output_var_names

  ! ───────────────────────────────────────────────────────────────────────────
  ! Variable-info functions
  ! ───────────────────────────────────────────────────────────────────────────

  subroutine bmi_get_var_grid(name, grid_id, ierr)
    !! Return 0 for element-only `(ne)` variables; 1 for layered `(nl, ne)`.
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: grid_id
    integer,          intent(out) :: ierr
    error stop "bmi_get_var_grid: stub."
  end subroutine bmi_get_var_grid

  subroutine bmi_get_var_type(name, vtype, ierr)
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: vtype
    integer,          intent(out) :: ierr
    error stop "bmi_get_var_type: stub."
  end subroutine bmi_get_var_type

  subroutine bmi_get_var_units(name, units, ierr)
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: units
    integer,          intent(out) :: ierr
    error stop "bmi_get_var_units: stub."
  end subroutine bmi_get_var_units

  subroutine bmi_get_var_itemsize(name, itemsize, ierr)
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: itemsize
    integer,          intent(out) :: ierr
    error stop "bmi_get_var_itemsize: stub."
  end subroutine bmi_get_var_itemsize

  subroutine bmi_get_var_nbytes(name, nbytes, ierr)
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: nbytes
    integer,          intent(out) :: ierr
    error stop "bmi_get_var_nbytes: stub."
  end subroutine bmi_get_var_nbytes

  subroutine bmi_get_var_location(name, location, ierr)
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: location
    integer,          intent(out) :: ierr
    location = "node"
    ierr = 0
  end subroutine bmi_get_var_location

  ! ───────────────────────────────────────────────────────────────────────────
  ! Time functions
  ! ───────────────────────────────────────────────────────────────────────────

  subroutine bmi_get_current_time(t, ierr)
    real(dp), intent(out) :: t
    integer,  intent(out) :: ierr
    error stop "bmi_get_current_time: stub."
  end subroutine bmi_get_current_time

  subroutine bmi_get_start_time(t, ierr)
    real(dp), intent(out) :: t
    integer,  intent(out) :: ierr
    error stop "bmi_get_start_time: stub."
  end subroutine bmi_get_start_time

  subroutine bmi_get_end_time(t, ierr)
    real(dp), intent(out) :: t
    integer,  intent(out) :: ierr
    error stop "bmi_get_end_time: stub."
  end subroutine bmi_get_end_time

  subroutine bmi_get_time_units(units)
    character(len=*), intent(out) :: units
    units = "s"
  end subroutine bmi_get_time_units

  subroutine bmi_get_time_step(dt, ierr)
    real(dp), intent(out) :: dt
    integer,  intent(out) :: ierr
    error stop "bmi_get_time_step: stub."
  end subroutine bmi_get_time_step

  ! ───────────────────────────────────────────────────────────────────────────
  ! Getter / setter functions
  ! ───────────────────────────────────────────────────────────────────────────

  subroutine bmi_get_value(name, dest, ierr)
    character(len=*), intent(in)  :: name
    real(dp),         intent(out) :: dest(:)
    integer,          intent(out) :: ierr
    error stop "bmi_get_value: stub."
  end subroutine bmi_get_value

  subroutine bmi_get_value_ptr(name, ierr)
    !! Pointer-return semantics in Fortran require `c_loc`/`c_f_pointer` if
    !! the caller is C/Python. Production version will need an interface
    !! variant per language; this stub is a placeholder.
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: ierr
    error stop "bmi_get_value_ptr: stub."
  end subroutine bmi_get_value_ptr

  subroutine bmi_get_value_at_indices(name, dest, inds, ierr)
    character(len=*), intent(in)  :: name
    real(dp),         intent(out) :: dest(:)
    integer,          intent(in)  :: inds(:)
    integer,          intent(out) :: ierr
    error stop "bmi_get_value_at_indices: stub."
  end subroutine bmi_get_value_at_indices

  subroutine bmi_set_value(name, src, ierr)
    character(len=*), intent(in)  :: name
    real(dp),         intent(in)  :: src(:)
    integer,          intent(out) :: ierr
    error stop "bmi_set_value: stub."
  end subroutine bmi_set_value

  subroutine bmi_set_value_at_indices(name, inds, src, ierr)
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: inds(:)
    real(dp),         intent(in)  :: src(:)
    integer,          intent(out) :: ierr
    error stop "bmi_set_value_at_indices: stub."
  end subroutine bmi_set_value_at_indices

  ! ───────────────────────────────────────────────────────────────────────────
  ! Grid functions
  ! ───────────────────────────────────────────────────────────────────────────

  subroutine bmi_get_grid_rank(grid_id, rank, ierr)
    integer, intent(in)  :: grid_id
    integer, intent(out) :: rank
    integer, intent(out) :: ierr
    error stop "bmi_get_grid_rank: stub."
  end subroutine bmi_get_grid_rank

  subroutine bmi_get_grid_size(grid_id, size_out, ierr)
    integer, intent(in)  :: grid_id
    integer, intent(out) :: size_out
    integer, intent(out) :: ierr
    error stop "bmi_get_grid_size: stub."
  end subroutine bmi_get_grid_size

  subroutine bmi_get_grid_type(grid_id, gtype, ierr)
    integer,          intent(in)  :: grid_id
    character(len=*), intent(out) :: gtype
    integer,          intent(out) :: ierr
    gtype = "rectilinear"
    ierr  = 0
  end subroutine bmi_get_grid_type

  subroutine bmi_get_grid_shape(grid_id, shape, ierr)
    integer, intent(in)  :: grid_id
    integer, intent(out) :: shape(:)
    integer, intent(out) :: ierr
    error stop "bmi_get_grid_shape: stub."
  end subroutine bmi_get_grid_shape

  subroutine bmi_get_grid_spacing(grid_id, spacing, ierr)
    integer,  intent(in)  :: grid_id
    real(dp), intent(out) :: spacing(:)
    integer,  intent(out) :: ierr
    error stop "bmi_get_grid_spacing: stub."
  end subroutine bmi_get_grid_spacing

  subroutine bmi_get_grid_origin(grid_id, origin, ierr)
    integer,  intent(in)  :: grid_id
    real(dp), intent(out) :: origin(:)
    integer,  intent(out) :: ierr
    error stop "bmi_get_grid_origin: stub."
  end subroutine bmi_get_grid_origin

  subroutine bmi_get_grid_x(grid_id, x, ierr)
    integer,  intent(in)  :: grid_id
    real(dp), intent(out) :: x(:)
    integer,  intent(out) :: ierr
    error stop "bmi_get_grid_x: stub."
  end subroutine bmi_get_grid_x

  subroutine bmi_get_grid_y(grid_id, y, ierr)
    integer,  intent(in)  :: grid_id
    real(dp), intent(out) :: y(:)
    integer,  intent(out) :: ierr
    error stop "bmi_get_grid_y: stub."
  end subroutine bmi_get_grid_y

  subroutine bmi_get_grid_z(grid_id, z, ierr)
    integer,  intent(in)  :: grid_id
    real(dp), intent(out) :: z(:)
    integer,  intent(out) :: ierr
    error stop "bmi_get_grid_z: stub."
  end subroutine bmi_get_grid_z

  subroutine bmi_get_grid_node_count(grid_id, count, ierr)
    integer, intent(in)  :: grid_id
    integer, intent(out) :: count
    integer, intent(out) :: ierr
    error stop "bmi_get_grid_node_count: stub."
  end subroutine bmi_get_grid_node_count

  subroutine bmi_get_grid_edge_count(grid_id, count, ierr)
    integer, intent(in)  :: grid_id
    integer, intent(out) :: count
    integer, intent(out) :: ierr
    error stop "bmi_get_grid_edge_count: stub."
  end subroutine bmi_get_grid_edge_count

  subroutine bmi_get_grid_face_count(grid_id, count, ierr)
    integer, intent(in)  :: grid_id
    integer, intent(out) :: count
    integer, intent(out) :: ierr
    error stop "bmi_get_grid_face_count: stub."
  end subroutine bmi_get_grid_face_count

  subroutine bmi_get_grid_edge_nodes(grid_id, edge_nodes, ierr)
    integer, intent(in)  :: grid_id
    integer, intent(out) :: edge_nodes(:)
    integer, intent(out) :: ierr
    error stop "bmi_get_grid_edge_nodes: stub."
  end subroutine bmi_get_grid_edge_nodes

  subroutine bmi_get_grid_face_edges(grid_id, face_edges, ierr)
    integer, intent(in)  :: grid_id
    integer, intent(out) :: face_edges(:)
    integer, intent(out) :: ierr
    error stop "bmi_get_grid_face_edges: stub."
  end subroutine bmi_get_grid_face_edges

  subroutine bmi_get_grid_face_nodes(grid_id, face_nodes, ierr)
    integer, intent(in)  :: grid_id
    integer, intent(out) :: face_nodes(:)
    integer, intent(out) :: ierr
    error stop "bmi_get_grid_face_nodes: stub."
  end subroutine bmi_get_grid_face_nodes

  subroutine bmi_get_grid_nodes_per_face(grid_id, nodes_per_face, ierr)
    integer, intent(in)  :: grid_id
    integer, intent(out) :: nodes_per_face(:)
    integer, intent(out) :: ierr
    error stop "bmi_get_grid_nodes_per_face: stub."
  end subroutine bmi_get_grid_nodes_per_face

end module gwswex_kernel_bmi
