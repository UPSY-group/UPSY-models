# Fortran style guide

This is a short summary of how code in UPSY-models is written. It isn't a rulebook to be enforced line by line; it describes the conventions that the existing code already follows, so that new code looks like it belongs. When in doubt, find a recent module that does something similar to what you're building and copy its layout.

Some older parts of the codebase don't follow all of these conventions yet (you'll still find `IF`, `DO` and `CALL` in capitals here and there). Don't mix in the old style when writing new code, but also don't reformat old code in the same PR as a functional change (as this makes reading Git diffs during code review unnecessarily tedious).

## Formatting

- **Lowercase everything**: keywords, intrinsics, variable names. `if`, `do`, `call`, `end do`, `elseif`. Exceptions are variables that are uppercase in mathematical notation (e.g. `H` for ice thickness) and procedure names based on proper nouns (e.g. `apply_Lloyds_algorithm`)
- **Indent with 2 spaces**, never tabs (pro tip: [set this as your default in VS Code](https://code.visualstudio.com/docs/editing/codebasics))
- Write **`end do`, `end if`, `end subroutine name`** in full, and always repeat the procedure name after `end subroutine` / `end function`.
- Use `elseif` (one word).
- **Space after an opening bracket or comma** in calls and array indexing, but not before the closing one: `call init_routine( routine_name)`, `self%H( vi)`, `size( poly,1)`, `max( 0._dp, x)`.
- **Align related lines** vertically: declarations (`::`), `intent`s, opening brackets of procedure calls, assignments to a group of similar variables, keyword arguments.
- Keep lines reasonably short (most are well under 100 characters). Break long calls with `&` and indent the continuation by 2 spaces.
- Leave a blank line after the declarations, around loops and `if` blocks that do something non-trivial, and before `end subroutine`.

```fortran
LI_xdy   = LI_xdy   + line_integral_xdy  ( p, q, tol_dist)
LI_mxydx = LI_mxydx + line_integral_mxydx( p, q, tol_dist)
LI_xydy  = LI_xydy  + line_integral_xydy ( p, q, tol_dist)
```

## Modules

- One module per file, and the file has the same name as the module.
- Always put `implicit none` and `private` in the module, and list what's `public` explicitly.
- Always import with `use ..., only: ...`. That way, anyone reading the module can see where each name comes from (including, importantly, `fortls` "go to definition").
- Use submodules (`submodule( parent) parent_submod_whatever`) to split up modules defining large derived types, specifically its type-bound procedures.

```fortran
module line_integrals

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash

  implicit none

  private

  public :: line_integral_xdy, calc_line_integrals_around_polygon

contains
  ...
end module line_integrals
```

## Procedures

Every non-trivial subroutine follows the same layout:

```fortran
  subroutine calc_something( mesh, d_in, d_out)
    ! One- or two-line description of what this does

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: d_in
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: d_out

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_something'
    integer                     :: vi

    ! Add routine to call stack
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      d_out( vi) = 2._dp * d_in( vi)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calc_something
```

- Split the declarations into `! In/output variables:` and `! Local variables:`.
- **Always give an `intent`**, written padded so they line up: `intent(in   )`, `intent(inout)`, `intent(  out)`.
- Every subroutine calls **`init_routine` / `finalise_routine`**, with a `routine_name` parameter that matches the subroutine name exactly. This drives the call stack, the timing output and the error messages, so don't skip it. Excepted are functions/subroutines that are called on a per-grid-cell basis, as updating the call stack there would significantly slow down the code.
- Mark functions `pure` where possible, and use `result( ...)` for function results.
- Prefer subroutines with several `intent(out)` arguments over functions with side effects.

## Naming

- Names are `lower_snake_case` and **descriptive, even when long**: `calc_line_integrals_around_all_overlaps` beats `calc_li`. Procedure names start with a verb (`calc_`, `initialise_`, `allocate_`, `remap_`, `run_`, ...).
- Physical quantities keep their usual symbols, including capitals: `H`, `Hs`, `T2m`, `u_3D`, `dHdt`, `beta_sq`.
- Whenever quantities are defined on the b-grid (triangles) or c-grid (edges), this is reflected in the variable name: `u_vav_b`, `u_par_c`. Quantities that are defined on the a-grid only have the `_a` extension if mapped copies of that quantity on the b- or c-grid also exist.
- **Mesh loop indices** follow a fixed convention: `vi` for vertices, `ti` for triangles, `ei` for edges, `k` for vertical layers, `m` for months. Process-local ranges are `mesh%vi1:mesh%vi2`, `mesh%ti1:mesh%ti2`, `mesh%ei1:mesh%ei2`.
- Derived types are called `type_something`; abstract types are called `atype_something`.
- The passed object in type-bound procedures is always called `self`.
- Use British spelling, as in the rest of the codebase: `initialise`, `finalise`, `neighbour`.

## Numbers and precision

- All real numbers are `real(dp)`, with `dp` imported from `precisions`.
- Write **every real literal with a kind suffix**: `0._dp`, `1e-3_dp`, `2._dp * pi`. Without the suffix, a literal is single precision and quietly loses digits.
- Convert integers explicitly: `real( m,dp)`.
- Take constants like `pi` and `NaN` from the `parameters` module instead of defining your own.

## Object-oriented models

New model components (climate, SMB, ocean, momentum balance, ...) are built as a hierarchy of types; see the [Wiki](https://github.com/UPSY-group/UPSY-models/wiki/Abstract-models-and-fields) for an in-depth explanation.

- An `atype_<thing>_model_data` holds the fields that are shared between different implementations of that model, and which other models are allowed to access.
- An abstract `atype_<thing>_model` extends it and implements the shared `allocate`, `deallocate`, `initialise`, `run` and `remap` routines. Each one first does what is common to all models, then what is common to this kind of model, and then calls a `deferred` procedure for what is specific to the implementation.
- Each concrete implementation extends the abstract type and provides the deferred procedures, whose interfaces are declared in an `abstract interface` block.
- Model fields are created with `self%create_field( ...)`, always with `name`, `long_name`, `units` and `remap_method` as keyword arguments.
- Call `self%checksum()` at the end of `initialise`, `run` and `remap`, so that the test suite picks up any change in model output.
- A new variant of an existing solver or model gets its own type, rather than an `if` branch inside an existing one.

## Errors and safety checks

- Stop on unrecoverable problems with `call crash('message')` and print non-fatal issues with `call warning('message')`. Never use `stop`. If needed, use the optional arguments `int_01`, `int_02`, ..., `dp_01`, `dp_02`, ..., to print the values of scalars to the warning/error message (with `{int_01}`, `{dp_02}`, etc. to marker where in the error string the values should be printed)
- Every `select` and every `if`/`elseif` chain over options gets a final `case default` / `else` branch that crashes with a clear message. For `select` constructs, start with the default case.

  ```fortran
  select case (choice_remapping_method)
  case default
    call crash('unknown choice_remapping_method "' // trim( choice_remapping_method) // '"')
  case ('nearest_neighbour')
    ...
  end select
  ```

- Fortran won't always crash when something goes wrong (out-of-bounds writes into shared memory just silently corrupt data). Where that's a risk, add an explicit check, commented with `! Safety:`.
- Loops that could get stuck (e.g. flood-fill searches) get a counter that crashes if it grows too large.

## Parallelisation

- Use the `allocate_dist_shared`, `deallocate_dist_shared`, `gather_dist_shared_to_all`, etc. procedures from `src/UPSY/basic/mpi_parallelisation/mpi_distributed_shared_memory.f90` to manipulate shared memory.
- Shared-memory variables require a Fortran intrinsic pointer array to hold the data, as well as an `MPI_win` object. The MPI window has the same name as the variable, prefixed with `w`.
- The pointer arrays must be declared `contiguous`, and must be set to `null()` before shared memory can be allocated; preferably do this in the type definition.

  ```fortran
  type :: type_ice_model
    real(dp), dimension(:,:), contiguous, pointer :: Hi => null()
    real(dp), dimension(:,:), contiguous, pointer :: Hb => null()
    real(dp), dimension(:,:), contiguous, pointer :: Hs => null()
    type(MPI_win) :: wHi, wHb, wHs
  end type type_ice_model
  ```

- Loop over the range of vertices/triangles/edges that is owned by the MPI process (`do vi = mesh%vi1, mesh%vi2`). Be aware that processes can read from/write to each other's ownership ranges in shared-memory variables; avoiding this is your responsibility!
- Shared-memory arrays need explicit synchronisation (`sync`, `sync_node`) before other processes read what you've written. Whenever possible, use existing procedures from `src/UPSY/basic/mpi_parallelisation`, as these already have the required safeguards in place.

## Comments

- Comments explain **why**, or summarise what a block does; they don't repeat the code. A short comment above each logical block of a routine (`! Calculate total area of overlap for this dst face`) makes long procedures much easier to follow.
- The now outdated practice of using underlined headers to split up long procedures or modules implies that these should be split into smaller procedures/submodules

  ```fortran
  ! If you feel the urge to add a separator like this, split up your code!
  ! ======================================================================
  ```

- Use `! NOTE:` for things that would surprise the reader, such as a workaround or a known limitation.

## New features

- New functionality is **off by default** and switched on through a config parameter (`C%choice_...` or `C%do_...`).
- Add a unit test for anything you can test in isolation (e.g. a geometry routine), following the existing `test_...` routines (see e.g. `src/UPSY/validation/unit_tests/`).
- Don't leave commented-out code, debug `print` statements or stray semicolons behind when you open a PR.
