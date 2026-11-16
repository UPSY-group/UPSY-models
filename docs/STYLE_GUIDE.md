# Fortran style guide

This is a short summary of how code in UPSY-models is written. It isn't a rulebook to be enforced line by line; it describes the conventions that the existing code already follows, so that new code looks like it belongs. When in doubt, find a recent module that does something similar to what you're building and copy its layout.

Some older parts of the codebase don't follow all of these conventions yet (you'll still find `IF`, `DO` and `CALL` in capitals here and there). Don't mix in the old style when writing new code, but also don't reformat old code in the same PR as a functional change.

## Formatting

- **Lowercase everything**: keywords, intrinsics, variable names. `if`, `do`, `call`, `end do`, `elseif`.
- **Indent with 2 spaces**, never tabs.
- Write **`end do`, `end if`, `end subroutine name`** in full, and always repeat the procedure name after `end subroutine` / `end function`.
- Use `elseif` (one word).
- **Space after an opening bracket** in calls and array indexing, but not before the closing one: `call init_routine( routine_name)`, `self%H( vi)`, `size( poly,1)`, `max( 0._dp, x)`.
- **Align related lines** vertically: declarations (`::`), `intent`s, assignments to a group of similar variables, keyword arguments.
- Keep lines reasonably short (most are well under 100 characters). Break long calls with `&` and indent the continuation by 2 spaces.
- Leave a blank line after the declarations, around loops and `if` blocks that do something non-trivial, and before `end subroutine`.

```fortran
LI_xdy   = LI_xdy   + line_integral_xdy(   p, q, tol_dist)
LI_mxydx = LI_mxydx + line_integral_mxydx( p, q, tol_dist)
LI_xydy  = LI_xydy  + line_integral_xydy(  p, q, tol_dist)
```

## Modules

- One module per file, and the file has the same name as the module.
- Always put `implicit none` and `private` in the module, and list what's `public` explicitly.
- Always import with `use ..., only: ...`. That way, anyone reading the module can see where each name comes from.
- Use submodules (`submodule( parent) parent_submod_whatever`) to split up large modules, e.g. a type and all its type-bound procedures.

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
- Every subroutine that does real work calls **`init_routine` / `finalise_routine`**, with a `routine_name` parameter that matches the subroutine name exactly. This drives the call stack, the timing output and the error messages, so don't skip it. Small `pure` helper functions are exempt.
- Mark functions `pure` where possible, and use `result( name)` for function results.
- Prefer subroutines with several `intent(out)` arguments over functions with side effects.

## Naming

- Names are `lower_snake_case` and **descriptive, even when long**: `calc_line_integrals_around_all_overlaps` beats `calc_li`. Procedure names start with a verb (`calc_`, `initialise_`, `allocate_`, `remap_`, `run_`, ...).
- Physical quantities keep their usual symbols, including capitals: `H`, `Hs`, `T2m`, `u_3D`, `dHdt`, `beta_sq`.
- **Mesh loop indices** follow a fixed convention: `vi` for vertices, `ti` for triangles, `ei` for edges, `k` for vertical layers, `m` for months. Process-local ranges are `vi1:vi2`, `ti1:ti2`, `ei1:ei2`.
- Derived types are called `type_something`; abstract types are called `atype_something`.
- The passed object in type-bound procedures is always called `self`.
- Use British spelling, as in the rest of the codebase: `initialise`, `finalise`, `neighbour`.

## Numbers and precision

- All real numbers are `real(dp)`, with `dp` imported from `precisions`.
- Write **every real literal with a kind suffix**: `0._dp`, `1e-3_dp`, `2._dp * pi`. Without the suffix, a literal is single precision and quietly loses digits.
- Convert integers explicitly: `real( m,dp)`.
- Take constants like `pi` and `NaN` from the `parameters` module instead of defining your own.

## Object-oriented models

New model components (climate, SMB, ocean, momentum balance, ...) are built as a hierarchy of types. Look at [`src/UPSY/models/demo_model/`](../src/UPSY/models/demo_model/) for a complete, minimal example.

- An `atype_<thing>_model_data` holds the fields that other models are allowed to access.
- An abstract `atype_<thing>_model` extends it and implements the shared `allocate`, `deallocate`, `initialise`, `run` and `remap` routines. Each one first does what is common to all models, then what is common to this kind of model, and then calls a `deferred` procedure for what is specific to the implementation.
- Each concrete implementation extends the abstract type and provides the deferred procedures, whose interfaces are declared in an `abstract interface` block.
- Model fields are created with `self%create_field( ...)`, always with `name`, `long_name`, `units` and `remap_method` as keyword arguments.
- Call `self%checksum()` at the end of `initialise`, `run` and `remap`, so that the test suite picks up any change in model output.
- A new variant of an existing solver or model gets its own type, rather than an `if` branch inside an existing one.

## Errors and safety checks

- Stop on unrecoverable problems with `call crash('message')` and print non-fatal issues with `call warning('message')`. Never use `stop`.
- Every `select type` and every `if`/`elseif` chain over options gets a final `class default` / `else` branch that crashes with a clear message:

  ```fortran
  select case (choice_remapping_method)
  case ('nearest_neighbour')
    ...
  case default
    call crash('unknown choice_remapping_method "' // trim( choice_remapping_method) // '"')
  end select
  ```

- Fortran won't always crash when something goes wrong (out-of-bounds writes into shared memory just silently corrupt data). Where that's a risk, add an explicit check, commented with `! Safety:`.
- Loops that could get stuck (e.g. flood-fill searches) get a counter that crashes if it grows too large.

## Parallelisation

- Loop over the process-local range only (`do vi = mesh%vi1, mesh%vi2`). Only fall back to looping over the full domain when you really need all of it, and add a comment saying why.
- Shared-memory arrays need explicit synchronisation (`sync`, `sync_node`) before other processes read what you've written. When you add one, add a comment saying what it protects against.

## Comments

- Comments explain **why**, or summarise what a block does; they don't repeat the code. A short comment above each logical block of a routine (`! Calculate total area of overlap for this dst face`) makes long routines much easier to follow.
- Use underlined headers to split up long routines or modules:

  ```fortran
  ! Simple intersecting squares
  ! ===========================
  ```

- Use `! NOTE:` for things that would surprise the reader, such as a workaround or a known limitation.
- In type definitions, document the type itself with `!<` comments.

## New features

- New functionality is **off by default** and switched on through a config parameter (`C%choice_...` or `C%do_...`).
- Add a unit test for anything you can test in isolation (e.g. a geometry routine), following the existing `test_...` routines (see e.g. `src/UPSY/validation/unit_tests/`).
- Don't leave commented-out code, debug `print` statements or stray semicolons behind when you open a PR.
