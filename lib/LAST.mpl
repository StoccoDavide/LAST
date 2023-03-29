# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                           _        _    ____ _____                          #
#                          | |      / \  / ___|_   _|                         #
#                          | |     / _ \ \___ \ | |                           #
#                          | |___ / ___ \ ___) || |                           #
#                          |_____/_/   \_\____/ |_|                           #
#                       Linear Algebra Symbolic Toolbox                       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#  Davide Stocco     (University of Trento)
#  Matteo Larcher    (University of Trento)
#  Enrico Bertolazzi (University of Trento)
#
# Inspired by the work of:
#   Wenqin Zhou       (University of Western Ontario) - Former affiliation
#   David J. Jeffrey  (University of Western Ontario)
#   Jacques Carette   (McMaster University)
#   Robert M. Corless (University of Western Ontario)
#
# License: BSD 3-Clause License
#
# This is a module for the 'LAST' (Linear Algebra Symbolic Toolbox) package. It
# contains the functions to solve linear systems of equations with large
# symbolic expressions. The module uses symbolic full pivoting LU and QR
# decompositions to solve linear systems. The 'LEM' (Large Expressions
# Management) package is used to avoid expression swell.
#
# The following code is hopefully an improved version of the original version
# provided inthe following PhD thesis:
#
#   Wenqin Zhou, Symbolic Computation Techniques for Solveing Large Expressions
#   Problems from Mathematics and Engineering (2007), Faculty of Graduate Studies,
#   The University of Western Ontario London, Ontario, Canada.
#
# We would like to thank Jacques Carette for providing the original code that we
# have used to develop this module.

unprotect('LAST');
module LAST()

  description "Linear Algebra Symbolic Toolbox module.";

  option object;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local m_Verbose           := false;
  local m_TimeLimit         := 1;
  local m_MinDegreeStrategy := "product_1";
  local m_LEM               := NULL;
  local m_Results           := NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info::static := proc()

    description "Print 'LAST' module information.";

    printf(
      "+--------------------------------------------------------------------------+\n"
      "| 'LAST' module version 1.0 - BSD 3-Clause License - Copyright (c) 2023    |\n"
      "| Current version authors:                                                 |\n"
      "|   D. Stocco, M. Larcher and E. Bertolazzi.                               |\n"
      "| Inspired by the work of:                                                 |\n"
      "|   W. Zhou, D. J. Jeffrey, J. Carette and R. M. Corless.                  |\n"
      "+--------------------------------------------------------------------------+\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad::static := proc()

    description "'LAST' module load procedure.";

    local i, lib_base_path;

    lib_base_path := NULL;
    for i in [libname] do
      if (StringTools:-Search("LAST", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = NULL) then
      error "Cannot find 'LAST' module";
    end if;
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload::static := proc()
    description "'LAST' module unload procedure.";
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Code inspired by:
  # https://www.mapleprimes.com/questions/235996-Is-There-Any-Command-Or-Function-For

  export GetDegrees::static := proc(
    _self::LAST,
    A::Matrix,
    $)::Matrix(nonnegint), Matrix(nonnegint);

    description "Get the degree matrices of the matrix <A>.";

    local i, j, k, m, n, r, c, ro, co;

    m, n := LinearAlgebra:-Dimensions(A);
    ro := Vector[column](m, k -> 1);
    co := Vector[row](n, k -> 1);
    r  := Vector[column](
      [seq(rtable_scanblock(A, [i,..], ':-NonZeros'), i = 1..m)]
    );
    c  := Vector[row](
      [seq(rtable_scanblock(A, [..,j], ':-NonZeros'), j = 1..n)]
    );
    return r.co, ro.c;
  end proc: # GetDegrees

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export InitLEM::static := proc(
    _self::LAST,
    label::{symbol, string} := NULL,
    $)

    description "Initialize the 'LEM' object with veiling label <label>.";

    _self:-m_LEM := Object(LEM);
    if type(label, symbol) then
      _self:-m_LEM:-SetVeilingLabel(_self:-m_LEM, label);
    elif type(label, string) then
      _self:-m_LEM:-SetVeilingLabel(_self:-m_LEM, parse(label));
    end if;
    return NULL;
  end proc: # InitLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ClearLEM::static := proc(
    _self::LAST,
    $)

    description "Clear the 'LEM' object.";

    _self:-m_LEM := NULL;
  end proc: # ClearLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetLEM::static := proc(
    _self::LAST,
    obj::LEM,
    $)

    description "Set the 'LEM' object <obj>.";

    _self:-m_LEM := obj;
    return NULL;
  end proc: # SetLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetLEM::static := proc(
    _self::LAST,
    $)::LEM;

    description "Get the 'LEM' object.";

    return _self:-m_LEM;
  end proc: # GetLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Spy::static := proc(
    _self::LAST,
    A::Matrix,
    $)::anything;

    description "Plot of non-zero values of the matrix <A>.";

    return plots:-sparsematrixplot(A, 'matrixview');
  end proc: # Spy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SpyLU::static := proc(
    _self::LAST,
    A::Matrix,
    L::Matrix,
    U::Matrix,
    $)::anything;

    description "Plot of non-zero values of the matrices <A>, <L> and <U> "
      "with fill-in values.";

    local mat_0, mat_1, mat_2;
    mat_0 := map(x -> `if`(x = 0, 0, 1), A);
    mat_1 := map(x -> `if`(x = 0, 0, 1), L + U);
    mat_2 := map(x -> `if`(x = 1, 1, 0), mat_0 + mat_1);

    return [plots:-sparsematrixplot(mat_0, 'matrixview', color = "Black"),
            plots:-sparsematrixplot(mat_2, 'matrixview', color = "Red")];
  end proc: # SpyLU

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export PermutationMatrices::static := proc(
    _self::LAST,
    r::Vector(nonnegint),
    c::Vector(nonnegint),
    $)::Matrix(nonnegint), Matrix(nonnegint);

    description "Compute the LU decomposition premutation matrices provided "
      "the rows pivot vector <r> and the columns pivot vector <c>.";

    local m, n, i, P, Q;

    m := LinearAlgebra:-RowDimension(r);
    n := LinearAlgebra:-RowDimension(c);
    P := Matrix(m, m);
    Q := Matrix(n, n);
    for i from 1 to m by 1 do
      P[i, r[i]] := 1;
    end do;
    for i from 1 to n by 1 do
      Q[c[i], i] := 1;
    end do;
    return P, Q;
  end proc: # Permutatiomnatrices

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVerbosity::static := proc(
    _self::LAST,
    x::boolean,
    $)

    description "Set the verbosity of the package to <x>.";

    _self:-m_Verbose := x;
    return NULL;
  end proc: # SetVerbosity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableVerbosity::static := proc(
    _self::LAST,
    $)

    description "Enable the verbosity of the package.";

    _self:-m_Verbose := true;
    return NULL;
  end proc: # EnableVerbosity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableVerbosity::static := proc(
    _self::LAST,
    $)

    description "Disable the verbosity of the package.";

    _self:-m_Verbose := false;
    return NULL;
  end proc: # DisableVerbosity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetTimeLimit::static := proc(
    _self::LAST,
    x::numeric,
    $)

    description "Set the time limit of the package to <x>.";

    if (x < 0) then
      error "time limit must be a non-negative number.";
    end if;

    _self:-m_TimeLimit := x;
    return NULL;
  end proc: # SetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SolveLinearSystem::static := proc(
    _self::LAST,
    b::Vector,
    $)

    description "Solve the factorized linear system (LU)*x=b or (QR)*x=b.";

    if (_self:-m_Results["method"] = "LU") then
      return _self:-LUsolve(_self, b);
    elif (_self:-m_Results["method"] = "QR") then
      return _self:-QRsolve(_self, b);
    else
      error "wrong or not available decomposition (use LAST::LU() or LAST::QR() "
        "first).";
    end if;

    # Work in progress: FFLU and QR2 methods.
    # elif (T["method"] = "FFLU") then
    #   return _self:-FFLUsolve(_self, b);
    # elif (T["method"] = "QR2") then
    #   return _self:-QR2solve(_self, b);
  end proc: # SolveLinearSystem

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ClearResults::static := proc( $ )

    description "Clear the results of the last factorization.";

    _self:-m_Results := table([]);
    return NULL;
  end proc: # ClearResults

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetResults::static := proc(
    _self::LAST,
    field::string := "all",
    $)

    description "Get the results of the last factorization. If <field> is "
      "specified, only the field Results['field'] is returned.";

    if (field = "all") then
      return _self:-m_Results;
    else
      return _self:-m_Results[field];
    end if;
  end proc: # GetResults

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$include "./lib/LAST_Pivoting.mpl"
$include "./lib/LAST_LU.mpl"
$include "./lib/LAST_QR.mpl"

# Work in progress: FFLU and QR2 methods.
# $include "./lib/LAST_FFLU.mpl"
# $include "./lib/LAST_QR2.mpl"

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # LAST

# That's all folks!
