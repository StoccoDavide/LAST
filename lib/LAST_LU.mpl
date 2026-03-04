# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   _    _   _
#  | |  | | | |
#  | |  | | | |
#  | |__| |_| |
#  |_____\___/
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LU::static := proc(
  _self::LAST,
  A::Matrix,
  $)

  description "Compute the LU decomposition of a square matrix <A> and check "
    "if the veiling symbol is already present in the matrix coefficients.";

  local V, M, L, U, pivot, pivot_list, m, n, mn, i, k, rnk, r, c, P, Q, tmp;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Get the veiling label
  V := _self:-m_LEM:-GetVeilingLabel(_self:-m_LEM);

  # Get matrix dimensions
  m, n := LinearAlgebra:-Dimensions(A):
  if _self:-m_VerboseMode then
    printf("LAST:-LU(...): %d x %d matrix detected.\n", m, n);
  end if;

  # Create pivot vector and matrix M
  M := copy(A);
  r := Vector(m, k -> k);
  c := Vector(n, k -> k);
  mn := min(m, n);
  rnk := mn;
  pivot_list := [];

  # Perform Gaussian elimination
  i := 1;
  for k from i to mn do
    if _self:-m_VerboseMode then
      printf(
        "LAST:-LU(...): processing %d-th row, veilings = %d.\n",
        k,  nops(_self:-m_LEM:-VeilList(_self:-m_LEM))
      );
    end if;

    pivot := _self:-Pivoting(_self, k, M, r, c);

    # Check if the pivot is zero
    if pivot["is_zero"] then
      rnk := k - 1;
      if _self:-m_VerboseMode then
        WARNING("LAST:-LU(...): the matrix appears to be not full rank.");
      end if;
      break;
    else
      pivot_list := [op(pivot_list), pivot["value"]];
    end if;

    if _self:-m_VerboseMode then
      printf(
        "LAST:-LU(...): M[%d,%d] = %a, cost = %d, degree_r = %d, degree_c = %d.\n",
        pivot["i"], pivot["j"], pivot["value"], pivot["cost"], pivot["degree_r"],
        pivot["degree_c"]
      );
      printf("LAST:-LU(...): performing Gaussian elimination...");
    end if;

    # Gaussian elimination
    tmp         := [k+1..-1];
    M[k,   k]   := _self:-m_LEM:-Veil(_self:-m_LEM,  pivot["value"]);
    M[tmp, k]   := _self:-m_LEM:-Veil~(_self:-m_LEM, M[tmp, k]/pivot["value"]);
    M[k,   tmp] := _self:-m_LEM:-Veil~(_self:-m_LEM, M[k, tmp]);
    M[tmp, tmp] := _self:-m_LEM:-Veil~(_self:-m_LEM, M[tmp, tmp]-M[tmp, k].M[k, tmp]);

    if _self:-m_VerboseMode then
      printf(" DONE\n");
    end if;
  end do;

  L := Matrix(M[1..m, 1..m], shape = triangular[lower, unit]);
  U := Matrix(M, shape = triangular[upper]);

  # Store the LU decomposition
  _self:-m_Results := table([
    "method" = "LU",
    "L"      = L,
    "U"      = U,
    "V"      = V,
    "r"      = r,
    "c"      = c,
    "rank"   = rnk,
    "pivots" = pivot_list,
    "L_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, L),
    "U_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, U),
    "V_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, _self:-m_LEM:-VeilList(_self:-m_LEM)),
    "A_cost" = _self:-m_LEM:-ExpressionCost(_self:-m_LEM, A),
    "L_nnz"  = nops(op(2, L)),
    "U_nnz"  = nops(op(2, U)),
    "A_nnz"  = nops(op(2, A))
  ]);
  return NULL;
end proc: # LU

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LUsolve::static := proc(
  _self::LAST,
  b::Vector,
  $)::Vector;

  description "Solve the linear system Ax=b using LU decomposition provided "
  "the vector <b>.";

  local L, U, r, c, m, n, p, q, x, y, i, j, s, rnk;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Check if the results are available
  _self:-CheckResults(_self);

  # Check if the LU decomposition is available
  if not (_self:-m_Results["method"] = "LU") then
    error("wrong or not available LU decomposition (use 'LAST:-LU()' first).");
  end if;

  # Extract the LU decomposition
  L   := _self:-m_Results["L"];
  U   := _self:-m_Results["U"];
  r   := _self:-m_Results["r"];
  c   := _self:-m_Results["c"];
  rnk := _self:-m_Results["rank"];

  # Get linear system dimension
  m, n := LinearAlgebra:-Dimensions(L);
  p, q := LinearAlgebra:-Dimensions(U);

  # Check if the linear system is consistent
  # sanity check
  if not ((m = n) and (p = q)) then
    error("only square system can be solved, got L = %d x %d, and U = %d x %d.",
      m, n, p, q);
  end if;

  # Check if the linear system is consistent
  if not (n = rnk) then
    error("only full rank linear system can be solved (got rank = %1, expected "
      "rank = %2).", rnk, n);
  end if;

  # Apply permutation P
  x := Vector(b[convert(r, list)], datatype = anything);

  # Perform forward substitution to solve Ly=b[r]
  x[1] := _self:-m_LEM:-Veil(_self:-m_LEM, x[1]);
  for i from 2 to m do
    if _self:-m_VerboseMode then
      printf("LAST:-LUsolve(...): forward substitution of %d-th row, veilings = %d.\n",
        i, nops(_self:-m_LEM:-VeilList(_self:-m_LEM)));
    end if;
    x[i] := _self:-m_LEM:-Veil(_self:-m_LEM, x[i] - add(_self:-m_LEM:-Veil(_self:-m_LEM, L[i, j]*x[j]), j=1..i-1));
  end do;

  # Perform backward substitution to solve Ux[c]=y
  if _self:-m_VerboseMode then
    printf("LAST:-LUsolve(...): division by U[%d,%d].\n", n, n);
  end if;
  x[n] := _self:-m_LEM:-Veil(_self:-m_LEM, x[n]/U[n, n]);
  for i from n-1 to 1 by -1 do
    if _self:-m_VerboseMode then
      printf("LAST:-LUsolve(...): backward substitution of %d-th column, veilings = %d.\n",
        i, nops(_self:-m_LEM:-VeilList(_self:-m_LEM)));
    end if;
    s := _self:-m_LEM:-Veil(_self:-m_LEM, x[i] - add(_self:-m_LEM:-Veil(_self:-m_LEM, U[i, j]*x[j]), j=i+1..n));
    x[i] := _self:-m_LEM:-Veil(_self:-m_LEM, s/U[i, i]);
  end do;

  # Apply inverse permutation Q
  y := Vector[column](n);
  for i from 1 to n do
    y[c[i]] := x[i];
  end do;

  # Return outputs
  return y;
end proc: # LUsolve

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LUapplyLP::static := proc(
  _self::LAST,
  b::Vector,
  $)::Vector;

  description "Apply L^(-1)*P to the vector <b>.";

  local L, r, x, m, i, j;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Check if the results are available
  _self:-CheckResults(_self);

  # Check if the LU decomposition is available
  if not (_self:-m_Results["method"] = "LU") then
    error("wrong or not available LU decomposition (use 'LAST:-LU()' first).");
  end if;

  # Extract the LU decomposition
  L := _self:-m_Results["L"];
  r := _self:-m_Results["r"];
  m := LinearAlgebra:-RowDimension(L);

  # Apply permutation P
  x := Vector(b[convert(r, list)], datatype = anything);

  # Perform forward substitution to solve Ly=b[r]
  x[1] := _self:-m_LEM:-Veil(_self:-m_LEM, x[1]);
  for i from 2 to m do
    if _self:-m_VerboseMode then
      printf("LAST:-LUapplyLP(...): forward substitution of %d-th row, veilings = %d.\n",
        i, nops(_self:-m_LEM:-VeilList(_self:-m_LEM)));
    end if;
    x[i] := _self:-m_LEM:-Veil(_self:-m_LEM, x[i] - add(L[i, j]*x[j], j=1..i-1));
  end do;

  # Return outputs
  return x;
end proc: # LUapplyLP

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LUgetUQT::static := proc(
  _self::LAST,
  $)::Matrix;

  description "Return the matrix U*Q^T.";

  local U, UT, c, n, i, rnk;

  # Check if the LEM object is initialized
  _self:-CheckInit(_self);

  # Check if the results are available
  _self:-CheckResults(_self);

  # Check if the LU decomposition is available
  if not (_self:-m_Results["method"] = "LU") then
    error("wrong or not available LU decomposition (use 'LAST:-LU()' first).");
  end if;

  # Extract the LU decomposition
  U   := _self:-m_Results["U"];
  c   := _self:-m_Results["c"];
  rnk := _self:-m_Results["rank"];
  n   := LinearAlgebra:-ColumnDimension(U);
  UT  := Matrix(rnk, n);

  # Apply inverse permutation Q
  for i from 1 to n do
    UT[1..rnk, c[i]] := U[1..rnk, i];
  end do;

  # Return outputs
  return UT;
end proc: # LUgetUQT

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
