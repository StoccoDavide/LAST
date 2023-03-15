# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

QR := proc( A::{Matrix}, V::{symbol}, $)::{table};

  description "Compute the Givens QR decomposition of a square matrix <A> "
              "using the veiling strategy <VeilingStrategy> and the veiling symbol <V>.";

  local m, n, Q, R, DG, k, j, a, b, z1, z2, r, Rk, Rj, apply_veil;

  LEM:-VeilForget(V);

  # Check if to veil or not
  apply_veil := (z) -> `if`( LULEM:-VeilingStrategy(z), LEM:-Veil[V](z), z );

  # Extract the dimensions of the matrix A
  m, n := LinearAlgebra[Dimensions](A):

  # Check if the matrix A valid
  assert(
    m >= n,
    "LULEM::QR(...): invalid matrix A(m,n) detected, got (m >= n)."
  );

  # Initialize some variables
  Q  := [];      # Orthogonal transformation as a list of Given rotations
  R  := copy(A); # Transformed matrix so far
  DG := Vector[column](m,k->1);

  # Compute the Householder QR decomposition with veiling
  for k from 1 to m-1 do
    if LULEM:-Verbose then
      printf(
        "LULEM::QR(...): processing %d-th colum. Length %d\n",
        k, length(LEM:-VeilList(V))
      );
    end;
    for j from k+1 to n do
      a := R[k,k];
      b := R[j,k];
      if not b = 0 then
        z1 := DG[k];
        z2 := DG[j];
        if a = 0 then
          # simple case do a swap
          Q          := [ op(Q), [ k, j, 0, 0, 0, 0, 0 ] ];
          Rk         := R[k,k..-1];
          Rj         := R[j,k..-1];
          R[k,k..-1] := Rj;
          R[j,k..-1] := Rk;
          DG[k]      := z2;
          DG[j]      := z1;
        else
          r  := apply_veil(Normalizer~(z2*a^2 + z1*b^2));
          Q  := [ op(Q), [ k, j, a, b, z1, z2, r ] ];
          Rk := R[k,k+1..-1];
          Rj := R[j,k+1..-1];

          R[k,k+1..-1] := apply_veil~(Normalizer~((a*z2)*Rk+(b*z1)*Rj));
          R[j,k+1..-1] := apply_veil~(Normalizer~(a*Rj-b*Rk));

          R[k,k] := r;
          R[j,k] := 0;
          #DG[k]  := apply_veil(Normalizer~(z1*z2*r));
          DG[k]  := z1*z2*r;
          DG[j]  := r;
        end if;
      end if;
    end do;
    if R[k,k] = 0 then
      error "R[%a,%a] = 0", k, k;
    end if;
  end do;

  # Return the QR decomposition
  return table([
    "method"   = "QR",
    "D"        = DG,
    "Q"        = Q,
    "R"        = R,
    "V"        = V,
    "Q_length" = length(convert(Q,list)),
    "D_length" = length(convert(DG,list)),
    "R_length" = length(convert(R,list)),
    "V_length" = length(LEM:-VeilList(V))
  ]);
end proc: # QR

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

QRsolve := proc( T::{table}, xb::{Vector}, V::{symbol, function}, $)

  description "Solve the linear system Ax=b using QR decomposition <T>, "
              "provided the vector <b>, the veiling symbol <V> and the veiling strategy "
              "<VeilingStrategy>.";

  local Q, R, DG, Di, Dj, DD, DDD, m, n, i, j, k, z, c, s, a, b, d, x, z1, z2, r, apply_veil;

  # Check if to veil or not
  apply_veil := (z) -> `if`( LULEM:-VeilingStrategy(z), LEM:-Veil[V](z), z );

  # apply Q^T a rhs
  Q  := T["Q"];
  R  := T["R"];
  DG := T["D"];

  # Extract the dimensions of the matrix R
  m, n := LinearAlgebra[Dimensions](R):
  x    := xb;

  for i from 1 to nops(Q) do
    k  := Q[i][1];
    j  := Q[i][2];
    a  := Q[i][3];
    if a = 0 then
      b    := x[k];
      x[k] := x[j];
      x[j] := b;
    else
      b  := Q[i][4];
      z1 := Q[i][5];
      z2 := Q[i][6];
      c  := apply_veil~(Normalizer~((a*z2)*x[k]+(b*z1)*x[j]));
      s  := apply_veil~(Normalizer~(a*x[j]-b*x[k]));
      x[k] := c;
      x[j] := s;
    end if;
  end do;
  # solve R^(-1)

  x[n] := apply_veil(x[n]/R[n,n]);
  for i from n-1 to 1 by -1 do
    if LULEM:-Verbose then
      printf("LULEM:-QRsolve, backward %d\n",i);
    end if;
    s    := apply_veil( x[i] - add(R[i,i+1..n] *~ x[i+1..n]));
    x[i] := apply_veil( s / R[i,i] );
  end do;

  # Return outputs
  return x;
end proc: # QRsolve

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -