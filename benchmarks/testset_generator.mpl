# Testset generator for the supplementary material of the paper

# Reference:
# J. Carette, W. Zhou, D. J. Jeffrey, M. B. Monagan, Linear algebra using Maple’s
# LargeExpressions package, in: Proceedings of Maple Conference, 2006, pp. 14-25.

size_list := [10, 25, 50]:
dens_list := [100, 50, 25]:
test_list := ["real", "symbol", "poly1", "poly2"]:

# Check initial data
if not evalb(dens_list[1] = 100) then
  error "The first density must be 100%";
end if:

# Total number of tests
printf("Total number of tests = %d\n",
  nops(size_list)* nops(dens_list)*nops(test_list)
):

# Generate the testset for the dense matrices
printf("Generating dense matrices and vectors for...\n"):
A_list := []:
for n in size_list do

  printf("\tsize = %d, density = %d\n", n, dens_list[1]):

  # Name of the generated matrices
  A_base := convert(cat("A_", n, "_", dens_list[1]), string):
  A_list := [op(A_list),
    convert(cat(A_base, "_", test_list[1]), string),
    convert(cat(A_base, "_", test_list[2]), string),
    convert(cat(A_base, "_", test_list[3]), string),
    convert(cat(A_base, "_", test_list[4]), string)
  ];
  d_tmp := evalf(dens_list[1]/100):

  # Random matrix with real entries
  assign(convert(A_list[-4], symbol) = LinearAlgebra:-RandomMatrix(n, n, generator = -1.0e12..1.0e12, density = d_tmp)):

  # Fully symbolic matrix
  assign(convert(A_list[-3], symbol) = Matrix(n, n, symbol = a));

  # Random matrices with univariate entries of degree 10 and 3 terms
  assign(convert(A_list[-2], symbol) = LinearAlgebra:-RandomMatrix(n, n, generator = (() -> randpoly([x], degree = 10, terms = 3)), density = d_tmp)):

  # Random matrices with bivariate entries of degree 10 and 6 terms
  assign(convert(A_list[-1], symbol) = LinearAlgebra:-RandomMatrix(n, n, generator = (() -> randpoly([x,y], degree = 10, terms = 6)), density = d_tmp)):

end do:
printf("DONE\n"):

# Generate the testset for the sparse matrices
printf("Generating sparse matrices and vectors for...\n"):
for n in size_list do
  for d in dens_list[2..-1] do

    # Name of the generated matrices
    A_base := convert(cat("A_", n, "_", d), string):
    A_list := [op(A_list),
      convert(cat(A_base, "_", test_list[1]), string),
      convert(cat(A_base, "_", test_list[2]), string),
      convert(cat(A_base, "_", test_list[3]), string),
      convert(cat(A_base, "_", test_list[4]), string)
    ];
    d_tmp := evalf(d/100):

    # Retrieve the full matrices
    A_base_full := convert(cat("A_", n, "_", dens_list[1]), string):
    A_list_full := [
      convert(cat(A_base_full, "_", test_list[1]), string),
      convert(cat(A_base_full, "_", test_list[2]), string),
      convert(cat(A_base_full, "_", test_list[3]), string),
      convert(cat(A_base_full, "_", test_list[4]), string)
    ];

    # Generate density patterns
    A_dens := LinearAlgebra:-RandomMatrix(n, n, generator = 1..1, density = d_tmp):

    # Random matrix with real entries
    assign(convert(A_list[-4], symbol) = eval(convert(A_list_full[1], symbol)) *~ A_dens):

    # Fully symbolic matrix
    assign(convert(A_list[-3], symbol) = eval(convert(A_list_full[2], symbol)) *~ A_dens):

    # Random matrices with univariate entries of degree 10 and 3 terms
    assign(convert(A_list[-2], symbol) = eval(convert(A_list_full[3], symbol)) *~ A_dens):

    # Random matrices with bivariate entries of degree 10 and 6 terms
    assign(convert(A_list[-1], symbol) = eval(convert(A_list_full[4], symbol)) *~ A_dens):

    printf("\tsize = %d, density = %d\n", n, d):
  end do:
end do:
printf("DONE\n"):

# Save the testset
printf("Saving the testset... "):
parse(cat(
  "save size_list, dens_list, test_list,",
  "     A_list, ", sprintf("%q", op(parse~(A_list))), ","
  "    `testset.mpl`"), statement);
printf("DONE\n"):
