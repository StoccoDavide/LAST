# Benchmark: Compare the performance of LULEM and LAST with Maple's LU decomposition

# Timing settings
time_limit := 7200:
expr_limit := 50:
generate_maple := true:
generate_lulem := true:
generate_last  := true:

# LAST initialization and settings
LAST_obj := Object(LAST);
LAST_obj:-InitLEM(LAST_obj, "V");
LAST_obj:-DisableVerboseMode(LAST_obj);
LAST_obj:-DisableWarningMode(LAST_obj);
LAST_obj:-DisableUnveiling(LAST_obj);
LAST_obj:-EnableFastPivoting(LAST_obj);
LAST_obj:-SetTimeLimit(LAST_obj, 0.5);
LEM_obj := LAST_obj:-GetLEM(LAST_obj);
LEM_obj:-SetExprMaxCost(LEM_obj, expr_limit);

# LULEM initialization and settings
Vs := LULEM:-VeilingStrategy_LB:
Zs := LULEM:-ZeroStrategy_length:

# Open file to append the results
fname_maple := "maple.csv":
fname_lulem := "lulem.csv":
fname_last  := "last.csv":
fremove(fname_maple):
fremove(fname_lulem):
fremove(fname_last):
fstream_maple := fopen(fname_maple, APPEND):
fstream_lulem := fopen(fname_lulem, APPEND):
fstream_last  := fopen(fname_last, APPEND):

fprintf(fstream_maple, "size, density, type, factime, rank, Alen, LUlen, Ands, LUnds, Annz, LUnnz, Vnum, Vnds, Vlen\n");
fprintf(fstream_lulem, "size, density, type, factime, rank, Alen, LUlen, Ands, LUnds, Annz, LUnnz, Vnum, Vnds, Vlen\n");
fprintf(fstream_last,  "size, density, type, factime, rank, Alen, LUlen, Ands, LUnds, Annz, LUnnz, Vnum, Vnds, Vlen\n");

# Total number of tests
test_tmp := ["real", "symbol", "poly1", "poly2"];
size_tmp := [10, 25, 50];
dens_tmp := [25, 50, 100];
total_tests := nops(cmpl_tmp)*nops(test_tmp)*nops(size_tmp)*nops(dens_tmp):
printf("Total number of tests = %d\n", total_tests):

printf("Benchmarking the algorithms...\n"):
test_num := 0:
for t in test_tmp do
  for n in size_tmp do
    for d in dens_tmp do
      test_num := test_num + 1:
      printf("\ttest = %d, entry = %s, size = %d, density = %d\n", test_num, t, n, d):

      # Retrieve the matrices and vectors
      A_tmp := parse(cat("A_", n, "_", d, "_", t)):
      b_tmp := parse(cat("b_", n, "_", d, "_", t)):

      A_nds := add(LEM_obj:-ExpressionCost~(LEM_obj, A_tmp));
      A_nnz := nops(op(2, A_tmp));
      A_len := length(A_tmp);

      # Maple - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if generate_maple then
        try
          st := time():
          P, L, U, rnk := timelimit(time_limit, LinearAlgebra:-LUDecomposition(A_tmp, output=['P','L','U','rank'])):
          fac_time := time() - st:
          LU_nds := add(LEM_obj:-ExpressionCost~(LEM_obj, L)) + add(LEM_obj:-ExpressionCost~(LEM_obj, U)):
          LU_nnz := nops(op(2, L)) + nops(op(2, U)):
          LU_len := length(L) + length(U):
          V_num  := 0;
          V_nds  := 0;
          V_len  := 0;
        catch "time expired":
          fac_time := 0;
          LU_nds   := -2;
          LU_nnz   := 0:
          LU_len   := 0;
          V_num    := 0;
          V_nds    := 0;
          V_len    := 0;
        catch "division by zero":
          fac_time := 0;
          LU_nds   := -3;
          LU_nnz   := 0:
          LU_len   := 0;
          V_num    := 0;
          V_nds    := 0;
          V_len    := 0;
        catch:
          fac_time := 0;
          LU_nds   := -4;
          LU_nnz   := 0:
          LU_len   := 0;
          V_num    := 0;
          V_nds    := 0;
          V_len    := 0;
        end try:
      else
        fac_time := 0;
        LU_nds   := 0;
        LU_nnz   := 0;
        LU_len   := 0;
        V_num    := 0;
        V_nds    := 0;
        V_len    := 0;
        rnk      := 0;
      end if:

      printf("\t\tMaple: fac_time = %.3f, rnk = %d, A_len = %d, LU_len = %d, A_nds = %d, LU_nds = %d, A_nnz = %d, LU_nnz = %d, V_num = %d, V_nds = %d, V_len = %d\n",
        fac_time, rnk, A_len, LU_len, A_nds, LU_nds, A_nnz, LU_nnz, V_num, V_nds, V_len);
      fprintf(fstream_maple,
        "%d, %d, %s, %.3f, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
        n, d, t, fac_time, rnk, A_len, LU_len, A_nds, LU_nds, A_nnz, LU_nnz, V_num, V_nds, V_len
      ):

      # LULEM - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if t = "real" then
        Ps := LULEM:-PivotStrategy_numeric:
      else
        Ps := LULEM:-PivotStrategy_Llength:
      end if;

      if generate_lulem then
        try
          st := time():
          P, L, U, r, rnk := timelimit(time_limit, LULEM:-LU(A_tmp, V, Vs, Ps, Zs)):
          fac_time := time() - st:
          LU_nds := add(LEM_obj:-ExpressionCost~(LEM_obj, L)) + add(LEM_obj:-ExpressionCost~(LEM_obj, U)):
          LU_nnz := nops(op(2, L)) + nops(op(2, U)):
          LU_len := length(L) + length(U):
          V_tmp  := LULEM:-ListVeil(V):
          V_num  := nops(V_tmp):
          V_nds  := add(LEM_obj:-ExpressionCost~(LEM_obj, V_tmp)):
          V_len  := add(length~(V_tmp)):
        catch "time expired":
          fac_time := 0;
          LU_nds   := -2;
          LU_nnz   := 0:
          LU_len   := 0;
          V_num    := 0;
          V_nds    := 0;
          V_len    := 0;
        catch "division by zero":
          fac_time := 0;
          LU_nds   := -3;
          LU_nnz   := 0:
          LU_len   := 0;
          V_num    := 0;
          V_nds    := 0;
          V_len    := 0;
        catch:
          fac_time := 0;
          LU_nds   := -4;
          LU_nnz   := 0:
          LU_len   := 0;
          V_num    := 0;
          V_nds    := 0;
          V_len    := 0;
        end try:
      else
        fac_time := 0;
        LU_nds   := 0;
        LU_nnz   := 0:
        LU_len   := 0;
        V_num    := 0;
        V_nds    := 0;
        V_len    := 0;
      end if:

      printf("\t\tLULEM: fac_time = %.3f, rnk = %d, A_len = %d, LU_len = %d, A_nds = %d, LU_nds = %d, A_nnz = %d, LU_nnz = %d, V_num = %d, V_nds = %d, V_len = %d\n",
        fac_time, rnk, A_len, LU_len, A_nds, LU_nds, A_nnz, LU_nnz, V_num, V_nds, V_len);
      fprintf(fstream_lulem,
        "%d, %d, %s, %.3f, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
        n, d, t, fac_time, rnk, A_len, LU_len, A_nds, LU_nds, A_nnz, LU_nnz, V_num, V_nds, V_len
      ):

      # LAST  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if generate_last then
        try
          st := time():
          timelimit(time_limit, LAST_obj:-LU(LAST_obj, A_tmp)):
          fac_time := time() - st:
          res  := LAST_obj:-GetResults(LAST_obj):
          rnk := res["rank"];
          LU_nds := add(LEM_obj:-ExpressionCost~(LEM_obj, res["L"])) + add(LEM_obj:-ExpressionCost~(LEM_obj, res["U"])):
          LU_nnz := nops(op(2, res["L"])) + nops(op(2, res["U"])):
          LU_len := length(res["L"]) + length(res["U"]):
          V_tmp  := LEM_obj:-VeilList(LEM_obj):
          V_num  := nops(V_tmp):
          V_nds  := add(LEM_obj:-ExpressionCost~(LEM_obj, V_tmp)):
          V_len  := add(length~(V_tmp)):
        catch "time expired":
          fac_time := 0;
          LU_nds   := -2;
          LU_nnz   := 0:
          LU_len   := 0;
          V_num    := 0;
          V_nds    := 0;
          V_len    := 0;
        catch "division by zero":
          fac_time := 0;
          LU_nds   := -3;
          LU_nnz   := 0:
          LU_len   := 0;
          V_num    := 0;
          V_nds    := 0;
          V_len    := 0;
        catch:
          fac_time := 0;
          LU_nds   := -4;
          LU_nnz   := 0:
          LU_len   := 0;
          V_num    := 0;
          V_nds    := 0;
          V_len    := 0;
        end try:
      else
        fac_time := 0;
        LU_nds   := 0;
        LU_nnz   := 0:
        LU_len   := 0;
        V_num    := 0;
        V_nds    := 0;
        V_len    := 0;
      end if:

      printf("\t\tLAST: fac_time = %.3f, rnk = %d, A_len = %d, LU_len = %d, A_nds = %d, LU_nds = %d, A_nnz = %d, LU_nnz = %d, V_num = %d, V_nds = %d, V_len = %d\n",
        fac_time, rnk, A_len, LU_len, A_nds, LU_nds, A_nnz, LU_nnz, V_num, V_nds, V_len);
      fprintf(fstream_last,
        "%d, %d, %s, %.3f, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n",
        n, d, t, fac_time, rnk, A_len, LU_len, A_nds, LU_nds, A_nnz, LU_nnz, V_num, V_nds, V_len
      ):

      # Clear the results - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LULEM:-ForgetVeil(V);
      LEM_obj:-ForgetVeil(LEM_obj);
      LAST_obj:-ClearResults(LAST_obj);

    end do:
  end do:
end do:
printf("DONE\n"):

fclose(fstream_maple):
fclose(fstream_lulem):
fclose(fstream_last):
