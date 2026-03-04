import time
import signal
import sys
import pandas as pd
import sympy as sp
from sympy.core.function import count_ops
from sympy import Matrix, SparseMatrix, symbols, IndexedBase

a = IndexedBase('a')

# Utility functions
def dag_nodes(expr):
    return len(expr.atoms()) + count_ops(expr)

def matrix_nodes(M):
    return sum(dag_nodes(e) for e in M)

def matrix_length(M):
    """Total symbolic size"""
    return sum(len(str(e)) for e in M)

def matrix_nnz(M):
    """Number of nonzeros"""
    return sum(1 for e in M if e != 0)

def timeout_handler(signum, frame):
    raise TimeoutError

def load_testset(filename):
    """
    Loads a Python-formatted SymPy testset file into
    a controlled namespace and returns the dictionary.
    """

    # Controlled execution environment
    namespace = {
        "sp": sp,
        "Matrix": Matrix,
        "SparseMatrix": SparseMatrix,
        "symbols": symbols,
        "IndexedBase": IndexedBase,
        # define indexed symbols used in file
        "a": IndexedBase("a"),
        "x": symbols("x"),
        "y": symbols("y"),
    }

    with open(filename, "r") as f:
        code = f.read()

    exec(code, namespace)

    return namespace

# Load testset
data = load_testset("testset_python.txt")

# Timing settings
time_limit = 7200

fname = "sympy.csv"

# Test configuration
test_tmp = ["real", "symbol", "poly1", "poly2"];
size_tmp = [10, 25, 50];
dens_tmp = [100, 50, 25];

total_tests = len(test_tmp) * len(size_tmp) * len(dens_tmp)
print(f"Total number of tests = {total_tests}")
print("Benchmarking the algorithms...")

results = []
test_num = 0

# Main loop
for t in test_tmp:
    for n in size_tmp:
        for d in dens_tmp:

            test_num += 1
            print(f"\ttest = {test_num}, entry = {t}, size = {n}, density = {d}")

            A_name = f"A_{n}_{d}_{t}"
            A_tmp = data[A_name]

            A_nds = matrix_nodes(A_tmp)
            A_nnz = matrix_nnz(A_tmp)
            A_len = matrix_length(A_tmp)

            # Defaults
            fac_time = 0
            LU_nds   = 0
            LU_nnz   = 0
            LU_len   = 0
            V_num    = 0
            V_nds    = 0
            V_len    = 0
            rnk      = 0

            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(time_limit)

            try:
                st = time.time()

                L, U, perm = A_tmp.LUdecomposition()
                rnk = A_tmp.rank()

                fac_time = time.time() - st

                # Remove identity from L
                for i in range(min(L.rows, L.cols)):
                    L[i, i] -= 1

                LU_nds = matrix_nodes(L) + matrix_nodes(U)
                LU_nnz = matrix_nnz(L) + matrix_nnz(U)
                LU_len = matrix_length(L) + matrix_length(U)

            except TimeoutError:
                LU_nds = -2

            except ZeroDivisionError:
                LU_nds = -3

            except Exception:
                LU_nds = -4

            finally:
                signal.alarm(0)

            print(
                f"\t\tSymPy: fac_time = {fac_time:.3f}, "
                f"rnk = {rnk}, A_len = {A_len}, LU_len = {LU_len}, "
                f"A_nds = {A_nds}, LU_nds = {LU_nds}, "
                f"A_nnz = {A_nnz}, LU_nnz = {LU_nnz}, "
            )

            results.append([
                n, d, t, fac_time, rnk,
                A_len, LU_len,
                A_nds, LU_nds,
                A_nnz, LU_nnz,
                V_num, V_nds, V_len
            ])

print("DONE")