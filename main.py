import math

def norma_vector(x, n):
    max_val = 0
    for i in range(1, n + 1):
        if math.fabs(x[i]) > max_val:
            max_val = math.fabs(x[i])
    return max_val

def LU_find(L, U, A, n):
    for k in range(1, n + 1):
        for i in range(k, n + 1):
            sum_val = 0.0
            for j in range(1, k):
                sum_val += L[i][j] * U[j][k]
            L[i][k] = A[i][k] - sum_val
        for i in range(k + 1, n + 1):
            sum_val = 0.0
            for j in range(1, k):
                sum_val += L[k][j] * U[j][i]
            U[k][i] = (A[k][i] - sum_val) / L[k][k]

def LU_solve(L, U, B, X, n):
    Z = [0.0] * (n + 1)
    Z[1] = B[1] / L[1][1]
    for k in range(2, n + 1):
        sum_val = 0.0
        for j in range(1, k):
            sum_val += L[k][j] * Z[j]
        Z[k] = (B[k] - sum_val) / L[k][k]
    
    X[n] = Z[n]
    for k in range(n - 1, 0, -1):
        sum_val = 0.0
        for j in range(k + 1, n + 1):
            sum_val += U[k][j] * X[j]
        X[k] = Z[k] - sum_val

def new_B(A, X, B0, n):
    for i in range(1, n + 1):
        sum_val = 0.0
        for j in range(1, n + 1):
            sum_val += A[i][j] * X[j]
        B0[i] = sum_val

def main():
    with open("matrix_A.txt", "r") as matrix_A, \
         open("matrix_B.txt", "w") as matrix_B, \
         open("matrix_X.txt", "w") as matrix_X, \
         open("matrix_L.txt", "w") as matrix_L, \
         open("matrix_U.txt", "w") as matrix_U:

        n = 100
        N = n + 1
        x_0 = 0.51
        eps = 1e-12
        kmax = int(1e+6)

        B = [0.0] * N
        B0 = [0.0] * N
        X = [0.0] * N
        R = [0.0] * N
        dX = [0.0] * N
        A = [[0.0] * N for _ in range(N)]
        L = [[0.0] * N for _ in range(N)]
        U = [[0.0] * N for _ in range(N)]

        # Initialize U
        for i in range(1, n + 1):
            U[i][i] = 1.0

        # Read matrix A from file
        for i in range(1, n + 1):
            line = matrix_A.readline().split()
            for j in range(1, n + 1):
                A[i][j] = float(line[j - 1])

        # Compute B[i] and write to matrix_B.txt
        for i in range(1, n + 1):
            B[i] = sum(A[i][j] for j in range(1, n + 1)) * x_0
            matrix_B.write(f"{B[i]:.10e}\n")

        # Perform LU decomposition
        LU_find(L, U, A, n)

        # Write matrices L and U
        for i in range(1, n + 1):
            matrix_L.write("\t".join(f"{L[i][j]:.10e}" for j in range(1, n + 1)) + f"\t{L[i][n]:.10e}\n")
            matrix_U.write("\t".join(f"{U[i][j]:.10e}" for j in range(1, n + 1)) + f"\t{U[i][n]:.10e}\n")

        # Solve initial system for X
        LU_solve(L, U, B, X, n)

        # Calculate initial error norm
        for i in range(1, n + 1):
            dX[i] = X[i] - x_0
        r_eps = norma_vector(dX, n)
        print(f"Start solution = {r_eps:.10e}")

        # Iterative process for refining solution
        k = 0
        while True:
            new_B(A, X, B0, n)
            for i in range(1, n + 1):
                R[i] = B[i] - B0[i]

            LU_solve(L, U, R, dX, n)

            for i in range(1, n + 1):
                X[i] += dX[i]

            if norma_vector(dX, n) < eps and norma_vector(R, n) < eps:
                break
            if k >= kmax:
                print("Max iterations reached")
                break

            k += 1

        print(f"Number of iterations = {k}")

        # Write result to matrix_X.txt
        for i in range(1, n + 1):
            matrix_X.write(f"{X[i]:.10e}\n")

# Call main function
if __name__ == "__main__":
    main()
