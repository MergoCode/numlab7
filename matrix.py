import random

def matrix_rand():
    n = 100
    a = 1.0
    b = 9.0

    with open("matrix_A.txt", "w") as matrix_A:
        for i in range(1, n + 1):
            for j in range(1, n):
                r_numb = a + b * random.random()
                matrix_A.write("{:e}\t".format(r_numb))
            r_numb = a + b * random.random()
            matrix_A.write("{:e}\n".format(r_numb))

matrix_rand()
