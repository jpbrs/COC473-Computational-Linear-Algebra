import copy
import math

example_matrix = [[5,-4,1,0],[-4,6,-4,1],[1,-4,6,-4],[0,1,-4,5]]

def Is_Squared(i,j):
    if i == j:
        return True
    else:
        return False

def Auxiliar(matrix, i : int):
    copy_matrix = copy.deepcopy(matrix)

    for row in range(0,len(copy_matrix)):
        copy_matrix[row] = copy_matrix[row][1:]

    return copy_matrix[:i] + copy_matrix[i+1:]

def Determinant(matrix):
    columns = len(matrix[0])
    rows = len(matrix)

    assert Is_Squared(columns,rows), "A matriz não é quadrada e portanto não é possível calcular seu determinante"

    if rows == 1:
        return matrix[0][0]
    else:
        determinant = 0
        for i in range(0,rows):
            determinant += matrix[i][0]*((-1)**i) * Determinant(Auxiliar(matrix, i))
        return determinant

def LU(matrix):
    columns = len(matrix[0])
    rows = len(matrix)

    assert Is_Squared(columns,rows), "A matriz não é quadrada e portanto não é possível calcular seu determinante"

    result = copy.deepcopy(matrix)

    for k in range(rows):
        for i in range(k+1, rows):
            result[i][k] = float(result[i][k]/result[k][k])

        for j in range(k+1, columns):
            for i in range(k+1, columns):
                result[i][j] = float(result[i][j]-result[i][k]*result[k][j])

    return result

def Is_Symmetric(matrix):
    columns = len(matrix[0])
    rows = len(matrix)

    assert Is_Squared(columns,rows), "A matriz não é quadrada e portanto não é possível calcular seu determinante"

    for i in range(0,len(matrix)):
        for j in range(0,len(matrix[i])):
            if(matrix[i][j] != matrix[j][i]):
                return False

    return True

def Positive_Definite(matrix):
    for i in range(0,len(matrix)):
        m2 = [[matrix[row][column] for row in range(0,i+1)] for column in range(0,i+1)]
        if Determinant(m2) <= 0 : return False
    return True

def Cholesky(matrix):
    columns = len(matrix[0])
    rows = len(matrix)

    assert Is_Squared(columns,rows), "A matriz não é quadrada e portanto não é possível calcular seu determinante"
    assert Is_Symmetric(matrix), "A matriz não é simétrica"
    assert Positive_Definite(matrix), "A matriz não é positivamente definida"

    result = [[0.0] * len(matrix) for _ in range(len(matrix))]

    for i in range(len(matrix)):
        for j in range(i + 1):

            if(i == j):
                summation = sum(result[i][k]**2 for k in range(i))
                result[i][i] = (matrix[i][i]-summation)**0.5
                continue

            summation = sum(result[i][k]*result[j][k] for k in range(i))
            result[i][j] = (1.0/result[j][j])*(matrix[i][j]-summation)

    return result

print("Decomposição LU")
for i in LU(example_matrix):
    print(i)

print("#########################")

print("Decomposição Cholesky")
for i in Cholesky(example_matrix):
    print(i)


