import math
from sys import exit
import argparse
from argparse import RawTextHelpFormatter
import copy


example_matrix = [[3,2,0],[2,3,-1],[0,-1,3]]

#region Exceptions

class NotSquaredException(Exception):
    pass

class NotSymmetricException(Exception):
    pass

#endregion

#region MatrixCalculations

def Determinant(matrix):
    columns = len(matrix[0])
    rows = len(matrix)

    if not Is_Squared(columns,rows): raise NotSquaredException("A matriz não é quadrada")

    if rows == 1:
        return matrix[0][0]
    else:
        determinant = 0
        for i in range(0,rows):
            determinant += matrix[i][0]*((-1)**i) * Determinant(Auxiliar(matrix, i))
        return determinant

#endregion

#region Conditions

def Is_Squared(i,j):
    if i == j:
        return True
    else:
        return False

def Is_Symmetric(matrix):
    columns = len(matrix[0])
    rows = len(matrix)

    if not Is_Squared(columns,rows): raise NotSquaredException("A matriz não é quadrada")

    for i in range(0,len(matrix)):
        for j in range(0,len(matrix[i])):
            if(matrix[i][j] != matrix[j][i]):
                return False

    return True

#endregion

#region MatrixManipulations

def Auxiliar(matrix, i : int):
    copy_matrix = copy.deepcopy(matrix)

    for row in range(0,len(copy_matrix)):
        copy_matrix[row] = copy_matrix[row][1:]

    return copy_matrix[:i] + copy_matrix[i+1:]

def Transposed_Matrix(matrix):
    answer = [[0.0]*len(matrix) for i in range(len(matrix[0]))]

    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            answer[j][i] = matrix[i][j]

    return answer

def Highest_Element(matrix):
    highest_value = -float("inf")
    rows = len(matrix)

    for i in range(rows):
        for j in range(rows):
            if highest_value < math.fabs(matrix[i][j]) and i != j:
                highest_value = math.fabs(matrix[i][j])
                index = (i, j)
    return index

def MatVect_Mul(matrix_a, vector):
    number_of_columns_of_a = len(matrix_a)
    number_of_columns_of_vector = range(len(vector))
    result = [0.0 for _ in range(number_of_columns_of_a)]

    for j in range(number_of_columns_of_a):
        summation = 0

        for i in number_of_columns_of_vector:
            summation += matrix_a[j][i]*vector[i]

        result[j] = summation

    return result

def MatMul(m1, m2):
    m1_rows = len(m1)
    m2_rows = len(m2)
    m2_columns = len(m2[0])

    result = [[0.0 for _ in range(m1_rows)]
              for _ in range(m1_rows)]
    for i in range(0, m1_rows):
        for j in range(0, m2_columns):
            for k in range(0, m2_rows):
                result[i][j] += m1[i][k] * m2[k][j]
    return result

def P_Matrix(matrix, indexes):
    rows = len(matrix)
    phi = 0

    p_matrix = [[0.0 for _ in range(rows)] for _ in range(rows)]
    for i in range(rows):
        p_matrix[i][i] = 1.0

    denominator = (matrix[indexes[0]][indexes[0]] -
                   matrix[indexes[1]][indexes[1]])

    if(matrix[indexes[0]][indexes[0]] == matrix[indexes[1]][indexes[1]]):
        phi = math.pi/4
    else:
        phi = math.atan(2*matrix[indexes[0]][indexes[1]]/denominator)/2

    p_matrix[indexes[0]][indexes[0]] = math.cos(phi)
    p_matrix[indexes[1]][indexes[1]] = math.cos(phi)
    p_matrix[indexes[0]][indexes[1]] = -math.sin(phi)
    p_matrix[indexes[1]][indexes[0]] = math.sin(phi)

    return p_matrix

#endregion

#region EigenValueMethods

def Metodo_Potencia(matrix, tolm = 10**-5):
    rows = len(matrix)
    autovetor = [1.0 for _ in range(0, rows)]

    y_vector = MatVect_Mul(matrix, autovetor)
    autovalor = y_vector[0]

    for i in range(0, rows):
        y_vector[i] = y_vector[i]/autovalor
        autovetor = y_vector

    first = 1
    residue = math.fabs(autovalor - first)/autovalor

    iterations = 1
    while (residue >= tolm):
        first = autovalor
        y_vector = MatVect_Mul(matrix, autovetor)
        autovalor = y_vector[0]

        for i in range(rows):
            y_vector[i] = y_vector[i]/autovalor

        autovetor = y_vector
        residue = math.fabs(autovalor-first)/autovalor

        iterations += 1

    print("Autovalor: " + str(autovalor))
    print("Autovetor: " + str(autovetor))
    print("Numero de Iterações: " + str(iterations))


def Metodo_Jacobi(matrix, tolm = 10**(-5)):

    try:
        if not Is_Symmetric(matrix): raise NotSymmetricException("A matriz inserida não é simétrica")
    except Exception as ex:
        print(ex)
        exit()


    highest = Highest_Element(matrix)

    rows = len(matrix)
    autovetor = [[float(i == j) for j in range(rows)] for i in range(rows)]

    iterations = 1
    result = []

    while (math.fabs(matrix[highest[0]][highest[1]]) > tolm):
        p_matrix = P_Matrix(matrix, highest)
        p_matrix_transposed = Transposed_Matrix(p_matrix)
        matrix = MatMul(p_matrix_transposed, MatMul(matrix, p_matrix))
        autovetor = MatMul(autovetor, p_matrix)
        highest = Highest_Element(matrix)

        iterations += 1

    for i in range(rows):
        result.append((i+1, matrix[i][i], autovetor[i]))
        temp_autovalor = matrix[i][i]
        print("Autovalor Número {}: {}".format(i+1,temp_autovalor))
        print("Autovetor: {}".format(autovetor[i])+"\n")

    print("Numero de Iterações: {}".format(iterations))

    return result

#endregion

#region Argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""COC 473 - Computational Linear Algebra Second Program \n
        Example of Usage: ALC2 matrixfile.txt --ICOD 1 --IDET 1
    """,formatter_class=RawTextHelpFormatter)
    parser.add_argument('matrix', help = "Matrix for Eigenvalues/vectors calculation" )
    parser.add_argument('ICOD',type=int, help="Select Resolution Method's Number : 1 for Power, 2 for Jacobi")
    parser.add_argument('--IDET','-d',type=int, help="Select a number greater than zero for return the Determinant of the Matrix")
    parser.add_argument('--TOLm','-t',type=float, help="Tolerance Number of Jacobi and Power Methods. Default : 10**(-5)")

    args = parser.parse_args()

    with open(args.matrix, 'r') as f:
        list_m = [[float(num) for num in line.split(',')] for line in f]

    print("Cálculo dos Autovalores e Autovetores")

    print("Matriz Fornecida :")
    for item in list_m:
        print(item)
    
    ########################## Resolution Methods #############################
    if args.ICOD == 1:
        print("\n Método da Potência")
        if args.TOLm:
            Metodo_Potencia(list_m, args.TOLm)
        else:
            Metodo_Potencia(list_m)
    if args.ICOD == 2:
        print("\nMétodo de Jacobi")
        if args.TOLm:
            Metodo_Jacobi(list_m, args.TOLm)
        else:
            Metodo_Jacobi(list_m)

    if args.IDET:
        print("Determinant :",Determinant(list_m))

    print("")

#endregion

# Metodo_Potencia(example_matrix)
# Metodo_Jacobi(example_matrix)




