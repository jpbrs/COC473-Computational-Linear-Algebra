import copy
import math

example_matrix = [[5,-4,1,0],[-4,6,-4,1],[1,-4,6,-4],[0,1,-4,5]]
example_vector = [-1,2,1,3]

#region Exceptions

class NotSquaredException(Exception):
    pass

class NotSymmetricException(Exception):
    pass

class NotPositiveDefiniteException(Exception):
    pass

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

def Positive_Definite(matrix):
    for i in range(0,len(matrix)):
        m2 = [[matrix[row][column] for row in range(0,i+1)] for column in range(0,i+1)]
        if Determinant(m2) <= 0 : return False
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

#region SubstitutionMethods

def Backward_sub(matrix, vector):
    rows = len(matrix)
    result_matrix = [0 for i in range(rows)]
    result_matrix[rows-1] = vector[rows-1] / matrix[rows-1][rows-1]

    for i in range(rows-2, -1, -1):
        total = vector[i]
        for j in range(i+1, rows):
            total -= matrix[i][j]*result_matrix[j]

        result_matrix[i] = total/float(matrix[i][i])

    return result_matrix

def Forward_Sub(m1, m2, control=False):
    rows = len(m1)
    vector = [0 for i in range(rows)]

    if control :
        vector[0] = m2[0]/m1[0][0]
    else:
        vector[0] = m2[0]

    for i in range(1, rows):
        total = m2[i]
        for j in range(i):
            total -= m1[i][j]*vector[j]

        if(not control):
            vector[i] = total
        else:
            vector[i] = total/m1[i][i]

    return vector

#endregion

#region SystemsSolvingMethods

def LU(matrix, vector):
    columns = len(matrix[0])
    rows = len(matrix)

    if not Is_Squared(columns,rows): raise NotSquaredException("A matriz não é quadrada")

    decomp_matrix = copy.deepcopy(matrix)

    for k in range(rows):
        for i in range(k+1, rows):
            decomp_matrix[i][k] = float(decomp_matrix[i][k]/decomp_matrix[k][k])

        for j in range(k+1, columns):
            for i in range(k+1, columns):
                decomp_matrix[i][j] = float(decomp_matrix[i][j]-decomp_matrix[i][k]*decomp_matrix[k][j])

    y = Forward_Sub(decomp_matrix, vector)
    return Backward_sub(decomp_matrix, y) 

def Cholesky(matrix, vector):
    columns = len(matrix[0])
    rows = len(matrix)

    if not Is_Squared(columns,rows): raise NotSquaredException("A matriz não é quadrada")
    if not Is_Symmetric(matrix): raise NotSymmetricException("A matriz não é simétrica")
    if not Positive_Definite(matrix): raise NotPositiveDefiniteException("A matriz não é positivamente definida")

    decomp_matrix = [[0.0] * len(matrix) for _ in range(len(matrix))]

    for i in range(len(matrix)):
        for j in range(i + 1):
            if(i == j):
                total = sum(decomp_matrix[i][k]**2 for k in range(i))
                decomp_matrix[i][i] = (matrix[i][i]-total)**0.5
            else:
                total = sum(decomp_matrix[i][k]*decomp_matrix[j][k] for k in range(i))
                decomp_matrix[i][j] = (1.0/decomp_matrix[j][j])*(matrix[i][j]-total)

    y = Forward_Sub(decomp_matrix, vector, True)
    return Backward_sub(Transposed_Matrix(decomp_matrix), y)

#endregion

#region Tests

if __name__ == '__main__':
    print("Decomposição LU")
    print(LU(example_matrix,example_vector))
    print("#########################")
    print(Cholesky(example_matrix,example_vector))

#endregion



