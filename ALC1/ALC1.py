#region Imports
import copy
import math
import argparse
from argparse import RawTextHelpFormatter
#endregion

#region ExampleMatrixes
example_matrix = [[5,-4,1,0],[-4,6,-4,1],[1,-4,6,-4],[0,1,-4,5]]
example_vector = [-1,2,1,3]
example_jacobi_matrix = [[6,1,1,1,1],[1,7,1,1,1],[1,1,8,1,1],[1,1,1,9,1],[1,1,1,1,10]]
example_jacobi_vector = [-10,-6,0,8,18]
#endregion

#region Exceptions

class NotSquaredException(Exception):
    pass

class NotSymmetricException(Exception):
    pass

class NotPositiveDefiniteException(Exception):
    pass

class DiagonallyDominantException(Exception):
    pass

class GeneralConvergenceException(Exception):
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

def Diagonally_Dominant(matrix):
    n = len(matrix)
    for i in range(0,n):
        diagonal = matrix[i][i]
        sigma_l = 0
        sigma_c = 0
        for j in range(0,n):
            if (i != j):
                sigma_l += math.fabs(matrix[i][j])
                sigma_c += math.fabs(matrix[j][i])

        if(sigma_l > diagonal or sigma_c > diagonal):
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

def Vect_Mul(vector_a, vector_b):
    result = 0

    for i in range(len(vector_a)):
        for j in range(len(vector_b)):
            if (i == j):
                result += vector_a[i]*vector_b[i]

    return result

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

def LU(matrix, vector, isRounded = True):
    columns = len(matrix[0])
    rows = len(matrix)

    try:
        if not Is_Squared(columns,rows): raise NotSquaredException("A matriz não é quadrada")
    except Exception as ex:
        print(ex)
        exit()

    decomp_matrix = copy.deepcopy(matrix)

    for k in range(rows):
        for i in range(k+1, rows):
            decomp_matrix[i][k] = float(decomp_matrix[i][k]/decomp_matrix[k][k])

        for j in range(k+1, columns):
            for i in range(k+1, columns):
                decomp_matrix[i][j] = float(decomp_matrix[i][j]-decomp_matrix[i][k]*decomp_matrix[k][j])

    y = Forward_Sub(decomp_matrix, vector)
    solution = Backward_sub(decomp_matrix, y)

    if isRounded:
        for i in range(0,len(solution)):
            solution[i] = round(solution[i],3)

    print("Solution: ",solution)

def Cholesky(matrix, vector, isRounded =  True):
    columns = len(matrix[0])
    rows = len(matrix)

    try:
        if not Is_Squared(columns,rows): raise NotSquaredException("A matriz não é quadrada")
        if not Is_Symmetric(matrix): raise NotSymmetricException("A matriz não é simétrica")
        if not Positive_Definite(matrix): raise NotPositiveDefiniteException("A matriz não é positivamente definida")
    except Exception as ex:
        print(ex)
        exit()

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
    solution = Backward_sub(Transposed_Matrix(decomp_matrix), y)

    if isRounded:
        for i in range(0,len(solution)):
            solution[i] = round(solution[i],3)

    print("Solution: ",solution)

def Iterative_Jacobi(matrix, vector, isRounded = True, tol = 10**(-5)):
    try:
        if not Diagonally_Dominant(matrix): raise DiagonallyDominantException("Warning : A matriz não é Diagonal Dominante")
    except Exception as ex:
        print(ex)
        exit()

    n = len(matrix)
    result = [0.0 for i in range(n)]
    zero = [0.0 for i in range(n)]
    residue = 1
    iteration = 0

    while (residue > tol):
        numerator = 0
        denominator = 0

        for j in range(n):
            result[j] = vector[j]
            for k in range(n):
                if (j != k):
                    result[j] += (-1)*(matrix[j][k] * zero[k])
            result[j] /= matrix[j][j]

        for z in range(n):
            numerator += (result[z]-zero[z])**2
            denominator += result[z]**2

        residue = float(numerator**0.5)/(denominator**0.5)

        for i in range(len(result)):
            zero[i] = result[i]
        iteration += 1

    if isRounded:
        for i in range(0,len(result)):
            result[i] = round(result[i],3)
        residue = round(residue,3)

    print("Solution: ", result)
    print("Residue: ", residue)
    print("Number of iterations: ", iteration)

def Gauss_Sidel(matrix, vector, isRounded = True, tol = 10**(-5)):
    try:
        if not ( Diagonally_Dominant(matrix) or ( Is_Symmetric(matrix) and Positive_Definite(matrix) ) ): raise GeneralConvergenceException("Warning : O algoritmo não converge para essa matriz")
    except Exception as ex:
        print(ex)
        exit()

    n = len(matrix)

    solution_zero = [1.0 for i in range(n)]
    solution = [0.0 for i in range(n)]
    residue = 1
    iteration = 0

    while (residue > tol):
        numerator = 0
        denominator = 0
        second_summation = 0
        second_summation = 0

        for j in range(n):
            first_summation = Vect_Mul(matrix[j][:j], solution[:j])
            second_summation = Vect_Mul(matrix[j][j+1:], solution_zero[j+1:])
            solution[j] = (vector[j] - first_summation - second_summation)/matrix[j][j]

        for z in range(n):
            numerator += (solution[z] - solution_zero[z])**2
            denominator += solution[z]**2

        residue = float(numerator**0.5)/(denominator**0.5)

        for i in range(len(solution)):
            solution_zero[i] = solution[i]

        iteration += 1

    if isRounded:
        for i in range(0,len(solution)):
            solution[i] = round(solution[i],3)
        residue = round(residue,3)

    print("Solution: ", solution)
    print("Residue: ", residue)
    print("Number of Iterations: ", iteration)
#endregion

#region ArgParse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""COC 473 - Computational Linear Algebra First Program \n
        Example of Usage: monolith.py matrixfile.txt vectorfile.txt --ICOD 1 --IDET 1
    """,formatter_class=RawTextHelpFormatter)
    parser.add_argument('matrix', help = "Linear System Matrix that must be solved" )
    parser.add_argument('vector', help = "Result Vector of the Linear System" )
    parser.add_argument('ICOD',type=int, help="Select Resolution Method's Number : 1 for LU, 2 for Cholesky, 3 for Jacobi, 4 for Gauss-Sidel")
    parser.add_argument('--IDET','-d',type=int, help="Select a number greater than zero for return the Determinant of the Matrix")
    parser.add_argument('--TOLm','-t',type=float, help="Tolerance Number of Jacobi and Gauss-Sidel Methods. Default : 10**(-5)")
    parser.add_argument('--long','-l',help="Return the long number in scientific notation. Default: Rounded Number",action='store_true')

    args = parser.parse_args()

    with open(args.matrix, 'r') as f:
        list_m = [[float(num) for num in line.split(',')] for line in f]
    with open(args.vector, 'r') as f:
        list_v = [[float(num) for num in line.split(',')] for line in f][0]

    print("Matriz Fornecida :")
    for item in list_m:
        print(item)

    print("\nVetor Fornecido :")
    print(list_v)
    
    ########################## Resolution Methods #############################
    if args.ICOD == 1:
        print("\nDecomposição LU")
        LU(list_m,list_v,not args.long)
    if args.ICOD == 2:
        print("\nDecomposição de Cholesky")
        Cholesky(list_m,list_v,not args.long)
    if args.ICOD == 3:
        print("\nMétodo Iterativo de Jacobi")
        if args.TOLm:
            Iterative_Jacobi(list_m,list_v, not args.long, args.TOLm)
        else:
            Iterative_Jacobi(list_m,list_v, not args.long)
    if args.ICOD == 4:
        print("\nMétodo de Gauss Seidel")
        if args.TOLm:
            Gauss_Sidel(list_m,list_v, not args.long, args.TOLm)
        else:
            Gauss_Sidel(list_m,list_v, not args.long)

    if args.IDET:
        print("Determinant :",Determinant(list_m))

    print("")

#endregion



