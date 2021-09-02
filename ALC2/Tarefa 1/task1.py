import math
import numpy as np
import argparse
from argparse import RawTextHelpFormatter


tolM = (10**(-7))
maxIter = 100

theta1 = 0.0
theta2 = 3.0

#region Utilitarios
def derivative(funcao, dv):
    lim_inf = 10**(-10)
    lim_sup = funcao(dv + h) - funcao(dv)
    return lim_sup/lim_inf

def parcialDerivative(multiVariableFunction, initialSolutionVector, dvIndex):
    lim_inf = 10**(-10)
    auxVec = initialSolutionVector[:]
    auxVec[dvIndex] += lim_inf
    auxFun = multiVariableFunction(initialSolutionVector)
    lim_sup = multiVariableFunction(auxVec) - auxFun
    return lim_sup/lim_inf

def norma(vetor):
    norma = 0
    for i in range(len(vetor)): norma += vetor[i]*vetor[i]
    return norma**(0.5)
#endregion

#region GetMatrixes
def GetJacobianMatrix(functionsList, initialSolutionVector):
    dimensionOfFunctions = len(functionsList)
    dimensionOfVector = len(initialSolutionVector)
    jacob = [[0 for i in range(dimensionOfVector)] for j in range(dimensionOfFunctions)]
    for i in range(dimensionOfFunctions):
        for j in range(dimensionOfVector):
            jacob[i][j] = parcialDerivative(functionsList[i], initialSolutionVector, j)
    return jacob

def GetFVector(functionsList, initialSolutionVector):
    dimensionOfFunctions = len(functionsList)
    vectorF = [0 for i in range(len(functionsList))]
    for i in range(len(functionsList)): vectorF[i] = functionsList[i](initialSolutionVector)
    return vectorF
#endregion

#region MatrixOperations
def MatVetMul(matrix, vector):
    rows = len(matrix)
    result = [0]*rows
    irange = range(len(vector))
    sum = 0
    for j in range(rows):
        r = matrix[j]
        for i in irange:
            sum += r[i]*vector[i]
        result[j],sum = sum,0
    return result

# Soma dois vetores
def VectorSum(v1, v2):
    sum = [0 for i in range(len(v1))]
    for i in range(len(v1)): sum[i] = v1[i] + v2[i]
    return sum

def VectorMul(v1, v2):
    prod = 0
    for i in range(len(v1)):
        for j in range(len(v2)):
            if (i == j):
                prod += v1[i]*v2[i]
    return prod

def MatrixSum(m1, m2):
    dimension = len(m1)
    result = [[0 for i in range(dimension)] for j in range(dimension)]

    for i in range(dimension):
        for j in range(dimension):
            result[i][j] = m1[i][j] + m2[i][j]

    return result
#endregion


#region Metodos
def NonLinearNewton(functionsList, initialSolutionVector):
    xi = initialSolutionVector
    for i in range(maxIter):
        matrizjacob = GetJacobianMatrix(functionsList, xi)
        vectorF = GetFVector(functionsList, xi)
        inversajacob = np.linalg.inv(matrizjacob)
        deltaX = MatVetMul(inversajacob, vectorF)
        deltaX = [x*(-1) for x in deltaX]
        xi = VectorSum(xi, deltaX)
        tolMk = norma(deltaX)/norma(xi)
        if (tolMk < tolM):
            break
    if (i == maxIter-1):
        print("Algoritmo não convergiu para o número máximo de iterações sendo {}".format(maxIter))

    print("Numero total de iteracoes: {}".format(i))
    print("Solução para as constantes: {}".format(xi))

def Broyden(functionsList, initialSolutionVector):
    dimension = len(initialSolutionVector)
    xi = initialSolutionVector

    for iteration in range(maxIter):

        matrizjacob1 = GetJacobianMatrix(functionsList, xi)
        matrizjacob2 = GetJacobianMatrix(functionsList, xi)

        for k in range(dimension):
            for j in range(dimension):
                matrizjacob2[k][j] = matrizjacob1[k][j]

        vector_f = GetFVector(functionsList, xi)

        inversamatrizjacob2 = np.linalg.inv(matrizjacob2)
        deltaX = MatVetMul(inversamatrizjacob2, vector_f)
        deltaX = [xi*(-1) for xi in deltaX]
        xi = VectorSum(xi, deltaX)

        vector_f2 = GetFVector(functionsList, xi)
        yK = VectorSum(vector_f2, [xi*(-1) for x in vector_f])
        residue = norma(deltaX)/norma(xi)
        if (residue < tolM):     print("Numero total de iteracoes: {}".format(iteration)); print("Solução para as constantes: {}".format(xi)); return
        matmul = MatVetMul(matrizjacob1, deltaX)
        matmul = [xi*(-1) for xi in matmul]

        sup = VectorSum(yK, matmul)
        sumDimension = len(sup)
        numerator = [[0 for i in range(sumDimension)] for j in range(sumDimension)]
        denominator = VectorMul(deltaX, deltaX)

        for i in range(sumDimension):
            for j in range(sumDimension):
                numerator[i][j] = sup[i]*deltaX[j]
        for k in range(len(matrizjacob1)):
            for j in range(len(matrizjacob1)):
                matrizjacob1[k][j] = numerator[k][j]/denominator

        matrizjacob1 = MatrixSum(matrizjacob2, matrizjacob1)

    if (iteration == maxIter-1):
        print("Algoritmo não convergiu para o número máximo de iterações sendo {}".format(maxIter))

#endregion

#region Funcoes
def funcao1(x):
    return 2*x[1]*x[1] + x[0]*x[0] + 6*x[2]*x[2] - 1

def funcao2(x):
    return 8*x[1]**3 + 6*x[1]*x[0]**2 + 36*x[1]*x[0]*x[2]+108*x[1]*x[2]**2 - theta1

def funcao3(x):
    return 60*(x[1]**4) + 60*(x[1]**2)*(x[0]**2)+576*(x[1]**2)*x[0]*x[2]+2232*(x[1]**2)*x[2]**2+252*(x[2]**2)*x[0]**2+1296*(x[2]**2)*x[0]+3348*(x[2]**4)*24*x[0]**3*(x[2])+3*x[0] - theta2
#endregion

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""COC 473 - Computational Linear Algebra - Second Exam - First Task \n
        Example of Usage: task1 ICOD theta1 theta2 --tolM
    """,formatter_class=RawTextHelpFormatter)
    parser.add_argument('ICOD',type=int, help="Select Resolution Method's Number : 1 for Non-Linear Newton Method, 2 for Broyden Method")
    parser.add_argument('Theta1', type=float, help = "Final result of Second Equation." )
    parser.add_argument('Theta2', type=float, help = "Final result of Third Equation." )
    parser.add_argument('--TOLm','-t',type=float, help = "Maximum Residue Tolerance Number. Default : 10^-7" )

    args = parser.parse_args()
    
    ########################## Resolution Methods #############################

    if(args.TOLm): tolM = args.TOLm
    theta1 = args.Theta1
    theta2 = args.Theta2

    if args.ICOD == 1:
        print("\nNon Linear Newton Method")
        NonLinearNewton([funcao1, funcao2, funcao3], [0.6874,0.3980,0.1908])
    if args.ICOD == 2:
        print("\nBroyden Method")
        Broyden([funcao1, funcao2, funcao3], [2,3,4])
    print("")

