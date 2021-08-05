  
from sys import exit
import argparse
from argparse import RawTextHelpFormatter

class InvalidLengthException(Exception):
    pass

def RegressaoLinear(x, y, n, p):
    nx = len(x)
    ny = len(y)

    try:
        if (nx != ny):
            raise InvalidLengthException("Erro, os arquivos de pontos X e Y possuem diferentes dimensões.")
        elif  (nx != n):
            raise InvalidLengthException("Erro, a dimensão N inserida não corresponde à dimensão dos arquivos")
        else:
            pass
    except Exception as ex:
        print(ex)
        exit()
        
    m11, m12, m21, m22 = float(n), 0.0, 0.0, 0.0

    for i in range(0, n):
        m12 += x[i]
        m21 += x[i]
        m22 += x[i]**2
        
    mDet = (m11*m22) - (m12*m21)
    mInv = [[m22/mDet, -m12/mDet],[-m21/mDet, m11/mDet]]

    v1 = 0
    v2 = 0

    for i in range(n):
        v1 += y[i]
        v2 += x[i]*y[i]

    V = [v1, v2]
    coef = MatVetMul(mInv, V)
    print("Coeficientes : {} , {}".format(round(coef[0],3), round(coef[1], 3)))
    reta = coef[1]*p + coef[0]
    
    print("Valor aproximado de Y no ponto x={} é {}".format(p, round(reta,3)))

    
def MatVetMul(matrix, vector):
    rows = len(matrix)
    columns = len(vector)
    sum = 0
    result = [0]*rows

    for j in range(rows):
        r = matrix[j]
        for i in range(columns):
            sum += r[i]*vector[i]
        result[j],sum = sum,0

    return result

def Lagrange(x, y, n, xp):

    nx = len(x)
    ny = len(y)

    try:
        if (nx != ny):
            raise InvalidLengthException("Erro, os arquivos de pontos X e Y possuem diferentes dimensões.")
        elif  (len(x) != n):
            raise InvalidLengthException("Erro, a dimensão N inserida não corresponde à dimensão dos arquivos")
        else:
            pass
    except Exception as ex:
        print(ex)
        exit()

    yp = 0.0

    for i in range(ny):
        
        p = 1.0
        
        for j in range(ny):
            if i != j:
                p = p * (xp - x[j])/(x[i] - x[j])
        
        yp = yp + p * y[i]    

    print("Valor aproximado de Y no ponto x={} é {}".format(xp, round(yp,3)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""COC 473 - Computational Linear Algebra Third Program \n
        Example of Usage: ALC2 Xfile.txt Yfile.txt ICOD N X
    """,formatter_class=RawTextHelpFormatter)
    parser.add_argument('xfile', help = "File with the list of X values" )
    parser.add_argument('yfile', help = "File with the list of Y values" )
    parser.add_argument('ICOD',type=int, help="Select Resolution Method's Number : 1 for Linear Regression, 2 for Lagrange Interpolation")
    parser.add_argument('N', type=int, help = "Dimension of the Dataset or Number of Points" )
    parser.add_argument('X', type=float, help = "X value of the point that we want to predict aproximately" )

    args = parser.parse_args()

    with open(args.xfile, 'r') as f:
        xfile = [[float(num) for num in line.split(',')] for line in f][0]
    with open(args.yfile, 'r') as f:
        yfile = [[float(num) for num in line.split(',')] for line in f][0]

    print("Cálculo da Função Aproximada")

    print("Conjunto de Pontos X fornercido : {}".format(xfile))
    print("Conjunto de Pontos Y fornercido : {}".format(yfile))

    
    ########################## Resolution Methods #############################
    if args.ICOD == 1:
        print("\nRegressão Linear")
        RegressaoLinear(xfile, yfile, args.N, args.X)
    if args.ICOD == 2:
        print("\nInterpolação de Lagrange")
        Lagrange(xfile, yfile, args.N, args.X)


    print("")