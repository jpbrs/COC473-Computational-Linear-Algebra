import math
import argparse
from argparse import RawTextHelpFormatter

TOLm = 10**-7
c1,c2,c3,c4 = 1,1,1,1

#region Exceptions
class BissectionMethodException(Exception):
    pass

class NumberOfPointsException(Exception):
    pass
#endregion

#region Weigths Quadratura de Gauss
W = {2:{"points": [-0.5773502691896257,0.5773502691896257],"weigths": [1.0,1.0]},
        3:{"points": [0.0,-0.7745966692414834,0.7745966692414834],"weigths": [0.8888888888888888,0.5555555555555556,0.5555555555555556]},
        4:{"points": [-0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526],"weigths": [0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538]},
        5:{"points": [0.0,-0.5384693101056831,0.5384693101056831,-0.9061798459386640,0.9061798459386640],"weigths": [0.5688888888888889,0.4786286704993665,0.4786286704993665,0.2369268850561891,0.2369268850561891]},
        6:{"points": [0.6612093864662645,-0.6612093864662645,-0.2386191860831969,0.2386191860831969,-0.9324695142031521,0.9324695142031521],"weigths": [0.3607615730481386,0.3607615730481386,0.4679139345726910,0.4679139345726910,0.1713244923791704,0.1713244923791704]},
        7:{"points": [0.0,0.4058451513773972,-0.4058451513773972,-0.7415311855993945,0.7415311855993945,-0.9491079123427585,0.9491079123427585],"weigths": [0.4179591836734694,0.3818300505051189,0.3818300505051189,0.2797053914892766,0.2797053914892766,0.1294849661688697,0.1294849661688697]},
        8:{"points": [-0.1834346424956498,0.1834346424956498,-0.5255324099163290,0.5255324099163290,-0.7966664774136267,0.7966664774136267,-0.9602898564975363,0.9602898564975363],"weigths": [0.3626837833783620,0.3626837833783620,0.3137066458778873,0.3137066458778873,0.2223810344533745,0.2223810344533745,0.1012285362903763,0.1012285362903763]},
        9:{"points": [0.0,-0.8360311073266358,0.8360311073266358,-0.9681602395076261,0.9681602395076261,-0.3242534234038089,0.3242534234038089,-0.6133714327005904,0.6133714327005904],"weigths": [0.3302393550012598,0.1806481606948574,0.1806481606948574,0.0812743883615744,0.0812743883615744,0.3123470770400029,0.3123470770400029,0.2606106964029354,0.2606106964029354]},
        10:{"points": [-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,-0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,-0.9739065285171717,0.9739065285171717],"weigths": [0.2955242247147529,0.2955242247147529,0.2692667193099963,0.2692667193099963,0.2190863625159820,0.2190863625159820,0.1494513491505806,0.1494513491505806,0.0666713443086881,0.0666713443086881]}}
#endregion

#region Raiz
def BissectionMethod(function, p1, p2):
    if not function(p1)*function(p2)<=0: raise BissectionMethodException("Não é possível calcular a raiz pelo método de bisseção")

    result = p1
    while((p2 - p1) >= TOLm):
        result = (p1+p2)/2

        if(function(result) == 0):
            break
        
        if (function(result)*function(p1) < 0):
            p2 = result
        else:
            p1 = result

    print("Raíz da Função : {}".format(result))
#endregion Raiz

#region Integral
def GaussQuadratureMethod(function, p1, p2, n_points):
    if not(n_points<=10 and n_points >=2): raise NumberOfPointsException("{} não é um número de pontos válido. Apenas valores entre 2 e 10 serão aceitos".format(n_points))
    weigths = W.get(n_points).get("weigths")
    points = W.get(n_points).get("points")

    L = p2-p1
    x = []
    for i in range(n_points):
        x.append((p1+p2+points[i]*L)/2)
    result = 0
    for i in range(n_points):
        result += function(x[i])*weigths[i]
    print("Valor da Integral: {}".format(result*L/2))

#endregion Integral

#region Derivada DF
def DFCentral(function, x, deltaX):
    lim_sup = function(x + deltaX) - function(x - deltaX)
    lim_inf = 2*deltaX
    return lim_sup/lim_inf

def DFStepForward(function, x, deltaX):
    lim_sup = function(x + deltaX) - function(x)
    lim_inf = deltaX
    return lim_sup/lim_inf

def DFStepBackward(function, x, deltaX):
    lim_sup = function(x) - function(x - deltaX)
    lim_inf = deltaX
    return lim_sup/lim_inf

#endregion Derivada DF

#region Derivada Interpolação de Richard
def RichardsonExtrapolationCentral(function, x, deltaX, p=1):
    d1 = DFCentral(function, x, deltaX)
    deltaX2 = deltaX/2
    d2 = DFCentral(function, x, deltaX2)
    q = deltaX/deltaX2
    result = d1 + (d1 - d2)/(q**(-p)-1)
    return result

def RichardsonExtrapolationStepForward(function, x, deltaX, p=1):
    d1 = DFStepForward(function, x, deltaX)
    deltaX2 = deltaX/2
    d2 = DFStepForward(function, x, deltaX2)
    q = deltaX/deltaX2
    result = d1 + (d1 - d2)/(q**(-p)-1)
    return result

def RichardsonExtrapolationStepBackward(function, x, deltaX, p=1):
    d1 = DFStepBackward(function, x, deltaX)
    deltaX2 = deltaX/2
    d2 = DFStepBackward(function, x, deltaX2)
    q = deltaX/deltaX2
    result = d1 + (d1 - d2)/(q**(-p)-1)
    return result
#endregion Derivada Interpolação de Richardson

def function(x):
    return c1**(c2*x)+c3*(x**(c4))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""COC 473 - Computational Linear Algebra - Second Exam - Second Task \n
        General example of Usage: task2 ICOD specific-parameters --TOLm

        Example for Bissection Method: 
        task2 1 c1 c2 c3 c4 --firstPoint/-fp --secondPoint/-sp --TOLm 0.0001(default : 0.0000001)
        
        Example for Gauss Quadrature Method: 
        task2 2 c1 c2 c3 c4 --firstPoint/-fp --secondPoint/-sp --numberOfPoints/-np(between 2 and 10)

        Example for Finite Difference Derivative:
        task2 3 c1 c2 c3 c4 --xPoint/-x --xDelta/-d --method/-m(1 for Central, 2 for Stepforward, 3 for Backforward)
        
        Example for Richard Extrapolation Derivative Method:
        task2 4 c1 c2 c3 c4 --xPoint/-x --xDelta/-d --method/-m(1 for Central, 2 for Stepforward, 3 for Backforward)
    
    """,formatter_class=RawTextHelpFormatter)
    parser.add_argument('ICOD',type=int, help="Select the functionality wanted: 1 for root, 2 for integral, 3 for DF derivative, 4 for RE derivative")
    parser.add_argument('c1',type=int, help="Value of C1")
    parser.add_argument('c2',type=int, help="Value of C2")
    parser.add_argument('c3',type=int, help="Value of C3")
    parser.add_argument('c4',type=int, help="Value of C4")

    parser.add_argument('--firstPoint','-fp', type=float, help = "First point for Bissection and Gauss Quadrature Methods." )
    parser.add_argument('--secondPoint','-sp', type=float, help = "Second point for Bissection and Gauss Quadrature Methods." )
    parser.add_argument('--numberOfPoints','-np', type=int, help = "Number of Points for Gauss Quadrature Method." )
    parser.add_argument('--xPoint','-x', type=float, help = "X point value for Derivative Methods." )
    parser.add_argument('--xDelta','-d', type=float, help = "Delta X value for Derivative Methods." )
    parser.add_argument('--method','-m', type=int, help = "Derivative Method: 1 for Central, 2 for Stepforward, 3 for Backforward" )
    parser.add_argument('--TOLm','-t',type=float, help = "Maximum Residue Tolerance Number. Default : 10^-7" )

    args = parser.parse_args()
    
    ########################## Resolution Methods #############################

    if(args.TOLm): tolM = args.TOLm
    c1 = args.c1
    c2 = args.c2
    c3 = args.c3
    c4 = args.c4

    if args.ICOD == 1:
        if(args.firstPoint and args.secondPoint):
            print("\nRoot with Bissection Method")
            BissectionMethod(function,args.firstPoint,args.secondPoint)
        else:
            print("First and Second Point arguments required")
    elif args.ICOD == 2:
        if(args.firstPoint and args.secondPoint and args.numberOfPoints):
            print("\nGauss Quadrature Method")
            GaussQuadratureMethod(function, args.firstPoint, args.secondPoint, args.numberOfPoints)
        else:
            print("First and Second Point arguments required with also the number of points")
    elif args.ICOD == 3:
        if(args.xPoint and args.xDelta and args.method):
            if(args.method == 1):
                print("\nFinite Difference Derivative Method with Central Difference")
                print("Resultado: {}".format(DFCentral(function, args.xPoint, args.xDelta)))
            elif(args.method == 2):
                print("\nFinite Difference Derivative Method with StepForward")
                print("Resultado: {}".format(DFStepForward(function, args.xPoint, args.xDelta)))
            elif(args.method ==3):
                print("\nFinite Difference Derivative Method with StepBackward")
                print("Resultado: {}".format(DFStepBackward(function, args.xPoint, args.xDelta)))
            else:
                print("Method must be a number between 1 and 3")
        else:
            print("xPoint, xDelta and Method number are required")
    elif args.ICOD == 4:
        if(args.xPoint and args.xDelta and args.method):
            if(args.method == 1):
                print("\nRichardson Extrapolation Derivative Method with Central Difference")
                print("Resultado: {}".format(RichardsonExtrapolationCentral(function, args.xPoint, args.xDelta)))
            elif(args.method == 2):
                print("\nRichardson Extrapolation Derivative Method with StepForward")
                print("Resultado: {}".format(RichardsonExtrapolationStepForward(function, args.xPoint, args.xDelta)))
            elif(args.method ==3):
                print("\nRichardson Extrapolation Derivative Method with StepBackward")
                print("Resultado: {}".format(RichardsonExtrapolationStepBackward(function, args.xPoint, args.xDelta)))
            else:
                print("Method must be a number between 1 and 3")
        else:
            print("xPoint, xDelta and Method number are required")
    else:
        print("ICOD method must be a number between 1 and 4")
        
    print("")

# BissectionMethod(function,-100,1000)
# GaussQuadratureMethod(function, 1, 3, 2)
# print(DFStepForward(function, 6, 0.5))
# print(DFStepBackward(function, 6, 0.5))
# print(DFCentral(function, 6, 0.5))
# print(RichardsonExtrapolationStepForward(function, 6, 0.5, 2))
# print(RichardsonExtrapolationStepBackward(function, 6, 0.5, 2))
# print(RichardsonExtrapolationCentral(function, 6, 0.5, 2))



