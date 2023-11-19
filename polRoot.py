import argparse
import sys
import numpy as np

parser = argparse.ArgumentParser()

# Arg Flags 
def check_positive(value):
    iValue = int(value)
    if iValue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return iValue

parser.add_argument('-newt', dest='isNewtEnabled', action='store_const',
                    const=True, default=False,
                    help='Uses Newtons method to find the root')
parser.add_argument('-sec', dest='isSecEnabled', action='store_const',
                    const=True, default=False,
                    help='Uses Secant method to find the root')
parser.add_argument('-hybrid', dest='isHybridEnabled', action='store_const',
                    const=True, default=False,
                    help='Uses mybrid between Newton and Second method to find the root')
parser.add_argument('-maxIt', dest='n', type=check_positive, default=10000,
                    help='Maximum number of iterations for finding root')

parser.add_argument('points', type=float, nargs='+')
parser.add_argument('file', type=argparse.FileType('r')) 

args = parser.parse_args()

# Error handling for multiple flags
numFlags = 0 
if args.isNewtEnabled :
    numFlags += 1
if args.isSecEnabled :
    numFlags += 1
if args.isHybridEnabled :
    numFlags += 1
if numFlags>1: 
    parser.error("Multiple methods detected, please choose one")

# Setting parameters as needed based on chosen method
maxIterations = args.n
initP1 = 0.0
initP2 = 0.0
filename = ""

if args.isNewtEnabled :
    if len(args.points) != 1:
        parser.error("Invalid argument number")
    else :
        initP1 = float(args.points[0])
        filename = sys.argv[1] 
else :
    if len(args.points) != 2:
        parser.error("Invalid argument number")
    else :
        initP1 = float(args.points[0])
        initP2 = float(args.points[1])
        filename = sys.argv[2] 

# Extract Input from file 
power = int(args.file.readline())
input = []

line =args.file.readline()
for value in line.split():
    input.append(float(value))

# Establish create function method
def function(x):
    sum = 0
    expo=power
    for value in input:
        sum += value * x ** expo
        expo -= 1
    return sum

def derivativeFunction(x):
    var = np.poly1d(input) 
    derivative = var.deriv() 
    return derivative(x)

# Declaring functions for root finding methods
def bisectionMethod(a, b, tol):
    fOfA = function(a)
    fOfB = function(b)

    if np.sign(fOfA) == np.sign(fOfB):
        raise Exception(
         "The scalars a and b do not bound a root")
        
    # get midpoint
    m = (a + b)/2
    fOfM = function(m)

    if abs(fOfA - fOfB) < 1e-14:
        return m
    elif tol == 0:
        # stopping condition, report m as root
        return m
    elif np.sign(fOfA) == np.sign(fOfM):
        # case where m is an improvement on a. 
        # Make recursive call with a = m
        return bisectionMethod(m, b, tol-1)
    elif np.sign(fOfB) == np.sign(fOfM):
        # case where m is an improvement on b. 
        # Make recursive call with b = m
        return bisectionMethod(a, m, tol-1)

def newtonsMethod(a, tol):
    derivativeValue = derivativeFunction(a)
    if derivativeValue == 0.0:
        return None
    
    value = a - function(a)/derivativeFunction(a)

    if abs(value - a) < 1e-14:
        return a
    elif tol == 0:
        return a
    else:
        return newtonsMethod(value, tol-1)

def secantMethod(a, b, tol):
    if function(b) - function(a) == 0:
        return a
    value = a - function(a) * (b-a) / (function(b) - function(a))

    if abs(b - a) < 1e-14:
        return a
    elif tol == 0:
        return a
    else:
        return secantMethod(value, a, tol-1)

def hybridMethod(a, b, tol):
    if tol < 2:
        return bisectionMethod(a, b, tol)
    value = bisectionMethod(a, b, tol/2 + tol/2 % 2)
    return newtonsMethod(value, tol/2)

# Main
solution = 0

if args.isNewtEnabled:
    solution = newtonsMethod(initP1, maxIterations)
elif args.isSecEnabled:
    solution = secantMethod(initP1, initP2, maxIterations)
elif args.isHybridEnabled:
    solution = hybridMethod(initP1, initP2, maxIterations)
else:
    solution = bisectionMethod(initP1, initP2, maxIterations)

print(solution)

# Write solution to file
output_file = open('./fun1.sol','a')

if solution == None:
    output_file.write(str(solution))
    output_file.write(" ")
    output_file.write(str(maxIterations))
    output_file.write(" fail\n")
else:
    output_file.write(str(solution))
    output_file.write(" ")
    output_file.write(str(maxIterations))
    output_file.write(" success\n")

output_file.close() 