# Adi buchris
#Esty fridmen
#Mlky rotshild
import pandas as pd
import numpy as np
def neville(datax, datay, x):
    """
    Finds an interpolated value using Neville's algorithm.
    Input
      datax: input x's in a list of size n
      datay: input y's in a list of size n
      x: the x value used for interpolation
    Output
      p[0]: the polynomial of degree n
    """
    lbl_2 = 1
    n = len(datax)
    p = n*[0]
    for k in range(n):
        lbl_2=k
        for i in range(n-k):
            if k == 0:
                p[i] = datay[i]
            else:
                p[i] = ((x-datax[i+k])*p[i]+ (datax[i]-x)*p[i+1])/(datax[i]-datax[i+k])
                print("p",i,lbl_2,"= ",p[i])
                lbl_2 += 1
        print("----------")

    return p[0]

datax=[1.2,1.3,1.4,1.5,1.6]
datay=[1.5095,1.6984,1.9043,2.1293,2.3756]


def jacobi(A, b, x0, tol, n_iterations=300):
    """
    Performs Jacobi iterations to solve the line system of
    equations, Ax=b, starting from an initial guess, ``x0``.

    Returns:
    x, the estimated solution
    """

    n = A.shape[0]
    x = x0.copy()
    x_prev = x0.copy()
    counter = 0
    x_diff = tol + 1

    while (x_diff > tol) and (counter < n_iterations):  # iteration level
        for i in range(0, n):  # element wise level for x
            s = 0
            for j in range(0, n):  # summation for i !=j
                if i != j:
                    s += A[i, j] * x_prev[j]

            x[i] = (b[i] - s) / A[i, i]

        # update values
        counter += 1
        x_diff = (np.sum((x - x_prev) ** 2)) ** 0.5
        x_prev = x.copy()  # use new x for next iteration

    print("Number of Iterations: ", counter)
    print("Norm of Difference: ", x_diff)
    return x
def printMenu():
    print("Select the interpolation method to find the approximate solution-\n1.neville\n2.Cubic Spline\n3.End of program")
def cubic_spline(x, y, leftedge=0, rightedge=0, tol=1e-10):
            """
            Interpolate using  cubic splines.

            Generates a strictly diagonal dominant matrix then applies Jacobi's method.

            Returns coefficients:
            b, coefficient of x of degree 1
            c, coefficient of x of degree 2
            d, coefficient of x of degree 3
            """
            x = np.array(x)
            y = np.array(y)
            ### check if sorted
            if np.any(np.diff(x) < 0):
                idx = np.argsort(x)
                x = x[idx]
                y = y[idx]

            size = len(x)
            delta_x = np.diff(x)
            delta_y = np.diff(y)

            ### Get matrix A
            A = np.zeros(shape=(size, size))
            b = np.zeros(shape=(size, 1))
            b[0] = leftedge
            b[-1] = rightedge
            A[0, 0] = 1
            A[-1, -1] = 1

            for i in range(1, size - 1):
                A[i, i - 1] = (1 / 6) * delta_x[i - 1]
                A[i, i + 1] = (1 / 6) * delta_x[i]
                A[i, i] = (1 / 3) * (delta_x[i - 1] + delta_x[i])
                ### Get matrix b
                b[i, 0] = (delta_y[i] / delta_x[i] - delta_y[i - 1] / delta_x[i - 1])

            ### Solves for c in Ac = b
            print('Jacobi Method Output:')
            m = jacobi(A, b, np.zeros(len(A)), tol=tol, n_iterations=1000)

            return m.squeeze()

def find_x(xxx, m, x, y):
            size = len(x)
            i = 0
            while (i < size and xxx > x[i]):
                i = i + 1

            i = i - 1
            delta_x = np.diff(x)

            yyy = (y[i + 1] * (xxx - x[i])) / delta_x[i] - (y[i] * (xxx - x[i + 1])) / delta_x[i] + (m[i + 1] / 6) * (
                        ((xxx - x[i]) ** 3) / delta_x[i] - delta_x[i] * (xxx - x[i])) - (m[i] / 6) * (
                              ((xxx - x[i + 1]) ** 3) / delta_x[i] - delta_x[i] * (xxx - x[i + 1]))

            print("x =" + str(xxx) + " ,y= " + str(yyy))


while True:
    printMenu()
    x = int(input())
    if x==1:
        print("x = ", 1.37, ",y = ", neville(datax, datay, 1.37))
    elif x==2:

        x = [1.2, 1.3, 1.4, 1.5, 1.6]
        y = [1.5095, 1.6984, 1.9043, 2.1293, 2.3756]
        m = cubic_spline(x, y)
        find_x(1.37, m, x, y)

    elif x==3:
        exit(1)
    else:
        print("Error.")
        continue


