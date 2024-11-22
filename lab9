import matplotlib.pyplot as plt
import math

def f1(x):
    return x**3 - 2*(x**2)- 4*x

def tabulation(a, b, N):
    xl, yl = [], []
    h = (b - a)/N
    for i in range(N+1):
        x = a + i*h
        xl.append(x)
        yl.append(f1(x))
    return xl, yl

N = 100
k = 17
xa = math.pi/2
xn = 1.5
eps = 10**(-10)

xL, yL = tabulation(k-7, k+3, N)
plt.plot(xL, yL)
plt.legend(['f(x)'])
plt.show()

def F(x):
    return math.cos(x)

def f(x):
    return -math.sin(x)

def f2(x):
    return -math.cos(x)

def simpleIteration(x):
    x0 = x
    t = -1/f(x0)
    x1 = x0 + t*F(x0)
    k = 1
    while (abs(x1 - x0) > eps) and (abs(f(x1)) > eps):
        x0 = x1
        x1 = x0 + t*F(x0)
        k += 1
    return x1, k

def newton(x):
    x0 = x
    x1 = x0 - F(x0)/f(x0)
    k = 1
    while (abs(x1 - x0) > eps) and (abs(F(x1)) > eps):
        x0 = x1
        x1 = x0 - F(x0)/f(x0)
        k += 1
    return x1, k

def chebyshev(x):
    x0 = x
    x1 = x0 - F(x0)/f(x0) - ((F(x0)**2) * f2(x0))/(2*f(x0)**3)
    k = 1
    while (abs(x1 - x0) > eps) and (abs(F(x1)) > eps):
        x0 = x1
        x1 = x0 - F(x0)/f(x0) - ((F(x0)**2) * f2(x0))/(2*f(x0)**3)
        k += 1
    return x1, k

def hord(x):
    x0 = x
    t = -1/f(x)
    x1 = x0 + t*F(x)
    x2 = x1 - F(x1)*(x1 - x0)/(F(x1) - F(x0))
    k = 1
    while (abs(x1 - x0) > eps) and (abs(F(x1)) > eps):
        x0 = x1
        x1 = x2
        x2 = x1 - F(x1)*(x1 - x0)/(F(x1) - F(x0))
        k += 1
    return x2, k

def reverseIteration(x):
    x0 = x
    t = -1/f(x)
    x1 = x0 + t*F(x)
    x2 = -x0*F(x1)/(F(x0) - F(x1)) - x1*F(x0)/(F(x1) - F(x0))
    k = 1
    while (abs(x1 - x0) > eps) and (abs(F(x1)) > eps):
        x0 = x1
        x1 = x2
        x2 = -x0*F(x1)/(F(x0) - F(x1)) - x1*F(x0)/(F(x1) - F(x0))
        k += 1
    return x2, k

print(simpleIteration(xn))
print(newton(xn))
print(chebyshev(xn))
print(hord(xn))
print(reverseIteration(xn))

def horner(coefficients, x):
    """Схема Горнера для обчислення значення полінома та його похідної"""
    n = len(coefficients)
    b = coefficients[0]
    c = coefficients[0]
    for i in range(1, n-1):
        b = coefficients[i] + b * x
        c = coefficients[i] + c * x
    b = coefficients[-1] + b * x
    return b, c  

def newton_horner(coefficients, x0, eps):
    """Метод Ньютона з використанням схеми Горнера для дійсних коренів"""
    x = x0
    k = 0
    while True:
        p, dp = horner(coefficients, x)
        x_new = x - p/dp
        k += 1
        if abs(x_new - x) < eps:
            break
        x = x_new
    return x, k

coefficients = [1, -6, 11, -6]
x0 = 3.5  
eps = 1e-10
root, iterations = newton_horner(coefficients, x0, eps)
print(f"Знайдений корінь: {root}, число ітерацій: {iterations}")

def lin_method(coefficients, alpha0, beta0, eps):
    """Метод Ліна для знаходження комплексних коренів"""
    alpha, beta = alpha0, beta0
    k = 0
    max_iterations = 1000 

    while k < max_iterations:
        p = [0] * len(coefficients)
        q = [0] * len(coefficients)
        p[0] = coefficients[0]
        p[1] = coefficients[1] + alpha * p[0]

        for i in range(2, len(coefficients)):
            p[i] = coefficients[i] + alpha * p[i-1] + beta * p[i-2]

        q[0] = p[0]
        for i in range(1, len(coefficients)-1):
            q[i] = p[i] + alpha * q[i-1]

        if abs(q[-2]) < eps:
            print("Помилка: похідна дуже мала або нульова. Можливо, початкове наближення потребує зміни.")
            break

        delta_alpha = -p[-2] / q[-2]
        delta_beta = -p[-1] / q[-2]
        alpha += delta_alpha
        beta += delta_beta

        k += 1

        if abs(delta_alpha) < eps and abs(delta_beta) < eps:
            break

    root1 = alpha + math.sqrt(beta) * 1j
    root2 = alpha - math.sqrt(beta) * 1j
    return root1, root2, k

coefficients = [1, -6, 11, -6]
alpha0 = 3  
beta0 = 1  
eps = 1e-10
root1, root2, iterations = lin_method(coefficients, alpha0, beta0, eps)
print(f"Знайдені комплексні корені: {root1}, {root2}, число ітерацій: {iterations}")
