import math
import matplotlib.pyplot as plt

def step(xx, nn):
    return xx ** nn if nn != 0 else 1

 
def F(X):
    return step((1 - X[0]), 2) + 100.0 * step((X[1] - X[0] ** 2), 2)

def investigating_search_full(X0, deltaX, eps, q):
    X1 = X0[:]
    for i in range(len(X0)):
        while True:
            X1[i] = X0[i] + deltaX[i]
            if F(X1) < F(X0):
                break
            X1[i] = X0[i] - deltaX[i]
            if F(X1) < F(X0):
                break
            deltaX[i] /= q
            X1[i] = X0[i]
    return X1

def investigating_search_simple(X0, deltaX):
    X1 = X0[:]
    for i in range(len(X0)):
        X1[i] = X0[i] + deltaX[i]
        if F(X1) >= F(X0):
            X1[i] = X0[i] - deltaX[i]
            if F(X1) >= F(X0):
                X1[i] = X0[i]
    return X1

def difference(X, Y, eps):
    return max(abs(Y[i] - X[i]) for i in range(len(X))) < eps

def sample_search(X0, X1, p):
    return [X0[i] + p * (X1[i] - X0[i]) for i in range(len(X0))]

def plot_rosenbrock():
    x = []
    y = []
    z = []
    
    for i in range(-20, 21):
        for j in range(-20, 21):
            X = [i * 0.1, j * 0.1]
            x.append(X[0])
            y.append(X[1])
            z.append(F(X))
    
    X_unique = sorted(set(x))
    Y_unique = sorted(set(y))
    Z = [[0] * len(Y_unique) for _ in range(len(X_unique))]

    for i in range(len(x)):
        X_index = X_unique.index(x[i])
        Y_index = Y_unique.index(y[i])
        Z[X_index][Y_index] = z[i]

    plt.figure()
    plt.contourf(X_unique, Y_unique, Z, levels=50, cmap='viridis')
    plt.colorbar()
    plt.title("Графік функції Розенброка")
    plt.xlabel("X1")
    plt.ylabel("X2")
    plt.plot(1, 1, 'ro')  
    plt.show()

def plot_system_functions():
    x_vals = [i * 0.1 for i in range(-20, 21)]
    y_vals_f1 = [1 - x for x in x_vals] 
    y_vals_f2 = [(x ** 2) for x in x_vals] 

    plt.figure()

    plt.plot(x_vals, y_vals_f1, label="f1 = 1 - x", color="blue")

    plt.plot(x_vals, y_vals_f2, label="f2 = x^2", color="green")

    plt.title("Графіки системи рівнянь")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    output_file = "output.txt"
    with open(output_file, "w", encoding="utf-8") as output:
        k = 0
        kmax = 1000
        eps1 = 1e-10
        eps2 = 1e-10
        q = 2
        p = 2
        n = 1
        deltaX0 = [0.01] * (n + 1)
        X0 = [0.0] * (n + 1)
        X1 = X0[:]

        while k < kmax:
            k += 1
            output.write(f"\nk{k}\n\n")
            X0 = X1[:]
            deltaX = deltaX0[:]

            X1 = investigating_search_full(X0, deltaX, eps1, q)
            output.write("ISF\n")
            if difference(X0, X1, eps2):
                output.write("D1\n")
                break
            
            X2 = X1[:]
            while True:
                X0 = X1[:]
                X1 = X2[:]
                X2p = sample_search(X0, X1, p)
                output.write("SS\n")
                X2 = investigating_search_simple(X2p, deltaX0)
                output.write("ISS\n")
                
                if difference(X2p, X2, eps2):
                    output.write("D2\n")
                    break
                
                if F(X2) >= F(X1):
                    break

        if k >= kmax:
            output.write("Max Iterations\n")
        else:
            output.write(f"Number of iterations k={k}\n")
        
        output.write("Знайдений мінімум:\n")
        for value in X1:
            output.write(f"{value:e}\n")
        print("Знайдений мінімум:", X1)
        print("Кількість кроків:", k)
    plot_rosenbrock()
    plot_system_functions()

if __name__ == "__main__":
    main()

