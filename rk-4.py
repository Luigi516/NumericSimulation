import math
import matplotlib.pyplot as plt

x0 = 0
y0 = 1
AscissaFinale = 5
h = 0.5
h_prex = 0.01

def dydx(x, y):
    return math.cos(x) - y

def RungeKutt4(x_corrente, y_corrente, h, dydx):
    k1 = dydx(x_corrente, y_corrente)
    k2 = dydx(x_corrente + h/2, y_corrente + k1*h/2)
    k3 = dydx(x_corrente + h/2, y_corrente + k2*h/2)
    k4 = dydx(x_corrente + h, y_corrente + k3*h)
    return y_corrente + h/6*(k1 + 2*k2 + 2*k3 + k4)

def y_exact(x0, y_0, x):
    term = 0.5 * (math.cos(x0) + math.sin(x0))
    c = (y_0 - term) / math.exp(-x0)
    return 0.5 * (math.cos(x) + math.sin(x)) + c * math.exp(-x)


# Computo la curva corretta
x = x0

ListaOrdinateEffettive = []
ListaAscisseEffettive = []

while x < AscissaFinale:
    ListaAscisseEffettive.append(x)
    ListaOrdinateEffettive.append(y_exact(x0, y0, x))
    x += h_prex


# Computo la curva approssimata
ListaAscisseAppx = [x0]
ListaOrdinateAppx = [y0]
y = y0
x = x0
N = (int)(round((AscissaFinale - x) / h))
for i in range(N): 
    y = RungeKutt4(x, y, h, dydx)
    x += h
    ListaAscisseAppx.append(x)
    ListaOrdinateAppx.append(y)


plt.grid(True)
plt.plot(ListaAscisseAppx, ListaOrdinateAppx, color='red', marker='o', label='Simulato-via-RK4')
plt.plot(ListaAscisseEffettive, ListaOrdinateEffettive, color='blue', label='Reale')
plt.legend()
plt.show()


