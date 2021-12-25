import sys
import json
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def prepare_data(scheme_data):
    lx = scheme_data["lx"]
    rx = scheme_data["rx"]
    ly = scheme_data["ly"]
    ry = scheme_data["ry"]
    T = scheme_data["T"]
    h1 = scheme_data["h1"]
    h2 = scheme_data["h2"]
    tau = scheme_data["tau"]
    return lx, rx, ly, ry, T, h1, h2, tau


def draw_error(test_scheme_data, analytical_data, name_test_scheme):
    lx, rx, ly, ry, T, h1, h2, tau = prepare_data(test_scheme_data)
    u_test_scheme = test_scheme_data["u"]
    u_analytical = analytical_data["u"]
    t = np.arange(0, T + tau, tau)
    y = np.arange(ly, ry + h2, h2)
    x = np.arange(lx, rx + h1, h1)
    error = np.zeros(len(t))
    z_test_scheme = np.array(u_test_scheme)
    z_analytical = np.array(u_analytical)
    for i in range(len(t)):
        error[i] = np.max(np.abs(z_test_scheme[i] - z_analytical[i]))
    plt.figure(figsize=(12, 7))
    plt.plot(y, error, label='Ошибка')
    plt.title(f'График изменения ошибки для метода - {name_test_scheme}')
    plt.xlabel('t')
    plt.ylabel('Ошибка')
    plt.grid(True)
    plt.show()


def draw_comparison_solutions(data1, data2, data3):
    lx, rx, ly, ry, T, h1, h2, tau = prepare_data(data1)
    u1 = data1["u"]
    u2 = data2["u"]
    u3 = data3["u"]
    x = np.arange(lx, rx + h1, h1)
    y = np.arange(ly, ry + h2, h2)
    t = np.arange(0, T + tau, tau)
    z1 = np.array(u1)
    z2 = np.array(u2)
    z3 = np.array(u3)
    x_idx = np.linspace(0, x.shape[0] - 1, 7, dtype=np.int32)
    y_idx = np.linspace(0, y.shape[0] - 1, 7, dtype=np.int32)
    fig, ax = plt.subplots(3, 2)
    fig.suptitle('Сравнение решений')
    fig.set_figheight(15)
    fig.set_figwidth(16)
    k = 0
    for i in range(3):
        for j in range(2):
            xx_idx = x_idx[k]
            yy_idx = y_idx[k]
            ax[i][j].plot(t, [u1[q][xx_idx][yy_idx] for q in range(len(t))], label='Схема переменных направлений')
            ax[i][j].plot(t, [u2[q][xx_idx][yy_idx] for q in range(len(t))], label='Схема дробных шагов')
            ax[i][j].plot(t, [u3[q][xx_idx][yy_idx] for q in range(len(t))], label='Аналитическое решение')
            ax[i][j].grid(True)
            ax[i][j].set_xlabel('t')
            ax[i][j].set_ylabel('u')
            ax[i][j].set_title(f'Решение при x = {x[xx_idx]}; y = {y[yy_idx]}')
            k += 1
    plt.legend()
    plt.show()


with open("../results/alternatingDirections.json") as f1, open("../results/fractionalSteps.json") as f2, \
        open("../results/analyticalSolution.json") as f3:
    data1 = json.load(f1)
    data2 = json.load(f2)
    data3 = json.load(f3)
    draw_error(data1, data3, "Схема переменных направлений")
    draw_error(data2, data3, "Схема дробных шагов")
    draw_comparison_solutions(data1, data2, data3)