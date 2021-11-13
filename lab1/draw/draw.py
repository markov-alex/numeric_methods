import sys
import json
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def prepare_data(scheme_data):
    left_bound = scheme_data["leftBound"]
    right_bound = scheme_data["rightBound"]
    tau = scheme_data["tau"]
    h = scheme_data["h"]
    T = scheme_data["T"]
    return left_bound, right_bound, tau, h, T


def draw3d(explicit_scheme_data, implicit_scheme_data, crank_nicolson_data, analytical_data):
    left_bound, right_bound, tau, h, T = prepare_data(explicit_scheme_data)
    u_explicit = explicit_scheme_data["u"]
    u_implicit = implicit_scheme_data["u"]
    u_crank_nicolson = crank_nicolson_data["u"]
    u_analytical = analytical_data["u"]
    fig = plt.figure(figsize=plt.figaspect(0.3))
    x = np.arange(left_bound, right_bound + h, h)
    t = np.arange(0, T + tau, tau)
    x, t = np.meshgrid(x, t)
    z_explicit = np.array(u_explicit)
    z_implicit = np.array(u_implicit)
    z_crank_nicolson = np.array(u_crank_nicolson)
    z_analytical = np.array(u_analytical)

    ax = fig.add_subplot(2, 2, 1, projection='3d')
    plt.title('Явная схема')
    ax.set_xlabel('x', fontsize=20)
    ax.set_ylabel('t', fontsize=20)
    ax.set_zlabel('u', fontsize=20)
    ax.plot_surface(x, t, z_explicit, cmap=cm.coolwarm,
                    linewidth=0, antialiased=True)

    ax = fig.add_subplot(2, 2, 2, projection='3d')
    plt.title('Неявная схема')
    ax.set_xlabel('x', fontsize=20)
    ax.set_ylabel('t', fontsize=20)
    ax.set_zlabel('u', fontsize=20)
    ax.plot_surface(x, t, z_implicit, cmap=cm.coolwarm,
                    linewidth=0, antialiased=True)

    ax = fig.add_subplot(2, 2, 3, projection='3d')
    plt.title('Схема Кранка-Николсона')
    ax.set_xlabel('x', fontsize=20)
    ax.set_ylabel('t', fontsize=20)
    ax.set_zlabel('u', fontsize=20)
    ax.plot_surface(x, t, z_crank_nicolson, cmap=cm.coolwarm,
                    linewidth=0, antialiased=True)

    ax = fig.add_subplot(2, 2, 4, projection='3d')
    plt.title('Аналитическое решение')
    ax.set_xlabel('x', fontsize=20)
    ax.set_ylabel('t', fontsize=20)
    ax.set_zlabel('u', fontsize=20)
    ax.plot_surface(x, t, z_analytical, cmap=cm.coolwarm,
                    linewidth=0, antialiased=True)

    surf = ax.plot_surface(x, t, z_analytical, cmap=cm.coolwarm,
                           linewidth=0, antialiased=True)

    fig.colorbar(surf, shrink=0.5, aspect=15)
    plt.show()


def draw_error(test_scheme_data, analytical_data, name_test_scheme):
    left_bound, right_bound, tau, h, T = prepare_data(test_scheme_data)
    u_test_scheme = test_scheme_data["u"]
    u_analytical = analytical_data["u"]
    t = np.arange(0, T + tau, tau)
    x = np.arange(left_bound, right_bound + h, h)
    t_idx = np.linspace(0, t.shape[0] - 1, 6, dtype=np.int32)
    error = np.zeros(len(t))
    for i in range(len(t)):
        error[i] = np.max(np.abs(np.array(u_test_scheme[i]) - np.array(u_analytical[i])))
    plt.figure(figsize=(12, 7))
    plt.plot(t, error, label='Ошибка')
    plt.title(f'График изменения ошибки для схемы - {name_test_scheme}')
    plt.xlabel('t')
    plt.ylabel('Ошибка')
    plt.grid(True)
    plt.show()


def draw_comparison_solutions(explicit_scheme_data, implicit_scheme_data, crank_nicolson_data, analytical_data):
    left_bound, right_bound, tau, h, T = prepare_data(explicit_scheme_data)
    u_explicit = explicit_scheme_data["u"]
    u_implicit = implicit_scheme_data["u"]
    u_crank_nicolson = crank_nicolson_data["u"]
    u_analytical = analytical_data["u"]
    t = np.arange(0, T + tau, tau)
    x = np.arange(left_bound, right_bound + h, h)
    t_idx = np.linspace(0, t.shape[0] - 1, 6, dtype=np.int32)
    fig, ax = plt.subplots(3, 2)
    fig.suptitle('Сравнение решений')
    fig.set_figheight(15)
    fig.set_figwidth(16)
    k = 0
    for i in range(3):
        for j in range(2):
            time_idx = t_idx[k]
            ax[i][j].plot(x, u_explicit[time_idx], label='Явная схема')
            ax[i][j].plot(x, u_implicit[time_idx], label='Неявная схема')
            ax[i][j].plot(x, u_crank_nicolson[time_idx], label='Кранка-Николсона')
            ax[i][j].plot(x, u_analytical[time_idx], label='Аналитическое решение')
            ax[i][j].grid(True)
            ax[i][j].set_xlabel('x')
            ax[i][j].set_ylabel('u')
            ax[i][j].set_title(f'Решения при t = {t[time_idx]}')
            k += 1
    plt.legend(bbox_to_anchor=(1.05, 2), loc='upper left', borderaxespad=0.)
    plt.show()


with open("../results/explicitScheme.json") as f1, open("../results/implicitScheme.json") as f2, \
        open("../results/crankNicolson.json") as f3, open("../results/analyticalSolution.json") as f4:
    data1 = json.load(f1)
    data2 = json.load(f2)
    data3 = json.load(f3)
    data4 = json.load(f4)
    draw3d(data1, data2, data3, data4)
    draw_error(data1, data4, "Явная")
    draw_error(data2, data4, "Неявная")
    draw_error(data3, data4, "Кранка-Николсона")
    draw_comparison_solutions(data1, data2, data3, data4)
