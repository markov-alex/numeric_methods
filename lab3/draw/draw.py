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
    h1 = scheme_data["h1"]
    h2 = scheme_data["h2"]
    return lx, rx, ly, ry, h1, h2


def draw3d(liebmann_data, seidel_data, relaxation_data, analytical_data):
    lx, rx, ly, ry, h1, h2 = prepare_data(liebmann_data)
    u_liebmann = liebmann_data["u"]
    u_seidel = seidel_data["u"]
    u_relaxation = relaxation_data["u"]
    u_analytical = analytical_data["u"]
    fig = plt.figure(figsize=plt.figaspect(0.3))
    x = np.arange(lx, rx + h1, h1)
    y = np.arange(ly, ry + h2, h2)
    x, y = np.meshgrid(x, y)
    z_liebmann = np.transpose(np.array(u_liebmann))
    z_seidel = np.transpose(np.array(u_seidel))
    z_relaxation = np.transpose(np.array(u_relaxation))
    z_analytical = np.transpose(np.array(u_analytical))

    ax = fig.add_subplot(2, 2, 1, projection='3d')
    plt.title('Метод Либмана')
    ax.set_xlabel('x', fontsize=20)
    ax.set_ylabel('y', fontsize=20)
    ax.set_zlabel('u', fontsize=20)
    ax.plot_surface(x, y, z_liebmann, cmap=cm.coolwarm,
                    linewidth=0, antialiased=True)

    ax = fig.add_subplot(2, 2, 2, projection='3d')
    plt.title('Метод Зейделя')
    ax.set_xlabel('x', fontsize=20)
    ax.set_ylabel('y', fontsize=20)
    ax.set_zlabel('u', fontsize=20)
    ax.plot_surface(x, y, z_seidel, cmap=cm.coolwarm,
                    linewidth=0, antialiased=True)

    ax = fig.add_subplot(2, 2, 3, projection='3d')
    plt.title('Верхняя релаксация')
    ax.set_xlabel('x', fontsize=20)
    ax.set_ylabel('y', fontsize=20)
    ax.set_zlabel('u', fontsize=20)
    ax.plot_surface(x, y, z_relaxation, cmap=cm.coolwarm,
                    linewidth=0, antialiased=True)

    ax = fig.add_subplot(2, 2, 4, projection='3d')
    plt.title('Аналитическое решение')
    ax.set_xlabel('x', fontsize=20)
    ax.set_ylabel('y', fontsize=20)
    ax.set_zlabel('u', fontsize=20)
    ax.plot_surface(x, y, z_analytical, cmap=cm.coolwarm,
                    linewidth=0, antialiased=True)

    surf = ax.plot_surface(x, y, z_analytical, cmap=cm.coolwarm,
                               linewidth=0, antialiased=True)

    fig.colorbar(surf, shrink=0.5, aspect=15)
    plt.show()


def draw_error(test_scheme_data, analytical_data, name_test_scheme):
    lx, rx, ly, ry, h1, h2 = prepare_data(test_scheme_data)
    u_test_scheme = test_scheme_data["u"]
    u_analytical = analytical_data["u"]
    y = np.arange(ly, ry + h2, h2)
    x = np.arange(lx, rx + h1, h1)
    y_idx = np.linspace(0, y.shape[0] - 1, 6, dtype=np.int32)
    error = np.zeros(len(y))
    z_test_scheme = np.transpose(np.array(u_test_scheme))
    z_analytical = np.transpose(np.array(u_analytical))
    for i in range(len(y)):
        error[i] = np.max(np.abs(np.array(z_test_scheme[i]) - np.array(z_analytical[i])))
    plt.figure(figsize=(12, 7))
    plt.plot(y, error, label='Ошибка')
    plt.title(f'График изменения ошибки для метода - {name_test_scheme}')
    plt.xlabel('y')
    plt.ylabel('Ошибка')
    plt.grid(True)
    plt.show()


def draw_comparison_solutions(liebmann_data, seidel_data, relaxation_data, analytical_data):
    lx, rx, ly, ry, h1, h2 = prepare_data(liebmann_data)
    u_liebmann = liebmann_data["u"]
    u_seidel = seidel_data["u"]
    u_relaxation = relaxation_data["u"]
    u_analytical = analytical_data["u"]
    z_liebmann = np.transpose(np.array(u_liebmann))
    z_seidel = np.transpose(np.array(u_seidel))
    z_relaxation = np.transpose(np.array(u_relaxation))
    z_analytical = np.transpose(np.array(u_analytical))
    x = np.arange(lx, rx + h1, h1)
    y = np.arange(ly, ry + h2, h2)
    y_idx = np.linspace(0, y.shape[0] - 1, 6, dtype=np.int32)
    fig, ax = plt.subplots(3, 2)
    fig.suptitle('Сравнение решений')
    fig.set_figheight(15)
    fig.set_figwidth(16)
    k = 0
    for i in range(3):
        for j in range(2):
            yy_idx = y_idx[k]
            ax[i][j].plot(x, z_liebmann[yy_idx], label='Метод Либмана')
            ax[i][j].plot(x, z_seidel[yy_idx], label='Метод Зейделя')
            ax[i][j].plot(x, z_relaxation[yy_idx], label='Верхняя релаксация')
            ax[i][j].plot(x, z_analytical[yy_idx], label='Аналитическое решение')
            ax[i][j].grid(True)
            ax[i][j].set_xlabel('x')
            ax[i][j].set_ylabel('u')
            ax[i][j].set_title(f'Решения при y = {y[yy_idx]}')
            k += 1
    plt.legend(bbox_to_anchor=(1.05, 2), loc='upper left', borderaxespad=0.)
    plt.show()


with open("../results/liebmann.json") as f1, open("../results/seidel.json") as f2, \
        open("../results/relaxation.json") as f3, open("../results/analyticalSolution.json") as f4:
    data1 = json.load(f1)
    data2 = json.load(f2)
    data3 = json.load(f3)
    data4 = json.load(f4)
    draw3d(data1, data2, data3, data4)
    draw_error(data1, data4, "Метод Либмана")
    draw_error(data2, data4, "Метод Зейделя")
    draw_error(data3, data4, "Верхняя релаксация")
    draw_comparison_solutions(data1, data2, data3, data4)