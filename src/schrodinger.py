# ------- импорт библиотек -------
import numpy as np
import json
from scipy.special import sph_harm, factorial
from math import sqrt, exp
import scipy.special as sp
import os

a0 = 1.0

# Радиальная часть R_{n,l}(r) для водорода (нормированная)
def R(n, l, r):
    rho = 2.0 * r / (n * a0)
    norm = sqrt((2.0 / (n * a0))**3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    L = sp.assoc_laguerre(rho, n - l - 1, 2 * l + 1)
    return norm * np.exp(-rho / 2.0) * rho**l * L

# Полная комплексная орбиталь ψ(r,θ,φ)
def psi(n, l, m, r, theta, phi):
    Y = sph_harm(m, l, phi, theta)  # комплексная сферическая гармоника
    return R(n, l, r) * Y

# Плотность вероятности = |ψ|² для выборки с отклонением
def probability(n, l, m, r, theta, phi):
    return abs(psi(n, l, m, r, theta, phi))**2

# Выборка точек в пространстве
def sample_points(n, l, m, N=50000):
    points = []
    max_prob = None

    for _ in range(N):
        r = np.random.exponential(scale=n**2)
        theta = np.random.uniform(0, np.pi)
        phi = np.random.uniform(0, 2 * np.pi)

        p = probability(n, l, m, r, theta, phi)
        if max_prob is None or p > max_prob:
            max_prob = p

        if np.random.rand() < p / max_prob:
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)

            psi_val = psi(n, l, m, r, theta, phi)
            points.append({
                "x": x,
                "y": y,
                "z": z,
                "psi_re": psi_val.real,
                "psi_im": psi_val.imag
            })
    return points

def main():
    print("\n=== Генератор водородных орбиталей (комплексный) ===\n")
    n = int(input("Введите n (1..7): "))
    l = int(input(f"Введите l (0..{n-1}): "))
    m = int(input(f"Введите m (-{l}..{l}): "))
    N = int(input("Сколько образцов частиц создать? (пример: 50000): "))

    print("\nСоздаются образцы... подождите...")
    pts = sample_points(n, l, m, N)

    os.makedirs("orbitals", exist_ok=True)
    filename = f"orbitals/orbital_n{n}_l{l}_m{m}_complex.json"
    with open(filename, "w") as f:
        json.dump({
            "n": n,
            "l": l,
            "m": m,
            "points": pts
        }, f)

    print(f"\nСохранено {len(pts)} образцов в файл: {filename}")
    print("\nЗагрузите этот JSON в C++-просмотрщик для визуализации.\n")

if __name__ == "__main__":
    main()
