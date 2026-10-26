import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

import dataclasses


FILE_PATH = r"./build/result.txt" 

def check_inside_region(x: float, y: float) -> int:
    if (np.abs(x) >= 2) or (np.abs(y) >= 2) or (x*x + (y-2)*(y-2) < 1):
        return -1
    if np.abs(x*x + (y-2)*(y-2) - 1) < 1e-6:
        return 0
    return 1

@dataclasses.dataclass
class SolutionData:
    M : int
    N : int
    x_min : float
    x_max : float
    y_min : float 
    y_max : float
    data : np.ndarray

def read_data(path: str) -> SolutionData:
    with open(path) as data_file:
        m, n, x_min, x_max, y_min, y_max = map(float, data_file.readline().split())
        M, N = int(m), int(n)
        data = np.empty((M+1, N+1),  dtype=float)
        i = 0
        while line := data_file.readline():
            data[i, :] = list(map(float, line.split()))
            i += 1

    return SolutionData(M, N, x_min, x_max, y_min, y_max, data)

def get_masked_solution_data(solution: SolutionData) -> np.ndarray:
    h_x = (solution.x_max - solution.x_min) / solution.M
    h_y = (solution.y_max - solution.y_min) / solution.N

    data_filter = np.array([[ bool(check_inside_region(solution.x_min + i*h_x, solution.y_min + j*h_y) >= 0) for j in range(solution.N+1)] for i in range(solution.M+1)])
    return np.where(data_filter, solution.data, 0)

def plot_solution(solution: SolutionData) -> None:
    _, (ax1, ax2) = plt.subplots(1, 2)
    ax1 : Axes
    ax2 : Axes

    ax1.imshow(solution.data.T, origin='lower', cmap="RdYlBu", extent=(solution.x_min, solution.x_max, solution.y_min, solution.y_max))
    ax1.set_aspect('equal')
    ax1.set_title("Unmasked")

    ax2.imshow(get_masked_solution_data(solution).T, origin='lower', cmap="RdYlBu", extent=(solution.x_min, solution.x_max, solution.y_min, solution.y_max))
    ax2.set_aspect('equal')
    ax2.set_title("Masked")

    plt.show()

if __name__ == "__main__":
    plot_solution(read_data(FILE_PATH))