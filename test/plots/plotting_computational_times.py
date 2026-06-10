import os
import matplotlib.pyplot as plt
import numpy as np


def plot_computational_cost(N, times, title, output_dir):
    N_logN = N * np.log(N)
    N_2 = N**2
    N_3 = N**3
    N_4 = N**4

    plt.figure(figsize=(7, 5))

    plt.loglog(N, N / N[0], color='yellow', label='N')
    plt.loglog(N, N_logN / N_logN[0], color='green', label='N log N')
    plt.loglog(N, N_2 / N_2[0], color='blue', label='N²')
    plt.loglog(N, N_3 / N_3[0], color='orange', label='N³')
    plt.loglog(N, N_4 / N_4[0], color='red', label='N⁴')

    plt.loglog(N, times / times[0], color='black', marker='o', label='True times')

    plt.xlabel("Number of points")
    plt.ylabel("Normalized computational time")
    plt.title(title)
    plt.grid(True, which="both")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir, dpi=300)
    plt.close()

def main():
	COMPUTATIONAL_COSTS_PATH = "test/data/computational_times"
	N = np.array([16, 32, 64, 128])
	N = N**2

	unclipped_times = np.loadtxt(os.path.join(COMPUTATIONAL_COSTS_PATH, "voronoi.txt"), dtype=np.int64)

	plot_computational_cost(N, unclipped_times, "Computational cost - Voronoi Unclipped", "test/plots/voronoi_unclipped.png")

	clipped_times = np.loadtxt(os.path.join(COMPUTATIONAL_COSTS_PATH, "voronoi_clipped.txt"), dtype=np.int64)

	plot_computational_cost(N, clipped_times, "Computational cost - Voronoi Clipped", "test/plots/voronoi_clipped.png")
 
 
if __name__ == "__main__":
    main()