import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def read_data(data_path):
    df = pd.read_csv(data_path)

    print(f"Reading {data_path}")

    u = np.array(df['u'])
    v = np.array(df["v"])

    u_min = u.min()
    u_max = u.max()

    v_min = v.min()
    v_max = v.max()

    C = np.zeros((u_max - u_min, v_max - v_min))

    for u_index in range(u_min, u_max):
        for v_index in range(v_min, v_max):
            data = float(df[(df["u"] == u_index) & (df["v"] == v_index)]["data"])

            if data > 0.1:
                C[u_index, v_index] = data

    return C


def main():
    data_paths = ["./debug/slice_data.csv", "./debug/kernel1.csv", "./debug/kernel2.csv", "./debug/kerneldifference.csv"]
    data_labels = ["Original Slice", "Blur σ=3", "Blur σ=5", "Difference of Gaussian"]
    fig, axs = plt.subplots(2,2)

    for ax, path, label in zip(axs.flat, data_paths, data_labels):
        C = read_data(path)
        # ax.set_title(path.split('/')[-1].split(".")[0])
        ax.set_title(label)
        ax.pcolor(C)

    plt.tight_layout()
    # plt.show()
    plt.savefig('./figure/slice_dog_visualisation.png')

if __name__ == "__main__":
    main()
