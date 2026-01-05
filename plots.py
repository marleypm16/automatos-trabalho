# plots.py
import numpy as np
import matplotlib.pyplot as plt

def plot_mass_history(filename="history_mass.npy"):
    h = np.load(filename)
    plt.figure(figsize=(8,4))
    plt.plot(h)
    plt.xlabel("Passos de tempo")
    plt.ylabel("Massa total (unidades)")
    plt.title("Evolução da massa total de poluente")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("mass_history.png", dpi=150)
    plt.show()

def plot_concentration_snapshot(filename="final_concentration.npy"):
    C = np.load(filename)
    plt.figure(figsize=(8,5))
    im = plt.imshow(C.T, origin='lower', aspect='auto')  # transpor para orientação visual
    plt.colorbar(im, label="Concentração")
    plt.title("Concentração final do poluente")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.tight_layout()
    plt.savefig("concentration_final.png", dpi=150)
    plt.show()

if __name__ == "__main__":
    plot_mass_history()
    plot_concentration_snapshot()
