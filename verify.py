"""
Compare some snapshots from preprocessing EC-Earth with
the output of the original code.
"""

from scipy.io import loadmat
import matplotlib.pyplot as plt

original_data = loadmat('../original_output_data/2002-01-01fluxes_storages.mat')
new_data = loadmat('../output_data/2002-01-01fluxes_storages.mat')

t = 0
for variable in ["Fa_E_top", "Fa_N_top", "Fa_E_down", "Fa_N_down", "E", "P", "W_top", "W_down", "Fa_Vert"]:

    original = original_data[variable][t, ::-1, :]
    new = new_data[variable][t, ::-1, :]
    diff = new - original

    fig, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(10, 4))
    pcm1 = ax1.pcolormesh(original)
    pcm2 = ax2.pcolormesh(new)
    pcm3 = ax3.pcolormesh(diff)

    ax1.set_title("original")
    ax2.set_title("new")
    ax3.set_title("diff")

    fig.colorbar(pcm1, ax=ax1)
    fig.colorbar(pcm2, ax=ax2)
    fig.colorbar(pcm3, ax=ax3)

    fig.suptitle(variable)
    fig.savefig(f'validation/{variable}.png')
