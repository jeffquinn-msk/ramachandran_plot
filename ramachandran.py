import os.path

from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.vectors import calc_dihedral
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import argparse


def get_args():
    parser = argparse.ArgumentParser(
        description="Plot ramachandran plot of a PDB file."
    )
    parser.add_argument(
        "--pdb",
        required=True,
        help="Path to the PDB file.",
    )
    parser.add_argument("--output", required=True, help="Path to the output file.")

    return parser.parse_args()


def plot_ramachandran(pdb_file, output_fn):
    parser = PDBParser()
    structure = parser.get_structure('id', pdb_file)
    angles = []
    indices = []
    for model in structure:
        for chain in model:
            polypeptides = PPBuilder().build_peptides(chain)
            for poly_index, poly in enumerate(polypeptides):
                phi_psi = poly.get_phi_psi_list()
                for res_index, residue in enumerate(poly):
                    phi, psi = phi_psi[res_index]
                    angles.append([phi, psi])
                    indices.append(poly_index)
    angles = pd.DataFrame(np.rad2deg(np.array(angles).astype(float)), columns=["phi", "psi"])
    angles["Polypeptide Id"] = indices

    angles = angles.dropna()

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.set_title(os.path.basename(pdb_file).replace(".pdb", ""))
    ax.set_xlabel("φ")
    ax.set_ylabel("ψ")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)

    sns.kdeplot(
        data=angles,
        x="phi",
        y="psi",
        cmap="Spectral_r",
        fill=True,
        ax=ax
    )

    sns.scatterplot(data=angles,
                    x="phi",
                    y="psi",
                    marker="x",
                    color="black",
                    ax=ax)

    plt.tight_layout()
    plt.savefig(output_fn)


def main():
    args = get_args()
    plot_ramachandran(args.pdb, args.output)


if __name__ == "__main__":
    main()
