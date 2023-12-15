import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker


class Cluster:
    def __init__(self, name: str, x: float, y: float, r: float):
        self.name = name
        self.x = x
        self.y = y
        self.r = r


class Point:
    def __init__(self, x: float, y: float, nseq: str, cseq: str, amino_acid: str):
        self.x = x
        self.y = y
        self.nseq = nseq
        self.cseq = cseq
        self.amino_acid = amino_acid


def read_point(filename: str) -> list[Point]:
    points: list[Point] = []
    with open(filename) as fs:
        hs = fs.readline().rstrip().split('\t')
        x_index = hs.index("x")
        y_index = hs.index("y")
        nseq_index = hs.index("nseq")
        cseq_index = hs.index("cseq")
        amino_acid_index = hs.index("amino_acid")
        for line in fs:
            spl = line.rstrip().split('\t')
            if '_' in spl[nseq_index] or '_' in spl[cseq_index]:
                continue
            points.append(
                Point(
                    float(spl[x_index]), float(spl[y_index]), spl[nseq_index], spl[cseq_index], spl[amino_acid_index]
                )
            )

    return points


def distance(a: (float, float), b: (float, float)):
    return math.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]))


def plot_logo(
        points: list[Point],
        clusters: list[Cluster],
        amino_acid_frequencies:
        dict[str, float],
        plot_type: str = "probability"
) -> None:
    for c in clusters:
        ps = [p for p in points if distance((p.x, p.y), (c.x, c.y)) <= c.r]
        ss: list[str] = [p.nseq + p.amino_acid + p.cseq for p in ps]
        n: int = len(ps[0].nseq)
        if plot_type == "information":
            out_df = logomaker.alignment_to_matrix(ss,
                                                   to_type="information",
                                                   background=[amino_acid_frequencies[aa] for aa in
                                                               "ACDEFGHIKLMNPQRSTVWY"])
            ss_logo = logomaker.Logo(out_df,
                                     width=.8,
                                     vpad=.05,
                                     fade_probabilities=True,
                                     stack_order='small_on_top')
        else:
            out_df = logomaker.alignment_to_matrix(ss, to_type="probability")
            ss_logo = logomaker.Logo(out_df,
                                     width=.8,
                                     vpad=.05,
                                     fade_probabilities=True,
                                     stack_order='small_on_top')
        ss_logo.style_spines(spines=['left', 'right'], visible=False)
        ss_logo.ax.set_xticks(range(len(out_df)))
        ss_logo.ax.set_title(f"{c.name}_{len(ss)}")
        ss_logo.ax.set_xticklabels('%+d' % x for x in range(-n, n + 1))
        ss_logo.ax.set_yticks([0, .5, 1])
        plt.savefig(f"all_wCentralAA_5_all_tsne\\{c.name}.pdf")
        # ss_logo.ax.axvline(2.5, color='k', linewidth=1, linestyle=':')


def process(points_filename: str, fasta_filename: str) -> None:
    points: list[Point] = read_point(points_filename)
    amino_acid_frequencies: dict[str, float] = get_amino_acid_frequencies(fasta_filename)
    clusters: list[Cluster] = [
        Cluster("c1", -2.0084, 2.6155, 1.0),
        Cluster("c2", -1.0084, 3.3747, 1.0),
        Cluster("c3", -0.0084, 3.3747, 1.0),
        Cluster("c4", 1.2481, 2.7285, 1.0),
        Cluster("c5", -1.3214, 0.7961, 1.0),
        Cluster("c6", 0.0087, 1.4348, 1.0),
        Cluster("c7", -0.4745, -0.2278, 1.0),
        Cluster("c8", 1.6727, 0.0, 1.0),
        Cluster("c9", -2.2023, -2.1306, 1.0),
        Cluster("c10", -0.0084, -3.3191, 1.0),
        Cluster("c11", 1.6727, -3.0139, 1.0),
        Cluster("c12", 0.7481, -4.3191, 1.0),
    ]
    plot_logo(points, clusters, amino_acid_frequencies, "probability")


def get_amino_acid_frequencies(fasta_filename: str) -> dict[str, float]:
    amino_acid: list[str] = list("ACDEFGHIKLMNPQRSTVWY")
    amino_acid_count = {aa: 0 for aa in amino_acid}
    with open(fasta_filename) as fs:
        for line in fs:
            if line.startswith('>'):
                continue
            for aa in line.rstrip():
                if aa not in amino_acid_count:
                    continue
                amino_acid_count[aa] += 1
    n: int = sum([value for aa, value in amino_acid_count.items()])
    return {aa: (value / n) for aa, value in amino_acid_count.items()}


if __name__ == "__main__":
    process(
        "all_wCentralAA_5_all_tsne\\file.out.txt",
        "Mouse_UP000000589_wIsoform.fasta"
    )
