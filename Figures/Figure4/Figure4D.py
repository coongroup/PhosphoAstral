import gzip
import os
import tarfile

import math
import os
import matplotlib.pyplot as plt

from enum import Enum

max_length = 1400

pLDDT_row = "ATOM"
pLDDT_atom = "CA "
pLDDT_atom_column = 2
pLDDT_amino_acid_column = 3
pLDDT_position_column = 5
pLDDT_amino_acids = {"SER", "THR", "TYR"}
pLDDT_value_column = 10


class Occurrence(Enum):
    ONE = 1
    REST = 2
    ALL = 3


class PhosphoSite:
    def __init__(self, protein: str, position: int, amino_acid: str, occurrence: Occurrence):
        self.protein = protein
        self.position = position
        self.amino_acid = amino_acid
        self.amino_acid = amino_acid
        self.occurrence = occurrence


def read_phospho_sites(filename: str) -> list[PhosphoSite]:
    phospho_sites: list[PhosphoSite] = []
    with open(filename) as fs:
        header = fs.readline().rstrip().split('\t')
        major_protein_index = header.index("Major Protein")
        major_position_index = header.index("Major Position")
        amino_acid_index = header.index("Amino Acid")
        expression_index = header.index("Flanking Region") + 1
        for line in fs:
            spl = line.rstrip().split('\t')
            position = int(spl[major_position_index]) - 1
            if position < max_length:
                occurrence = Occurrence.REST
                cnt = 0
                for i in range(expression_index, len(spl)):
                    if float(spl[i]) > 0:
                        cnt += 1
                if cnt == 1:
                    occurrence = Occurrence.ONE
                elif cnt == (len(spl) - expression_index):
                    occurrence = Occurrence.ALL
                phospho_sites.append(
                    PhosphoSite(
                        spl[major_protein_index],
                        position,
                        spl[amino_acid_index],
                        occurrence
                    )
                )
    return phospho_sites


def get_suitable_proteins(sites: list[PhosphoSite], af_archive_file: str) -> list[str]:
    with tarfile.open(af_archive_file, 'r') as tarf:
        all_files = [i.split('-')[1] for i in tarf.getnames() if i.endswith(".pdb.gz")]
    proteins: set[str] = set()
    for site in sites:
        if site.position < max_length and site.protein in all_files:
            proteins.add(site.protein)
    return list(proteins)


class pLDDTDistanceCount:
    def __init__(self, d: float = 12, n: int = 120):
        self.d: float = d
        self.n: int = n
        self.all_sty: list[list[float]] = [[] for i in range(n)]
        self.phospho_sty: list[list[float]] = [[] for i in range(n)]


def get_ddistance(x: tuple[float, float, float], y: tuple[float, float, float]) -> float:
    return sum([(x[i] - y[i]) * (x[i] - y[i]) for i in range(3)])


def get_plddt_distance_count(sites: list[PhosphoSite], af_archive_file: str,
                             proteins: list[str], d: float = 12, n: int = 120) -> pLDDTDistanceCount:
    cnt = pLDDTDistanceCount(d, n)
    protein2sites: dict[str, list[PhosphoSite]] = {}
    for site in sites:
        if site.protein in protein2sites:
            protein2sites[site.protein].append(site)
        else:
            protein2sites[site.protein] = [site]
    dd = d * d
    with tarfile.open(af_archive_file, 'r') as tarf:
        nn = 0
        for protein in proteins:
            file = f"AF-{protein}-F1-model_v4.pdb"
            file_gz = f"{file}.gz"
            tarf.extract(file_gz, "")
            values: list[float] = []
            positions: list[tuple[float, float, float]] = []
            amino_acids: list[bool] = []
            with gzip.open(file_gz, 'rt') as fs:
                for line in fs:
                    if line.startswith(pLDDT_row):
                        # spl = line.rstrip().split()
                        atom = line[13:16]
                        if atom == pLDDT_atom:
                            values.append(float(line[61:66]))
                            amino_acids.append(line[17:20] in pLDDT_amino_acids)
                            positions.append((float(line[31:38]), float(line[39:46]), float(line[47:54])))
            m = len(amino_acids)
            for site in protein2sites[protein]:
                if site.position < len(values) and amino_acids[site.position]:
                    for j in range(m):
                        if site.position == j:
                            continue
                        dd0 = get_ddistance(positions[site.position], positions[j])
                        if dd0 < dd:
                            ii = int(math.sqrt(dd0) * n / d)
                            cnt.phospho_sty[ii].append(values[j])
                    amino_acids[site.position] = False
            for i in range(len(amino_acids)):
                if amino_acids[i]:
                    for j in range(m):
                        if i == j:
                            continue
                        dd0 = get_ddistance(positions[i], positions[j])
                        if dd0 < dd:
                            ii = int(math.sqrt(dd0) * n / d)
                            cnt.all_sty[ii].append(values[j])
            os.remove(file_gz)
            if nn % 1000 == 0:
                print(nn)
            nn += 1
    return cnt


def plot_plddt_distance_distribution(sites: list[PhosphoSite], af_archive_file: str, proteins: list[str]) -> None:
    d: float = 20
    n: int = 200
    k: int = 30
    cnt: pLDDTDistanceCount = get_plddt_distance_count(sites, af_archive_file, proteins, d, n)
    fig, ax = plt.subplots()
    xs = [i for i in range(k, n)]
    m = len(cnt.phospho_sty)
    for i in range(m):
        cnt.phospho_sty[i].sort()
        cnt.all_sty[i].sort()

    ax.plot(xs, [cnt.phospho_sty[i][int(0.25 * len(cnt.phospho_sty[i]))] for i in range(k, m)], color="blue")
    ax.plot(xs, [cnt.phospho_sty[i][int(0.50 * len(cnt.phospho_sty[i]))] for i in range(k, m)], color="blue")
    ax.plot(xs, [cnt.phospho_sty[i][int(0.75 * len(cnt.phospho_sty[i]))] for i in range(k, m)], color="blue")

    ax.plot(xs, [cnt.all_sty[i][int(0.25 * len(cnt.all_sty[i]))] for i in range(k, m)], color="red")
    ax.plot(xs, [cnt.all_sty[i][int(0.50 * len(cnt.all_sty[i]))] for i in range(k, m)], color="red")
    ax.plot(xs, [cnt.all_sty[i][int(0.75 * len(cnt.all_sty[i]))] for i in range(k, m)], color="red")

    ax.set_xlim([0, n])
    ax.set_ylim([20, 100])
    # ax.legend()
    plt.show()


if __name__ == "__main__":
    fasta_file = "UP000000589.2023_08_25.canonical.fasta"
    peptide_file = "Protein_Expression.canonical_swiss\\peptides.txt"

    af_archive_file = "UP000000589_10090_MOUSE_v4.tar"
    phospho_swiss_file = "Phosphorylation.canonical.imputed.ver1.txt"
    process(af_archive_file, phospho_swiss_file)
