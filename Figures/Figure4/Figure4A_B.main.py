from enum import Enum

import logging
import random
import sys
import numpy
from sklearn.manifold import MDS, TSNE
from matplotlib import pyplot as plt
import seaborn as sns

M: int = 15
N: int = 5
K: int = -1


class Occurrence(Enum):
    ONE = 1
    REST = 2
    ALL = 3


tissues: list[str] = []


class PhosphoSite:
    def __init__(self, protein: str, position: int, amino_acid: str, flanking_region: str, occurrence: Occurrence,
                 expression: list[bool], only_tissue: str):
        self.protein = protein
        self.position = position
        self.amino_acid = amino_acid
        self.flanking_region = flanking_region
        self.occurrence = occurrence
        self.only_tissue = only_tissue
        self.expression = expression


def read_phospho_sites(filename: str) -> list[PhosphoSite]:
    phospho_sites: list[PhosphoSite] = []
    with open(filename) as fs:
        header = fs.readline().rstrip().split('\t')
        major_protein_index = header.index("Major Protein")
        major_position_index = header.index("Major Position")
        amino_acid_index = header.index("Amino Acid")
        flanking_region_index = header.index("Flanking Region")
        expression_index = flanking_region_index + 1
        global tissues
        tissues = [header[i] for i in range(expression_index, len(header))]
        for line in fs:
            spl = line.rstrip().split('\t')
            position = int(spl[major_position_index]) - 1
            occurrence = Occurrence.REST
            cnt = 0
            expression: list[bool] = []
            for i in range(expression_index, len(spl)):
                if float(spl[i]) > 0:
                    cnt += 1
                    expression.append(True)
                else:
                    expression.append(False)
            only_tissue = "None"
            if cnt == 1:
                only_tissue = tissues[[i for i in range(len(tissues)) if expression[i]][0]]
                occurrence = Occurrence.ONE
            elif cnt == (len(spl) - expression_index):
                occurrence = Occurrence.ALL
            phospho_sites.append(
                PhosphoSite(
                    spl[major_protein_index],
                    position,
                    spl[amino_acid_index],
                    spl[flanking_region_index],
                    occurrence,
                    expression,
                    only_tissue
                )
            )
    return phospho_sites


def read_blosum(blosum_file: str) -> dict[str, dict[str, int]]:
    table: dict[str, dict[str, int]] = {}
    with open(blosum_file) as fs:
        aas = fs.readline().rstrip().split()
        for line in fs:
            spl = line.rstrip().split()
            aa0 = spl[0]
            table[aa0] = {}
            for aa, score in zip(aas, spl[1:]):
                table[aa0][aa] = int(score)
    return table


def distance(a: str, b: str, sim: dict[str, dict[str, int]]) -> float:
    values: list[float] = []
    for i in range(center - N, center + N + 1):
        a0 = a[i]
        b0 = b[i]
        if a0 == '_':
            a0 = 'X'
        if b0 == '_':
            b0 = 'X'
        ab = sim[a0][b0] + 4
        aa = sim[a0][a0] + 4
        bb = sim[b0][b0] + 4
        values.append(2 * ab / (aa + bb))
    values.sort(reverse=True)
    return 1 - sum(values[:4]) / 4


def get_dissimilarity_array(sequences: list[str], blosum: dict[str, dict[str, int]]) -> numpy.array:
    n: int = len(sequences)
    dissimilarity_matrix: numpy.array = numpy.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            score = distance(sequences[i], sequences[j], blosum)
            dissimilarity_matrix[i, j] = score
            dissimilarity_matrix[j, i] = score
    return dissimilarity_matrix


class CollapseMethod(Enum):
    TSNE = 1,
    MDS = 2


def process_collapse(phospho_sites: list[PhosphoSite], blosum: dict[str, dict[str, int]],
                     method: CollapseMethod) -> numpy.array:
    global embedding
    sequences: list[str] = [s.flanking_region for s in phospho_sites]
    logging.info(f"Calculate dissimilarity matrix (n={len(sequences)})")
    dissimilarity_matrix: numpy.array = get_dissimilarity_array(sequences, blosum)
    logging.info(f"Calculate {method.name}")
    if method == CollapseMethod.TSNE:
        embedding = TSNE(
            n_components=2,
            metric="precomputed",
            init="random",
            n_iter=500,
            n_iter_without_progress=150,
            n_jobs=-1,
            random_state=42
        )
    elif method == CollapseMethod.MDS:
        embedding = MDS(
            n_components=2,
            dissimilarity='precomputed',
            metric=False,
            random_state=42,
            normalized_stress='auto',
            n_init=1,
            max_iter=120,
            n_jobs=-1
        )
    else:
        logging.error(f"No method like {method.name}")
        exit(-1)
    x = embedding.fit_transform(dissimilarity_matrix)
    numpy.savetxt("x.txt", x, delimiter="\t")
    return x


def plot(mds_array: numpy.array) -> None:
    plt.scatter(mds_array[:, 0], mds_array[:, 1], alpha=0.5, edgecolors=None)
    plt.show()


def write(sites: list[PhosphoSite], collapse_array: numpy.array, filename: str,
          kinase_substrates: dict[(str, int), str]) -> None:
    with open(filename, "w") as fs:
        fs.write("\t".join(
            ["x", "y", "nseq", "cseq", "amino_acid", "protein", "position", "kinase", "occurence",
             "only_tissue"] + tissues) + "\n")
        for i in range(len(sites)):
            kinase = ""
            if (sites[i].protein, sites[i].position) in kinase_substrates:
                kinase = kinase_substrates[(sites[i].protein, sites[i].position)]
            fs.write("\t".join(
                [str(collapse_array[i, 0]), str(collapse_array[i, 1]),
                 sites[i].flanking_region[(M // 2 - N): (M // 2)],
                 sites[i].flanking_region[(M // 2 + 1): (M // 2 + N + 1)],
                 sites[i].amino_acid,
                 sites[i].protein, str(sites[i].position), kinase, sites[i].occurrence.name, sites[i].only_tissue] + [
                    str(j) for j in sites[
                        i].expression]) + "\n")


def read_kinase_substrates(filename: str) -> dict[(str, int), str]:
    result: dict[(str, int), str] = {}
    with open(filename) as fs:
        spl = fs.readline().rstrip().split('\t')
        kin_index = spl.index("KINASE")
        kin_org_index = spl.index("KIN_ORGANISM")
        sub_acc_index = spl.index("SUB_ACC_ID")
        sub_org_index = spl.index("SUB_ORGANISM")
        sub_mod_rsd_index = spl.index("SUB_MOD_RSD")
        for line in fs:
            spl = line.rstrip().split('\t')
            if spl[kin_org_index] == "mouse" and spl[sub_org_index] == "mouse":
                key = (spl[sub_acc_index], int(spl[sub_mod_rsd_index][1:]) - 1)
                val = spl[kin_index]
                result[key] = val
    return result


def process(phospho_file: str, blosum_file: str, kinase_substrate_file: str):
    phospho_sites: list[PhosphoSite] = read_phospho_sites(phospho_file)
    blosum: dict[str, dict[str, int]] = read_blosum(blosum_file)
    kinase_substrates: dict[(str, int), str] = read_kinase_substrates(kinase_substrate_file)
    # ss = [s for s in phospho_sites if s.amino_acid == "T" or s.amino_acid == "Y"]
    # ss = [s for s in phospho_sites if s.amino_acid == "Y"]
    random.seed(42)
    if K == -1:
        ss = phospho_sites
    else:
        ss = random.sample(phospho_sites, K)
    collapse_array: numpy.array = process_collapse(ss, blosum, CollapseMethod.TSNE)
    write(ss, collapse_array, "file.out.txt", kinase_substrates)
    plot(collapse_array)


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    phospho_swiss_file = "Phosphorylation.canonical.imputed.ver1.txt"
    blosum_file = "blosum62.txt"
    kinase_substrate_file = "Kinase_Substrate_Dataset"
    process(phospho_swiss_file, blosum_file, kinase_substrate_file)
