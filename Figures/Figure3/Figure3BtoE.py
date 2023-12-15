import math

import matplotlib.pyplot as plt
import numpy


class SiteRecord:
    def __init__(self, protein_groups: list[str], positions: list[int], genes: list[str], amino_acid: str,
                 flanking_region: str, expression: dict[str, float]):
        self.protein_groups = protein_groups
        self.positions = positions
        self.genes = genes
        self.amino_acid = amino_acid
        self.flanking_region = flanking_region
        self.expression = expression


def read_site_records(filename: str) -> list[SiteRecord]:
    sites: list[SiteRecord] = []
    with open(filename) as fs:
        header = fs.readline().rstrip().split('\t')
        indexes = {i: header.index(i) for i in
                   ["Protein Groups", "Positions", "Genes", "Amino Acid", "Flanking Region"]}
        conditions = {header[i]: i for i in range(indexes["Flanking Region"] + 1, len(header))}
        for line in fs:
            spl = line.rstrip().split('\t')
            sites.append(
                SiteRecord(
                    spl[indexes["Protein Groups"]].split(';'),
                    [int(i) for i in spl[indexes["Positions"]].split(';')],
                    spl[indexes["Genes"]].split(';'),
                    spl[indexes["Amino Acid"]],
                    spl[indexes["Flanking Region"]],
                    {condition: float(spl[i]) for condition, i in conditions.items()}
                )
            )
    return sites


def plot_count_dist(records: list[SiteRecord], maxvalue: int):
    n = len(records[0].expression)
    counts = [0 for i in range(n)]
    values = []
    for record in records:
        counts[sum([1 for condition, intensity in record.expression.items() if intensity > 0.0]) - 1] += 1
        values += [math.log2(intensity) for condition, intensity in record.expression.items() if intensity > 0.0]
    print(counts)

    # fig, ax = plt.subplots()
    # ax.bar(range(1, n + 1), counts, align='center')
    # ax.set_ylim((0, maxvalue))
    # plt.show()

    mean = numpy.mean(values)
    std = numpy.std(values)
    fig, ax = plt.subplots()
    violin_values = [[] for _ in range(1, n + 1)]
    for record in records:
        x = sum([1 for condition, intensity in record.expression.items() if intensity > 0.0]) - 1
        for condition, intensity in record.expression.items():
            if intensity > 0.0:
                violin_values[x].append((math.log2(intensity) - mean) / std)

    ax.violinplot(violin_values, range(1, n + 1), vert=True, widths=0.75, showextrema=False,
                  points=1000, quantiles=[[0.1, 0.25, 0.5, 0.75, 0.9] for _ in range(1, n + 1)])
    ax.set_ylim((-4, 4))
    plt.show()

    cdf_counts = [sum(counts[:(i + 1)]) for i in range(n)]
    xvalues = [i + 1 for i in range(n)]
    print(cdf_counts, xvalues)

    fig, ax = plt.subplots()
    ax.plot(xvalues, cdf_counts)
    ax.set_ylim((0, maxvalue * 4))
    plt.show()


def plot_common_unique_counts(left_records: list[SiteRecord], right_records: list[SiteRecord], maxvalue: int):
    right = set(right_records[0].expression.keys())
    count_sets = {tissue: [set(), set()] for tissue in right}
    for record in left_records:
        for tissue, value in record.expression.items():
            if value > 0.0:
                count_sets[tissue][0].add(record.flanking_region)
    for record in right_records:
        for tissue, value in record.expression.items():
            if value > 0.0:
                count_sets[tissue][1].add(record.flanking_region)
    counts = {tissue: [len(cs[0] - cs[1]), len(cs[0] & cs[1]), len(cs[1] - cs[0])] for tissue, cs in count_sets.items()}

    fig, ax = plt.subplots()
    ys = [k for k, v in sorted(counts.items(), key=lambda item: sum(item[1]), reverse=True)]

    a0 = [counts[y][0] for y in ys]
    a1 = [counts[y][1] for y in ys]
    a01 = [a0[i] + a1[i] for i in range(len(ys))]
    a2 = [counts[y][2] for y in ys]
    ax.barh(range(len(ys)), a0, align='center')
    ax.barh(range(len(ys)), a1, left=a0, align='center')
    ax.barh(range(len(ys)), a2, left=a01, align='center')
    for i in reversed(range(len(ys))):
        print(ys[i], a0[i] + a1[i], a1[i] + a2[i])

    ax.set_yticks(range(len(ys)), labels=ys)
    plt.show()


def plot_sty_counts(records: list[SiteRecord], maxvalue: int):
    aas = {"S": 0, "T": 1, "Y": 2}
    counts = {tissue: [0, 0, 0] for tissue, cs in records[0].expression.items()}
    counts['all'] = [0, 0, 0]
    for record in records:
        trigger = False
        for tissue, cs in record.expression.items():
            if cs > 0.0:
                trigger = True
                counts[tissue][aas[record.amino_acid]] += 1
        if trigger:
            counts["all"][aas[record.amino_acid]] += 1

    fig, ax = plt.subplots()
    ys = [k for k, v in sorted(counts.items(), key=lambda item: sum(item[1]), reverse=True)]

    a0 = [counts[y][0] for y in ys]
    a1 = [counts[y][1] for y in ys]
    a01 = [a0[i] + a1[i] for i in range(len(ys))]
    a2 = [counts[y][2] for y in ys]
    ax.barh(range(len(ys)), a0, align='center')
    ax.barh(range(len(ys)), a1, left=a0, align='center')
    ax.barh(range(len(ys)), a2, left=a01, align='center')
    for i in reversed(range(len(ys))):
        print(ys[i], a0[i], a1[i], a2[i], sum(counts[ys[i]]))

    ax.set_yticks(range(len(ys)), labels=ys)
    ax.set_xlim((0, maxvalue))
    plt.show()


if __name__ == "__main__":
    huttlin_phospho_data = "C:\\Data\\Astral_Phospho\\data\\Phosphorylation.huttlin2010.ver1.txt"
    huttlin_phospho_records = read_site_records(huttlin_phospho_data)
    print(len(huttlin_phospho_records))
    plot_count_dist(huttlin_phospho_records, 25000)
    our_phospho_data = "C:\\Data\\Astral_Phospho\\data\\Phosphorylation.ver4.txt"
    our_phospho_records = read_site_records(our_phospho_data)
    print(len(our_phospho_records))
    plot_count_dist(our_phospho_records, 25000)

    # plot_common_unique_counts(huttlin_phospho_records, our_phospho_records, 25000)
    # plot_sty_counts(our_phospho_records, 100000)
    # print()
    # plot_sty_counts(huttlin_phospho_records, 100000)
