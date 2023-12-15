import matplotlib.pyplot as plt
import math

label = "Phospho (STY)"


class Precursor:
    def __init__(self, precursor: str, qvalue: float, is_human: bool, is_phospho: bool):
        self.precursor = precursor
        self.qvalue = qvalue
        self.is_human = is_human
        self.is_phospho = is_phospho
        # self.is_decoy = is_decoy


def read(filename: str) -> list[Precursor]:
    precursors: list[Precursor] = []
    with open(filename, 'r') as fs:
        header = fs.readline().rstrip().split('\t')
        precursor_index = header.index("EG.PrecursorId")
        qvalue_index = header.index("EG.Qvalue")
        proteins_index = header.index("PG.ProteinNames")
        # is_decoy_index = header.index("EG.IsDecoy")
        for line in fs:
            spl = line.rstrip().split('\t')
            proteins = spl[proteins_index].split(';')
            cnt_h = sum([1 for p in proteins if p.endswith("_HUMAN")])
            cnt_m = sum([1 for p in proteins if p.endswith("_MAIZE")])
            if (cnt_m == 0 and cnt_h == 0) or (cnt_m != 0 and cnt_h != 0):
                continue
            precursors.append(
                Precursor(
                    spl[precursor_index],
                    float(spl[qvalue_index]),
                    cnt_h > 0,
                    label in spl[precursor_index],
                    # spl[is_decoy_index] == "True"
                )
            )
    return precursors


def plot_qq(precursors: list[Precursor]) -> None:
    fig, ax = plt.subplots()
    maxv = 0.0
    minv = 1.0
    for color, label in zip(["red", "blue"], ["all", "phospho"]):
        xs = []
        ys = []
        cnt_h = 1.0
        cnt_m = 1.0
        for precursor in precursors:
            if label == "phospho" and not precursor.is_phospho:  # all
                continue
            if precursor.is_human:
                cnt_h += 1
            else:
                cnt_m += 1
            xs.append(precursor.qvalue)
            ys.append(cnt_m / (cnt_m + cnt_h))
        minv = min(minv, min(min(xs), min(ys)))
        # maxv = max(maxv, max(xs[-1], ys[-1]))
        ax.plot(xs, ys, linestyle="-", marker="", c=color, label=label, linewidth=0.8)
    ax.plot((0, 1), (0, 1), linestyle="-", marker="", c="grey")
    maxv = 1.0
    print(minv, maxv)
    ax.set_xlim([minv, maxv])
    ax.set_xlabel("Internal FDR")
    ax.set_xscale('log')
    ax.set_ylim([minv, maxv])
    ax.set_ylabel("External FDR")
    ax.set_yscale('log')
    # ax.set_aspect('equal', adjustable='box')
    # for aa, probs in cnt.items():
    #     n = len(probs)
    #     xs = []
    #     ys = []
    #     probs0 = sorted(probs)
    #     for i in range(n):
    #         xs.append(((i + 1) / n) * 100)
    #         ys.append(probs0[i])
    #     ax.plot(xs, ys, linestyle="-", marker="", c=aas_colors[aa], label=f"p{aa}", linewidth=3)
    #     # break
    ax.legend()
    plt.show()


def plot_count(precursors: list[Precursor]) -> None:
    fig, ax = plt.subplots()
    maxv = 0.0
    for color, label in zip(["red", "blue"], ["all", "phospho"]):
        xs0 = []
        ys0 = []
        xs1 = []
        ys1 = []
        cnt_h = 0.0
        cnt_m = 0.0
        for precursor in precursors:
            # if precursor.is_decoy:
            #     continue
            if label == "phospho" and not precursor.is_phospho:  # all
                continue
            if precursor.is_human:
                cnt_h += 1
            else:
                cnt_m += 1
            ys0.append(cnt_m + cnt_h)
            xs0.append(cnt_m / (cnt_m + cnt_h))
            ys1.append(cnt_m + cnt_h)
            xs1.append(precursor.qvalue)
        maxv = max(maxv, max(ys0[-1], ys1[-1]))
        # maxv = max(maxv, max(xs[-1], ys[-1]))
        ax.plot(xs0, ys0, linestyle="-", marker="", c=color, label=f"{label} external", linewidth=0.8)
        ax.plot(xs1, ys1, linestyle="--", marker="", c=color, label=f"{label} internal", linewidth=0.8)
        # break
    # ax.plot((0, 1), (0, 1), linestyle="-", marker="", c="grey")
    print(maxv)
    ax.set_xlim([-0.0001, 0.0015])
    ax.set_xlabel("Estimated FDR")
    ax.set_ylim([1000, maxv * 1.1])
    ax.set_ylabel("Precursor Count")
    ax.set_yscale('log')
    # ax.set_aspect('equal', adjustable='box')
    # for aa, probs in cnt.items():
    #     n = len(probs)
    #     xs = []
    #     ys = []
    #     probs0 = sorted(probs)
    #     for i in range(n):
    #         xs.append(((i + 1) / n) * 100)
    #         ys.append(probs0[i])
    #     ax.plot(xs, ys, linestyle="-", marker="", c=aas_colors[aa], label=f"p{aa}", linewidth=3)
    #     # break
    ax.legend()
    plt.show()


if __name__ == "__main__":
    filename = "zea_mays.PSM100_PRT100\\20231005_181451_PSM100_PRT100\\PSM100_PRT100_Report_Peptide Quant (Normal).tsv"
    precursors = read(filename)
    precursors.sort(key=lambda x: x.qvalue)
    plot_count(precursors)
