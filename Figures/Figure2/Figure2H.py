import matplotlib.pyplot as plt

filenames = {
    "STYP.w_filter": "C:\\Users\\psinitcyn\\Desktop\\Experiments\\STYP.w_filter\\20230927_152846_STYP_HEK293_w\\STYP_HEK293_w_-PTMSiteReport.tsv",
}

aas = "YTSP"
mod_title = "Phospho (STYP)"
aas_colors = {
    "Y": "#A8A9AD",
    "T": "#6D6D70",
    "S": "#3F3F40",
    "P": "#0E6638"
}


class ModifiedSite:
    def __init__(self, group, amino_acid, site_prob):
        self.group = group
        self.amino_acid = amino_acid
        self.site_prob = site_prob


def read_mod_sites(filename: str) -> dict[str, ModifiedSite]:
    sites: dict[str, ModifiedSite] = {}
    with open(filename, 'r') as fs_input:
        header = fs_input.readline().rstrip().split('\t')
        file_index = header.index("R.FileName")
        group_index = header.index("PTM.Group")
        title_index = header.index("PTM.ModificationTitle")
        amino_acid_index = header.index("PTM.SiteAA")
        site_prob_index = header.index("PTM.SiteProbability")
        for line in fs_input:
            spl = line.rstrip().split('\t')
            if spl[title_index] != mod_title or spl[group_index] == "":
                continue
            # files.add(spl[file_index])
            # if spl[file_index] not in experiment_list[2:]:
            #     continue
            if spl[group_index] not in sites:  # and float(spl[site_prob_index]) > 0.2
                sites[spl[group_index]] = ModifiedSite(
                    spl[group_index],
                    spl[amino_acid_index],
                    float(spl[site_prob_index])
                )
    return sites


def count(sites: dict[str, ModifiedSite]) -> None:
    # cnt0 = {aa: 0 for aa in aas}
    # for group, site in sites.items():
    #     cnt0[site.amino_acid] += 1
    # print(cnt0)
    # # #
    # cnt1 = {aa: 0 for aa in aas}
    # for group, site in sites.items():
    #     if site.site_prob > 0.05:
    #         continue
    #     cnt1[site.amino_acid] += 1
    # print(cnt1)
    # #
    cnt2 = {aa: 0 for aa in aas}
    for group, site in sites.items():
        if site.site_prob < 0.75:
            continue
        cnt2[site.amino_acid] += 1
    print(cnt2)

    cnt3 = {aa: 0 for aa in aas}
    for group, site in sites.items():
        if site.site_prob < 0.99:
            continue
        cnt3[site.amino_acid] += 1
    print(cnt3)


def plot(sites: dict[str, ModifiedSite], title: str) -> None:
    cnt = {aa: [] for aa in aas}
    for group, site in sites.items():
        cnt[site.amino_acid].append(site.site_prob)

    fig, ax = plt.subplots()
    for aa, probs in cnt.items():
        n = len(probs)
        xs = []
        ys = []
        probs0 = sorted(probs)
        for i in range(n):
            xs.append(((i + 1) / n) * 100)
            ys.append(probs0[i])
        ax.plot(xs, ys, linestyle="-", marker="", c=aas_colors[aa], label=f"p{aa}", linewidth=3)
        # break
    ax.set_title(title)
    ax.set_xlim([0, 100])
    ax.set_xlabel("% Identified Phosphosites")
    ax.set_ylim([0, 1])
    ax.set_ylabel("Localization Probability")
    ax.legend()
    plt.show()


if __name__ == "__main__":
    for exp_name, file_name in filenames.items():
        sites = read_mod_sites(file_name)
        count(sites)
        plot(sites, exp_name)
