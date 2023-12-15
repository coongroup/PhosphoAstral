import json
import random
import math
from random import sample
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm
from matplotlib.colors import Normalize

from typing import List, Dict, Set, Tuple
from scipy.stats import gaussian_kde


class Experiment:
    def __init__(self, name: str, spectra_file: str, msms_file: str, psite_file: str, id: int, xlim: Tuple[int, int],
                 ylim: Tuple[int, int]):
        self.name = name
        self.spectra_file = spectra_file
        self.msms_file = msms_file
        self.psite_file = psite_file
        self.id = id
        self.xlim = xlim
        self.ylim = ylim


class PhosphoSite:
    def __init__(self, proteins: List[str], positions: List[int], scan_number: int, loc_prob: float):
        self.proteins = proteins
        self.positions = positions
        self.scan_number = scan_number
        self.msms = None
        self.spectrum = None
        self.loc_prob = loc_prob


def read_phospho_sites(filename: str) -> List[PhosphoSite]:
    sites: List[PhosphoSite] = []
    with open(filename, 'r') as fs:
        spl = fs.readline().rstrip().split('\t')
        reverse_index = spl.index("Reverse")
        contaminant_index = spl.index("Potential contaminant")
        localization_index = spl.index("Localization prob")
        protein_index = spl.index("ï»¿Proteins")
        position_index = spl.index("Positions within proteins")
        scan_index = spl.index("Best localization scan number")
        # msms_index = spl.index("Best localization MS/MS ID")
        for line in fs:
            spl = line.rstrip().split('\t')
            prob = float(spl[localization_index])
            if spl[reverse_index] == '+' or spl[contaminant_index] == '+' or prob < 0.75:
                continue
            sites.append(
                PhosphoSite(
                    spl[protein_index].split(';'),
                    [int(i) for i in spl[position_index].split(';')],
                    int(spl[scan_index]),
                    prob
                )
            )
    return sites


class Fragment:
    def __init__(self, match: str, mass: float, intensity: float, mass_error: float):
        self.match = match
        self.mass = mass
        self.intensity = intensity
        self.mass_error = mass_error
        self.spectrum_fragment = None


class Msms:
    def __init__(self, scan_number: int, fragments: List[Fragment]):
        self.scan_number = scan_number
        self.fragments = sorted(fragments, key=lambda x: x.mass)


def read_msms(filename: str, phospho_sites: List[PhosphoSite]) -> None:
    phospho_sites_dict: Dict[int, PhosphoSite] = {site.scan_number: site for site in phospho_sites}
    with open(filename, 'r') as fs:
        spl = fs.readline().rstrip().split('\t')
        scan_index = spl.index("Scan number")
        match_index = spl.index("Matches")
        intensity_index = spl.index("Intensities")
        mass_error_index = spl.index("Mass deviations [ppm]")
        mass_index = spl.index("Masses")
        for line in fs:
            spl = line.rstrip().split('\t')
            scan_number = int(spl[scan_index])
            if scan_number not in phospho_sites_dict:
                continue
            phospho_sites_dict[scan_number].msms = Msms(
                scan_number,
                [Fragment(match, float(mass), float(intensity), float(mass_error)) for
                 match, mass, intensity, mass_error in
                 zip(spl[match_index].split(';'), spl[mass_index].split(';'), spl[intensity_index].split(';'),
                     spl[mass_error_index].split(';'))]
            )


class FragmentShort:
    def __init__(self, mass: float, intensity: float, resolution: float):
        self.mass = mass
        self.intensity = intensity
        self.resolution = resolution


class Spectrum:
    def __init__(self, scan_number: int, fragments: List[FragmentShort]):
        self.scan_number = scan_number
        self.fragments = sorted(fragments, key=lambda x: x.mass)


def merge(f0: List[Fragment], f1: List[FragmentShort]):
    j = 0
    for i in range(len(f0)):
        while j < len(f1) and abs(f0[i].mass - f1[j].mass) > 0.00001 and \
                f1[j].mass - f0[i].mass < 0.00001:
            j += 1
        if j == len(f1):
            continue
        if f1[j].mass - f0[i].mass >= 0.00001:
            continue
        if abs(f0[i].intensity - f1[j].intensity) < 1.0:
            f0[i].spectrum_fragment = f1[j]


def find_closest(arr, target):
    n = len(arr)
    if target <= arr[0].mass:
        return 0
    if target >= arr[n - 1].mass:
        return n - 1
    i = 0
    j = n
    mid = 0
    while i < j:
        mid = (i + j) // 2
        if arr[mid].mass == target:
            return mid
        if target < arr[mid].mass:
            if mid > 0 and target > arr[mid - 1].mass:
                return get_closest(arr, mid - 1, mid, target)
            j = mid
        else:
            if mid < n - 1 and target < arr[mid + 1].mass:
                return get_closest(arr, mid, mid + 1, target)

            i = mid + 1
    return mid


def get_closest(a, i1, i2, target):
    v1 = a[i1].mass
    v2 = a[i2].mass
    if target - v1 >= v2 - target:
        return i2
    else:
        return i1


def merge0(f0: List[Fragment], f1: List[FragmentShort]):
    for i in range(len(f0)):
        j = find_closest(f1, f0[i].mass)
        f0[i].spectrum_fragment = f1[j]


def read_spectra(filename: str, phospho_sites: List[PhosphoSite]) -> None:
    phospho_sites_dict: Dict[int, PhosphoSite] = {site.scan_number: site for site in phospho_sites}
    with open(filename, 'r') as fs:
        spl = fs.readline().rstrip().split('\t')
        scan_index = spl.index("Id")
        mass_index = spl.index("Mass")
        intensity_index = spl.index("Intensity")
        res_index = spl.index("Resolution")
        for line in fs:
            spl = line.rstrip().split('\t')
            scan_number = int(spl[scan_index])
            if scan_number not in phospho_sites_dict:
                continue
            phospho_sites_dict[scan_number].spectrum = Spectrum(
                scan_number,
                [FragmentShort(float(mass), float(intensity), int(resolution)) for
                 mass, intensity, resolution in
                 zip(spl[mass_index].split(';'), spl[intensity_index].split(';'), spl[res_index].split(';'))]
            )
            merge0(phospho_sites_dict[scan_number].msms.fragments, phospho_sites_dict[scan_number].spectrum.fragments)


def print_stats(phospho_sites: List[PhosphoSite]):
    print(f"N={len(phospho_sites)}")
    a = 0
    b = 0
    for site in phospho_sites:
        if site.msms is None:
            a += 1
        if site.spectrum is None:
            b += 1
    print(f"MSMS: {a}; Spectra: {b}")
    c = 0
    n = 0
    for site in phospho_sites:
        if site.msms is None:
            continue
        for fragment in site.msms.fragments:
            if fragment.spectrum_fragment is None:
                c += 1
            n += 1
    print(f"{c}/{n}")


def make_colours(vals):
    norm = Normalize(vmin=vals.min(), vmax=vals.max())
    colours = [cm.ScalarMappable(norm=norm, cmap='jet').to_rgba(val) for val in vals]
    return colours


def plot_resolution(sites, xlim, ylim, n=30000, fig_size=5, dot_size=2):
    fig, ax = plt.subplots(1, 1, figsize=(fig_size, fig_size))
    xys = []
    for site in sites:
        if site.msms is None or site.spectrum is None:
            continue
        for fragment in site.msms.fragments:
            if fragment.spectrum_fragment is None:
                continue
            xys.append((fragment.mass, fragment.spectrum_fragment.resolution))
    xys1 = sample(xys, n)
    xs = [x for x, y in xys1]
    ys = [y for x, y in xys1]
    values = np.vstack([xs, ys])
    z = gaussian_kde(values)(values)
    idx = z.argsort()
    ax.scatter(np.array(xs)[idx], np.array(ys)[idx], color=make_colours(z[idx]), s=dot_size)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.show()


def plot_resolution2(sites, xlim, ylim, bin=100, n=100000, fig_size=5, line_width=2):
    p = 0.05
    n = min(n, len(sites))
    fig, ax = plt.subplots(1, 1, figsize=(fig_size, fig_size))
    xys = []
    for site in sites:
        if site.msms is None or site.spectrum is None:
            continue
        for fragment in site.msms.fragments:
            if fragment.spectrum_fragment is None:
                continue
            xys.append((fragment.mass, fragment.spectrum_fragment.resolution))
    xys1 = sample(xys, n)
    xys1.sort()

    xs = []
    ys = []
    ys0 = []
    ys1 = []
    for i in range(0, (n // bin) * bin, bin):
        j = i + bin // 2
        xs.append(xys1[j][0])
        y = sorted([xys1[i0][1] for i0 in range(i, min(n, i + bin))])
        ys.append(y[bin // 2])
        ys0.append(y[int(bin * p)])
        ys1.append(y[int(bin * (1 - p))])

    ax.plot(xs, ys, '-', linewidth=line_width, color="black")
    ax.plot(xs, ys0, '-', linewidth=line_width, color="grey")
    ax.plot(xs, ys1, '-', linewidth=line_width, color="grey")

    yss = [120000 * math.sqrt(200) * math.sqrt(1 / mz) for mz in range(100, 1000, 10)]
    xss = [mz for mz in range(100, 1000, 10)]
    ax.plot(xss, yss, '-', linewidth=line_width, color="blue")
    ax.set_yscale('log', base=10)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.show()


def plot_mass_error(sites):
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    xs = []
    ys = []
    ss = 0.0
    for site in sites:
        if site.msms is None or site.spectrum is None:
            continue
        for fragment in site.msms.fragments:
            if fragment.spectrum_fragment is None:
                continue
            xs.append(fragment.mass)
            ys.append(fragment.mass_error)
            ss += fragment.mass_error
    ss = ss / len(ys)
    # values = np.vstack([xs[:10000], ys[:10000]])
    # z = gaussian_kde(values)(values)
    # idx = z.argsort()
    # ax.scatter(np.array(xs)[idx], np.array(ys)[idx], color=make_colours(z[idx]), s=1)
    ax.set_xlim((-15, 15))
    ax.set_ylim((0, 0.32))
    ax.hist([y - ss for y in ys], bins=100, density=True)
    print(f"Length: {len(ys)}; -1...+1: {len([1 for y in ys if abs(y - ss) <= 5])}")
    plt.show()
    print(np.std(ys))


def plot(parameters: Dict):
    experiments = []
    for name, description in parameters["input"]["experiments"].items():
        experiments.append(
            Experiment(
                name,
                description["spectra"],
                description["msms"],
                description["psite"],
                int(description["id"]),
                (description["xlim"]["min"], description["xlim"]["max"]),
                (description["ylim"]["min"], description["ylim"]["max"])
            )
        )
    # print(experiments[0].id)
    for experiment in experiments:
        sites = read_phospho_sites(experiment.psite_file)
        read_msms(experiment.msms_file, sites)
        read_spectra(experiment.spectra_file, sites)
        print_stats(sites)
        plot_resolution2(sites, experiment.xlim, experiment.ylim, bin=parameters["input"]["bin_size"])
        plot_mass_error(sites)


if __name__ == "__main__":
    with open("parameters.json", 'r') as parameters_fs:
        plot(json.loads(parameters_fs.read()))