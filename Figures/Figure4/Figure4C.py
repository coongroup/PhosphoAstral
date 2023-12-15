import gzip
import os
import tarfile

import math
import os
import matplotlib.pyplot as plt

max_length = 1400

pLDDT_row = "ATOM"
pLDDT_atom = "CA "
pLDDT_atom_column = 2
pLDDT_amino_acid_column = 3
pLDDT_position_column = 5
pLDDT_amino_acids = {"SER", "THR", "TYR"}
pLDDT_value_column = 10


class PhosphoSite:
    def __init__(self, protein: str, position: int, amino_acid: str, is_in_all_tissues: bool):
        self.protein = protein
        self.position = position
        self.amino_acid = amino_acid
        self.amino_acid = amino_acid
        self.is_in_all_tissues = is_in_all_tissues


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
                cnt = [1 for i in range(expression_index, len(spl)) if float(spl[i]) > 0]
                phospho_sites.append(
                    PhosphoSite(
                        spl[major_protein_index],
                        position,
                        spl[amino_acid_index],
                        sum(cnt) == (len(spl) - expression_index)
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


def read_fasta(fasta_file: str) -> dict[str, str]:
    fasta_records: dict[str, str] = {}
    with open(fasta_file) as fs:
        name = ""
        for line in fs:
            if line.startswith('>'):
                name = line.split('|')[1]
                fasta_records[name] = ""
            else:
                fasta_records[name] += line.rstrip()
    return fasta_records


def read_peptides(peptide_file: str) -> dict[str, list[str]]:
    peptides: dict[str, list[str]] = {}
    with open(peptide_file) as fs:
        header = fs.readline().rstrip().split('\t')
        contaminant_index = header.index("Potential contaminant")
        reverse_index = header.index("Reverse")
        peptide_index = 0
        protein_index = header.index("Proteins")
        for line in fs:
            spl = line.rstrip().split('\t')
            if spl[reverse_index] == '+' and spl[contaminant_index] == '+':
                continue
            peptides[spl[peptide_index]] = [p for p in spl[protein_index].split(';')]
    return peptides


def read_protein_coverage(fasta_file: str, peptide_file: str) -> dict[str, list[bool]]:
    fasta_records: dict[str, str] = read_fasta(fasta_file)
    peptides: dict[str, list[str]] = read_peptides(peptide_file)
    protein_coverage: dict[str, list[bool]] = {pname: [False for _ in pseq] for pname, pseq in fasta_records.items()}
    for peptide, proteins in peptides.items():
        for protein in proteins:
            if protein not in protein_coverage:
                continue
            index = fasta_records[protein].index(peptide)
            if index == -1:
                continue
            for i in range(index, index + len(peptide)):
                protein_coverage[protein][i] = True
    return protein_coverage


class Count:
    def __init__(self):
        self.all_sty: list[float] = []
        self.phospho_sty: list[float] = []
        self.phospho_sty_all_tissues: list[float] = []
        self.all: list[float] = []
        self.protein_covered: list[float] = []


def get_count(sites: list[PhosphoSite], af_archive_file: str, proteins: list[str],
              protein_coverage: dict[str, list[bool]]) -> Count:
    cnt = Count()
    protein2sites: dict[str, list[PhosphoSite]] = {}
    for site in sites:
        if site.protein not in protein_coverage:
            continue
        if site.protein in protein2sites:
            protein2sites[site.protein].append(site)
        else:
            protein2sites[site.protein] = [site]

    with tarfile.open(af_archive_file, 'r') as tarf:
        n = 0
        for protein in proteins:
            file = f"AF-{protein}-F1-model_v4.pdb"
            file_gz = f"{file}.gz"
            tarf.extract(file_gz, "")
            values: list[float] = []
            amino_acids: list[str] = []
            with gzip.open(file_gz, 'rt') as fs:
                for line in fs:
                    if line.startswith(pLDDT_row):
                        # spl = line.rstrip().split()
                        atom = line[13:16]
                        if atom == pLDDT_atom:
                            values.append(float(line[61:66]))
                            amino_acids.append(line[17:20])

            for site in protein2sites[protein]:
                if site.position < len(values) and amino_acids[site.position] in pLDDT_amino_acids:
                    cnt.phospho_sty.append(values[site.position])
                    if site.is_in_all_tissues:
                        cnt.phospho_sty_all_tissues.append(values[site.position])
            for i in range(len(amino_acids)):
                cnt.all.append(values[i])
                if i < len(protein_coverage[protein]) and protein_coverage[protein][i]:
                    cnt.protein_covered.append(values[i])
                if amino_acids[i] in pLDDT_amino_acids:
                    cnt.all_sty.append(values[i])
            os.remove(file_gz)
            if n % 1000 == 0:
                print(n)
            n += 1
    return cnt


def plot_density(sites: list[PhosphoSite], af_archive_file: str, proteins: list[str],
              protein_coverage: dict[str, list[bool]]) -> None:
    cnt: Count = get_count(sites, af_archive_file, proteins, protein_coverage)
    fig, ax = plt.subplots()
    for a, lbl, clr in zip([cnt.phospho_sty, cnt.all, cnt.protein_covered, cnt.all_sty, cnt.phospho_sty_all_tissues],
                           ["phospho_sty", "all", "protein_covered", "all_sty", "phospho_sty_all_tissues"],
                           ["red", "green", "blue", "yellow", "green"]):
        ax.hist(a, bins=[i for i in range(0, 100)], histtype="step", label=lbl, density=True)
    ax.set_xlim([0, 100])
    ax.legend()
    plt.show()

def process(fasta_file: str, peptide_file: str, af_archive_file: str, phospho_swiss_file: str):
    protein_coverage: dict[str, list[bool]] = read_protein_coverage(fasta_file, peptide_file)
    sites: list[PhosphoSite] = read_phospho_sites(phospho_swiss_file)
    proteins: list[str] = get_suitable_proteins(sites, af_archive_file)
    plot_density(sites, af_archive_file, proteins, protein_coverage)


if __name__ == "__main__":
    fasta_file = 'UP000000589.2023_08_25.canonical.fasta'
    peptide_file = "Protein_Expression.canonical_swiss/peptides.txt"
    af_archive_file = 'UP000000589_10090_MOUSE_v4.tar'
    phospho_swiss_file = "Phosphorylation.canonical.imputed.ver1.txt"
    process(fasta_file, peptide_file, af_archive_file, phospho_swiss_file)
