#################################################################
# Copyright(c) 2004 Whitehead Institute for Biomedical Research.
#              All Right Reserve
#
# Created:     1/23/2004
# author:      Bingbing Yuan
#################################################################
 
# ============================================================================
# Objetivo: calcular a energia termodinâmica da extremidade 3' de uma fita
# dangling (em balanço)
# Requer 2 arquivos: dangling_3_energy.txt e nearest_neighbor_energy.txt
# ============================================================================
from pathlib import Path
from Bio.Seq import Seq


class Thermodynamics:
    def __init__(self, seq, nearest_neighbor_file='sirna_seeker/perl/nearest_neighbor_energy.txt', dangling_file="sirna_seeker/perl/dangling_3_energy.txt"):
        self.seq = seq.upper()
        self.base = self.seq[5]
        self.nearest_neighbor_file = nearest_neighbor_file
        self.dangling_file = dangling_file
        self.energy = None
        self.antisense_energy = None
        self.energy_diff = None

    def cal_energy(self):
        nb_hash = self.file_to_hash(self.nearest_neighbor_file, key_pos=0)
        dangling_hash = self.file_to_hash(self.dangling_file, key_pos=1)

        ns_seq = self.seq[:5]
        sum_nb = self.cal_nb(ns_seq, nb_hash)

        dangling_seq = self.seq[4:6]
        dangling_score = self.cal_dangling(dangling_seq, self.base, dangling_hash)

        self.energy = sum_nb + dangling_score

    def cal_duplex_energy_diff(self):
        # Gerar a fita antisenso: complemento reverso da fita senso
        antisense_seq = Seq(self.seq).reverse_complement()
        antisense_base = antisense_seq[5]

        nb_hash = self.file_to_hash(self.nearest_neighbor_file, key_pos=0)
        dangling_hash = self.file_to_hash(self.dangling_file, key_pos=1)

        ns_seq = antisense_seq[:5]
        sum_nb = self.cal_nb(ns_seq, nb_hash)

        dangling_seq = antisense_seq[4:6]
        dangling_score = self.cal_dangling(dangling_seq, antisense_base, dangling_hash)

        self.antisense_energy = sum_nb + dangling_score
        if self.energy is None:
            self.cal_energy()
        self.energy_diff = round((self.energy - self.antisense_energy), 2)

    def cal_nb(self, seq, nb_hash):
        return sum(float(nb_hash.get(seq[i:i+2], 0)) for i in range(len(seq) - 1))


    def cal_dangling(self, seq, base, dangling_hash):
        values = dangling_hash.get(seq, [])
        for i in range(0, len(values), 2):
            if values[i] == base:
                return float(values[i + 1])
        return 0

    def file_to_hash(self, file_path, key_pos=0):
        hash_data = {}
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {file_path}")
        with path.open('r', encoding='utf-8') as file:
            for line in file:
                if line.strip() == "":
                    continue
                fields = line.strip().split('\t')
                key = fields[key_pos]
                values = [v for i, v in enumerate(fields) if i != key_pos]
                if len(fields) == 2:
                    hash_data[key] = values[0]
                else:
                    hash_data[key] = values
        return hash_data

    def reverse_complement(self, seq):
        comp = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(comp.get(nt, 'N') for nt in reversed(seq))

    def complement(self, nt):
        comp = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        return comp.get(nt, 'N')


thermo = Thermodynamics(
    seq="AGCTGAC"
)


