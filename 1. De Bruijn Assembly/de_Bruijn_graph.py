import networkx as nx
from tqdm import tqdm
# import matplotlib.pyplot as plt


class Kmer():
    def __init__(self, s):
        self.kmer = s
        self.k = len(s)
        self.prefix = s[:-1]
        self.suffix = s[1:]


class De_Bruijn_graph():
    def __init__(self, k=10):
        self.G = nx.MultiDiGraph()
        self.k = k
        self.components = []
        self.num_components = 0

    def get_kmers(self, seq):
        return [Kmer(seq[i:i + self.k]) for i in range(len(seq) - self.k + 1)]

    def update_graph(self, fq):
        for i in tqdm(range(len(fq))):
            for kmer in self.get_kmers(fq[i].seq):
                self.G.add_nodes_from([kmer.suffix, kmer.prefix])
                self.G.add_edge(kmer.suffix, kmer.prefix, label=kmer.kmer)
        self.get_components()

    def get_components(self):
        comp_nodes = nx.weakly_connected_components(self.G)
        self.components = [self.G.subgraph(nodes) for nodes in comp_nodes]
        self.num_components = len(self.components)

    @staticmethod
    def make_contig_of_path(path):
        contig = ''
        for u, v in path:
            if len(contig) == 0:
                contig += u + v[-1]
            else:
                contig += v[-1]
        return contig

    def get_contigs(self):
        contigs = []
        for component in tqdm(self.components):
            if nx.has_eulerian_path(component):
                contig = self.make_contig_of_path(nx.eulerian_path(component))
                contigs.append(contig)
        return contigs


