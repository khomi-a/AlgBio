import networkx as nx
import matplotlib.pyplot as plt

class Kmer():
    def __init__(self, s):
        self.kmer = s
        self.k = len(s)
        self.prefix = s[:-1]
        self.suffix = s[1:]


class De_Bruijn_graph():
    def __int__(self, k=10):
        self.G = nx.Graph()
        self.k = k

    def get_kmers(self, seq):
        return [Kmer(seq[i:i+self.k]) for i in range(len(seq) - self.k)]


    def build_graph(self, fq):
        for read in fq:
            for kmer in self.get_kmers(read.seq):
                self.G.add_nodes_from([kmer.suffix, kmer.prefix)

    def build_graph_from_strings(self, s):
        for read in s:
            for kmer in self.get_kmers(read):
                self.G.add_nodes_from([kmer.suffix, kmer.prefix)

def visualize_graph(graph):
    """
    A naive method to plot graph
    :param G:
    :return:
    """
    nx.draw(graph, with_labels=True, cmap=plt.get_cmap('jet'))
    plt.show()




