from Bio.Seq import Seq


def chromosome_to_cycle(chromosome):
    '''Converts genome made of synteny blocks to oriented cycle'''

    if type(chromosome) == str:
        chromosome = [int(val) for val in chromosome.split()]

    nodes = [0] * 2 * len(chromosome)

    for j in range(0, len(chromosome)):
        i = chromosome[j]
        if i > 0:
            nodes[2 * j] = 2 * i - 1
            nodes[2 * j + 1] = 2 * i
        else:
            nodes[2 * j] = -2 * i
            nodes[2 * j + 1] = -2 * i - 1
    return nodes


def cycle_to_chromosome(nodes, output = False):
    chromosome = []
    for j in range(1,len(nodes),2):
        if nodes[j-1] < nodes[j]:
            chromosome.append(nodes[j] // 2)
        else:
            chromosome.append(-nodes[j-1] // 2)
    if output:
        print(' '.join(('+' if i > 0 else '') + str(i) for i in chromosome))
    else:
        return chromosome

def preprocess(genome):
  genome = [[int(val) for val in chrom.split()] for chrom in genome[1:-1].split(')(')]
  return genome

def formation(p):
    for j in range(len(p)):
        p[j] = ('(' + ' '.join(('+' if i > 0 else '') + str(i) for i in p[j]) + ')')
    print(''.join(p))


def graph_to_genome(genome_graph):
    '''Converts genome graph to genome'''

    genome = []
    cycles = []

    if len(genome_graph) == 0:
        return -1

    while len(genome_graph) != 0:
        cycle = []
        startnode = genome_graph[0]
        while startnode != -1:
            cycle += startnode
            genome_graph.remove(startnode)
            startnode = move(startnode, genome_graph)
        cycles.append(cycle)

    for cycle in cycles:
        chromosome = cycle_to_chromosome([cycle[-1]] + cycle[:-1])
        genome.append(chromosome)
    return genome


def move(startnode, edges, i=0):
    if len(edges) == 0:
        return -1
    elif len(edges) != 0:
        while not (startnode[1] + 1 == edges[i][0] or startnode[1] - 1 == edges[i][0]):
            i += 1
            if i == len(edges):
                return -1
    return edges[i]


def colored_edges(genome):
    """Outputs edges between genomes' synteny blocks"""

    edges = []
    for chromosome in genome:
        nodes = chromosome_to_cycle(chromosome)
        for j in range(1, abs(len(nodes)), 2):
            if j != len(nodes) - 1:
                edges.append([nodes[j], nodes[j + 1]])
            else:
                edges.append([nodes[j], nodes[0]])
    return edges


def break_2_on_genome_graph(graph, i1, i2, i3, i4):
    if type(graph) == str:
        graph = graph.split(', (')
        graph = [(el.split(', ')[0], el.split(', ')[1][:-1]) for el in graph]
        graph = [(int(val[0]), int(val[1])) for val in graph]

    for pair in ((i1, i2), (i3, i4)):
        if pair in graph:
            graph.remove(pair)
        elif pair[::-1] in graph:
            graph.remove(pair[::-1])
    for pair in ((i1, i3), (i2, i4)):
        graph.append((pair[0], pair[1]))

    return graph


def break_2_on_genome(p, i1, i2, i3, i4):
    genome_graph = colored_edges(p)
    genome_graph = break_2_on_genome_graph(genome_graph, i1, i2, i3, i4)
    p = graph_to_genome(genome_graph)
    return p

def dist_2_break(p,q):
    '''
    The 2-break distance between genomes P and Q is equal to Blocks(P, Q) - Cycles(P, Q).
    Blocks(P, Q) = Cycles(Q, Q)
    '''
    colored_edges_P = colored_edges(p)
    colored_edges_Q = colored_edges(q)
    colored_edges_all = colored_edges_P + colored_edges_Q
    colored_edges_all = [list(val) for val in colored_edges_all]
    n_cycle = 0
    while colored_edges_all:
        n_cycle += 1
        current_cycle = colored_edges_all[0]
        colored_edges_all.pop(0)
        for number in current_cycle:
            for edge in colored_edges_all:
                if number in edge:
                    edge.remove(number)
                    current_cycle.extend(edge)
                    colored_edges_all.remove(edge)

    return len(colored_edges_Q) - n_cycle

def shared_kmers(k, dna1, dna2):
    dnadict = {}
    for i in range(len(dna1) - k + 1):
        if dna1[i:i+k] in dnadict:
            dnadict[dna1[i:i+k]].append(i)
        else:
            dnadict[dna1[i:i+k]] = [i]
    dn = []
    reverse_comp = str(Seq(dna2).reverse_complement())

    for j in range(len(dna2) - k + 1):
        if dna2[j:j+k] in dnadict:
            for seq in dnadict[dna2[j:j+k]]:
                dn.append([seq, j])
        elif str(Seq(dna2[j:j+k]).reverse_complement()) in dnadict:
            for seq in dnadict[str(Seq(dna2[j:j+k]).reverse_complement())]:
                dn.append([seq, j])

    return dn

def read_dna(file_path):
    with open(file_path) as file:
        for line in file:
            return line

if __name__ == "__main__":
    k = 15
    dna1 = read_dna('data/dna1.txt')
    dna2 = read_dna('data/dna2.txt')
    sharedkmers = shared_kmers(k, dna1, dna2)
    for lines in sharedkmers:
        print('(' + ', '.join(map(str, lines)) + ')')
