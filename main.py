from de_Bruijn_graph import De_Bruijn_graph as DBG
import pyfastx


def save_contigs(contigs):
    """
    Save contigs in ~ .fa format
    :param contigs: list of strings
    :return: None
    """
    with open('contigs/output_contigs.fa', 'w') as file:
        for i, c in enumerate(contigs):
            file.write(f'> contig {i + 1}\n' + c + '\n')


def invert_sequence(s):
    """
    Makes complimentary inversion of the string
    :param s: input sequence
    :return: inverted sequence
    """
    d_inv = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
    return ''.join([d_inv[s[i]] for i in range(len(s) - 1, 0, -1)])


if __name__ == "__main__":
    dbg = DBG(k=92)

    lfq = pyfastx.Fastq('data/Carsonella_ruddii_reads_paired_reads_left.fastq')
    rfq = pyfastx.Fastq('data/Carsonella_ruddii_reads_paired_reads_right.fastq')

    print('Graph is being created:\n')
    dbg.update_graph(lfq)
    print('\nAssembling contigs:\n')
    contigs = dbg.get_contigs()
    save_contigs(contigs)

