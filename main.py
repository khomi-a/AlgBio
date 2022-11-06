def construct_dB_grap():
    pass


def save_contigs():
    pass


def invert_sequence(s):
    """
    Makes complimentary inversion of the string
    :param s: input sequence
    :return: inverted sequence
    """
    d_inv = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
    return ''.join([d_inv[s[i]] for i in range(len(s) - 1, 0, -1)])


if __name__ == "__main__":
    construct_dB_grap()
    save_contigs()
