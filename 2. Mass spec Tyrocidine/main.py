import numpy as np

spectrum_unknown =np.array([371.5, 375.4, 390.4, 392.2, 409.0, 420.2, 427.2, 443.3, 446.4, 461.3,
    471.4, 477.4, 491.3, 505.3, 506.4, 519.2, 536.1, 546.5, 553.3, 562.3,
    588.2, 600.3, 616.2, 617.4, 618.3, 633.4, 634.4, 636.2, 651.5, 652.4,
    702.5, 703.4, 712.5, 718.3,721.0, 730.3, 749.4, 762.6, 763.4, 764.4,
    779.6, 780.4,781.4, 782.4, 797.3, 862.4, 876.4, 877.4, 878.6, 879.4,
    893.4, 894.4, 895.4, 896.5, 927.4, 944.4, 975.5, 976.5, 977.4, 979.4,
    1005.5, 1007.5, 1022.5, 1023.7, 1024.5, 1039.5, 1040.3, 1042.5, 1043.4, 1057.5,
    1119.6, 1120.6, 1137.6, 1138.6, 1139.5, 1156.5, 1157.6, 1168.6, 1171.6, 1185.4,
    1220.6, 1222.5, 1223.6, 1239.6, 1240.6, 1250.5, 1256.5, 1266.5, 1267.5, 1268.6])


aminoacids = {'Ala': 71.0, 'Cys':103.0, 'Asp': 115.0, 'Glu': 129.0, 'Phe': 147.1, 'Gly': 57.0, 'His': 137.1, 'Ile': 113.1,
    'Lys':128.1, 'Leu':113.1, 'Met': 131.0, 'Asn': 114.0, 'Pro': 97.1, 'Gln': 128.1, 'Arg': 156.1, 'Ser': 87.0,
    'Thr': 101.1, 'Val':99.1, 'Trp': 186.1, 'Tyr': 163.1}


def get_diffs(spectrum, M=20):
    '''
    :param spectrum: -
    :param M: desired number of aa
    :return: convolution for all positive differences in spectrum
    '''
    spectrum = [0] + sorted(spectrum)
    diffs = {}
    n = len(spectrum)
    for i in range(n-1):
        for j in range(i+1, n):
            diff = round(spectrum[j] - spectrum[i],0)
            if 57 <= diff <= 200:
                if diff in diffs:
                    diffs[diff] += 1
                else:
                    diffs[diff] = 1
    diffs = [p[0] for p in sorted(diffs.items(), key=lambda item: item[1], reverse=True)[:M]]
    return diffs


def get_possible_aa(spectrum, m=20):
    '''
    :param spectrum: -
    :param m: -
    :return: most likely aminoacids included in protein
    '''
    possible_aa = []
    diffs = get_diffs(spectrum, m)
    for aa in aminoacids.keys():
        for mass in diffs:
            if abs(aminoacids[aa] - mass) <= 1:
                possible_aa.append(aa)
    return np.unique(sorted(possible_aa))


def get_spectrum(peptide, aa_masses = aminoacids):
    '''
    :param peptide: -
    :param aa_masses: -
    :return: theoretical spectrum of a given peptide
    '''
    prefix_mass = [0]
    prefix_mass_last = 0
    for aa in peptide:
        prefix_mass_last += aa_masses[aa]
        prefix_mass.append(prefix_mass_last)
    # peptide_mass = prefix_mass[-1]

    n = len(prefix_mass)
    cyclic_spectrum = []
    for i in range(n-1):
        for j in range(i+1,n):
            cyclic_spectrum.append(round(prefix_mass[j] - prefix_mass[i],0))
            if i > 0 and j < n:
                cyclic_spectrum.append(round(prefix_mass[-1] -(prefix_mass[j] - prefix_mass[i]),0))
    return sorted(np.unique(cyclic_spectrum))


def score(spectrum, peptide):
    '''
    :param spectrum: experimental spectrum
    :param peptide: theoretical peptide
    :return: scoring result
    '''
    theor = set(get_spectrum(peptide))
    exp = set(round(el, 0) for el in spectrum)
    return len(set.intersection(theor,exp))

def expand_leaderboard(board, possible_aa):
    new_board = []
    for peptide in board:
        for aa in possible_aa:
            new_board.append(peptide + [aa])
    return new_board

def trim_leaderboard(board, spectrum, N):
    return sorted(board, key=lambda item: score(spectrum,item), reverse=True)[:N]

def get_mass(peptide, masses = aminoacids):
    '''
    :param peptide: -
    :param masses: -
    :return: mass of a given peptide
    '''
    m = 0
    for aa in peptide:
        m += masses[aa]
    return m

def leaderboard_cyclopeptide_sequencing(spectrum, N=1000, M = 20,):
    '''
    Basic algorithm that finds peptide for a given spectrum
    :param spectrum: -
    :param N: number of holded possible peptides in leaderboard
    :param M: number of amino acids to choose from spectrum
    :return: resulting peptide
    '''
    aa_in_use = get_possible_aa(spectrum, m=M)
    leaderboard = [[]]
    leader_peptide = []
    k = 0
    while len(leaderboard) > 0:
        leaderboard = expand_leaderboard(leaderboard, aa_in_use)
        for peptide in leaderboard[:]:
            if abs(get_mass(peptide) - max(spectrum)) < 0.5:
                if score(spectrum, peptide) > score(spectrum, leader_peptide):
                    leader_peptide = peptide
            if get_mass(peptide) > max(spectrum):
                leaderboard.remove(peptide)

        leaderboard = trim_leaderboard(leaderboard, spectrum, N)

    return concatenate(leader_peptide) #, score(spectrum, leader_peptide) / len(spectrum)

def concatenate(arr_aa):
    return '-'.join(arr_aa)


if __name__ == "__main__":
    print(leaderboard_cyclopeptide_sequencing(spectrum_unknown, 1000, 12))
    # Phe-Pro-Gln-Val-Tyr-Glu-Asn-Met-Met-Glu
    
