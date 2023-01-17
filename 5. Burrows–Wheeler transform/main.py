from tqdm import tqdm


def BWT(a):
    words = list(a)
    bwt = []

    for i in tqdm(range(len(words))):
        word = a[-1] + a[:-1]
        new = ''.join(word)
        a = new
        bwt.append(new)
        i += 1

    sort = sorted(bwt)

    output = []
    for i in tqdm(range(len(words))):
        element = sort[i]
        last = element[-1]
        i = i + 1
        output.append(last)
    return "".join(output)


def inverseBWT(bwt):
    l = len(bwt)
    output = [''] * l
    count = dict()
    shortcut = [0] * l
    byteStart = dict()
    for i in range(l):
        lastChar = bwt[i]
        currCount = count.get(lastChar, 0)
        shortcut[i] = currCount
        count[lastChar] = currCount + 1

    currIndex = 0
    firstCol = []
    for char, currCount in sorted(count.items(), key=lambda x: x[0]):
        firstCol += [char] * currCount
        byteStart[char] = currIndex
        currIndex += currCount

    currIndex = 0
    for i in range(l):
        output[l - i - 1] = firstCol[currIndex]
        currIndex = byteStart[bwt[currIndex]] + shortcut[currIndex]

    return ''.join(output)


def read_dna(file_path):
    s = ''
    with open(file_path) as file:
        for line in file:
            s += line
    return s


if __name__ == "__main__":
    # text = 'GCGTGCCTGGTCA$'
    # print(BWT(text)) - > ACTGGCT$TGCGGC
    # print(inverseBWT(BWT(text))) -> GCGTGCCTGGTCA$
    text = read_dna('Mycoplasma pneumoniae.txt')
    with open('BWT(mycoplasma pneumoniae).txt', 'w') as file:
        file.write(BWT(text))
    # print(BWT(text) == inverseBWT(BWT(text)))


##########################################
# EXTRA: suffix array bwt implementation #
##########################################
# def suffix_array(string):
#     return (list(sorted(range(len(string)), key=lambda i: string[i:])))
#
#
# def bwt_from_suffix(string, s_array=None):
#     if s_array is None:
#         s_array = suffix_array(string)
#     return ("".join(string[idx - 1] for idx in s_array))
#
#
# def generate_all(input_string, s_array=None, eos="$"):
#
#     input_string = "".join([input_string, eos])
#     if s_array is None:
#         s_array = suffix_array(input_string)
#     bwt = bwt_from_suffix(input_string, s_array)
#
#     return bwt
##########################################
