import numpy as np

dict = {
        'N': {'N': 0, 'A': 1, 'T': 2, 'G': 1, 'C': 3},
        'A': {'N': 1, 'A': 0, 'T': 1, 'G': 5, 'C': 1},
        'T': {'N': 2, 'A': 1, 'T': 0, 'G': 9, 'C': 1},
        'G': {'N': 1, 'A': 5, 'T': 9, 'G': 0, 'C': 1},
        'C': {'N': 3, 'A': 1, 'T': 1, 'G': 1, 'C': 0}
    }

"""
The base substitution cost table. 
dict['A']['C'] represents the cost of replacing an A base with a C base.
dict['N']'['C'] represents the cost of adding a C base
"""


def edit_distance(s1, s2):            # The edit distance of two sequence.
    global dict
    s1_match = []
    s2_match = []
    line = []

    distance = np.zeros([len(s1) + 1, len(s2) + 1])
    path = np.zeros((len(s1) + 1, len(s2) + 1))           # Record the path to output the alignment scheme
    distance[0][0] = 0
    for i in range(1, len(s1) + 1):
        distance[i][0] = distance[i-1][0] + dict['N'][s1[i-1]]

    for j in range(1, len(s2) + 1):
        distance[0][j] = distance[0][j-1] + dict['N'][s2[j-1]]

    for i in range(1, len(s1) + 1):
        path[i][0] = 3

    for j in range(1, len(s2) + 1):
        path[0][j] = 2

    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            distance[i][j] = min(distance[i-1][j-1] + dict[s1[i-1]][s2[j-1]],
                                 distance[i][j-1] + dict['N'][s2[j-1]],
                                 distance[i-1][j] + dict['N'][s1[i-1]])
            if distance[i][j] == distance[i-1][j-1] + dict[s1[i-1]][s2[j-1]]:
                path[i][j] = 1
                continue
            if distance[i][j] == distance[i][j-1] + dict['N'][s2[j-1]]:
                path[i][j] = 2
                continue
            if distance[i][j] == distance[i-1][j] + dict['N'][s1[i-1]]:
                path[i][j] = 3
                continue
    i = len(s1)
    j = len(s2)
    while i is not 0 or j is not 0:
        if path[i][j] == 1:
            s1_match.insert(0, s1[i-1])
            s2_match.insert(0, s2[j-1])
            line.insert(0, "|")
            i = i - 1
            j = j - 1
            continue
        if path[i][j] == 2:
            s1_match.insert(0, '-')
            s2_match.insert(0, s2[j-1])
            line.insert(0, ' ')
            j = j - 1
            continue
        if path[i][j] == 3:
            s1_match.insert(0, s1[i-1])
            s2_match.insert(0, '-')
            line.insert(0, ' ')
            i = i - 1
            continue

    return distance[len(s1)][len(s2)], s1_match, s2_match, line


def align(s1, s2):                  # Scrolling array to calculate x coordinate
    global dict
    previous = [0] * (len(s1) + 1)
    current = [0] * (len(s1) + 1)
    previous[0] = 0
    for i in range(1, len(previous)):
        previous[i] = previous[i-1] + dict['N'][s1[i-1]]
    for j in range(1, len(s2) + 1):
        current[0] = previous[0] + dict['N'][s2[j-1]]
        for i in range(1, len(current)):
            current[i] = min(previous[i-1] + dict[s1[i-1]][s2[j-1]],
                             previous[i] + dict['N'][s2[j-1]],
                             current[i-1] + dict['N'][s1[i-1]])
        previous = current.copy()
    return current


def getx(current1, current2):
    i_min = 0
    min_distance = float("inf")
    k = len(current1)
    for i in range(k):
        distance = current1[i] + current2[k-1-i]
        if distance < min_distance:
            min_distance = distance
            i_min = i
    return i_min


def Hirschberg(s1, s2, h1, h2, w1, w2):
    if (w2 - w1) == 1:
        result, s1_m, s2_m, line = edit_distance(s1[h1: h2], s2[w1: w2])
        return result, s1_m, s2_m, line
    y = int((w1 + w2) / 2)
    current1 = align(s1[h1: h2], s2[w1: y])
    current2 = align(s1[h1: h2][:: -1], s2[y: w2][:: -1])
    x = getx(current1, current2)
    dist_l, s1_match_l, s2_match_l, line_l = Hirschberg(s1, s2, h1, h1+x, w1, y)
    dist_r, s1_match_r, s2_match_r, line_r = Hirschberg(s1, s2, h1+x, h2, y, w2)
    return dist_l + dist_r, s1_match_l + s1_match_r, s2_match_l + s2_match_r, line_l + line_r


def Hirschberg_anl(s1, s2):
    w2 = len(s2)
    h2 = len(s1)
    if w2 == 0 and h2 == 0:
        distance = 0
        s1_match = []
        s2_match = []
        line = []

    if w2 == 0:
        distance = 0
        s2_match = []
        s1_match = ''.join(s1)
        line = []
        for i in range(h2):
            distance += dict['N'][s1[i]]
            s2_match.append('-')
            line.append(' ')

    else:
        distance, s1_match, s2_match, line = Hirschberg(s1, s2, 0, h2, 0, w2)

    print("The optimal alignment costs " + str(distance))
    print(''.join(s1_match))
    print(''.join(line))
    print(''.join(s2_match))


if __name__ == "__main__":
    f1 = open("DNA1", "r")
    f2 = open("DNA2", "r")
    s1 = list(f1.read())
    s2 = list(f2.read())
    Hirschberg_anl(s1, s2)


