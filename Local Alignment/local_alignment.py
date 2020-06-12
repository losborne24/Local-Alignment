def dynprog(alphabet, matrix, sequence_1, sequence_2):      # basic dynamic programming
    s1_length = len(sequence_1)
    s2_length = len(sequence_2)
    max_value = 0
    scoring_matrix = [[Item(0, []) for _ in range(s1_length+1)] for _ in range(s2_length+1)]
    max_loc = [0, 0]
    for i in range(1, s2_length + 1):
        for j in range(1, s1_length + 1):
            # compute scores
            diagonal_score = scoring_matrix[i-1][j-1].get_value() + substitution_lookup(alphabet, matrix, sequence_1[j-1], sequence_2[i - 1])
            left_score = scoring_matrix[i][j-1].get_value() + substitution_lookup(alphabet, matrix, sequence_1[j-1], '_')
            up_score = scoring_matrix[i-1][j].get_value() + substitution_lookup(alphabet, matrix, '_', sequence_2[i - 1])
            # select maximum
            maximum = max([diagonal_score, left_score, up_score, 0])
            scoring_matrix[i][j].set_value(maximum)
            direction = [i for i, j in enumerate([diagonal_score, left_score, up_score]) if j == maximum]
            if direction:
                scoring_matrix[i][j].set_direction(direction)
            if maximum >= max_value:
                max_value = maximum
                max_loc = [i, j]

    # align sequences
    temp_1 = []
    temp_2 = []
    alignment_1 = []
    alignment_2 = []
    x = max_loc[0]
    y = max_loc[1]
    while scoring_matrix[x][y].get_direction():
        if 0 in scoring_matrix[x][y].get_direction():
            temp_1.insert(0, sequence_1[y-1])
            temp_2.insert(0, sequence_2[x-1])
            alignment_1.insert(0, y-1)
            alignment_2.insert(0, x-1)
            x -= 1
            y -= 1
        elif 1 in scoring_matrix[x][y].get_direction():
            temp_1.insert(0, sequence_1[y - 1])
            temp_2.insert(0, "_")
            y -= 1
        elif 2 in scoring_matrix[x][y].get_direction():
            temp_1.insert(0, "_")
            temp_2.insert(0, sequence_2[x - 1])
            x -= 1
        else:
            break

    return [max_value, alignment_1, alignment_2]


def dynproglin(alphabet, matrix, sequence_1, sequence_2):   # Dynamic programming in linear space
    s1_length = len(sequence_1)
    s2_length = len(sequence_2)
    forward_max_value = 0
    forward_max_loc = [0, 0]
    current_row = [0] * (s1_length + 1)
    # forward pass
    for i in range(1, s2_length+1):
        diagonal_score = substitution_lookup(alphabet, matrix, sequence_1[0], sequence_2[i-1])
        for j in range(1, s1_length+1):
            left_score = current_row[j-1] + substitution_lookup(alphabet, matrix, sequence_1[j-1], '_')
            up_score = current_row[j] + substitution_lookup(alphabet, matrix, '_', sequence_2[i - 1])
            # select maximum
            maximum = max([diagonal_score, left_score, up_score, 0])
            if j < s1_length:
                diagonal_score = current_row[j] + substitution_lookup(alphabet, matrix, sequence_1[j], sequence_2[i - 1])
            current_row[j] = maximum
            if maximum >= forward_max_value:
                forward_max_value = maximum
                forward_max_loc = [i, j]
    # backwards pass
    reverse_max_value = 0
    reverse_max_loc = [0, 0]
    current_row = [0] * (s1_length +1)
    for i in range(s2_length-1, -1, -1):
        diagonal_score = substitution_lookup(alphabet, matrix, sequence_1[s1_length-1], sequence_2[i])
        for j in range(s1_length - 1, -1, -1):
            right_score = current_row[j + 1] + substitution_lookup(alphabet, matrix, sequence_1[j], '_')
            down_score = current_row[j] + substitution_lookup(alphabet, matrix, '_', sequence_2[i])
            maximum = max([diagonal_score, right_score, down_score, 0])
            if j > 0:
                diagonal_score = current_row[j] + substitution_lookup(alphabet, matrix, sequence_1[j-1], sequence_2[i])
            current_row[j] = maximum
            if maximum >= reverse_max_value:
                reverse_max_value = maximum
                reverse_max_loc = [i, j]
    # align sequences
    subsequence_1 = sequence_1[reverse_max_loc[1]:forward_max_loc[1]]
    subsequence_2 = sequence_2[reverse_max_loc[0]:forward_max_loc[0]]
    alignment = global_alignment_recursion(alphabet, matrix, subsequence_1[::-1], subsequence_2[::-1], 0, len(subsequence_1), 0, len(subsequence_2), ["", ""])
    alignment[0]= alignment[0][::-1]
    alignment[1] = alignment[1][::-1]
    alignment_1 = [reverse_max_loc[1]]
    alignment_2 = [reverse_max_loc[0]]
    pointer_1 = reverse_max_loc[1] + 1
    pointer_2 = reverse_max_loc[0] + 1
    for i in range(1, len(alignment[0])):
        if alignment[0][i] == "_" and not alignment[1][i] == "_":
            pointer_1 += 1
        elif not alignment[0][i] == "_" and alignment[1][i] == "_":
            pointer_2 += 1
        else:
            alignment_1.append(pointer_1)
            alignment_2.append(pointer_2)
            pointer_1 += 1
            pointer_2 += 1

    return [forward_max_value, alignment_1, alignment_2]


def global_alignment_recursion(alphabet, matrix, sequence_1, subsequence_2, s1_start, s1_end, s2_start, s2_end, alignment):
    subsequence_1 = sequence_1[s1_start:s1_end]
    midpoint = s2_start + int(((s2_end-s2_start) + 1) / 2)
    s_2a = subsequence_2[s2_start:midpoint]
    s_2b = subsequence_2[midpoint: s2_end]
    if not subsequence_1:
        string_length = len(s_2a) + len(s_2b)
        blank_string = "_" * string_length
        return [str(alignment[0]) + s_2a + s_2b, str(alignment[1]) + blank_string]
    elif not s_2b and len(subsequence_1) == 1:
        return [str(alignment[0]) + s_2a + s_2b, str(alignment[1]) + subsequence_1]

    prefix = global_alignment(alphabet,matrix, subsequence_1, s_2a, "forward")
    suffix = global_alignment(alphabet,matrix, subsequence_1, s_2b, "reverse")
    length = [prefix[i] + suffix[i] for i in range(len(prefix))]
    max_loc = length.index(max(length)) + s1_start
    if not s_2b:
        blank_string_1 = "_" * (length.index(max(length)) - 1)
        blank_string_2 = "_" * (len(subsequence_1) - length.index(max(length)))
        return [str(alignment[0]) + blank_string_1 + s_2a + blank_string_2, str(alignment[1]) + subsequence_1]
    alignment = global_alignment_recursion(alphabet, matrix, sequence_1, subsequence_2, s1_start, max_loc, s2_start, midpoint, alignment)
    alignment = global_alignment_recursion(alphabet, matrix, sequence_1, subsequence_2, max_loc, s1_end, midpoint, s2_end, alignment)
    return alignment


def global_alignment(alphabet, matrix, sequence_1, sequence_2, direction):
    s1_length = len(sequence_1)
    s2_length = len(sequence_2)
    if direction == "forward":
        current_row = [0] * (s1_length + 1)
        for i in range(1, s1_length+1):
            current_row[i] = current_row[i - 1] + substitution_lookup(alphabet, matrix, sequence_1[i - 1], '_')
        for i in range(1, s2_length + 1):
            diagonal_score = current_row[0] + substitution_lookup(alphabet, matrix, sequence_1[0], sequence_2[i - 1])
            current_row[0] = current_row[0] + substitution_lookup(alphabet, matrix, '_', sequence_2[i - 1])
            for j in range(1, s1_length + 1):
                left_score = current_row[j - 1] + substitution_lookup(alphabet, matrix, sequence_1[j - 1], '_')
                up_score = current_row[j] + substitution_lookup(alphabet, matrix, '_', sequence_2[i - 1])
                maximum = max([diagonal_score, left_score, up_score])
                if j < s1_length:
                    diagonal_score = current_row[j] + substitution_lookup(alphabet, matrix, sequence_1[j],
                                                                          sequence_2[i - 1])
                current_row[j] = maximum
        return current_row
    else:
        current_row = [0] * (s1_length + 1)
        for i in range(s1_length - 1, -1, -1):
            current_row[i] = current_row[i + 1] + substitution_lookup(alphabet, matrix, sequence_1[i], '_')
        for i in range(s2_length - 1, -1, -1):
            diagonal_score = current_row[s1_length]+ substitution_lookup(alphabet, matrix, sequence_1[s1_length - 1], sequence_2[i])
            current_row[s1_length] = current_row[s1_length] + substitution_lookup(alphabet, matrix, '_', sequence_2[i])
            for j in range(s1_length - 1, -1, -1):
                right_score = current_row[j + 1] + substitution_lookup(alphabet, matrix, sequence_1[j], '_')
                down_score = current_row[j] + substitution_lookup(alphabet, matrix, '_', sequence_2[i])
                maximum = max([diagonal_score, right_score, down_score])
                if j > 0:
                    diagonal_score = current_row[j] + substitution_lookup(alphabet, matrix, sequence_1[j - 1],
                                                                          sequence_2[i])
                current_row[j] = maximum
        return current_row


def heuralign(alphabet, matrix, sequence_1, sequence_2):
    is_swap = False     # set sequence 1 to be the longer sequence
    if len(sequence_2) > len(sequence_1):
        is_swap = True
        swap = sequence_1
        sequence_1 = sequence_2
        sequence_2 = swap

    max_score = 0
    for i in range(0, len(matrix)):
        for j in range(i, len(matrix)):
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]

    found = False
    window = 2
    couples_1 = {}
    couples_2 = {}
    while not found:    # try window size two. If no matches, try window size one
        couples_1 = {}
        couples_2 = {}
        if len(sequence_1) > (window - 1) and len(sequence_2)>(window - 1):
            for i in range(0, len(sequence_1) - (window - 1)):
                if sequence_1[i:i+window] in couples_1:
                    item = couples_1.get(sequence_1[i:i+window])
                    item.append(i)
                    couples_1[sequence_1[i:i+window]] = item
                elif window == 2:
                    if substitution_lookup(alphabet, matrix, sequence_1[i], sequence_1[i]) + substitution_lookup(alphabet, matrix, sequence_1[i+1], sequence_1[i+1]) >= max_score:
                        couples_1[sequence_1[i:i + window]] = [i]
                else:
                    couples_1[sequence_1[i:i + window]] = [i]
            for i in range(0, len(sequence_2) - (window-1)):
                if sequence_2[i:i + window] in couples_1:
                    if sequence_2[i:i + window] in couples_2:
                        item = couples_2.get(sequence_2[i:i + window])
                        item.append(i)
                        couples_2[sequence_2[i:i + window]] = item
                    else:
                        couples_2[sequence_2[i:i + window]] = [i]
            window -= 1
        if couples_2:
            found = True
        elif window == 0:
            return None
    difference = {}
    width = int((len(sequence_1)) / 2) + 1
    for i in range((- len(sequence_2)) + (window - 1), len(sequence_1) + 1 - (window-1)):
        difference[i] = 0
    for key in couples_2:
        for i in range(0, len(couples_2.get(key))):
            for j in range(0, len(couples_1.get(key))):
                diff = couples_1.get(key)[j] -couples_2.get(key)[i]
                for k in range(0, window + 1):
                    difference[diff] += substitution_lookup(alphabet, matrix, key[k], key[k])
    dict_maximum = 0
    for i in range(1, width*2):
        dict_maximum += difference[(-len(sequence_2)) + i]
    dict_max_loc = -len(sequence_2)
    score = dict_maximum
    for i in range((-len(sequence_2)) + 2, len(sequence_1) - width):
        score = score + difference[i + width] - difference[i-1]
        if score > dict_maximum:
            dict_max_loc = i
            dict_maximum = score
    max_loc = [0, 0]
    max_value = 0
    scoring_matrix = [[Item(0, []) for _ in range((width * 2)+1)] for _ in range(len(sequence_2)+1)]
    for i in range(1, len(sequence_2) + 1):
        prev_score = 0
        if ((width * 2)+i+dict_max_loc -1) >= 1:
            for j in range(0, (width * 2) + 1):
                if not (j+i+ dict_max_loc - 1) < 1 and not (j+i+ dict_max_loc - 1) > len(sequence_1):
                    diagonal_score = scoring_matrix[i-1][j].get_value() + substitution_lookup(alphabet, matrix,sequence_1[j+i+ dict_max_loc - 2], sequence_2[i-1])
                    left_score = prev_score + substitution_lookup(alphabet, matrix,sequence_1[j+i+ dict_max_loc -2 ], '_')
                    if j == (width * 2):
                        up = 0
                    else:
                        up = scoring_matrix[i - 1][j + 1].get_value()
                    up_score = up + substitution_lookup(alphabet, matrix, '_', sequence_2[i-1])
                    maximum = max([diagonal_score, left_score, up_score, 0])
                    scoring_matrix[i][j].set_value(maximum)
                    prev_score = maximum
                    direction = [i for i, j in enumerate([diagonal_score, left_score, up_score]) if j == maximum]
                    if direction:
                        scoring_matrix[i][j].set_direction(direction)
                    if maximum >= max_value:
                        max_value = maximum
                        max_loc = [i, j]
    alignment_1 = []
    alignment_2 = []
    x = max_loc[0]
    y = max_loc[1]
    while scoring_matrix[x][y].get_direction():
        if 0 in scoring_matrix[x][y].get_direction():
            alignment_1.insert(0, y + x + dict_max_loc - 2)
            alignment_2.insert(0, x - 1)
            x -= 1
        elif 1 in scoring_matrix[x][y].get_direction():
            y -= 1
        elif 2 in scoring_matrix[x][y].get_direction():
            x -= 1
            y += 1
        else:
            break
    if is_swap:
        alignment = alignment_1
        alignment_1 = alignment_2
        alignment_2 = alignment

    return [max_value, alignment_1, alignment_2]


def substitution_lookup(alphabet, matrix, char_1, char_2):
    if char_1 == '_':
        x = len(alphabet)
    else:
        x = alphabet.find(char_1)
    if char_2 == '_':
        y = len(alphabet)
    else:
        y = alphabet.find(char_2)
    return matrix[x][y]


class Item:
    def __init__(self, value, direction):
        self.value = value
        self.direction = direction

    def set_value(self, value):
        self.value = value

    def get_value(self):
        return self.value

    def set_direction(self, direction):
        self.direction = direction

    def get_direction(self):
        return self.direction

