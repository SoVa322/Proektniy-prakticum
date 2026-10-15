#  Задача 3. Локальное выравнивание Смита–Ватермана
#  с аффинными гэпами (BLOSUM62, gap_open=11, gap_extend=1)

import numpy as np


# Парсинг файла blosum62.txt
def parse_blosum62(filepath: str) -> dict:
    
    with open(filepath, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
    header = lines[0].split()
    blosum62 = {}
    for line in lines[1:]:
        parts = line.split()
        row_aa = parts[0]
        values = [int(x) for x in parts[1:]]
        blosum62[row_aa] = dict(zip(header, values))
    return blosum62

# Инициализация трёх матриц

def init_matrices(n: int, m: int):
    M  = np.zeros((n + 1, m + 1), dtype=np.int32)
    Ix = np.zeros((n + 1, m + 1), dtype=np.int32)
    Iy = np.zeros((n + 1, m + 1), dtype=np.int32)
    return M, Ix, Iy

# Заполнение матриц
def fill_matrices(seq1, seq2, blosum62, gap_open=11, gap_extend=1):
    n, m = len(seq1), len(seq2)
    M, Ix, Iy = init_matrices(n, m)
    max_score = 0
    max_positions = []
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            aa1 = seq1[i - 1]
            aa2 = seq2[j - 1]
            s = blosum62[aa1][aa2]
            Ix[i][j] = max(
                M[i - 1][j] - gap_open,
                Ix[i - 1][j] - gap_extend
            )
            Iy[i][j] = max(
                M[i][j - 1] - gap_open,     
                Iy[i][j - 1] - gap_extend   
            )
            M[i][j] = max(
                0,
                M[i - 1][j - 1],
                Ix[i - 1][j - 1],
                Iy[i - 1][j - 1]
            ) + s

# Отслеживание максимума
            best_here = max(M[i][j], Ix[i][j], Iy[i][j])
            if best_here > max_score:
                max_score = best_here
                max_positions = [(i, j)]
            elif best_here == max_score and best_here > 0:
                max_positions.append((i, j))

    return M, Ix, Iy, max_score, max_positions

# Обратный ход
def traceback(M, Ix, Iy, seq1, seq2, blosum62,
              i, j, gap_open=11, gap_extend=1):
    aligned1, aligned2 = [], []
    ambiguous = False
    # Определяем из какой матрицы стартуем
    if M[i][j] >= Ix[i][j] and M[i][j] >= Iy[i][j]:
        state = "M"
    elif Ix[i][j] >= Iy[i][j]:
        state = "Ix"
    else:
        state = "Iy"

    while i > 0 and j > 0:
        if state == "M":
            if M[i][j] == 0:
                break
            s = blosum62[seq1[i - 1]][seq2[j - 1]]
            
            candidates = []
            if M[i][j] == M[i - 1][j - 1] + s:
                candidates.append("M")
            if M[i][j] == Ix[i - 1][j - 1] + s:
                candidates.append("Ix")
            if M[i][j] == Iy[i - 1][j - 1] + s:
                candidates.append("Iy")
            if M[i][j] == s and (M[i-1][j-1] <= 0 and
                                 Ix[i-1][j-1] <= 0 and
                                 Iy[i-1][j-1] <= 0):
                candidates.append("START")

            if len(candidates) > 1:
                ambiguous = True

            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])

            prev = candidates[0]
            i -= 1
            j -= 1
            if prev == "START":
                break
            state = prev
        elif state == "Ix":
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            if Ix[i][j] == M[i - 1][j] - gap_open:
                i -= 1
                state = "M"
            elif Ix[i][j] == Ix[i - 1][j] - gap_extend:
                i -= 1
            else:
                break
        elif state == "Iy":
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            if Iy[i][j] == M[i][j - 1] - gap_open:
                j -= 1
                state = "M"
            elif Iy[i][j] == Iy[i][j - 1] - gap_extend:
                j -= 1
            else:
                break
    aligned1 = "".join(reversed(aligned1))
    aligned2 = "".join(reversed(aligned2))
    start1 = i + 1
    start2 = j + 1

    return aligned1, aligned2, start1, start2, ambiguous


# Основная функция
def smith_waterman_affine(seq1, seq2, blosum62,
                          gap_open=11, gap_extend=1):
    M, Ix, Iy, max_score, max_positions = fill_matrices(
        seq1, seq2, blosum62, gap_open, gap_extend
    )
    end_i, end_j = max_positions[0]
    aligned1, aligned2, start1, start2, ambiguous = traceback(
        M, Ix, Iy, seq1, seq2, blosum62,
        end_i, end_j, gap_open, gap_extend
    )
    is_unique = (len(max_positions) == 1) and (not ambiguous)

    return (max_score,
            aligned1, aligned2,
            start1, end_i, start2, end_j,
            is_unique)

# Пример использования
if __name__ == "__main__":
    # Вариант A - blosum62 из файла
    # blosum62 = parse_blosum62("blosum62.txt")

    # Вариант B - словарь прямо в коде
    amino_acids = "ARNDCQEGHILKMFPSTWYV"
    blosum62 = {
        'A': [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
        'R': [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
        'N': [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
        'D': [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
        'C': [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
        'Q': [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
        'E': [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
        'G': [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
        'H': [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
        'I': [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
        'L': [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
        'K': [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
        'M': [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
        'F': [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
        'P': [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
        'S': [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
        'T': [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
        'W': [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
        'Y': [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
        'V': [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
    }
    blosum62_dict = {
        aa: dict(zip(amino_acids, blosum62[aa])) for aa in amino_acids
    }
    seq1 = "HEAGAWGHEE"
    seq2 = "PAWHEAE"
    (score, a1, a2, s1, e1, s2, e2, uniq) = smith_waterman_affine(
        seq1, seq2, blosum62_dict
    )

    print(f"alignment_score = {score}")
    print(f"aligned_seq1    = {a1}")
    print(f"aligned_seq2    = {a2}")
    print(f"start1, end1    = {s1}, {e1}")
    print(f"start2, end2    = {s2}, {e2}")
    print(f"is_unique       = {uniq}")