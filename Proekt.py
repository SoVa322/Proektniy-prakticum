import os
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import pandas as pd
import numpy as np

# Алфавит аминокислот
AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV"

# --- 1. ПАРСИНГ ФАЙЛОВ ---
def parse_scoring_matrix(filepath: str) -> dict:
    ext = os.path.splitext(filepath)[1].lower()
    if ext in ['.xls', '.xlsx']:
        df = pd.read_excel(filepath, index_col=0)
    else:
        df = pd.read_csv(filepath, sep=None, engine='python', index_col=0)
        
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    return df.to_dict(orient='index')

def load_sequences_from_file(filepath: str):
    ext = os.path.splitext(filepath)[1].lower()
    if ext == '.csv':
        df = pd.read_csv(filepath)
        return str(df.iloc[0, 0]).strip(), str(df.iloc[1, 0]).strip()
    elif ext in ['.xls', '.xlsx']:
        df = pd.read_excel(filepath)
        return str(df.iloc[0, 0]).strip(), str(df.iloc[1, 0]).strip()
    else:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read().strip()
        if content.startswith(">"):
            blocks = [b.strip() for b in content.split(">") if b.strip()]
            seqs = []
            for b in blocks[:2]:
                lines = b.splitlines()
                seqs.append("".join(line.strip() for line in lines[1:]))
            if len(seqs) >= 2:
                return seqs[0], seqs[1]
        lines = [line.strip() for line in content.splitlines() if line.strip()]
        if len(lines) >= 2:
            return lines[0], lines[1]
        raise ValueError("В файле должно быть минимум две последовательности.")

# --- 2. АЛГОРИТМ СМИТА-ВАТЕРМАНА ---
def init_matrices(n: int, m: int):
    M  = np.zeros((n + 1, m + 1), dtype=np.int32)
    Ix = np.zeros((n + 1, m + 1), dtype=np.int32)
    Iy = np.zeros((n + 1, m + 1), dtype=np.int32)
    return M, Ix, Iy

def fill_matrices(seq1, seq2, scoring_matrix, gap_open=11, gap_extend=1):
    n, m = len(seq1), len(seq2)
    M, Ix, Iy = init_matrices(n, m)
    max_score = 0
    max_positions = []
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            aa1 = seq1[i - 1]
            aa2 = seq2[j - 1]
            s = scoring_matrix[aa1][aa2]
            
            Ix[i][j] = max(M[i - 1][j] - gap_open, Ix[i - 1][j] - gap_extend)
            Iy[i][j] = max(M[i][j - 1] - gap_open, Iy[i][j - 1] - gap_extend)
            M[i][j] = max(0, M[i - 1][j - 1] + s, Ix[i - 1][j - 1] + s, Iy[i - 1][j - 1] + s)

            best_here = max(M[i][j], Ix[i][j], Iy[i][j])
            if best_here > max_score:
                max_score = best_here
                max_positions = [(i, j)]
            elif best_here == max_score and best_here > 0:
                max_positions.append((i, j))

    return M, Ix, Iy, max_score, max_positions

def traceback(M, Ix, Iy, seq1, seq2, scoring_matrix, i, j, gap_open=11, gap_extend=1):
    aligned1, aligned2 = [], []
    ambiguous = False
    
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
            s = scoring_matrix[seq1[i - 1]][seq2[j - 1]]
            candidates = []
            if M[i][j] == M[i - 1][j - 1] + s:
                candidates.append("M")
            if M[i][j] == Ix[i - 1][j - 1] + s:
                candidates.append("Ix")
            if M[i][j] == Iy[i - 1][j - 1] + s:
                candidates.append("Iy")
            if M[i][j] == s and (M[i-1][j-1] <= 0 and Ix[i-1][j-1] <= 0 and Iy[i-1][j-1] <= 0):
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

def smith_waterman_affine(seq1, seq2, scoring_matrix, gap_open=11, gap_extend=1):
    M, Ix, Iy, max_score, max_positions = fill_matrices(seq1, seq2, scoring_matrix, gap_open, gap_extend)
    if not max_positions:
        return 0, "", "", 0, 0, 0, 0, True
        
    end_i, end_j = max_positions[0]
    aligned1, aligned2, start1, start2, ambiguous = traceback(
        M, Ix, Iy, seq1, seq2, scoring_matrix, end_i, end_j, gap_open, gap_extend
    )
    is_unique = (len(max_positions) == 1) and (not ambiguous)
    return int(max_score), aligned1, aligned2, start1, end_i, start2, end_j, is_unique

# --- 3. ГРАФИЧЕСКАЯ ОБОЛОЧКА ---
class AlignmentApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Локальное выравнивание Смита–Ватермана")
        self.root.geometry("860x550")
        
        self.scoring_matrix = None

        top_frame = tk.Frame(root, pady=10, padx=10)
        top_frame.pack(fill=tk.X)

        self.matrix_lbl = tk.Label(top_frame, text="Матрица замен: Не загружена", fg="red", font=("Arial", 9, "bold"))
        self.matrix_lbl.pack(side=tk.LEFT, padx=5)

        tk.Button(top_frame, text="Загрузить матрицу замен (TXT/CSV/XLSX)", command=self.load_matrix).pack(side=tk.LEFT, padx=5)
        tk.Button(top_frame, text="Загрузить seq из файла", command=self.load_seqs_file).pack(side=tk.RIGHT, padx=5)

        seq_frame = tk.LabelFrame(root, text="Входные аминокислотные последовательности", padx=10, pady=8)
        seq_frame.pack(fill=tk.X, padx=10, pady=5)

        tk.Label(seq_frame, text="Последовательность 1 (seq1):").pack(anchor=tk.W)
        self.seq1_entry = tk.Text(seq_frame, height=3, font=("Courier", 10))
        self.seq1_entry.pack(fill=tk.X, pady=3)

        tk.Label(seq_frame, text="Последовательность 2 (seq2):").pack(anchor=tk.W)
        self.seq2_entry = tk.Text(seq_frame, height=3, font=("Courier", 10))
        self.seq2_entry.pack(fill=tk.X, pady=3)

        param_frame = tk.Frame(root, pady=8, padx=10)
        param_frame.pack(fill=tk.X)

        tk.Label(param_frame, text="Штраф за открытие гэпа (gap_open):").pack(side=tk.LEFT)
        self.gap_open_var = tk.StringVar(value="11")
        tk.Entry(param_frame, textvariable=self.gap_open_var, width=5).pack(side=tk.LEFT, padx=5)

        tk.Label(param_frame, text="Штраф за продолжение (gap_extend):").pack(side=tk.LEFT, padx=(15, 0))
        self.gap_ext_var = tk.StringVar(value="1")
        tk.Entry(param_frame, textvariable=self.gap_ext_var, width=5).pack(side=tk.LEFT, padx=5)

        tk.Button(param_frame, text="РАССЧИТАТЬ ВЫРАВНИВАНИЕ", bg="#0052cc", fg="white",
                  font=("Arial", 10, "bold"), padx=12, pady=4, command=self.run_alignment).pack(side=tk.RIGHT)

        res_frame = tk.LabelFrame(root, text="Результаты выравнивания", padx=10, pady=8)
        res_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        self.result_text = scrolledtext.ScrolledText(res_frame, font=("Courier", 10))
        self.result_text.pack(fill=tk.BOTH, expand=True)

    def load_matrix(self):
        filepath = filedialog.askopenfilename(
            title="Выберите файл матрицы замен",
            filetypes=[("Supported files", "*.txt *.csv *.xlsx *.xls"), ("All files", "*.*")]
        )
        if filepath:
            try:
                self.scoring_matrix = parse_scoring_matrix(filepath)
                self.matrix_lbl.config(text=f"Матрица замен: {os.path.basename(filepath)}", fg="green")
                messagebox.showinfo("Успех", "Матрица замен успешно загружена!")
            except Exception as e:
                messagebox.showerror("Ошибка чтения матрицы", f"Не удалось прочитать файл:\n{e}")

    def load_seqs_file(self):
        filepath = filedialog.askopenfilename(
            title="Выберите файл с последовательностями",
            filetypes=[("Supported files", "*.txt *.fasta *.csv *.xlsx *.xls"), ("All files", "*.*")]
        )
        if filepath:
            try:
                s1, s2 = load_sequences_from_file(filepath)
                self.seq1_entry.delete("1.0", tk.END)
                self.seq1_entry.insert(tk.END, s1)
                self.seq2_entry.delete("1.0", tk.END)
                self.seq2_entry.insert(tk.END, s2)
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось загрузить последовательности:\n{e}")

    def run_alignment(self):
        if self.scoring_matrix is None:
            messagebox.showwarning("Внимание", "Сначала загрузите файл матрицы замен!")
            return

        seq1 = "".join(self.seq1_entry.get("1.0", tk.END).split()).upper()
        seq2 = "".join(self.seq2_entry.get("1.0", tk.END).split()).upper()

        if not seq1 or not seq2:
            messagebox.showwarning("Внимание", "Введите обе последовательности!")
            return

        invalid1 = set(seq1) - set(AMINO_ACIDS)
        invalid2 = set(seq2) - set(AMINO_ACIDS)
        if invalid1 or invalid2:
            messagebox.showerror("Ошибка алфавита", f"Недопустимые символы: {invalid1 | invalid2}\nРазрешены только 20 стандартных аминокислот: {AMINO_ACIDS}")
            return

        try:
            gap_open = int(self.gap_open_var.get())
            gap_ext = int(self.gap_ext_var.get())
        except ValueError:
            messagebox.showerror("Ошибка", "Штрафы за гэпы должны быть целыми числами!")
            return

        try:
            score, a1, a2, s1, e1, s2, e2, uniq = smith_waterman_affine(
                seq1, seq2, self.scoring_matrix, gap_open, gap_ext
            )
        except KeyError as e:
            messagebox.showerror("Ошибка матрицы", f"Алгоритм не смог найти в таблице букву {e}.\nУбедитесь, что загрузили правильную матрицу.")
            return
        except Exception as e:
            messagebox.showerror("Системная ошибка", f"Сбой при расчете:\n{e}")
            return

        if score == 0:
            self.result_text.delete("1.0", tk.END)
            self.result_text.insert(tk.END, "Выравнивание не найдено. Последовательности не имеют значимых общих участков.")
            return

        match_line = "".join(
            "|" if c1 == c2 else (" " if c1 == "-" or c2 == "-" else ".")
            for c1, c2 in zip(a1, a2)
        )

        aln_length = len(a1)
        matches = match_line.count('|')
        gaps = match_line.count(' ')
        mismatches = match_line.count('.')
        identity_pct = (matches / aln_length) * 100 if aln_length > 0 else 0

        unique_str = "Алгоритм нашел единственный оптимальный вариант (is_unique = True)." if uniq else "Найдено несколько равноценных путей выравнивания, показан один из них (is_unique = False)."

        output = (
            f"=== ОСНОВНЫЕ ПАРАМЕТРЫ ===\n"
            f"Итоговый балл выравнивания     = {score}\n"
            f"Начало 1-й последовательности  = {s1}\n"
            f"Конец 1-й последовательности   = {e1} (длина участка: {e1 - s1 + 1})\n"
            f"Начало 2-й последовательности  = {s2}\n"
            f"Конец 2-й последовательности   = {e2} (длина участка: {e2 - s2 + 1})\n"
            f"Единственный вариант решения   = {'Да' if uniq else 'Нет (есть другие пути с таким же баллом)'}\n\n"
            f"=== СХЕМА ВЫРАВНИВАНИЯ ===\n"
            f"1-я посл.: {a1}\n"
            f"           {match_line}\n"
            f"2-я посл.: {a2}\n\n"
            f"=== БИОЛОГИЧЕСКИЙ ВЫВОД ===\n"
            f"Найдено локальное сходство длиной {aln_length} позиций. Из них:\n"
            f" • Точных совпадений: {matches} ({identity_pct:.1f}% идентичности).\n"
            f" • Точечных замен (мутаций): {mismatches}.\n"
            f" • Вставок/выпадений (гэпов): {gaps}.\n\n"
            f"Наличие протяженного участка с {identity_pct:.1f}% совпадений указывает на "
            f"{'высокую' if identity_pct > 50 else 'умеренную' if identity_pct > 30 else 'низкую'} "
            f"гомологию исследуемых фрагментов."
        )

        self.result_text.delete("1.0", tk.END)
        self.result_text.insert(tk.END, output)

if __name__ == "__main__":
    root = tk.Tk()
    app = AlignmentApp(root)
    root.mainloop()
