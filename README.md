# Реализация локального выравнивания белковых последовательностей по алгоритму Смита–Ватермана с аффинной моделью гэпов.
## Выполнили: Беседин, Кульков, Розинько

# Текст задания:
<img width="476" height="910" alt="image" src="https://github.com/user-attachments/assets/dc44a94f-8b2b-4320-b2c1-eb1ad9c3e804" />
<img width="501" height="328" alt="image" src="https://github.com/user-attachments/assets/66f22dcf-71d6-4b03-8e25-e5f0962ca482" />
Стек: Python.
Вход: - `seq1`, `seq2` - две строки (однобуквенный код аминокислот), `blosum62` - словарь весов замен (загружается из `blosum62.txt` через `parse_blosum62`).

Выход: `alignment_score` - int счёт локального выравнивания
`aligned_seq1` - str; выровненная seq1 (с гэпами `-`)
`aligned_seq2` - str; выровненная seq2 (с гэпами `-`)
`start1`, `end1` - int; координаты в seq1 (1-based, включительно)
`start2`, `end2` - int; координаты в seq2 (1-based, включительно)
`is_unique` - bool; `True`, если оптимум единственный
