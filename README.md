Assignment: Local alignment

Example execution: `$ python main.py dynprog "ABCD" "[[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]]" "AAAAACCDDCCDDAAAAACC" "CCAAADDAAAACCAAADDCCAAAA"`

- Execution mode: `dynprog`
- Alphabet: `ABCD`
- Substitution matrix: `[[1, -5, -5, -5, -1], [-5, 1, -5, -5, -1], [-5, -5, 5, -5, -4], [-5, -5, -5, 6, -4], [-1, -1, -4, -4, -9]]`
- Sequence 1: "AAAAACCDDCCDDAAAAACC"
- Sequence 2: "CCAAADDAAAACCAAADDCCAAAA"
