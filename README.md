## DNA Sequencing Using Binary Chips with Negative Errors

---

### **Problem Overview**

#### Binary Chip

- Consists of two parts using alphabets:
`{W, S, A, C, G, T}` and `{R, Y, A, C, G, T}`.
    - `W = {A, T}` – Weak nucleotides
    - `S = {C, G}` – Strong nucleotides
    - `R = {A, G}` – Purines
    - `Y = {C, T}` – Pyrimidines
- **Encoding:**
Each chip element has a length of `k+1` and uses two types of probes `{W,S}` and `{R,Y}` to improve error control during hybridization. This dual-path approach ensures sequence reconstruction by verifying correctness.


#### Binary SBH with Negative Errors

- Negative errors reduce elements in the spectrum parts. However, corresponding elements in both spectra are rarely lost simultaneously.
- Errors affect reconstruction and successor selection.

---

### **Exact Algorithm**

1. **Input Data:**
Input is processed from a URL. The algorithm uses a recursive function `reconstruct()` to build the DNA sequence.
2. **Steps:**
    - Generate candidate oligonucleotides based on overlaps.
    - Add candidates to the solution while ensuring constraints are met.
    - Backtrack if necessary to find the optimal sequence.
3. **Output:**
Combines both paths into a reconstructed DNA sequence.

---

### **Greedy Algorithm**

- Iteratively selects the pair with the largest overlap until no further addition is possible.

```python
def greedy(ws: WSRY, ry: WSRY, r: ReconstructionData) -> tuple[WSRY, WSRY]:
    """Greedy algorithm for DNA reconstruction."""
    while True:
        pair = first_nonzero_overlap_pair(ws, ry)
        if pair is None:
            break
        # Update paths and depths
    return ws, ry
```

---

### **Metaheuristic - Tabu Search**

1. **Evaluation Functions:**
    - Global grade: Maximizes the number of oligonucleotides in the solution.
    - Condensation grade: Maximizes overlap between oligonucleotides.
2. **Neighbor Generation:**
Five types of moves are used (e.g., adding/removing oligonucleotides or clusters).
3. **Algorithm Steps:**
    - Start with a greedy solution.
    - Generate neighbors while considering the Tabu list.
    - Select the best neighbor using evaluation functions.

---

### **Performance Analysis**

- Tested on datasets with parameters: `n={50, 60,...200}`, `k=10`, `sqne=12`.
- Results:
    - Exact algorithm performs well for small sequences but struggles with larger ones due to recursion limits.
    - Tabu search achieves near-optimal results within time constraints (e.g., 60 seconds).


#### Example Results:

| Length | Errors | Time (s) |
| :-- | :-- | :-- |
| 50 | 12 | 0.37 |
| 100 | 15 | 60.54 |
| 200 | 14 | 64.73 |

---

### **Execution**

Run the program from the appropriate directory (`metaheuristic` or `exact algorithm`) using:

```bash
python3 main.py path_to_input_file.xml
```

---

### **Conclusion**

Python's ease of XML parsing made it suitable for this project but introduced performance limitations for larger sequences due to recursion and memory usage. For larger datasets (`n > 100`), switching to a lower-level language and using JSON input could improve efficiency.

Despite these challenges, the approximate algorithm provided high-quality results under reasonable constraints.

---

### References

1. M. Radom, *Combinatorial Aspects of Non-Classical DNA Sequencing by Hybridization*, Poznań 2011.
2. Kamil Kwarciak, Piotr Formanowicz, *Tabu Search Algorithm for DNA Sequencing by Hybridization*, Poznań 2014.
3. J. Błażewicz et al., *Tabu Search for DNA Sequencing with False Negatives and Positives*, Poznań 1999.
