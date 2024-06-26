"""Main module of the program."""

import sys
import time
from tabu import Tabu
from wsry import WSRY
from greedy import greedy
from reconstruction_data import ReconstructionData
from common import (
    fetch_test_data,
    nucleotide_to_weak_strong,
    nucleotide_to_purine_pyrimidine,
)


def main() -> None:
    """Main function of the program."""
    if len(sys.argv) != 2:
        print("Usage: python main.py <xml_file>")
        return
    filename = sys.argv[1]
    r: ReconstructionData = fetch_test_data(filename)
    start = time.time()
    ws: WSRY = WSRY(nucleotide_to_weak_strong, r.start, r.ws_probe.cells)
    ry: WSRY = WSRY(nucleotide_to_purine_pyrimidine, r.start, r.ry_probe.cells)
    t1 = Tabu(3, 100, 100)
    # print(ws.not_used_oligos())
    # print(greedy(ws, ry, r))
    sol = t1.find_solution(ws, ry, r, greedy(ws, ry, r))
    end = time.time()
    print(Tabu.reconstruct_dna(sol[0], sol[1]))
    # print("find_solution: ", sol)
    # print(
    #     "Długość odtworzonego DNA: ",
    #     # Tabu.reconstruct_dna(sol[0], sol[1]),
    #     len(Tabu.reconstruct_dna(sol[0], sol[1])),
    # )
    # print(
    #     f"Suma wartości nieidelnych nałożeń (idealnie {r.sqne}): {Tabu.count_negative_errors(sol[0])}"
    # )
    # print("Czas wykonania: ", end - start)


if __name__ == "__main__":
    main()
