"""Main module of the program."""

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
    r: ReconstructionData = fetch_test_data(n=100, k=8, sqne=13)
    ws: WSRY = WSRY(nucleotide_to_weak_strong, r.start, r.ws_probe.cells)
    ry: WSRY = WSRY(nucleotide_to_purine_pyrimidine, r.start, r.ry_probe.cells)
    t1 = Tabu(4, 100, 100)
    # print(ws.not_used_oligos())
    # print(greedy(ws, ry, r))
    t1.find_solution(ws, ry, r, greedy(ws, ry, r))


if __name__ == "__main__":
    main()
