"""Main module of the program."""

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
    r: ReconstructionData = fetch_test_data(n=500, k=10, sqne=0)
    ws: WSRY = WSRY(nucleotide_to_weak_strong, r.start, r.ws_probe.cells)
    ry: WSRY = WSRY(nucleotide_to_purine_pyrimidine, r.start, r.ry_probe.cells)
    print(greedy(ws, ry, r))


if __name__ == "__main__":
    main()
