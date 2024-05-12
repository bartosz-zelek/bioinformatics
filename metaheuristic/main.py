from greedy import greedy
from reconstruction_data import ReconstructionData
from wsry import WSRY
from common import (
    fetch_test_data,
    nucleotide_to_weak_strong,
    nucleotide_to_purine_pyrimidine,
)


def main() -> None:
    """Main function of the program."""
    while True:
        r: ReconstructionData = fetch_test_data(sqne=4)
        ws: WSRY = WSRY(nucleotide_to_weak_strong, r.start, r.ws_probe.cells)
        ry: WSRY = WSRY(nucleotide_to_purine_pyrimidine, r.start, r.ry_probe.cells)
        g = greedy(ws, ry, r)
        if len(g[0].path) == 1:
            print(r)
            print(g)
            break


main()
