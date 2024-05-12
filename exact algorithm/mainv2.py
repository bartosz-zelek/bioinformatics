from reconstruction import WSRY, reconstruct
from wsry import nucleotide_to_weak_strong, nucleotide_to_purine_pyrimidine
from data import ReconstructionData, fetch_test_data


def main() -> None:
    """Main function of the program."""
    r: ReconstructionData = fetch_test_data(sqne=4)
    ws: WSRY = WSRY(nucleotide_to_weak_strong, r.start, r.ws_probe.cells)
    ry: WSRY = WSRY(nucleotide_to_purine_pyrimidine, r.start, r.ry_probe.cells)
    solutions: list[tuple[WSRY, WSRY, int]] = list()
    reconstruct(ws, ry, r, solutions)
    best_solution = sorted(solutions, key=lambda x: x[2], reverse=True)[0]
    reconstructed_dna = ""
    for ws_oligo, ry_oligo, depth in zip(
        best_solution[0].path, best_solution[1].path, best_solution[0].depth
    ):
        connected = WSRY.connect_ws_ry(ws_oligo, ry_oligo)
        reconstructed_dna += connected[depth - len(ws_oligo) :]
    print(reconstructed_dna)


if __name__ == "__main__":
    main()
