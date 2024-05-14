"""Main module of the program."""

import copy
import sys
from typing import Optional

import requests
import xmltodict
from reconstruction_data import ReconstructionData

sys.setrecursionlimit(10**6)

nucleotide_to_weak_strong = {
    "A": "W",
    "T": "W",
    "C": "S",
    "G": "S",
}

nucleotide_to_purine_pyrimidine = {
    "A": "R",
    "G": "R",
    "C": "Y",
    "T": "Y",
}


class WSRY:
    """Store data for WS or RY cells and provide methods to work with them."""

    def convert_oligo(self, oligo: str) -> str:
        """convert oligo to WS or RY according to the given dictionary without last nucleotide"""
        half = ""
        for i in range(len(oligo) - 1):
            half += self.dict_convertion[oligo[i]]
        half += oligo[-1]

        return half

    def get_tmp_length_solution(self) -> int:
        """Return length of the current solution."""
        return sum([len(self.path[0]) - depth for depth in self.depth])

    @staticmethod
    def connect_ws_ry(oligo_ws: str, oligo_ry: str) -> str:
        """Connects WS and RY oligos according to the rules."""
        connected = ""
        for nucleotide_ws, nucleotide_ry in zip(oligo_ws, oligo_ry):
            temp = nucleotide_ws + nucleotide_ry
            if temp == "SR":
                connected += "G"
            elif temp == "SY":
                connected += "C"
            elif temp == "WR":
                connected += "A"
            elif temp == "WY":
                connected += "T"
            elif nucleotide_ws == nucleotide_ry:
                connected += nucleotide_ws
            else:
                raise ValueError(
                    f"Invalid nucleotides: {nucleotide_ws}, {nucleotide_ry}"
                )

        return connected

    def __init__(self, dict_convertion: dict, oligo: str, cells: dict):
        self.dict_convertion = dict_convertion
        self.start_converted = self.convert_oligo(oligo)  # set_first - Z_set_first

        self.cells_dict = {}  # ols
        for cell in cells:
            self.cells_dict[cell] = False
        self.cells_dict[self.start_converted] = True

        self.path = [self.start_converted]
        self.depth = [0]

    def __repr__(self) -> str:
        return f"Start: {self.start_converted} Path: {self.path} Depth: {self.depth}"


def fetch_test_data(
    n: int = 16,
    k: int = 4,
    mode: str = "binary",
    intensity: int = 0,
    position: int = 0,
    sqpe: int = 0,
    sqne: int = 0,
    pose: int = 0,
) -> ReconstructionData:
    """Fetches test data from server"""
    content = requests.get(
        f"https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php?n={n}&k={k}&mode={mode}&intensity={intensity}&position={position}&sqpe={sqpe}&sqne={sqne}&pose={pose}",
        timeout=10,
    ).content
    if not content:
        raise requests.exceptions.RequestException("Failed to fetch test data")
    data = xmltodict.parse(content)
    return ReconstructionData(data)


def check_overlap(oligo1: str, oligo2: str, probe: int) -> int:
    """Return maximum overlap between two oligos."""
    for offset in range(probe - 1, 0, -1):
        if oligo1[probe - offset :] == oligo2[:offset]:
            return offset
    return 0


# Jeżeli przez liczbę błędów negatywnych obu części spektrum
# oznaczymy różnicę pomiędzy liczbą ich elementów a liczbą elementów w spektrum
# idealnym[...]


def add_ongoing_vertices_to_list(
    ws: WSRY, ry: WSRY
) -> tuple[tuple[str, str, int], ...]:
    """Add not used vertices to the list of candidates(ws,ry,overlap). Return sorted by overlap tuple of candidates."""
    ws = copy.deepcopy(ws)
    ry = copy.deepcopy(ry)

    candidates: list[tuple[str, str, int]] = list()
    last_added_path_ws = ws.path[-1]
    last_added_path_ry = ry.path[-1]
    # temporarly convert last nucleotide for easier comparison
    tmp_last_ws = (
        last_added_path_ws[:-1] + nucleotide_to_weak_strong[last_added_path_ry[-1]]
    )
    tmp_last_ry = (
        last_added_path_ry[:-1]
        + nucleotide_to_purine_pyrimidine[last_added_path_ws[-1]]
    )

    for vertex_ws in ws.cells_dict:  # for (VertexW S ← OverlapSet) do
        if not ws.cells_dict[vertex_ws]:
            for vertex_ry in ry.cells_dict:  # for (VertexRY ← OverlapSet) do
                if not ry.cells_dict[vertex_ry]:
                    if (
                        vertex_ws[-1] == vertex_ry[-1]
                    ):  # if (sameLastNucleotide(VertexWS, VertexRY) = TRUE) then
                        overlap_ws, overlap_ry = check_overlap(
                            tmp_last_ws, vertex_ws, len(vertex_ws)
                        ), check_overlap(tmp_last_ry, vertex_ry, len(vertex_ry))
                        if overlap_ws == overlap_ry and overlap_ws > 0:
                            candidates.append(
                                (vertex_ws, vertex_ry, overlap_ws)
                            )  # Candidates ← addPair(VertexWS, VertexRY);

    return tuple(sorted(candidates, key=lambda x: x[2], reverse=True))


def add_new_vertex_to_solution(
    candidates: tuple[tuple[str, str, int], ...],
    ws: WSRY,
    ry: WSRY,
    r: ReconstructionData,
) -> tuple[tuple[tuple[str, str, int], ...], WSRY, WSRY]:
    """return copies of ws and ry with added candidate to the solution. Remove candidate from the list of candidates. raise ValueError if no candidates to add."""
    ws = copy.deepcopy(ws)
    ry = copy.deepcopy(ry)

    # TODO: sprawdź czy rozszerzone ścieżki nie przekroczą maksymalnej dopuszczalnej długości oraz
    # czy w przypadku nałożenia mniejszego niż maksymalne, powiększony o odpowiednią
    # wartość licznik takich nałożeń wciąż znajduje się poniżej limitu ustalonego dla poszukiwanego rozwiązania
    for candidate in candidates:
        ws.path.append(candidate[0])
        ws.depth.append(candidate[2])
        ry.path.append(candidate[1])
        ry.depth.append(candidate[2])

        if ws.get_tmp_length_solution() > r.length:
            ws.path.pop()
            ws.depth.pop()
            ry.path.pop()
            ry.depth.pop()
            continue

        ws.cells_dict[candidate[0]] = True
        ry.cells_dict[candidate[1]] = True
        ret_candidates = list(candidates)
        ret_candidates.remove(candidate)
        return tuple(ret_candidates), ws, ry

    raise ValueError("No candidates to add to the solution.")


def reconstruct(
    ws: WSRY,
    ry: WSRY,
    r: ReconstructionData,
    solutions: list[tuple[WSRY, WSRY, int]],
    narrowed_candidates: Optional[tuple[tuple[str, str, int], ...]] = None,
):
    if narrowed_candidates is not None and len(narrowed_candidates) == 0:
        return

    if narrowed_candidates is None:
        candidates = add_ongoing_vertices_to_list(ws, ry)
    else:
        candidates = narrowed_candidates

    try:
        _, ws, ry = add_new_vertex_to_solution(candidates, ws, ry, r)
        reconstruct(ws, ry, r, solutions)
        return
    except ValueError:
        tmp_solution_length = ws.get_tmp_length_solution()
        if tmp_solution_length == r.length:
            solutions.append((ws, ry, tmp_solution_length))
            return
        else:
            if len(solutions) == 0:
                solutions.append((ws, ry, tmp_solution_length))
            else:
                if tmp_solution_length > solutions[0][2]:
                    solutions.clear()
                    solutions.append((ws, ry, tmp_solution_length))

    # reverse steps
    ws = copy.deepcopy(ws)
    ry = copy.deepcopy(ry)
    narr_candidates: list[tuple[str, str, int]] = list(candidates)
    try:
        # sometimes it throws a ValueError exception, because the candidate is not in the list
        # then solution got only start oligo
        # what to do in this case?
        narr_candidates.remove((ws.path[-1], ry.path[-1], ws.depth[-1]))
    except ValueError:
        pass
    ws.cells_dict[ws.path[-1]] = False
    ry.cells_dict[ry.path[-1]] = False
    ws.path.pop()
    ws.depth.pop()
    ry.path.pop()
    ry.depth.pop()
    reconstruct(ws, ry, r, solutions, tuple(narr_candidates))


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
    print(reconstructed_dna, end=" ")
    print(len(reconstructed_dna))


main()
