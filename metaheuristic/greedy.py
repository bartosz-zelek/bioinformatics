import copy
from typing import Optional
from wsry import WSRY
from reconstruction_data import ReconstructionData
from common import (
    check_overlap,
    nucleotide_to_weak_strong,
    nucleotide_to_purine_pyrimidine,
)


def first_nonzero_overlap_pair(ws: WSRY, ry: WSRY) -> Optional[tuple[str, str, int]]:
    """
    Find the biggest overlap between last added and not used oligos (according to the rules).
    Is it too slow for greedy?
    Probably want to use a different data structure, than list for candidates.
    """
    last_added_path_ws = ws.path[-1]
    last_added_path_ry = ry.path[-1]
    # temporarly convert last nucleotide for easier comparison
    tmp_last_ws = (
        last_added_path_ws[:-1] + nucleotide_to_weak_strong[last_added_path_ws[-1]]
    )
    tmp_last_ry = (
        last_added_path_ry[:-1]
        + nucleotide_to_purine_pyrimidine[last_added_path_ry[-1]]
    )

    candidates: list[tuple[str, str, int]] = list()

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
    try:
        return tuple(sorted(candidates, key=lambda x: x[2], reverse=True))[0]
    except IndexError:
        return None


def greedy(ws: WSRY, ry: WSRY, r: ReconstructionData) -> tuple[WSRY, WSRY]:
    """Greedy algorithm for DNA reconstruction."""
    ws = copy.deepcopy(ws)
    ry = copy.deepcopy(ry)

    reconstructed_dna_length = len(ws.start_converted)

    while True:
        pair = first_nonzero_overlap_pair(
            ws, ry
        )  # tuple[str, str, int] -> (VertexWS, VertexRY, overlap)

        if pair is None:
            break
        reconstructed_dna_length += len(pair[0]) - pair[2]
        if reconstructed_dna_length > r.length:
            break

        ws.cells_dict[pair[0]] = True
        ry.cells_dict[pair[1]] = True
        ws.path.append(pair[0])
        ry.path.append(pair[1])
        ws.depth.append(pair[2])
        ry.depth.append(pair[2])

    return ws, ry
