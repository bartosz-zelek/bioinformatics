import copy
from typing import Optional
from wsry import WSRY, nucleotide_to_weak_strong, nucleotide_to_purine_pyrimidine
from data import ReconstructionData


def check_overlap(oligo1: str, oligo2: str, probe: int) -> int:
    """Return maximum overlap between two oligos."""
    for offset in range(probe - 1, 0, -1):
        if oligo1[probe - offset :] == oligo2[:offset]:
            return offset
    return 0


def add_ongoing_vertices_to_list(
    ws: WSRY, ry: WSRY, no_more_errors: bool = False
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
                        if not no_more_errors:
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
        # Use depth to check if we go all errors already to set no_more_errors flag
        # if ws.depth[-1] == r.length:
        #     no_more_errors = True
        # else:
        #     no_more_errors = False
        no_more_errors = ws.depth[-1] == r.length
        candidates = add_ongoing_vertices_to_list(ws, ry, no_more_errors)
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
    narr_candidates.remove((ws.path[-1], ry.path[-1], ws.depth[-1]))
    ws.cells_dict[ws.path[-1]] = False
    ry.cells_dict[ry.path[-1]] = False
    ws.path.pop()
    ws.depth.pop()
    ry.path.pop()
    ry.depth.pop()
    reconstruct(ws, ry, r, solutions, tuple(narr_candidates))