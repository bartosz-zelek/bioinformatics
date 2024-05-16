"""Tabu search metaheuristic implementation."""

import copy
import pprint
import random
from enum import Enum
from typing import Optional
from wsry import WSRY
from common import (
    check_overlap,
    nucleotide_to_weak_strong,
    nucleotide_to_purine_pyrimidine,
)
from reconstruction_data import ReconstructionData


class Moves(Enum):
    """All possible moves in tabu search."""

    INSERT_OLIGO = 1
    DELETE_OLIGO = 2
    DELETE_CLUSTER = 3
    SHIFT_OLIGO = 4
    SHIFT_CLUSTER = 5


class Tabu:
    """Class for tabu search metaheuristic."""

    def __init__(self, tabu_size, number_of_iterations, number_of_neighbours):
        self.tabu_size = tabu_size
        self.tabu_list_ws = []
        self.tabu_list_ry = []
        self.number_of_iterations = number_of_iterations
        self.number_of_neighbours = number_of_neighbours

    def add(self, move_ws, move_ry):
        """Add move to tabu list and remove oldest move if list is full."""
        self.tabu_list_ws.append(move_ws)
        self.tabu_list_ry.append(move_ry)
        if len(self.tabu_list_ws) > self.tabu_size:
            self.tabu_list_ws.pop(0)
            self.tabu_list_ry.pop(0)

    def is_tabu(self, move):
        """Check if move is in tabu list."""
        return (move in self.tabu_list_ws) or (move in self.tabu_list_ry)

    def generate_neighbour_insert_oligo(
        self,
        ws: WSRY,
        ry: WSRY,
        not_used_not_tabu_oligos_ws,
        not_used_not_tabu_oligos_ry,
    ) -> tuple[WSRY, WSRY] | tuple[()]:
        """
        Generate neighbours for current solution by inserting oligo.
        WS and RY must have the same last nucleotide and have the same overlap with the oligo before the insertion point.
        """

        ### INDEX OUTSIDE CLUSTER - BEGIN
        biggest_possible_overlap = len(ws.start_converted) - 1
        possible_idxs = []
        for idx in range(1, len(ws.path) + 1):
            if (
                idx != len(ws.path)
                and ws.depth[idx - 1] == biggest_possible_overlap
                and ws.depth[idx] == biggest_possible_overlap
            ):
                continue
            possible_idxs.append(idx)
        ### INDEX OUTSIDE CLUSTER - END
        if possible_idxs:

            random.shuffle(not_used_not_tabu_oligos_ws)
            for tmp_oligo_ws in not_used_not_tabu_oligos_ws:
                for tmp_oligo_ry in not_used_not_tabu_oligos_ry:
                    if tmp_oligo_ws[-1] == tmp_oligo_ry[-1]:
                        for idx_to_insert in possible_idxs:
                            new_ws = copy.deepcopy(ws)
                            new_ry = copy.deepcopy(ry)
                            new_ws.path.insert(idx_to_insert, tmp_oligo_ws)
                            new_ry.path.insert(idx_to_insert, tmp_oligo_ry)
                            new_ws.update_depth()
                            new_ry.update_depth()

                            # if new_ws.depth != new_ry.depth or 0 in new_ws.depth[1:]:
                            if new_ws.depth != new_ry.depth:
                                continue
                            return (new_ws, new_ry)

        return ()

    def generate_neighbour_delete_oligo(
        self,
        ws: WSRY,
        ry: WSRY,
        skip_tabu_check=False,
    ) -> tuple[WSRY, WSRY] | tuple[()]:
        """Generate neighbours for current solution by deleting oligo."""

        ### INDEX OUTSIDE CLUSTER - BEGIN
        biggest_possible_overlap = len(ws.start_converted) - 1
        while True:
            random_delete_idx = random.randint(
                1, len(ws.path) - 1
            )  # randint or incrementing index?
            if (
                random_delete_idx != len(ws.path) - 1
                and ws.depth[random_delete_idx] == biggest_possible_overlap
                and ws.depth[random_delete_idx + 1] == biggest_possible_overlap
            ):
                continue
            break
        ### INDEX OUTSIDE CLUSTER - END

        if not self.is_tabu(ws.path[random_delete_idx]) or skip_tabu_check:

            ws = copy.deepcopy(ws)
            ry = copy.deepcopy(ry)
            ws.path.pop(random_delete_idx)
            ry.path.pop(random_delete_idx)
            ws.update_depth()
            ry.update_depth()

            return (ws, ry)
        return ()

    def generate_neighbour_shift_oligo(
        self, ws: WSRY, ry: WSRY
    ) -> tuple[WSRY, WSRY] | tuple[()]:
        """Generate neighbours for current solution by shifting oligo."""
        ### INDEX OUTSIDE CLUSTER - BEGIN - choose random oligo to shift
        biggest_possible_overlap = len(ws.start_converted) - 1
        while True:
            random_oligo_idx = random.randint(
                1, len(ws.path) - 1
            )  # randint or incrementing index?
            if (
                random_oligo_idx != len(ws.path) - 1
                and ws.depth[random_oligo_idx] == biggest_possible_overlap
                and ws.depth[random_oligo_idx + 1] == biggest_possible_overlap
            ):
                continue
            break
        ### INDEX OUTSIDE CLUSTER - END

        if not self.is_tabu(ws.path[random_oligo_idx]):

            ws = copy.deepcopy(ws)
            ry = copy.deepcopy(ry)
            oligo_ws = ws.path.pop(random_oligo_idx)
            oligo_ry = ry.path.pop(random_oligo_idx)
            ws.update_depth()
            ry.update_depth()

            idxs = []
            ### INDEX OUTSIDE CLUSTER - BEGIN - choose random index to insert
            for idx_to_insert in range(1, len(ws.path) + 1):
                if (
                    idx_to_insert != len(ws.path)
                    and ws.depth[idx_to_insert - 1] == biggest_possible_overlap
                    and ws.depth[idx_to_insert] == biggest_possible_overlap
                ):
                    continue
                idxs.append(idx_to_insert)
            ### INDEX OUTSIDE CLUSTER - END
            random.shuffle(idxs)
            for i in idxs:
                ws.path.insert(i, oligo_ws)
                ry.path.insert(i, oligo_ry)
                ws.update_depth()
                ry.update_depth()
                if ws.depth == ry.depth and 0 not in ws.depth[1:]:
                    return (ws, ry)
                ws.path.pop(i)
                ry.path.pop(i)
        return ()

    def generate_neighbours(self, ws: WSRY, ry: WSRY) -> tuple[tuple[WSRY, WSRY], ...]:
        """Generate neighbours for current solution."""
        neighbours: list[tuple[WSRY, WSRY]] = []
        ws_paths_set = {str(ws.path)}  # terrible solution, but it works
        not_used_not_tabu_oligos_ws = [
            oligo for oligo in ws.not_used_oligos() if not self.is_tabu(oligo)
        ]
        not_used_not_tabu_oligos_ry = [
            oligo for oligo in ry.not_used_oligos() if not self.is_tabu(oligo)
        ]

        for _ in range(self.number_of_neighbours):
            # chosen_move = random.choice(list(Moves))
            # chosen_move = Moves.INSERT_OLIGO
            # chosen_move = Moves.DELETE_OLIGO
            chosen_move = Moves.SHIFT_OLIGO
            match chosen_move:
                case Moves.INSERT_OLIGO:
                    neighbour = self.generate_neighbour_insert_oligo(
                        ws, ry, not_used_not_tabu_oligos_ws, not_used_not_tabu_oligos_ry
                    )
                case Moves.DELETE_OLIGO:
                    neighbour = self.generate_neighbour_delete_oligo(ws, ry)
                case Moves.DELETE_CLUSTER:
                    pass
                case Moves.SHIFT_OLIGO:
                    neighbour = self.generate_neighbour_shift_oligo(ws, ry)
                case Moves.SHIFT_CLUSTER:
                    pass

            if neighbour and str(neighbour[0].path) not in ws_paths_set:
                ws_paths_set.add(str(neighbour[0].path))
                neighbours.append(neighbour)

        return tuple(neighbours)

    def find_solution(
        self,
        ws: WSRY,
        ry: WSRY,
        r: ReconstructionData,
        greedy_solution: tuple[WSRY, WSRY],
    ):
        """
        Find solution using tabu search.
        ws, ry - WSRY objects with initial paths
        """
        best_solution = greedy_solution
        # best_solution = (ws, ry)
        print(best_solution)
        neighbours = self.generate_neighbours(best_solution[0], best_solution[1])
        for neighbour in neighbours:
            print(neighbour[0], neighbour[1])
