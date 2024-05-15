"""Tabu search metaheuristic implementation."""

import copy
import pprint
import random
from enum import Enum
from typing import Optional
from wsry import WSRY
from common import check_overlap
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
        self, ws: WSRY, ry: WSRY
    ) -> tuple[WSRY, WSRY] | tuple[()]:
        """Generate neighbours for current solution by inserting oligo."""

        # ws i ry muszą mieć ten sam ostatni nukleotyd oraz mieć identyczny overlap z oligo przed miejscem wstawienia
        not_used_not_tabu_oligos_ws = [
            oligo for oligo in ws.not_used_oligos() if not self.is_tabu(oligo)
        ]
        not_used_not_tabu_oligos_ry = [
            oligo for oligo in ry.not_used_oligos() if not self.is_tabu(oligo)
        ]

        # find index where to insert (take care of clusters). skip 0 index
        idx_to_insert = random.randint(1, len(ws.path))
        previous_oligo_ws = ws.path[idx_to_insert - 1]
        previous_oligo_ry = ry.path[idx_to_insert - 1]

        random.shuffle(not_used_not_tabu_oligos_ws)
        random.shuffle(not_used_not_tabu_oligos_ry)  # is it necessary?
        for tmp_oligo_ws in not_used_not_tabu_oligos_ws:
            for tmp_oligo_ry in not_used_not_tabu_oligos_ry:
                if tmp_oligo_ws[-1] == tmp_oligo_ry[-1]:
                    oligo_ws_overlap = check_overlap(
                        previous_oligo_ws, tmp_oligo_ws, len(tmp_oligo_ws)
                    )
                    oligo_ry_overlap = check_overlap(
                        previous_oligo_ry, tmp_oligo_ry, len(tmp_oligo_ry)
                    )
                    if oligo_ws_overlap == oligo_ry_overlap:
                        new_ws = copy.deepcopy(ws)
                        new_ry = copy.deepcopy(ry)
                        new_ws.path.insert(idx_to_insert, tmp_oligo_ws)
                        new_ry.path.insert(idx_to_insert, tmp_oligo_ry)
                        new_ws.depth.insert(idx_to_insert, oligo_ws_overlap)
                        new_ry.depth.insert(idx_to_insert, oligo_ry_overlap)
                        return (new_ws, new_ry)

        return ()

    def generate_neighbours(self, ws: WSRY, ry: WSRY) -> tuple[tuple[WSRY, WSRY], ...]:
        """Generate neighbours for current solution."""
        neighbours: list[tuple[WSRY, WSRY]] = []
        for _ in range(self.number_of_neighbours):
            # chosen_move = random.choice(list(Moves))
            chosen_move = Moves.INSERT_OLIGO
            match chosen_move:
                case Moves.INSERT_OLIGO:
                    neighbour = self.generate_neighbour_insert_oligo(ws, ry)
                    if neighbour:
                        neighbours.append(neighbour)
                case Moves.DELETE_OLIGO:
                    pass
                case Moves.DELETE_CLUSTER:
                    pass
                case Moves.SHIFT_OLIGO:
                    pass
                case Moves.SHIFT_CLUSTER:
                    pass
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
        neighbours = self.generate_neighbours(best_solution[0], best_solution[1])
        print(best_solution)
        for neighbour in neighbours:
            print(neighbour[0], neighbour[1])

        # for i in range(self.number_of_iterations):
        #     pass
