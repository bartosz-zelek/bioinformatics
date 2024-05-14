"""Tabu search metaheuristic implementation."""

import random
from enum import Enum
from wsry import WSRY
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

    def generate_neighbours(self):
        """Generate neighbours for current solution."""
        chosen_move = random.choice(list(Moves))
        match chosen_move:
            case Moves.INSERT_OLIGO:
                pass
            case Moves.DELETE_OLIGO:
                pass
            case Moves.DELETE_CLUSTER:
                pass
            case Moves.SHIFT_OLIGO:
                pass
            case Moves.SHIFT_CLUSTER:
                pass

    def find_solution(
        self,
        ws: WSRY,
        ry: WSRY,
        r: ReconstructionData,
        greedy_solution: tuple[WSRY, WSRY, int],
    ):
        """
        Find solution using tabu search.
        ws, ry - WSRY objects with initial paths
        """
        best_solution = greedy_solution
        for i in range(self.number_of_iterations):
            pass
