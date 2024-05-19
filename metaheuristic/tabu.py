"""Tabu search metaheuristic implementation."""

import copy
import random
from enum import Enum
from typing import Optional
from wsry import WSRY
from reconstruction_data import ReconstructionData


class Moves(Enum):
    """All possible moves in tabu search."""

    INSERT_OLIGO = 1
    DELETE_OLIGO = 2
    DELETE_CLUSTER = 3
    SHIFT_OLIGO = 4
    SHIFT_CLUSTER = 5


# TODO: cluster list can be generated once at iteration
class Tabu:
    """Class for tabu search metaheuristic."""

    @staticmethod
    def reconstruct_dna(ws: WSRY, ry: WSRY) -> str:
        """
        Reconstruct DNA from WS and RY paths.
        """
        reconstructed_dna = ""
        for ws_oligo, ry_oligo, depth in zip(ws.path, ry.path, ws.depth):
            connected = WSRY.connect_ws_ry(ws_oligo, ry_oligo)
            reconstructed_dna += connected[depth - len(ws_oligo) :]
        return reconstructed_dna

    @staticmethod
    def count_negative_errors(ws: WSRY) -> int:
        """Count negative errors (not perfect overlaps) in the solution."""
        cnt = 0
        for el in ws.depth[1:]:
            cnt += ws.perfect_overlap - el
        return cnt

    def __init__(self, tabu_size, number_of_iterations, number_of_neighbours):
        self.tabu_size = tabu_size
        self.tabu_list_ws = []
        self.tabu_list_ry = []
        self.number_of_iterations = number_of_iterations
        self.number_of_neighbours = number_of_neighbours
        self.ommit_tabu = False

    def add(self, move_ws, move_ry):
        """Add move to tabu list and remove oldest move if list is full."""
        self.tabu_list_ws.append(move_ws)
        self.tabu_list_ry.append(move_ry)
        # if len(self.tabu_list_ws) > self.tabu_size:
        #     self.tabu_list_ws.pop(0)
        #     self.tabu_list_ry.pop(0)

    def is_tabu(self, move):
        """Check if move is in tabu list."""
        if self.ommit_tabu:
            return False
        return (move in self.tabu_list_ws) or (move in self.tabu_list_ry)

    def generate_neighbour_insert_oligo(
        self,
        ws: WSRY,
        ry: WSRY,
        r: ReconstructionData,
        not_used_not_tabu_oligos_ws,
        not_used_not_tabu_oligos_ry,
    ) -> Optional[tuple[WSRY, WSRY, str, str]]:
        """
        Generate neighbours for current solution by inserting oligo.
        WS and RY must have the same last nucleotide and have the same overlap with the oligo before the insertion point.
        """

        ### INDEXES OUTSIDE CLUSTER - BEGIN
        possible_idxs = []
        for idx in range(1, len(ws.path) + 1):
            if idx != len(ws.path) and (
                ws.depth[idx - 1] == ws.perfect_overlap
                and ws.depth[idx] == ws.perfect_overlap
            ):
                continue
            possible_idxs.append(idx)
        ### INDEXES OUTSIDE CLUSTER - END
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
                            if (
                                new_ws.depth == new_ry.depth
                                and len(Tabu.reconstruct_dna(new_ws, new_ry))
                                <= r.length
                            ):
                                new_ws.cells_dict[tmp_oligo_ws] = True
                                new_ry.cells_dict[tmp_oligo_ry] = True
                                return (new_ws, new_ry, tmp_oligo_ws, tmp_oligo_ry)

        return None

    def generate_neighbour_delete_oligo(
        self,
        ws: WSRY,
        ry: WSRY,
        r: ReconstructionData,
        skip_tabu_check=False,
    ) -> Optional[tuple[WSRY, WSRY, str, str]]:
        """Generate neighbours for current solution by deleting oligo."""

        ### INDEX OUTSIDE CLUSTER - BEGIN
        possible_idxs = []
        for idx in range(1, len(ws.path)):
            if (
                idx != len(ws.path) - 1
                and (
                    ws.depth[idx] == ws.perfect_overlap
                    and ws.depth[idx + 1] == ws.perfect_overlap
                )
            ) or (self.is_tabu(ws.path[idx]) is True and skip_tabu_check is False):
                continue
            possible_idxs.append(idx)
        ### INDEX OUTSIDE CLUSTER - END
        random.shuffle(possible_idxs)
        org_ws = copy.deepcopy(ws)
        org_ry = copy.deepcopy(ry)
        for random_delete_idx in possible_idxs:
            ws = copy.deepcopy(org_ws)
            ry = copy.deepcopy(org_ry)
            deleted_ws = ws.path.pop(
                random_delete_idx
            )  # IndexError: pop index out of range
            deleted_ry = ry.path.pop(random_delete_idx)
            ws.update_depth()
            ry.update_depth()

            if ws.depth == ry.depth and len(Tabu.reconstruct_dna(ws, ry)) <= r.length:
                ws.cells_dict[deleted_ws] = False
                ry.cells_dict[deleted_ry] = False
                return (ws, ry, deleted_ws, deleted_ry)
        return None

    def generate_neighbour_shift_oligo(
        self, ws: WSRY, ry: WSRY, r: ReconstructionData
    ) -> Optional[tuple[WSRY, WSRY]]:
        """Generate neighbours for current solution by shifting oligo."""
        ### INDEX OUTSIDE CLUSTER - BEGIN - choose random oligo to shift
        possible_idxs: list[int] = []
        idxs: list[int] = []
        for idx in range(1, len(ws.path)):
            if (
                idx != len(ws.path) - 1
                and (
                    ws.depth[idx] == ws.perfect_overlap
                    and ws.depth[idx + 1] == ws.perfect_overlap
                )
            ) or (self.is_tabu(ws.path[idx])):
                continue
            possible_idxs.append(idx)
        ### INDEX OUTSIDE CLUSTER - END
        random.shuffle(possible_idxs)
        org_ws = copy.deepcopy(ws)
        org_ry = copy.deepcopy(ry)
        for random_oligo_idx in possible_idxs:
            ws = copy.deepcopy(org_ws)
            ry = copy.deepcopy(org_ry)
            oligo_ws = ws.path.pop(random_oligo_idx)
            oligo_ry = ry.path.pop(random_oligo_idx)
            ws.update_depth()
            ry.update_depth()

            if not idxs:
                ### INDEX OUTSIDE CLUSTER - BEGIN - choose random index to insert
                for idx_to_insert in range(1, len(ws.path) + 1):
                    if idx_to_insert == idx or (
                        idx_to_insert != len(ws.path)
                        and ws.depth[idx_to_insert - 1] == ws.perfect_overlap
                        and ws.depth[idx_to_insert] == ws.perfect_overlap
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
                # if ws.depth == ry.depth and 0 not in ws.depth[1:]: # do we allow 0 in depth?
                if (
                    random_oligo_idx != i
                    and ws.depth == ry.depth
                    and len(Tabu.reconstruct_dna(ws, ry)) <= r.length
                ):
                    # print(f"shifted: {random_oligo_idx} to {i}")
                    return (ws, ry)
                ws.path.pop(i)
                ry.path.pop(i)
        return None

    def generate_neigbour_delete_cluster(
        self, ws: WSRY, ry: WSRY, r: ReconstructionData
    ) -> Optional[tuple[WSRY, WSRY]]:
        """Generate neighbours for current solution by deleting cluster."""
        # find indexes of first elements of clusters
        ws = copy.deepcopy(ws)
        ry = copy.deepcopy(ry)

        # FIND CLUSTERS - BEGIN
        begin_of_cluster = False
        cluster_indexes = []
        for i in range(1, len(ws.path) - 1):  # last element is never begin of cluster
            if ws.depth[i] == ws.perfect_overlap:
                if begin_of_cluster is False and ws.depth[i + 1] == ws.perfect_overlap:
                    cluster_indexes.append(i)
                    begin_of_cluster = True
            else:
                begin_of_cluster = False
        # FIND CLUSTERS - END

        if cluster_indexes:  # just one check - empty neigbour if some oligo won't fit
            random_cluster_idx = random.choice(cluster_indexes)
            while (  # remove cluster
                random_cluster_idx < len(ws.path)
                and ws.depth[random_cluster_idx] == ws.perfect_overlap
            ):
                deleted_ws = ws.path.pop(random_cluster_idx)
                deleted_ry = ry.path.pop(random_cluster_idx)
                ws.depth.pop(random_cluster_idx)
                ry.depth.pop(random_cluster_idx)
                ws.cells_dict[deleted_ws] = False
                ry.cells_dict[deleted_ry] = False

            ws.update_depth()
            ry.update_depth()
            if ws.depth == ry.depth and len(Tabu.reconstruct_dna(ws, ry)) <= r.length:
                return (ws, ry)

        return None

    def generate_neighbour_shift_cluster(
        self, ws: WSRY, ry: WSRY, r: ReconstructionData
    ) -> Optional[tuple[WSRY, WSRY]]:
        """Generate neighbours for current solution by shifting cluster."""
        ws = copy.deepcopy(ws)
        ry = copy.deepcopy(ry)

        # FIND CLUSTERS - BEGIN
        begin_of_cluster = False
        cluster_indexes = []
        for i in range(1, len(ws.path) - 1):  # last element is never begin of cluster
            if ws.depth[i] == ws.perfect_overlap:
                if begin_of_cluster is False and ws.depth[i + 1] == ws.perfect_overlap:
                    cluster_indexes.append(i)
                    begin_of_cluster = True
            else:
                begin_of_cluster = False
        # FIND CLUSTERS - END

        cluster_to_shift_ws = []
        cluster_to_shift_ry = []

        if cluster_indexes:
            random_cluster_idx = random.choice(cluster_indexes)
            while (
                random_cluster_idx < len(ws.path)
                and ws.depth[random_cluster_idx] == ws.perfect_overlap
            ):
                cluster_to_shift_ws.append(ws.path.pop(random_cluster_idx))
                cluster_to_shift_ry.append(ry.path.pop(random_cluster_idx))
                ws.depth.pop(random_cluster_idx)
                ry.depth.pop(random_cluster_idx)

            ws.update_depth()
            ry.update_depth()
            ### INDEXES OUTSIDE CLUSTER - BEGIN
            possible_idxs = []
            for idx in range(1, len(ws.path) + 1):
                if (
                    idx != len(ws.path)
                    and ws.depth[idx - 1] == ws.perfect_overlap
                    and ws.depth[idx] == ws.perfect_overlap
                ) or idx == random_cluster_idx:
                    continue
                possible_idxs.append(idx)
            ### INDEXES OUTSIDE CLUSTER - END
            random.shuffle(possible_idxs)
            ws_org = copy.deepcopy(ws)
            ry_org = copy.deepcopy(ry)
            for random_idx_to_insert in possible_idxs:
                ws = copy.deepcopy(ws_org)
                ry = copy.deepcopy(ry_org)
                _cluster_to_shift_ws = cluster_to_shift_ws.copy()
                _cluster_to_shift_ry = cluster_to_shift_ry.copy()
                while _cluster_to_shift_ws:
                    ws.path.insert(random_idx_to_insert, _cluster_to_shift_ws.pop())
                    ry.path.insert(random_idx_to_insert, _cluster_to_shift_ry.pop())

                ws.update_depth()
                ry.update_depth()
                if (
                    ws.depth == ry.depth
                    and len(Tabu.reconstruct_dna(ws, ry)) <= r.length
                ):
                    return (ws, ry)
        return None

    def generate_neighbours(
        self, ws: WSRY, ry: WSRY, r: ReconstructionData
    ) -> tuple[tuple[WSRY, WSRY] | tuple[WSRY, WSRY, str, str], ...]:
        """Generate neighbours for current solution."""
        neighbours: list[tuple[WSRY, WSRY] | tuple[WSRY, WSRY, str, str]] = []
        ws_paths_set = {str(ws.path)}  # terrible solution, but it works
        not_used_not_tabu_oligos_ws = [
            oligo for oligo in ws.not_used_oligos() if not self.is_tabu(oligo)
        ]
        not_used_not_tabu_oligos_ry = [
            oligo for oligo in ry.not_used_oligos() if not self.is_tabu(oligo)
        ]

        for _ in range(self.number_of_neighbours):
            chosen_move = random.choice(list(Moves))
            # chosen_move = Moves.INSERT_OLIGO
            # chosen_move = Moves.DELETE_OLIGO
            # chosen_move = Moves.SHIFT_OLIGO
            # chosen_move = Moves.DELETE_CLUSTER
            # chosen_move = Moves.SHIFT_CLUSTER
            neighbour: Optional[tuple[WSRY, WSRY] | tuple[WSRY, WSRY, str, str]] = None
            match chosen_move:
                case Moves.INSERT_OLIGO:
                    # print("INSERT_OLIGO")
                    neighbour = self.generate_neighbour_insert_oligo(
                        ws,
                        ry,
                        r,
                        not_used_not_tabu_oligos_ws,
                        not_used_not_tabu_oligos_ry,
                    )  # seems ok
                case Moves.DELETE_OLIGO:
                    # print("DELETE_OLIGO")
                    neighbour = self.generate_neighbour_delete_oligo(
                        ws, ry, r
                    )  # seems ok
                case Moves.DELETE_CLUSTER:
                    # print("DELETE_CLUSTER")
                    neighbour = self.generate_neigbour_delete_cluster(ws, ry, r)
                case Moves.SHIFT_OLIGO:
                    # print("SHIFT_OLIGO")
                    neighbour = self.generate_neighbour_shift_oligo(
                        ws, ry, r
                    )  # seems ok
                case Moves.SHIFT_CLUSTER:
                    # print("SHIFT_CLUSTER")
                    neighbour = self.generate_neighbour_shift_cluster(
                        ws, ry, r
                    )  # can't get neigbour - is the possibility of it really low or is there a bug?

            if neighbour and str(neighbour[0].path) not in ws_paths_set:
                ws_paths_set.add(str(neighbour[0].path))
                neighbours.append(neighbour)

        return tuple(neighbours)

    def global_grade(self, ws: WSRY, ry: WSRY):
        """Calculate global grade of the solution - number of oligos in the path."""
        return len(ws.path)

    def condensation_grade(self, ws: WSRY, ry: WSRY):
        """Calculate condensation grade of the solution - number of clusters in the path."""
        return len(ws.path) / len(Tabu.reconstruct_dna(ws, ry))

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
        best_solution: tuple[WSRY, WSRY] | tuple[WSRY, WSRY, str] = copy.deepcopy(
            greedy_solution
        )
        solution_to_return: Optional[tuple[WSRY, WSRY, int, int]] = (
            None  # ws ry abs(errors) solution_length
        )

        print(best_solution)
        reconstructed_dna = Tabu.reconstruct_dna(best_solution[0], best_solution[1])
        print(reconstructed_dna, len(reconstructed_dna), f"max: {r.length}")
        for _ in range(self.number_of_iterations):
            max_grade_neighbour: Optional[
                tuple[int, tuple[WSRY, WSRY] | tuple[WSRY, WSRY, str, str]]
            ] = None

            neighbours = self.generate_neighbours(best_solution[0], best_solution[1], r)
            for neighbour in neighbours:
                if random.randint(0, 1000) % 2:
                    grade = self.global_grade(neighbour[0], neighbour[1])
                else:
                    grade = self.condensation_grade(neighbour[0], neighbour[1])

                if max_grade_neighbour is None or grade > max_grade_neighbour[0]:
                    max_grade_neighbour = (grade, neighbour)

            if max_grade_neighbour:
                reconstructed_dna = Tabu.reconstruct_dna(
                    max_grade_neighbour[1][0], max_grade_neighbour[1][1]
                )
                print(max_grade_neighbour[1][:2])
                print(reconstructed_dna, len(reconstructed_dna), f"max: {r.length}")
                print(self.count_negative_errors(max_grade_neighbour[1][0]))
                if (
                    len(reconstructed_dna) == r.length
                    and self.count_negative_errors(max_grade_neighbour[1][0]) == r.sqne
                ):
                    print("Found perfect solution")
                    return (max_grade_neighbour[1][0], max_grade_neighbour[1][1])

                if solution_to_return is None or (
                    len(reconstructed_dna) > solution_to_return[3]
                    or (
                        len(reconstructed_dna) >= solution_to_return[3]
                        and abs(
                            self.count_negative_errors(max_grade_neighbour[1][0])
                            - r.sqne
                        )
                        < solution_to_return[2]
                    )
                ):
                    solution_to_return = (
                        max_grade_neighbour[1][0],
                        max_grade_neighbour[1][1],
                        abs(
                            self.count_negative_errors(max_grade_neighbour[1][0])
                            - r.sqne
                        ),
                        len(reconstructed_dna),
                    )

                best_solution = max_grade_neighbour[1]
                if len(best_solution) == 4:
                    self.add(best_solution[2], best_solution[3])
                if self.ommit_tabu:
                    self.ommit_tabu = False
            else:
                print("No better solution found")
                self.ommit_tabu = True
            if len(self.tabu_list_ws) > self.tabu_size:
                self.tabu_list_ws.pop(0)
                self.tabu_list_ry.pop(0)
                # add move to tabu list
        return (solution_to_return[0], solution_to_return[1])
