"""WSRY class to store data for WS or RY cells and provide methods to work with them."""


class WSRY:
    """Store data for WS or RY cells and provide methods to work with them."""

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

    def not_used_oligos(self) -> tuple[str, ...]:
        """Return oligos that are not used in the solution."""
        return tuple([oligo for oligo in self.cells_dict if not self.cells_dict[oligo]])

    def __repr__(self) -> str:
        return f"Start: {self.start_converted} Path: {self.path} Depth: {self.depth}"
