import requests
import xmltodict
import copy
from reconstruction_data import ReconstructionData


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
    # convert oligo to WS or RY according to the given dictionary without last nucleotide
    def convert_oligo(self, oligo: str) -> str:
        half = ""
        for i in range(len(oligo) - 1):
            half += self.dict_convertion[oligo[i]]
        half += oligo[-1]

        return half

    def connect_WS_RY(oligo_WS: str, oligo_RY) -> str:
        connected = ""
        for nucleotide_WS, nucleotide_RY in zip(oligo_WS, oligo_RY):
            temp = nucleotide_WS + nucleotide_RY
            if temp == "SR":
                connected += "G"
            elif temp == "SY":
                connected += "C"
            elif temp == "WR":
                connected += "A"
            elif temp == "WY":
                connected += "T"
        return connected

    def __init__(self, dict_convertion: dict, oligo: str, cells: dict):
        self.dict_convertion = dict_convertion
        self.start_converted = self.convert_oligo(oligo)

        self.cells_dict = {}
        for cell in cells:
            self.cells_dict[cell] = False
        self.cells_dict[self.start_converted] = True

        self.path = [self.start_converted]
        self.depth = [0]

    def __repr__(self) -> str:
        return f"Start: {self.start_converted} Path: {self.path} Depth: {self.depth}"

    start_converted: str = None  # set_first - Z_set_first
    dict_convertion: dict = None
    cells_dict: dict = None  # ols
    path: list[str] = None
    depth: list[int] = None


def fetch_test_data(
    n: int = 500,
    k: int = 10,
    mode: str = "binary",
    intensity: int = 0,
    position: int = 0,
    sqpe: int = 0,
    sqne: int = 100,
    pose: int = 0,
):
    content = requests.get(
        f"https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php?n={n}&k={k}&mode={mode}&intensity={intensity}&position={position}&sqpe={sqpe}&sqne={sqne}&pose={pose}"
    ).content
    data = xmltodict.parse(content)
    return ReconstructionData(data)


def check_overlap(oligo1, oligo2, probe):
    max_overlap = 0
    for offset in range(probe - 1, 0, -1):
        if oligo1[probe - offset :] == oligo2[:offset]:
            return max_overlap


def search_overlapings(
    oligo: str,
    cells: list[str],
    cells_dict: dict[str, bool],
    overlapings: dict[str, int],
):
    for cell in cells:
        max_overlap = check_overlap(oligo, cell, len(cell))
        if not cells_dict[cell] and not max_overlap:
            overlapings[cell] = max_overlap


def reconstruct(
    is_reconstructed: bool,
    reconstructed_dna_length: int,
    reconstructed_dna: str,
    ws: WSRY,
    ry: WSRY,
    rd: ReconstructionData,
):
    _ws = copy.deepcopy(ws)
    _ry = copy.deepcopy(ry)
    _rd = copy.deepcopy(rd)
    if is_reconstructed:
        return (
            is_reconstructed,
            reconstructed_dna_length,
            reconstructed_dna,
            _ws,
            _ry,
            _rd,
        )
    if _rd.length == reconstructed_dna_length:
        is_reconstructed = True
        reconstructed_dna = WSRY.connect_WS_RY(_ws.start_converted, _ry.start_converted)
        return (
            is_reconstructed,
            reconstructed_dna_length,
            reconstructed_dna,
            _ws,
            _ry,
            _rd,
        )

    _ws.start_converted[-1] = nucleotide_to_weak_strong[
        _ws.start_converted[-1]
    ]  # change last nucleotide to WS
    _ry.start_converted[-1] = nucleotide_to_purine_pyrimidine[
        _ry.start_converted[-1]
    ]  # change last nucleotide to RY

    overlapings_ws = dict()
    overlapings_ry = dict()
    search_overlapings(
        _ws.start_converted, _rd.WS_probe.cells, _ws.cells_dict, overlapings_ws
    )
    search_overlapings(
        _ry.start_converted, _rd.RY_probe.cells, _ry.cells_dict, overlapings_ry
    )
    if not overlapings_ry or not overlapings_ws:
        # what to return here?
        return (
            is_reconstructed,
            reconstructed_dna_length,
            reconstructed_dna,
            _ws,
            _ry,
            _rd,
        )


def main():
    r = fetch_test_data()
    ws = WSRY(nucleotide_to_weak_strong, r.start, r.WS_probe.cells)
    ry = WSRY(nucleotide_to_purine_pyrimidine, r.start, r.RY_probe.cells)
    # print(ws)
    # print(ry)
    # print(r)
    print(reconstruct(True, 0, "", ws, ry, r))


main()
