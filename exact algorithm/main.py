import requests
import xmltodict
import pprint
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

    def connect_WS_RY(self, oligo_WS: str, oligo_RY) -> str:
        connected = ""
        for nucleotide_WS, nucleotide_RY in zip(oligo_WS, oligo_RY):
            temp = nucleotide_WS + nucleotide_RY
            if temp == "SR":
                connected += "A"
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

    start_converted: str = None
    dict_convertion: dict = None
    cells_dict: dict = None
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


def main():
    r = fetch_test_data()
    ws = WSRY(nucleotide_to_weak_strong, r.start, r.WS_probe.cells)
    ry = WSRY(nucleotide_to_purine_pyrimidine, r.start, r.RY_probe.cells)
    print(ws)
    print(ry)
    print(r)


main()
