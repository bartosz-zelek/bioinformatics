"""Common functions for metaheuristic algorithms."""

import requests
import xmltodict
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


def fetch_test_data(
    filename: str,
    # n: int = 16,
    # k: int = 4,
    # mode: str = "binary",
    # intensity: int = 0,
    # position: int = 0,
    # sqpe: int = 0,
    # sqne: int = 0,
    # pose: int = 0,
) -> ReconstructionData:
    """Fetches test data from server"""
    # content = requests.get(
    #     f"https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php?n={n}&k={k}&mode={mode}&intensity={intensity}&position={position}&sqpe={sqpe}&sqne={sqne}&pose={pose}",
    #     timeout=10,
    # ).content
    # file = open(f"../test_data/n{n}k{k}sqne{sqne}.xml", "r", encoding="utf-8")
    file = open(filename, "r", encoding="utf-8")
    content = file.read().encode("utf-8")
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
