import requests
import xmltodict
from typing import Optional
from reconstruction_data import ReconstructionData


def fetch_test_data(
    n: int = 16,
    k: int = 4,
    mode: str = "binary",
    intensity: int = 0,
    position: int = 0,
    sqpe: int = 0,
    sqne: int = 0,
    pose: int = 0,
) -> ReconstructionData:
    """Fetches test data from server"""
    content = requests.get(
        f"https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php?n={n}&k={k}&mode={mode}&intensity={intensity}&position={position}&sqpe={sqpe}&sqne={sqne}&pose={pose}",
        timeout=10,
    ).content
    if not content:
        raise requests.exceptions.RequestException("Failed to fetch test data")
    data = xmltodict.parse(content)
    return ReconstructionData(data)
