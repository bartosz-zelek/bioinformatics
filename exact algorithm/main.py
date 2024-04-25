import requests
import xmltodict
import pprint
from reconstruction_data import ReconstructionData


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
    print(r)


main()
