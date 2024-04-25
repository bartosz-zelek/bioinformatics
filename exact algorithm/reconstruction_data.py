class Probe:
    def __init__(self, json_xml: dict):
        self.pattern = json_xml["@pattern"]
        self.cells = json_xml["cell"]

    def __repr__(self):
        return f"Pattern: {self.pattern}, Cells: {self.cells[:3]}"

    pattern: str = None
    cells: list[str] = None


class ReconstructionData:
    def __init__(self, json: dict):
        self.length = json["dna"]["@length"]
        self.start = json["dna"]["@start"]
        self.probes = []

        for probe in json["dna"]["probe"]:
            self.probes.append(Probe(probe))

    def __repr__(self):
        return f"Length: {self.length}, Start: {self.start}, Probes: {self.probes}"

    length: int = None
    start: str = None
    probes: list[Probe] = None
