class Probe:
    def __init__(self, json_xml: dict):
        self.pattern = json_xml["@pattern"]
        self.cells = json_xml["cell"]

    def __repr__(self):
        return f"Pattern: {self.pattern}, Cells: {self.cells[:3]}"

    pattern: str = None
    cells: list[str] = None  # sx - s1, s2


class ReconstructionData:
    def __init__(self, json: dict):
        self.length = json["dna"]["@length"]
        self.start = json["dna"]["@start"]
        # self.probes = []

        # for probe in json["dna"]["probe"]:
        #     self.probes.append(Probe(probe))
        self.WS_probe = Probe(json["dna"]["probe"][0])
        self.RY_probe = Probe(json["dna"]["probe"][1])

    def __repr__(self):
        return f"Length: {self.length}, Start: {self.start}, Probes:[ {self.WS_probe}, {self.RY_probe} ]"

    length: int = None
    start: str = None
    # probes: list[Probe] = None
    WS_probe: Probe = None
    RY_probe: Probe = None
