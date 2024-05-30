"""Reconstruction data class"""


class Probe:
    """Probe class"""

    def __init__(self, json_xml: dict):
        self.pattern = json_xml["@pattern"]
        self.cells = json_xml["cell"]  # sx - s1, s2

    def __repr__(self):
        return f"Pattern: {self.pattern}, Cells: {self.cells[:3]}"

    # pattern: str = None
    # cells: list[str] = None  # sx - s1, s2


class ReconstructionData:
    """Reconstruction data class"""

    def __init__(self, json: dict):
        self.length = int(json["dna"]["@length"])
        self.start = json["dna"]["@start"]
        # self.probes = []

        # for probe in json["dna"]["probe"]:
        #     self.probes.append(Probe(probe))
        self.ws_probe = Probe(json["dna"]["probe"][0])
        self.ry_probe = Probe(json["dna"]["probe"][1])
        sis = self.length - len(self.ws_probe.cells[0]) + 1
        self.sqne = (sis - len(self.ws_probe.cells)) + (sis - len(self.ry_probe.cells))

    def __repr__(self):
        return f"Length: {self.length}, Start: {self.start}, Probes:[ {self.ws_probe}, {self.ry_probe} ]"

    # length: int = None
    # start: str = None
    # # probes: list[Probe] = None
    # WS_probe: Probe = None
    # RY_probe: Probe = None
