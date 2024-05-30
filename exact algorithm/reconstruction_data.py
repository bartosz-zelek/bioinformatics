"""Reconstruction data class"""


class Probe:
    """Probe class"""

    def __init__(self, json_xml: dict):
        self.pattern = json_xml["@pattern"]
        self.cells = json_xml["cell"]  # sx - s1, s2

    def __repr__(self):
        return f"Pattern: {self.pattern}, Cells: {self.cells}"


class ReconstructionData:
    """Reconstruction data class"""

    def __init__(self, json: dict):
        self.length = int(json["dna"]["@length"])
        self.start = json["dna"]["@start"]

        self.ws_probe = Probe(json["dna"]["probe"][0])
        self.ry_probe = Probe(json["dna"]["probe"][1])
        sis = self.length - len(self.ws_probe.cells[0]) + 1
        self.sqne = (sis - len(self.ws_probe.cells)) + (sis - len(self.ry_probe.cells))

    def __repr__(self):
        return f"Length: {self.length}, Start: {self.start}, Probes:[ {self.ws_probe},\n\t {self.ry_probe} ]"
