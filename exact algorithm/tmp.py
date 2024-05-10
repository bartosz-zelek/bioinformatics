#!/bin/python3


def find_overlaps(oligo1, oligos):
    overlaps = []
    for oligo2 in oligos:
        if oligo1[1:] == oligo2[:-1]:
            overlaps.append(oligo2)
    return overlaps


def construct_sequence(solution, oligos):
    if not oligos:
        yield solution
    overlaps = find_overlaps(solution[-1], oligos)
    if not overlaps:
        yield solution
    for oligo in overlaps:
        if oligo not in solution:
            solution_cpy = solution.copy()
            solution_cpy.append(oligo)

            oligos_cpy = oligos.copy()
            oligos_cpy.remove(oligo)
            yield from construct_sequence(solution_cpy, oligos_cpy)


if __name__ == "__main__":
    keyLength, _length = input().rstrip().split()
    keyStart, start = input().rstrip().split()
    keyProbe, _probe = input().rstrip().split()
    probe = int(_probe)
    length = int(_length)

    oligonucleotides = [input() for _ in range(length - probe + 1)]
    start_lst = [start]

    try:
        solution = sorted(
            list(construct_sequence(start_lst, oligonucleotides)),
            key=lambda x: len(x),
            reverse=True,
        )[0]
        dna = ""
        for oligo in solution:
            dna += oligo[0]
        dna += solution[-1][1 - probe :]
        print(dna)
    except ValueError as e:
        print("No")
        print(str(e))
