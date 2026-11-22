"""Verify coordinate-weight obstructions for the two frozen binary sign codes."""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def binary_span(generators):
    words = [0]
    for generator in generators:
        words.extend([word ^ generator for word in words])
    return words


def check(document):
    generators = document["generator"]
    if len(generators) != 22 or any(type(row) is not int or not 0 < row < 1 << 38
                                    for row in generators):
        raise ValueError("Invalid 38-coordinate generator")
    words = binary_span(generators)
    if len(words) != len(set(words)):
        raise ValueError("The generator lacks full rank")
    for case in document["shortenings"]:
        dimension = case["dimension"]
        load = case["coordinate_load"]
        if dimension not in (34, 35) or type(load) is not int or load <= 0:
            raise ValueError("Invalid shortening or load")
        short = {word for word in words if word < 1 << dimension}
        if len(short) != 1 << (dimension - 16):
            raise ValueError("Wrong shortening dimension")
        minimum = min(word.bit_count() for word in short if word)
        if minimum != 8:
            raise ValueError("Expected Hamming distance eight")
        number_of_contacts = sum(word.bit_count() == 8 for word in short)
        if number_of_contacts != case["nonzero_weight_eight"]:
            raise ValueError("Incorrect number of weight-eight words")

        column_loads = [0] * dimension
        coefficient_sum = 0
        for word, multiplicity in case["witness"]:
            if word not in short or word.bit_count() != 8:
                raise ValueError("Obstruction uses a word outside the weight-eight code")
            if type(multiplicity) is not int or multiplicity <= 0:
                raise ValueError("Multiplicities must be positive integers")
            coefficient_sum += multiplicity
            for axis in range(dimension):
                if word & (1 << axis):
                    column_loads[axis] += multiplicity
        if column_loads != [load] * dimension or dimension * load != 8 * coefficient_sum:
            raise ValueError("Weighted supports do not balance the coordinates")
        if 4 * load >= coefficient_sum:
            raise ValueError("This balance does not force a strict violation")
        print(json.dumps({
            "dimension": dimension,
            "binary_code_size": len(short),
            "weight_eight_words": number_of_contacts,
            "balanced_rows": len(case["witness"]),
            "maximum_possible_minimum_weighted_distance_as_fraction_of_total": f"8/{dimension}",
            "required_fraction": "1/4",
        }))


if __name__ == "__main__":
    check(json.loads((HERE / "weighted_sign_witness.json").read_text()))
