"""Exact two-block nonorthogonal union; unsupported implicit backends reject."""
from ..core import check_params
from ..compatibility import SUPPORTED_KINDS, certified_union, verify_glued_artifact


def union(inputs, params, context):
    check_params(params, ("cross_gram", "scales", "heads"))
    if len(inputs) != 2:
        raise ValueError("gluing.union requires two source artifacts")
    if "cross_gram" not in params:
        raise ValueError("gluing.union requires explicit unscaled cross_gram entries; unrelated source frames are not inferred")
    return certified_union(inputs[0], inputs[1], params["cross_gram"], context,
                           scales=params.get("scales"), heads=params.get("heads"))


verify_artifact = verify_glued_artifact
OPERATORS = {"gluing.union": union}
INPUT_KINDS = {"gluing.union": (SUPPORTED_KINDS, SUPPORTED_KINDS)}
