"""Empirical exact-input squared-error floor. Basic decomposition, not novelty.

Caller supplies exact full-input byte keys for fixed observed rows. No labels
are used to choose groups. Nothing here certifies the keys' scientific meaning.
"""
from collections import defaultdict
from math import fsum, isfinite
from numbers import Real


def collision_mse_floor(keys, labels):
    keys, labels = list(keys), list(labels)
    if not keys or len(keys) != len(labels):
        raise ValueError('nonempty equal-length keys and labels required')
    groups = defaultdict(list)
    for key, label in zip(keys, labels):
        if not isinstance(key, bytes):
            raise TypeError('keys must be exact full-input bytes, not inferred identifiers')
        if isinstance(label, bool) or not isinstance(label, Real):
            raise TypeError('labels must be finite real numbers, not booleans')
        label = float(label)
        if not isfinite(label):
            raise ValueError('non-finite label')
        groups[key].append(label)
    sums = []
    ambiguous = 0
    for values in groups.values():
        # Shift first to avoid summing many large, nearly equal labels.
        offset = values[0]
        try:
            mean = offset + fsum((v-offset)/len(values) for v in values)
            ss = fsum((v-mean)**2 for v in values)
        except (OverflowError, ValueError):
            raise ValueError('numerical range exceeds finite squared-error calculation') from None
        if not isfinite(ss):
            raise ValueError('non-finite squared error')
        sums.append(ss)
        ambiguous += len(set(values)) > 1
    try:
        floor = fsum(x/len(keys) for x in sums)
    except OverflowError:
        raise ValueError('non-finite total squared error') from None
    if not isfinite(floor):
        raise ValueError('non-finite floor')
    return {'n_rows':len(keys), 'n_exact_input_groups':len(groups),
            'n_conflicting_groups':ambiguous, 'empirical_mse_floor':floor,
            'scope':'Fixed observed rows; deterministic predictions constant within supplied exact full-input groups. Not a population bound or a benchmark gate.'}
