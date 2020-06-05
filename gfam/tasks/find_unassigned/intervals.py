def _substract_interval_from_interval(interval1, interval2):
    """Substract from an interval (left, right) another interval
    (leftr,rightr), returning a list of resulting intervals
    (which may be empty) equivalent to the set difference
    """
    (left, right) = interval1
    (leftr, rightr) = interval2

    if leftr > right or rightr < left:
        return [(left, right)]
    elif leftr <= left:
        if rightr >= right:
            return []
        else:
            return [(rightr+1, right)]
    elif rightr >= right:
        # leftr is > left
        return [(left, leftr-1)]
    else:
        # case leftr > left and rightr < right
        return [(left, leftr-1), (rightr+1, right)]


def substract(region, list_to_substract):
    (reg_id, left, right) = region

    # we obtain those regions which overlaps with our region
    # if our region is (left,right) and the set of iterated regions is
    # (leftr, rightr), the overlapping are those (leftr,rightr) where
    subtrahend = [(leftr, rightr) for (leftr, rightr)
                  in list_to_substract
                  if leftr <= right and rightr >= left]

    if not subtrahend:
        # no overlap, the region is returned "as is"
        return [region]
    minuend = [(left, right)]
    for sub in subtrahend:
        minuend.extend(
            _substract_interval_from_interval(minuend.pop(), sub))
        if not minuend:
            return []
    return [(reg_id, le, ri) for (le, ri) in minuend]
