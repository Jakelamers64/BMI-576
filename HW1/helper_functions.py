def overlap_length(left, right):
    """Returns the length of the longest suffix of left that is a prefix of right

    Args:
        left: a string
        right: a string
    Returns:
        An integer length of the longest overlap (0 if there is no overlap)
    """
    if left == right:
        return len(left)
    else:
        for i in range(0,len(left))[::-1]:
            if right.startswith(left[(len(left) - i - 1):(len(left))]):
                return len(left[(len(left) - i - 1):(len(left))])

        return 0
