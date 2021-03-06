import pysam

from collections import Counter, namedtuple

code2str = {
    0: "MATCH",  # M - alignment match (can be a sequence match or mismatch)
    1: "INS",  # I - insertion to the reference
    2: "DEL",  # D - deletion from the reference
    3: "SKIP",  # N - skipped region from the reference
    4: "SOFT_CLIP",  # S - soft clipping (clipped sequences present in SEQ)
    5: "HARD_CLIP",  # H - hard clipping (clipped sequences NOT present in SEQ)
    6: "PAD",  # P - padding (silent deletion from padded reference)
    7: "EQUAL",  # = - sequence match
    8: "DIFF",  # X - sequence mismatch
}


def cigar_profile(cigar_tuples):
    cigar_prof = Counter()
    for cigar_tuple in cigar_tuples:
        cigar_operation = code2str[cigar_tuple[0]]
        cigar_prof[cigar_operation] += cigar_tuple[1]
    return cigar_prof


def get_cigar_counts(alignment):
    cigartuples = alignment.cigartuples
    if cigartuples is None:
        return None

    return cigar_profile(cigartuples)


def get_MID_from(cigar_counts):
    return cigar_counts["INS"] + cigar_counts["DEL"] + cigar_counts["MATCH"]


def get_clip_size(cigar_counts):
    return cigar_counts["SOFT_CLIP"] + cigar_counts["HARD_CLIP"]


def get_edit_distance_from(alignment):
    return alignment.get_tag("NM")


def get_num_matches(alignment, cigar_counts):
    edit_dist = get_edit_distance_from(alignment)
    return edit_dist - cigar_counts["INS"] - cigar_counts["DEL"]


AlignmentStats = namedtuple("AlignmentStats",
                            "edit_distance MID query_start query_end "
                            "ref_start ref_end clipped")


def get_all_alignments_from(file_path):
    file = pysam.Samfile(file_path)

    res = []
    for alignment in file:
        cigar_counts = get_cigar_counts(alignment)
        if cigar_counts is None:
            continue
        MID = get_MID_from(cigar_counts)
        clipped = get_clip_size(cigar_counts)

        query_start = alignment.query_alignment_start
        query_end = alignment.query_alignment_end
        query_len = alignment.query_alignment_length
        ref_start = alignment.reference_start
        ref_end = alignment.reference_end
        ref_len = alignment.reference_length
        edit_distance = alignment.get_tag("NM")

        res.append(AlignmentStats(edit_distance, MID, query_start, query_end,
                                  ref_start, ref_end, clipped))

    return res


def get_best_alignment(file_path):
    """ The best alignment is the alignment with the greatest
    number of matches. """
    file = pysam.Samfile(file_path)

    best_alignment = None
    greatest_num_matches = 0
    for alignment in file:
        cigar_counts = get_cigar_counts(alignment)
        if cigar_counts is None:
            continue

        MID = get_MID_from(cigar_counts)
        edit_distance = get_edit_distance_from(alignment)
        clipped = get_clip_size(cigar_counts)
        matches = get_num_matches(alignment, cigar_counts)

        if best_alignment is None or greatest_num_matches < matches:
            best_alignment = alignment
            greatest_num_matches = matches

    return best_alignment
