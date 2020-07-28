from count_table.dataclasses import Junction

def get_variants_around_junction(vcf, junction, sample_id=None,
                                 overhang=(100, 100), variant_filter=None):
    junction = Junction.from_str(junction) if type(
        junction) == str else junction

    acceptor = junction.acceptor_region(overhang=overhang)
    donor = junction.donor_region(overhang=overhang)
    query = vcf.query_variants([donor, acceptor], sample_id)

    if variant_filter:
        query = query.filter(variant_filter)
    variant_intervals = list(query.variant_intervals)

    return (
        list(variant_intervals[0][0]),
        list(variant_intervals[1][0])
    )