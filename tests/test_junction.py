from splicemap.dataclasses import Junction


def test_Junction():
    junc = Junction.from_str('1:10100-10500:+')

    assert junc.donor == 10100
    assert junc.acceptor == 10500

    donor = junc.donor_region()
    assert donor.start == 10000
    assert donor.end == 10200

    acceptor = junc.acceptor_region()
    assert acceptor.start == 10400
    assert acceptor.end == 10600

    junc = Junction.from_str('1:10100-10500:-')

    assert junc.donor == 10500
    assert junc.acceptor == 10100

    acceptor = junc.acceptor_region()
    assert acceptor.start == 10000
    assert acceptor.end == 10200

    donor = junc.donor_region()
    assert donor.start == 10400
    assert donor.end == 10600


def test_Junction_dinucleotide_region():
    junc = Junction.from_str('1:10100-10500:+')
    dinucleotide = junc.dinucleotide_region()
    assert dinucleotide[0].chrom == '1'
    assert dinucleotide[0].start == 10100
    assert dinucleotide[0].end == 10102
    assert dinucleotide[0].strand == '.'

    assert dinucleotide[1].chrom == '1'
    assert dinucleotide[1].start == 10498
    assert dinucleotide[1].end == 10500
    assert dinucleotide[1].strand == '.'
