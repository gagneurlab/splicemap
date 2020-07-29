from typing import List
from kipoiseq import Interval


class Junction(Interval):

    @property
    def acceptor(self):
        return self.start if self.strand == '-' else self.end

    @property
    def donor(self):
        return self.end if self.strand == '-' else self.start

    def dinucleotide_region(self):
        return Interval(self.chrom, self.start, self.start + 2), \
            Interval(self.chrom, self.end - 2, self.end)

    def acceptor_region(self, overhang=(100, 100)):
        return Interval(self.chrom, self.acceptor,
                        self.acceptor, strand=self.strand) \
            .slop(upstream=overhang[0], downstream=overhang[1])

    def donor_region(self, overhang=(100, 100)):
        return Interval(self.chrom, self.donor,
                        self.donor, strand=self.strand) \
            .slop(upstream=overhang[0], downstream=overhang[1])


class Event:

    def __init__(self, junctions: List):
        self.junctions = junctions

    @classmethod
    def from_str(cls, event_str: str):
        return cls([
            Junction.from_str(junc)
            for junc in event_str.split(';')
        ])


class EventPSI5(Event):

    def __init__(self, junctions: List):
        super().__init__(junctions)
        assert all(
            i.donor == self.junctions[0].donor
            for i in self.junctions
        )

    @property
    def donor_str(self):
        j = self.junctions[0]
        return f'{j.chrom}:{j.donor}:{j.strand}'


class EventPSI3(Event):

    def __init__(self, junctions: List):
        super().__init__(junctions)
        assert all(
            i.acceptor == self.junctions[0].acceptor
            for i in self.junctions
        )

    @property
    def acceptor_str(self):
        j = self.junctions[0]
        return f'{j.chrom}:{j.acceptor}:{j.strand}'
