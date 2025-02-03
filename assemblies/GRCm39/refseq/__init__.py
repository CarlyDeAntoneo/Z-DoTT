from pathlib import Path
from typing import Literal

from attr import define, field
from biobit.toolkit import annotome as at

ROOT = Path(__file__).parent

gff = ROOT / "GRCm39.refseq.gff.gz"
index = ROOT / "GRCm39.refseq.annotome.pkl"

GeneSource = Literal[
    'Gnomon', 'RefSeq', 'BestRefSeq', 'BestRefSeq%2CGnomon', 'Curated Genomic',
    'tRNAscan-SE', 'BestRefSeq%2Ccmsearch', 'cmsearch'
]
GeneBiotype = Literal[
    'Y_RNA', 'J_segment_pseudogene', 'telomerase_RNA', 'ncRNA', 'snRNA', 'J_segment', 'V_segment', 'lncRNA',
    'antisense_RNA', 'snoRNA', 'ncRNA_pseudogene', 'RNase_MRP_RNA', 'RNase_P_RNA', 'misc_RNA', 'D_segment', 'miRNA',
    'other', 'C_region', 'tRNA', 'protein_coding', 'pseudogene', 'transcribed_pseudogene', 'rRNA',
    'D_segment_pseudogene', 'C_region_pseudogene', 'scRNA', 'V_segment_pseudogene'
]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrGene:
    source: GeneSource
    name: str
    description: str | None
    biotype: GeneBiotype
    partial: bool
    synonyms: frozenset[str] = field(converter=lambda x: frozenset(x))

    def __attrs_post_init__(self):
        if self.source not in GeneSource.__args__:
            raise ValueError(f"Invalid source: {self.source}")
        if self.biotype not in GeneBiotype.__args__:
            raise ValueError(f"Invalid biotype: {self.biotype}")


RNASource = Literal['Gnomon', 'RefSeq', 'BestRefSeq', 'Curated Genomic', 'tRNAscan-SE', 'cmsearch']
RNABiotype = Literal[
    'miRNA_primary_transcript', 'scaRNA', 'Y_RNA', 'telomerase_RNA', 'unknown', 'snRNA', 'antisense_RNA', 'snoRNA',
    'C_gene_segment_pseudogene', 'V_gene_segment', 'RNase_MRP_RNA', 'J_gene_segment_pseudogene', 'RNase_P_RNA',
    'lnc_RNA', 'D_gene_segment_pseudogene', 'mRNA', 'miRNA', 'D_gene_segment', 'V_gene_segment_pseudogene', 'tRNA',
    'pseudogene', 'rRNA', 'lnc_RNA_pseudogene', 'C_gene_segment', 'scRNA', 'J_gene_segment'
]
RNAExperiment = Literal[
    'COORDINATES: cap analysis [ECO:0007248]',
    'COORDINATES: polyA evidence [ECO:0006239]',
    'COORDINATES: cap analysis [ECO:0007248] and polyA evidence [ECO:0006239]',
]
RNATag = Literal['MANE Select', 'RefSeq Select', 'MANE Plus Clinical']


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrRNA:
    source: RNASource
    name: str | None
    product: str | None
    partial: bool
    biotype: RNABiotype
    tags: frozenset[RNATag]
    experiment: RNAExperiment | None

    def __attrs_post_init__(self):
        if self.source not in RNASource.__args__:
            raise ValueError(f"Invalid source: {self.source}")
        if self.biotype not in RNABiotype.__args__:
            raise ValueError(f"Invalid biotype: {self.biotype}")
        if any(x not in RNATag.__args__ for x in self.tags):
            raise ValueError(f"Invalid tags: {self.tags}")
        if self.experiment and self.experiment not in RNAExperiment.__args__:
            raise ValueError(f"Invalid experiment: {self.experiment}")


CDSSource = Literal[
    'Curated Genomic', 'Gnomon', 'BestRefSeq', 'RefSeq'
]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrCDS:
    source: CDSSource
    partial: bool
    product: str | None
    transcripts: frozenset[str] = field(converter=lambda x: frozenset(x))

    def __attrs_post_init__(self):
        if self.source not in CDSSource.__args__:
            raise ValueError(f"Invalid source: {self.source}")


def load() -> at.Annotome[AttrGene, AttrRNA, AttrCDS]:
    return at.read_pkl(index.as_posix())
