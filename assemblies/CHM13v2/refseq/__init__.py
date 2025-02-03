from pathlib import Path
from typing import Literal

from attr import define, field
from biobit.toolkit import annotome as at

ROOT = Path(__file__).parent

gff = ROOT / "CHM13v2.refseq.gff.gz"
index = ROOT / "CHM13v2.refseq.annotome.pkl"

GeneSource = Literal[
    'Gnomon', 'Curated Genomic', 'tRNAscan-SE', 'Curated Genomic%2Ccmsearch',
    'BestRefSeq', 'BestRefSeq%2CGnomon', 'cmsearch'
]
GeneBiotype = Literal[
    'rRNA', 'misc_RNA', 'D_segment', 'ncRNA_pseudogene', 'antisense_RNA', 'ncRNA', 'Y_RNA', 'C_region', 'snRNA',
    'lncRNA', 'miRNA', 'V_segment_pseudogene', 'tRNA', 'scRNA', 'other', 'transcribed_pseudogene', 'V_segment',
    'vault_RNA', 'snoRNA', 'C_region_pseudogene', 'J_segment', 'protein_coding', 'pseudogene', 'RNase_P_RNA',
    'J_segment_pseudogene', 'telomerase_RNA', 'RNase_MRP_RNA'
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


RNASource = Literal['Gnomon', 'Curated Genomic', 'tRNAscan-SE', 'BestRefSeq', 'cmsearch']
RNABiotype = Literal[
    'rRNA', 'antisense_RNA', 'Y_RNA', 'unknown', 'V_gene_segment', 'RNase_P_RNA',
    'snRNA', 'J_gene_segment', 'mRNA', 'miRNA', 'tRNA', 'miRNA_primary_transcript', 'scRNA', 'lnc_RNA_pseudogene',
    'scaRNA', 'C_gene_segment', 'vault_RNA', 'snoRNA', 'RNase_MRP_RNA', 'V_gene_segment_pseudogene', 'pseudogene',
    'lnc_RNA', 'C_gene_segment_pseudogene', 'telomerase_RNA', 'J_gene_segment_pseudogene', 'D_gene_segment',
]
RNAExperiment = Literal[
    'COORDINATES: polyA evidence [ECO:0006239]',
    'COORDINATES: cap analysis [ECO:0007248]',
    'COORDINATES: cap analysis [ECO:0007248] and polyA evidence [ECO:0006239]'
]
RNATag = Literal['RefSeq Select', 'RefSeq Plus Clinical']


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
    'Gnomon', 'Curated Genomic', 'BestRefSeq'
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
