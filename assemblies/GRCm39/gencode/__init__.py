from pathlib import Path
from typing import Literal

from attr import define, field
from biobit.toolkit import annotome as at

ROOT = Path(__file__).parent

gff3 = ROOT / "gencode.vM36.primary_assembly.annotation.gff3.gz"
gtf = ROOT / "gencode.vM36.primary_assembly.annotation.gtf.gz"
index = ROOT / "GRCm39.gencode.M36.annotome.pkl"

GeneSource = Literal['ENSEMBL', 'HAVANA']

GeneType = Literal[
    'TEC', 'transcribed_unprocessed_pseudogene', 'IG_D_gene', 'TR_V_gene', 'transcribed_unitary_pseudogene',
    'IG_C_pseudogene', 'IG_V_gene', 'sRNA', 'IG_LV_gene', 'Mt_tRNA', 'IG_pseudogene', 'snoRNA', 'Mt_rRNA', 'rRNA',
    'unprocessed_pseudogene', 'protein_coding', 'lncRNA', 'TR_C_gene', 'scRNA', 'translated_unprocessed_pseudogene',
    'transcribed_processed_pseudogene', 'TR_J_gene', 'IG_J_gene', 'TR_V_pseudogene', 'ribozyme', 'snRNA',
    'IG_V_pseudogene', 'IG_D_pseudogene', 'pseudogene', 'misc_RNA', 'TR_J_pseudogene', 'unitary_pseudogene',
    'processed_pseudogene', 'IG_C_gene', 'miRNA', 'scaRNA', 'TR_D_gene'
]
GeneLevel = Literal[1, 2, 3]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrGene:
    source: GeneSource
    level: GeneLevel
    name: str
    type: GeneType

    def __attrs_post_init__(self):
        if self.type not in GeneType.__args__:
            raise ValueError(f"Invalid biotype: {self.type}")
        if self.source not in GeneSource.__args__:
            raise ValueError(f"Invalid source: {self.source}")


RNASource = Literal['ENSEMBL', 'HAVANA']

RNAType = Literal[
    'IG_C_pseudogene', 'TR_V_gene', 'TR_D_gene', 'Mt_rRNA', 'nonsense_mediated_decay', 'pseudogene',
    'unitary_pseudogene', 'non_stop_decay', 'TR_V_pseudogene', 'IG_V_pseudogene', 'transcribed_unprocessed_pseudogene',
    'rRNA', 'protein_coding_LoF', 'TR_J_gene', 'IG_D_gene', 'unprocessed_pseudogene', 'sRNA', 'misc_RNA',
    'protein_coding_CDS_not_defined', 'protein_coding', 'miRNA', 'IG_V_gene', 'transcribed_unitary_pseudogene',
    'IG_LV_gene', 'IG_C_gene', 'retained_intron', 'IG_D_pseudogene', 'IG_J_gene', 'lncRNA', 'TEC', 'scaRNA',
    'processed_pseudogene', 'IG_pseudogene', 'TR_J_pseudogene', 'translated_unprocessed_pseudogene', 'snoRNA',
    'transcribed_processed_pseudogene', 'ribozyme', 'processed_transcript', 'Mt_tRNA', 'scRNA', 'TR_C_gene', 'snRNA'
]
RNATag = Literal[
    'Ensembl canonical', 'GENCODE basic',
    'cds_end_NF', 'downstream_ATG', 'upstream_ATG', 'retained_intron_final', 'retained_intron_first',
    'dotter_confirmed', 'not_best_in_genome_evidence', 'alternative_3_UTR', 'cds_start_NF',
    '5_standard_supported_extension', 'CAGE_supported_TSS', 'appris_alternative_1', 'overlapping_uORF',
    'non_canonical_U12', 'stop_codon_readthrough', 'retained_intron_CDS', 'alternative_5_UTR',
    'RNA_Seq_supported_partial', '3_standard_supported_extension', 'inferred_exon_combination', 'RP_supported_TIS',
    'non_canonical_polymorphism', '5_nested_supported_extension', 'non_canonical_conserved',
    'non_submitted_evidence', 'NAGNAG_splice_site', 'appris_principal_3', 'appris_principal_2', 'non_ATG_start',
    'RNA_Seq_supported_only', 'not_organism_supported', 'appris_principal_4', 'readthrough_transcript',
    'seleno', 'NMD_likely_if_extended', 'appris_alternative_2', 'non_canonical_genome_sequence_error', 'TAGENE',
    'inferred_transcript_model', 'sequence_error', 'mRNA_end_NF', 'appris_principal_1', 'non_canonical_other',
    'low_sequence_quality', '3_nested_supported_extension', 'upstream_uORF', 'NMD_exception', 'mRNA_start_NF',
    'exp_conf', 'appris_principal_5', 'non_canonical_TEC', 'CCDS'
]
RNATSL = Literal[1, 2, 3, 4, 5, None]
RNALevel = Literal[1, 2, 3]


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrRNA:
    source: RNASource
    level: RNALevel
    name: str
    type: RNAType
    tags: frozenset[RNATag]
    TSL: RNATSL
    CDS: frozenset[str]

    def __attrs_post_init__(self):
        if self.source not in RNASource.__args__:
            raise ValueError(f"Invalid source: {self.source}")
        if self.type not in RNAType.__args__:
            raise ValueError(f"Invalid biotype: {self.type}")
        if any(x not in RNATag.__args__ for x in self.tags):
            raise ValueError(f"Invalid tags: {self.tags}")
        if self.TSL not in RNATSL.__args__:
            raise ValueError(f"Invalid TSL: {self.TSL}")


CDSSource = Literal['ENSEMBL', 'HAVANA']


@define(hash=True, slots=True, frozen=True, eq=True, order=True, repr=True, str=True)
class AttrCDS:
    source: CDSSource
    transcripts: frozenset[str] = field(converter=lambda x: frozenset(x))

    def __attrs_post_init__(self):
        if self.source not in CDSSource.__args__:
            raise ValueError(f"Invalid source: {self.source}")


def load() -> at.Annotome[AttrGene, AttrRNA, AttrCDS]:
    return at.read_pkl(index.as_posix())
