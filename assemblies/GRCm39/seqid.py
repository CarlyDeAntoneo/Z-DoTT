def from_refseq(refseq: str) -> str:
    return {
        'NC_000067.7': 'chr1',
        'NC_000068.8': 'chr2',
        'NC_000069.7': 'chr3',
        'NC_000070.7': 'chr4',
        'NC_000071.7': 'chr5',
        'NC_000072.7': 'chr6',
        'NC_000073.7': 'chr7',
        'NC_000074.7': 'chr8',
        'NC_000075.7': 'chr9',
        'NC_000076.7': 'chr10',
        'NC_000077.7': 'chr11',
        'NC_000078.7': 'chr12',
        'NC_000079.7': 'chr13',
        'NC_000080.7': 'chr14',
        'NC_000081.7': 'chr15',
        'NC_000082.7': 'chr16',
        'NC_000083.7': 'chr17',
        'NC_000084.7': 'chr18',
        'NC_000085.7': 'chr19',
        'NC_000086.8': 'chrX',
        'NC_000087.8': 'chrY',
        'NT_166280.1': 'GL456210.1',
        'NT_166281.1': 'GL456211.1',
        'NT_166282.1': 'GL456212.1',
        'NT_162750.1': 'GL456221.1',
        'NT_166338.1': 'GL456239.1',
        'NW_023337852.1': 'MU069434.1',
        'NT_187055.1': 'JH584295.1',
        'NT_166438.1': 'GL456354.1',
        'NT_187056.1': 'JH584296.1',
        'NT_187057.1': 'JH584297.1',
        'NT_187058.1': 'JH584298.1',
        'NT_187059.1': 'JH584299.1',
        'NT_166307.1': 'GL456219.1',
        'NT_165789.3': 'GL456233.2',
        'NT_187060.1': 'JH584300.1',
        'NT_187061.1': 'JH584301.1',
        'NT_187062.1': 'JH584302.1',
        'NT_187063.1': 'JH584303.1',
        'NT_166443.1': 'GL456359.1',
        'NT_166444.1': 'GL456360.1',
        'NT_166450.1': 'GL456366.1',
        'NT_166451.1': 'GL456367.1',
        'NT_166452.1': 'GL456368.1',
        'NT_166454.1': 'GL456370.1',
        'NT_166456.1': 'GL456372.1',
        'NT_166462.1': 'GL456378.1',
        'NT_166463.1': 'GL456379.1',
        'NT_166465.1': 'GL456381.1',
        'NT_166466.1': 'GL456382.1',
        'NT_166467.1': 'GL456383.1',
        'NT_166469.1': 'GL456385.1',
        'NT_166471.1': 'GL456387.1',
        'NT_166473.1': 'GL456389.1',
        'NT_166474.1': 'GL456390.1',
        'NT_166476.1': 'GL456392.1',
        'NT_166478.1': 'GL456394.1',
        'NT_166480.1': 'GL456396.1',
        'NT_187064.1': 'JH584304.1',
        'NW_023337853.1': 'MU069435.1',
        'NC_005089.1': 'chrM'
    }[refseq]


def from_ensembl(ensembl: str) -> str:
    return {
        '1': 'chr1',
        '2': 'chr2',
        '3': 'chr3',
        '4': 'chr4',
        '5': 'chr5',
        '6': 'chr6',
        '7': 'chr7',
        '8': 'chr8',
        '9': 'chr9',
        '10': 'chr10',
        '11': 'chr11',
        '12': 'chr12',
        '13': 'chr13',
        '14': 'chr14',
        '15': 'chr15',
        '16': 'chr16',
        '17': 'chr17',
        '18': 'chr18',
        '19': 'chr19',
        'X': 'chrX',
        'Y': 'chrY',
        'MT': 'chrM',
        'GL456210.1': 'GL456210.1',
        'GL456211.1': 'GL456211.1',
        'GL456212.1': 'GL456212.1',
        'GL456219.1': 'GL456219.1',
        'GL456221.1': 'GL456221.1',
        'GL456233.2': 'GL456233.2',
        'GL456239.1': 'GL456239.1',
        'GL456354.1': 'GL456354.1',
        'GL456359.1': 'GL456359.1',
        'GL456360.1': 'GL456360.1',
        'GL456366.1': 'GL456366.1',
        'GL456367.1': 'GL456367.1',
        'GL456368.1': 'GL456368.1',
        'GL456370.1': 'GL456370.1',
        'GL456372.1': 'GL456372.1',
        'GL456378.1': 'GL456378.1',
        'GL456379.1': 'GL456379.1',
        'GL456381.1': 'GL456381.1',
        'GL456382.1': 'GL456382.1',
        'GL456383.1': 'GL456383.1',
        'GL456385.1': 'GL456385.1',
        'GL456387.1': 'GL456387.1',
        'GL456389.1': 'GL456389.1',
        'GL456390.1': 'GL456390.1',
        'GL456392.1': 'GL456392.1',
        'GL456394.1': 'GL456394.1',
        'GL456396.1': 'GL456396.1',
        'JH584295.1': 'JH584295.1',
        'JH584296.1': 'JH584296.1',
        'JH584297.1': 'JH584297.1',
        'JH584298.1': 'JH584298.1',
        'JH584299.1': 'JH584299.1',
        'JH584300.1': 'JH584300.1',
        'JH584301.1': 'JH584301.1',
        'JH584302.1': 'JH584302.1',
        'JH584303.1': 'JH584303.1',
        'JH584304.1': 'JH584304.1',
        'MU069434.1': 'MU069434.1',
        'MU069435.1': 'MU069435.1',
    }[ensembl]


def all() -> tuple[str, ...]:
    return (
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
        'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY', 'chrM', 'GL456210.1', 'GL456211.1',
        'GL456212.1', 'GL456219.1', 'GL456221.1', 'GL456233.2', 'GL456239.1', 'GL456354.1', 'GL456359.1', 'GL456360.1',
        'GL456366.1', 'GL456367.1', 'GL456368.1', 'GL456370.1', 'GL456372.1', 'GL456378.1', 'GL456379.1', 'GL456381.1',
        'GL456382.1', 'GL456383.1', 'GL456385.1', 'GL456387.1', 'GL456389.1', 'GL456390.1', 'GL456392.1', 'GL456394.1',
        'GL456396.1', 'JH584295.1', 'JH584296.1', 'JH584297.1', 'JH584298.1', 'JH584299.1', 'JH584300.1', 'JH584301.1',
        'JH584302.1', 'JH584303.1', 'JH584304.1', 'MU069434.1', 'MU069435.1'
    )


def sizes() -> dict[str, int]:
    return {
        'chr1': 195154279,
        'chr2': 181755017,
        'chr3': 159745316,
        'chr4': 156860686,
        'chr5': 151758149,
        'chr6': 149588044,
        'chr7': 144995196,
        'chr8': 130127694,
        'chr9': 124359700,
        'chr10': 130530862,
        'chr11': 121973369,
        'chr12': 120092757,
        'chr13': 120883175,
        'chr14': 125139656,
        'chr15': 104073951,
        'chr16': 98008968,
        'chr17': 95294699,
        'chr18': 90720763,
        'chr19': 61420004,
        'chrX': 169476592,
        'chrY': 91455967,
        'chrM': 16299,
        'GL456210.1': 169725,
        'GL456211.1': 241735,
        'GL456212.1': 153618,
        'GL456219.1': 175968,
        'GL456221.1': 206961,
        'GL456233.2': 559103,
        'GL456239.1': 40056,
        'GL456354.1': 195993,
        'GL456359.1': 22974,
        'GL456360.1': 31704,
        'GL456366.1': 47073,
        'GL456367.1': 42057,
        'GL456368.1': 20208,
        'GL456370.1': 26764,
        'GL456372.1': 28664,
        'GL456378.1': 31602,
        'GL456379.1': 72385,
        'GL456381.1': 25871,
        'GL456382.1': 23158,
        'GL456383.1': 38659,
        'GL456385.1': 35240,
        'GL456387.1': 24685,
        'GL456389.1': 28772,
        'GL456390.1': 24668,
        'GL456392.1': 23629,
        'GL456394.1': 24323,
        'GL456396.1': 21240,
        'JH584295.1': 1976,
        'JH584296.1': 199368,
        'JH584297.1': 205776,
        'JH584298.1': 184189,
        'JH584299.1': 953012,
        'JH584300.1': 182347,
        'JH584301.1': 259875,
        'JH584302.1': 155838,
        'JH584303.1': 158099,
        'JH584304.1': 114452,
        'MU069434.1': 8412,
        'MU069435.1': 31129,
    }
