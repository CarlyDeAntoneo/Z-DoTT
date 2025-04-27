from pathlib import Path

import pandas as pd
from biobit.toolkit import seqproj

from utils.seqproj import JCCSeq

ROOT = Path(__file__).parent
FASTQ = ROOT / "fastq"


def sample_builder(ind: str, _: pd.DataFrame) -> seqproj.Sample:
    match ind:
        case 'WT_Input_Mock_A+WT_Z22_Mock_A+KO_Input_Mock_A+KO_Z22_Mock_A':
            condition = 'Mock'
            replica = 'A'
        case 'WT_Input_Mock_B+WT_Z22_Mock_B+KO_Input_Mock_B+KO_Z22_Mock_B':
            condition = 'Mock'
            replica = 'B'
        case 'WT_Input_Mock_C+WT_Z22_Mock_C+KO_Input_Mock_C+KO_Z22_Mock_C':
            condition = 'Mock'
            replica = 'C'
        case 'WT_Input_IFN_A+WT_Z22_IFN_A+KO_Input_IFN_A+KO_Z22_IFN_A':
            condition = 'IFN'
            replica = 'A'
        case 'WT_Input_IFN_B+WT_Z22_IFN_B+KO_Input_IFN_B+KO_Z22_IFN_B':
            condition = 'IFN'
            replica = 'B'
        case 'WT_Input_IFN_C+WT_Z22_IFN_C+KO_Input_IFN_C+KO_Z22_IFN_C':
            condition = 'IFN'
            replica = 'C'
        case _:
            raise ValueError(f"Unknown sample: {ind}")


    organism = {"Mus musculus"}

    attributes = {"condition": condition, "replica": replica}
    return seqproj.Sample(ind=ind, organism=organism, attributes=attributes)


def library_builder(_: str, __: str, data: pd.DataFrame) -> seqproj.Library:
    tags = set(data['Tags'])
    assert len(tags) == 1, f"Expected 1 unique Tags set, got: {tags}"

    tag_tuple = tags.pop()
    well = "_".join(tag_tuple)

    if well in {'WT_Input_Mock_A', 'WT_Input_Mock_B', 'WT_Input_Mock_C', 'WT_Input_IFN_A', 'WT_Input_IFN_B', 'WT_Input_IFN_C'}:
        selection = {"Total RNA", "rRNA depletion"}
        RIP = 'Input'
        Cell = 'WT'
    elif well in {'KO_Input_Mock_A', 'KO_Input_Mock_B', 'KO_Input_Mock_C', 'KO_Input_IFN_A', 'KO_Input_IFN_B', 'KO_Input_IFN_C'}:
        selection = {"Total RNA", "rRNA depletion"}
        RIP = 'Input'
        Cell = 'KO'
    elif well in {'WT_Z22_Mock_A', 'WT_Z22_Mock_B', 'WT_Z22_Mock_C', 'WT_Z22_IFN_A', 'WT_Z22_IFN_B', 'WT_Z22_IFN_C'}:
        selection = {"Total RNA", "rRNA depletion", "Z22 RIP"}
        RIP = 'Z22'
        Cell = 'WT'
    elif well in {'KO_Z22_Mock_A', 'KO_Z22_Mock_B', 'KO_Z22_Mock_C', 'KO_Z22_IFN_A', 'KO_Z22_IFN_B', 'KO_Z22_IFN_C'}:
        selection = {"Total RNA", "rRNA depletion", "Z22 RIP"}
        RIP = 'Z22'
        Cell = 'KO'
    else:
        raise ValueError(f"Unknown well: {well}")
        
    return seqproj.Library({"RNA"}, selection, seqproj.Strandedness.Reverse, attributes={"RIP": RIP, "Cell": Cell})


def experiment_builder(
        ind: str, sample: seqproj.Sample, library: seqproj.Library, runs: tuple[seqproj.Run, ...], _: pd.DataFrame
) -> seqproj.Experiment:
    selection = library.attributes['RIP']
    title = f"{library.attributes['Cell']}_{selection}_{sample.attributes['condition']}_{sample.attributes['replica']}"
    return seqproj.Experiment(ind=ind, sample=sample, library=library, runs=runs, attributes={"title": title})


def project_builder(
        inds: tuple[str, ...], experiments: tuple[seqproj.Experiment, ...], samples: tuple[seqproj.Sample, ...]
):
    assert len(inds) == 1
    return seqproj.Project(f"JCC400-0822", experiments, samples)


data = JCCSeq.initialize(FASTQ.glob("**/*.fastq.gz"), ROOT)

# Remap input & RIP samples
remapping = {
'WT_Input_Mock_A':'WT_Input_Mock_A+WT_Z22_Mock_A+KO_Input_Mock_A+KO_Z22_Mock_A', 'WT_Z22_Mock_A':'WT_Input_Mock_A+WT_Z22_Mock_A+KO_Input_Mock_A+KO_Z22_Mock_A', 
'KO_Input_Mock_A':'WT_Input_Mock_A+WT_Z22_Mock_A+KO_Input_Mock_A+KO_Z22_Mock_A', 'KO_Z22_Mock_A':'WT_Input_Mock_A+WT_Z22_Mock_A+KO_Input_Mock_A+KO_Z22_Mock_A', 
'WT_Input_Mock_B':'WT_Input_Mock_B+WT_Z22_Mock_B+KO_Input_Mock_B+KO_Z22_Mock_B', 'WT_Z22_Mock_B':'WT_Input_Mock_B+WT_Z22_Mock_B+KO_Input_Mock_B+KO_Z22_Mock_B', 
'KO_Input_Mock_B':'WT_Input_Mock_B+WT_Z22_Mock_B+KO_Input_Mock_B+KO_Z22_Mock_B', 'KO_Z22_Mock_B':'WT_Input_Mock_B+WT_Z22_Mock_B+KO_Input_Mock_B+KO_Z22_Mock_B',
'WT_Input_Mock_C':'WT_Input_Mock_C+WT_Z22_Mock_C+KO_Input_Mock_C+KO_Z22_Mock_C', 'WT_Z22_Mock_C':'WT_Input_Mock_C+WT_Z22_Mock_C+KO_Input_Mock_C+KO_Z22_Mock_C', 
'KO_Input_Mock_C':'WT_Input_Mock_C+WT_Z22_Mock_C+KO_Input_Mock_C+KO_Z22_Mock_C', 'KO_Z22_Mock_C':'WT_Input_Mock_C+WT_Z22_Mock_C+KO_Input_Mock_C+KO_Z22_Mock_C', 
'WT_Input_IFN_A':'WT_Input_IFN_A+WT_Z22_IFN_A+KO_Input_IFN_A+KO_Z22_IFN_A', 'WT_Z22_IFN_A':'WT_Input_IFN_A+WT_Z22_IFN_A+KO_Input_IFN_A+KO_Z22_IFN_A', 
'KO_Input_IFN_A':'WT_Input_IFN_A+WT_Z22_IFN_A+KO_Input_IFN_A+KO_Z22_IFN_A', 'KO_Z22_IFN_A':'WT_Input_IFN_A+WT_Z22_IFN_A+KO_Input_IFN_A+KO_Z22_IFN_A',
'WT_Input_IFN_B':'WT_Input_IFN_B+WT_Z22_IFN_B+KO_Input_IFN_B+KO_Z22_IFN_B', 'WT_Z22_IFN_B':'WT_Input_IFN_B+WT_Z22_IFN_B+KO_Input_IFN_B+KO_Z22_IFN_B', 
'KO_Input_IFN_B':'WT_Input_IFN_B+WT_Z22_IFN_B+KO_Input_IFN_B+KO_Z22_IFN_B', 'KO_Z22_IFN_B':'WT_Input_IFN_B+WT_Z22_IFN_B+KO_Input_IFN_B+KO_Z22_IFN_B', 
'WT_Input_IFN_C':'WT_Input_IFN_C+WT_Z22_IFN_C+KO_Input_IFN_C+KO_Z22_IFN_C', 'WT_Z22_IFN_C':'WT_Input_IFN_C+WT_Z22_IFN_C+KO_Input_IFN_C+KO_Z22_IFN_C', 
'KO_Input_IFN_C':'WT_Input_IFN_C+WT_Z22_IFN_C+KO_Input_IFN_C+KO_Z22_IFN_C', 'KO_Z22_IFN_C':'WT_Input_IFN_C+WT_Z22_IFN_C+KO_Input_IFN_C+KO_Z22_IFN_C'
}


data["Lane"] = [f"{s}-{l}" for s, l in data[["Sample", "Lane"]].itertuples(index=False)]
data["Sample"] = data["Tags"].apply(lambda x: "_".join(x)).map(remapping)
sample_names = data["Tags"].apply(lambda x: "_".join(x))
missing = sample_names[~sample_names.isin(remapping)]
print("Missing tag strings (not in remapping):")
print(missing.unique())

assert data.isnull().sum().sum() == 0, data.isnull().sum()

project = JCCSeq.parse(
    data,
    seq_machine="Illumina NovaSeq 6000",
    sample_builder=sample_builder,
    library_builder=library_builder,
    experiment_builder=experiment_builder,
    project_builder=project_builder
)
seqproj.adapter.yaml.dump(project, ROOT / "seq-project.yaml")
