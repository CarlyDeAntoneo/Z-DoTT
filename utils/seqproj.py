from collections import defaultdict
from os import PathLike
from pathlib import Path
from typing import Iterable, Callable

import pandas as pd
from biobit.toolkit import nfcore, seqproj


def initialize_rnaseq(initpy: PathLike[str], fname: str = "seq-project.yaml") -> seqproj.Project:
    root = Path(initpy).parent
    project = seqproj.adapter.yaml.load(root / fname)
    return nfcore.rnaseq.parse.into_seqproj(
        project, root / "results",
        seqexp2descriptor=lambda exp: nfcore.rnaseq.descriptor.from_seqexp(exp, title_builder="title")
    )


class JCCSeq:
    @staticmethod
    def initialize(fastq: Iterable[PathLike[str]], reroot: PathLike[str] | None = None) -> pd.DataFrame:
        path = [Path(fq) for fq in fastq]
        if reroot:
            path = [x.relative_to(reroot) for x in path]
        fname = [x.name for x in path]
        df = pd.DataFrame({"fastq": path, "fname": fname})

        exp, prj, tags, smpl, lane, rmate, postfix = [], [], [], [], [], [], []
        for fname in df["fname"]:
            parts = str(fname).split("_")
            postfix.append(parts.pop())
            rmate.append(parts.pop())
            lane.append(parts.pop())
            smpl.append(parts.pop())
            exp.append(parts.pop(0))
            prj.append(parts.pop(0))
            
            # Exclude 'Feb2025' from tags
            parts = [p for p in parts if p != "Feb2025"]
            tags.append(tuple(parts))

        df["Experiment"], df["Project"], df["Tags"], df["Sample"], df["Lane"], df["Read Mate"], df["Postfix"] = \
            exp, prj, tags, smpl, lane, rmate, postfix

        return df

    @staticmethod
    def parse(
            data: pd.DataFrame,
            seq_machine: str,
            sample_builder: Callable[[str, pd.DataFrame], seqproj.Sample],
            library_builder: Callable[[str, str, pd.DataFrame], seqproj.Library],
            experiment_builder: Callable[[str, seqproj.Sample, seqproj.Library,
                                          tuple[seqproj.Run, ...], pd.DataFrame], seqproj.Experiment],
            project_builder: Callable[[tuple[str, ...], tuple[seqproj.Experiment, ...],
                                       tuple[seqproj.Sample, ...]], seqproj.Project]
    ) -> seqproj.Project:
        # 1. Build samples
        samples = {}
        for smplid, subdf in data.groupby("Sample"):
            samples[smplid] = sample_builder(smplid, subdf)

        # 2. Build seqruns
        seqruns = defaultdict(list)
        for (expid, sample, seqlane), subdf in data.groupby(["Experiment", "Sample", "Lane"]):
            if len(subdf) != 2:
                raise ValueError(f"Expected 2 files for {expid} {sample} {seqlane}: {subdf}")
            if set(subdf["Read Mate"]) != {"R1", "R2"}:
                raise ValueError(f"Expected R1 & R2 for {expid} {sample} {seqlane}: {subdf}")

            r1 = subdf[subdf["Read Mate"] == "R1"]
            r2 = subdf[subdf["Read Mate"] == "R2"]

            layout = seqproj.Layout.Paired(seqproj.MatesOrientation.Inward, (r1["fastq"].iloc[0], r2["fastq"].iloc[0]))
            seqruns[(expid, sample)].append(seqproj.Run(
                f"{expid}-{sample}-{seqlane}", layout, seq_machine,
            ))

        # 3. Build experiments
        experiments = {}
        for expid, subdf in data.groupby("Experiment"):
            if subdf["Sample"].nunique() != 1:
                raise ValueError(f"Multiple samples for an experiment are not allowed: {subdf['Sample'].unique()}")
            smplid = subdf["Sample"].iloc[0]

            library = library_builder(expid, smplid, subdf)

            assert expid not in experiments
            experiments[expid] = experiment_builder(
                expid, samples[smplid], library, tuple(seqruns[(expid, smplid)]), subdf
            )

        prjind = tuple(data["Project"].unique())
        return project_builder(prjind, tuple(experiments.values()), tuple(samples.values()))