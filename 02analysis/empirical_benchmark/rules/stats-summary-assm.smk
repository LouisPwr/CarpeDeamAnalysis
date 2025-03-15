from itertools import chain
import pathlib
import os


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def get_stats_assm(wildcards, file_paths):
    files = expand(
        file_paths,
        smp=sample_label_read,
    )
    return files


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def remove_suffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[: -len(suffix)]
    return input_string


def read_assm_tables(file):
    if is_non_zero_file(file):
        df = pd.read_csv(file, sep="\t", header=None)
        df.columns = ["label", "contig", "length"]
        return df


rule stats_summary_assm:
    input:
        stats_assm_refined=lambda wc: get_stats_assm(
            wc, file_paths=config["rdir"] + "/stats/{smp}.assembly-refined.txt"
        ),
        stats_assm_combined=lambda wc: get_stats_assm(
            wc, file_paths=config["rdir"] + "/stats/{smp}.assembly-combined.txt"
        ),
        stats_assm_extended=lambda wc: get_stats_assm(
            wc, file_paths=config["rdir"] + "/stats/{smp}.assembly-extended.txt"
        ),
    output:
        stats_assm_refined_summary=(
            config["rdir"] + "/stats/all.stats-assm-refined-summary.tsv.gz"
        ),
        stats_assm_combined_summary=(
            config["rdir"] + "/stats/all.stats-assm-combined-summary.tsv.gz"
        ),
        stats_assm_extended_summary=(
            config["rdir"] + "/stats/all.stats-assm-extended-summary.tsv.gz"
        ),
    threads: config["threads_16"]
    log:
        config["rdir"] + "/logs/stats/stats-summary-assm.log",
    benchmark:
        config["rdir"] + "/benchmarks/stats/stats-summary-assm.bmk"
    message:
        """--- Summarize stats"""
    run:
        df = pd.concat(
            map(
                lambda file: read_assm_tables(file),
                fast_flatten(list({input.stats_assm_refined})),
            )
        )
        df["label"] = df["label"].astype(str).map(sample_label_dict_assm)
        df.to_csv(
            output.stats_assm_refined_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )

        df = pd.concat(
            map(
                lambda file: read_assm_tables(file),
                fast_flatten(list({input.stats_assm_combined})),
            )
        )
        df["label"] = df["label"].astype(str).map(sample_label_dict_assm)
        df.to_csv(
            output.stats_assm_combined_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )

        df = pd.concat(
            map(
                lambda file: read_assm_tables(file),
                fast_flatten(list({input.stats_assm_extended})),
            )
        )
        df["label"] = df["label"].astype(str).map(sample_label_dict_assm)
        df.to_csv(
            output.stats_assm_extended_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )
