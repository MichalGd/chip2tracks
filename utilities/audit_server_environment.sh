#!/usr/bin/env bash
set -u

output="chip2tracks_server_audit_$(date -u +%Y%m%dT%H%M%SZ).txt"
host_index=""
composite_index=""
reference_root=""

usage() {
    echo "Usage: audit_server_environment.sh [--output FILE] [--host-index PREFIX] [--composite-index PREFIX] [--reference-root DIR]"
}
while (( $# )); do
    case "$1" in
        --output) output="${2:?missing output path}"; shift 2 ;;
        --host-index) host_index="${2:?missing Bowtie2 prefix}"; shift 2 ;;
        --composite-index) composite_index="${2:?missing composite Bowtie2 prefix}"; shift 2 ;;
        --reference-root) reference_root="${2:?missing reference root}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

tools=(bash python3 Rscript fastqc multiqc trim_galore cutadapt bowtie2 samtools bedtools picard \
       bamCoverage multiBamSummary plotCorrelation plotPCA plotFingerprint computeMatrix plotProfile \
       macs3 epic2 idr run_spp.R preseq bedGraphToBigWig bedClip)

{
    printf 'chip2tracks server audit\ncreated_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '\n[host]\n'
    uname -a 2>&1 || true
    id 2>&1 || true
    printf 'shell\t%s\npath\t%s\n' "${SHELL:-unknown}" "$PATH"
    printf '\n[resources]\n'
    command -v nproc >/dev/null 2>&1 && nproc || true
    command -v free >/dev/null 2>&1 && free -h || true
    df -h . 2>&1 || true
    ulimit -a 2>&1 || true

    printf '\n[commands]\ncommand\tstatus\tpath\tversion\n'
    for tool in "${tools[@]}"; do
        path="$(command -v "$tool" 2>/dev/null || true)"
        if [[ -z "$path" ]]; then
            printf '%s\tMISSING\t.\t.\n' "$tool"
            continue
        fi
        version="$($tool --version 2>&1 | head -n 1 || $tool -v 2>&1 | head -n 1 || true)"
        version="${version//$'\t'/ }"
        printf '%s\tFOUND\t%s\t%s\n' "$tool" "$path" "$version"
    done

    printf '\n[conda]\n'
    if command -v conda >/dev/null 2>&1; then
        conda info 2>&1 || true
        conda env list 2>&1 || true
        conda list 2>&1 || true
    else
        echo 'conda: MISSING'
    fi

    printf '\n[R_packages]\n'
    if command -v Rscript >/dev/null 2>&1; then
        Rscript -e 'p<-c("DESeq2","DiffBind","GenomicAlignments","GenomicRanges","Rsamtools","rtracklayer","BiocParallel","spp"); for(x in p) cat(x,"\t",requireNamespace(x,quietly=TRUE),"\t",if(requireNamespace(x,quietly=TRUE)) as.character(packageVersion(x)) else ".","\n",sep="")' 2>&1 || true
    fi

    printf '\n[python_packages]\n'
    python3 -m pip list --format=freeze 2>&1 | grep -Ei '^(multiqc|cutadapt|macs3|epic2|idr|pybigwig|pandas|matplotlib)==' || true

    printf '\n[references]\n'
    if [[ -n "$reference_root" ]]; then
        printf 'reference_root\t%s\n' "$reference_root"
        find "$reference_root" -maxdepth 4 -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.chrom.sizes' -o -name '*blacklist*.bed' -o -name '*.gtf' \) -printf '%s\t%p\n' 2>&1 | sort || true
    else
        echo 'reference_root: NOT_PROVIDED'
    fi
    for label_prefix in "host:$host_index" "host_dm6_composite:$composite_index"; do
        label="${label_prefix%%:*}"; prefix="${label_prefix#*:}"
        [[ -n "$prefix" ]] || { printf '%s\tNOT_PROVIDED\n' "$label"; continue; }
        if [[ -f "${prefix}.1.bt2" || -f "${prefix}.1.bt2l" ]]; then
            printf '%s\tFOUND\t%s\n' "$label" "$prefix"
            bowtie2-inspect -n "$prefix" 2>&1 | head -n 20 || true
        else
            printf '%s\tMISSING\t%s\n' "$label" "$prefix"
        fi
    done
} > "$output"

printf 'Audit written to %s\n' "$output"
printf 'No software or references were installed or modified.\n'
