"""Content fingerprints and exact installed versions; no third-party dependencies.

Called by write_run_info.sh. Paths that cannot be resolved are reported explicitly.
Environment specifications are installation recipes; this records what actually ran.
"""
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def fingerprint(path):
    path = Path(path)
    result = {"path": str(path.absolute())}
    try:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        result.update(sha256=digest.hexdigest(), bytes=path.stat().st_size)
    except OSError as error:
        result["error"] = str(error)
    return result


def reference_paths(config):
    paths = set()
    for key in ("sample_table", "genome_fasta", "chrom_sizes"):
        if config.get(key):
            paths.add(config[key])
    homer = config["program"]["homer"]
    if homer.get("tss", {}).get("gtf"):
        paths.add(homer["tss"]["gtf"])
    for key in ("mirnas", "trnas", "snrnas"):
        if config.get("qc", {}).get(key):
            paths.add(config["qc"][key])
    index = Path(config["program"]["genome_index"])
    if index.is_dir():
        paths.update(str(p) for p in index.rglob("*") if p.is_file())
    else:
        paths.update(str(p) for p in index.parent.glob(index.name + ".*") if p.is_file())
    if Path(homer["genome"]).is_file():
        paths.add(homer["genome"])
    else:
        executable = shutil.which(homer.get("tss", {}).get("program", "findcsRNATSS.pl"))
        resolved = False
        if executable:
            root = Path(executable).resolve().parent.parent
            conf = root / "config.txt"
            if conf.is_file():
                paths.add(str(conf))
                section = ""
                for line in conf.read_text().splitlines():
                    fields = line.split("\t")
                    if len(fields) == 1 and fields[0].isupper():
                        section = fields[0]
                    if section == "GENOMES" and fields[0] == homer["genome"] and len(fields) >= 5:
                        genome_dir = root / fields[4]
                        files = [p for p in genome_dir.rglob("*") if p.is_file()]
                        paths.update(str(p) for p in files)
                        resolved = bool(files)
        if not resolved:
            paths.add("UNRESOLVED_HOMER_GENOME/" + homer["genome"])
    return sorted(paths)


def collect(config, repo):
    tools = ("bwa", "STAR", "bowtie2", "hisat2", "samtools", "bedtools", "pigz",
             "bfqutils", "bam2td", config["program"]["homer"]["tss"].get("program", "findcsRNATSS.pl"),
             "findPeaks", "annotatePeaks.pl", "makeTagDirectory", "makeUCSCfile")
    executables = {name: fingerprint(shutil.which(name)) if shutil.which(name)
                   else {"error": "not on PATH"} for name in tools}
    python_packages = sorted({(d.metadata["Name"], d.version)
                              for d in importlib.metadata.distributions() if d.metadata["Name"]})
    # Conda's installed records include build strings and source URLs. Avoid exporting
    # channel credentials: retain name/version/build/subdir only.
    conda_packages = []
    for record in sorted((Path(sys.prefix) / "conda-meta").glob("*.json")):
        data = json.loads(record.read_text())
        conda_packages.append({k: data.get(k) for k in ("name", "version", "build", "subdir")})
    try:
        result = subprocess.run(["Rscript", "--vanilla", "-e", 'p <- installed.packages(); write.table(p[,c("Package","Version","Built")], row.names=FALSE, sep="\\t", quote=FALSE); sessionInfo()'],
                                capture_output=True, text=True, timeout=60)
        r_packages = result.stdout if result.returncode == 0 else result.stderr
    except (OSError, subprocess.TimeoutExpired) as error:
        r_packages = str(error)
    sources = [Path(repo) / "workflow" / "Snakefile"]
    for folder in ("scripts", "schemas"):
        sources.extend(p for p in (Path(repo) / "workflow" / folder).iterdir() if p.is_file())
    sources.extend(Path(repo) / "workflow/envs" / name for name in ("minimal_env.yaml", "full_env.yaml"))
    homer_installation = {}
    homer_executable = shutil.which(config["program"]["homer"]["tss"].get("program", "findcsRNATSS.pl"))
    if homer_executable:
        homer_config = Path(homer_executable).resolve().parent.parent / "config.txt"
        homer_installation["config"] = fingerprint(homer_config)
        if homer_config.is_file():
            for line in homer_config.read_text().splitlines():
                fields = line.split("\t")
                if len(fields) >= 2 and fields[0] == "homer":
                    homer_installation["version"] = fields[1]
                    break
    return {"homer_installation": homer_installation, "references": [fingerprint(p) for p in reference_paths(config)],
            "workflow_sources": [fingerprint(p) for p in sorted(sources)],
            "executables": executables, "python": sys.version,
            "python_packages": python_packages, "conda_packages": conda_packages,
            "R_packages_and_session": r_packages,
            "container": os.environ.get("APPTAINER_CONTAINER") or os.environ.get("SINGULARITY_CONTAINER")}


if __name__ == "__main__":
    with open(sys.argv[1]) as handle:
        config = json.load(handle)
    print(json.dumps(collect(config, sys.argv[2]), indent=2))
