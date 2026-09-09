"""Check actual HOMER region coordinates and the run's content fingerprints."""
import csv
import hashlib
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])


def regions(name):
    with (root / name).open() as handle:
        rows = list(csv.reader(handle, delimiter="\t"))[1:]
    return {row[0]: row[1:5] for row in rows}


assert regions("tss.consensus.homer.raw.txt") == regions("tss.consensus.in.homer.raw.txt")
manifest = (root / "run_info.txt").read_text()
text = manifest.split("content fingerprints and installed packages\n", 1)[1]
details, _ = json.JSONDecoder().raw_decode(text)
references = details["references"]
assert any(item["path"].endswith("genome.fa.bwt") for item in references)
for item in references:
    assert "error" not in item, item
    assert hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest() == item["sha256"]
assert details["python_packages"] and details["R_packages_and_session"]
assert all("sha256" in item for item in details["workflow_sources"])
print("Consensus/input coordinates and run reference fingerprints verified")
