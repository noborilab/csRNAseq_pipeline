"""Regression checks for content hashes and installed HOMER reference resolution."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[2] / "workflow/scripts/provenance.py"
spec = importlib.util.spec_from_file_location("provenance", MODULE)
provenance = importlib.util.module_from_spec(spec)
spec.loader.exec_module(provenance)


class ProvenanceTests(unittest.TestCase):
    def test_content_changes_and_missing_files(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference with spaces.fa"
            path.write_text(">chr1\nACGT\n")
            before = provenance.fingerprint(path)
            path.write_text(">chr1\nACGA\n")
            after = provenance.fingerprint(path)
            self.assertNotEqual(before["sha256"], after["sha256"])
            self.assertEqual(before["bytes"], after["bytes"])
            self.assertIn("error", provenance.fingerprint(path.with_suffix(".missing")))

    def test_homer_alias_and_all_index_components(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "bin").mkdir()
            (root / "data/genomes/test").mkdir(parents=True)
            annotation = root / "data/genomes/test/annotation.gtf"
            annotation.write_text("annotation")
            (root / "config.txt").write_text("GENOMES\ntest\tv1\tdescription\turl\tdata/genomes/test/\n")
            for suffix in (".bwt", ".ann", ".sa"):
                (root / ("index" + suffix)).write_text("index")
            config = {"chrom_sizes": str(root / "chrom.sizes"), "program": {
                "genome_index": str(root / "index"), "homer": {"genome": "test"}}}
            with patch.object(provenance.shutil, "which", return_value=str(root / "bin/findcsRNATSS.pl")):
                paths = provenance.reference_paths(config)
            self.assertIn(str(annotation.resolve()), paths)
            self.assertTrue(all(str(root / ("index" + suffix)) in paths
                                for suffix in (".bwt", ".ann", ".sa")))
            with patch.object(provenance.shutil, "which", return_value=None):
                self.assertIn("UNRESOLVED_HOMER_GENOME/test", provenance.reference_paths(config))


if __name__ == "__main__":
    unittest.main()
