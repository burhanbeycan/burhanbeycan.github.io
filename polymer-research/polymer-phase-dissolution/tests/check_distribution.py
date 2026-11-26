"""Validate internal assets, replay chapters and scored PS cloud-point records."""
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
import csv

here = Path(__file__).resolve().parent.parent
dist = here / "dist" if (here / "dist").is_dir() else here

class References(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if "id" in attrs:
            assert attrs["id"] not in self.ids, (page.name, attrs["id"])
            self.ids.add(attrs["id"])
        for key in ("href", "src", "poster"):
            raw = attrs.get(key, "")
            link = urlsplit(raw)
            if raw and not link.scheme and not link.netloc and link.path:
                self.links.append(unquote(link.path))

for page in sorted(dist.rglob("*.html")):
    links = References()
    links.feed(page.read_text())
    for target in links.links:
        assert (page.parent / target).is_file(), (page.name, target)

with (dist / "ps-cloud-points.csv").open() as handle:
    rows = list(csv.DictReader(handle))
assert len(rows) == 97
assert Counter(r["split"] for r in rows) == {"calibration": 73, "holdout": 24}
for coefficient, expected in (("T_cloud_FH0_K", 2.45), ("T_cloud_FH1_K", 0.774)):
    residuals = [
        float(row[coefficient]) - float(row["T_cloud_observed_K"])
        for row in rows if row["split"] == "holdout"
    ]
    rmse = (sum(r * r for r in residuals) / len(residuals)) ** 0.5
    assert abs(rmse - expected) < 0.002, (coefficient, rmse)
print("Static HTML references, unique IDs, and 97 experimental rows verified.")
