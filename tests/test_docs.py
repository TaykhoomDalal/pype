from html.parser import HTMLParser
from pathlib import Path


ROOT = Path(__file__).parents[1]


class SiteParser(HTMLParser):
	def __init__(self):
		super().__init__()
		self.ids = set()
		self.links = []

	def handle_starttag(self, tag, attrs):
		attrs = dict(attrs)
		if "id" in attrs:
			assert attrs["id"] not in self.ids
			self.ids.add(attrs["id"])
		for name in ("href", "src"):
			if name in attrs:
				self.links.append(attrs[name])


def validate_page(path):
	parser = SiteParser()
	parser.feed(path.read_text())

	for link in parser.links:
		if link.startswith("#"):
			assert link[1:] in parser.ids
		elif link.startswith(("http://", "https://", "mailto:")):
			continue
		elif link.startswith("assets/"):
			assert (ROOT / link).is_file()
		else:
			target = path.parent / link
			if target.is_dir():
				target /= "index.html"
			assert target.is_file()


def test_documentation_links_and_anchors_exist():
	for page in (ROOT / "docs").glob("*.html"):
		validate_page(page)


if __name__ == "__main__":
	test_documentation_links_and_anchors_exist()
