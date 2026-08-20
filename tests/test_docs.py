from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlsplit
from xml.etree import ElementTree


ROOT = Path(__file__).parents[1]


class SiteParser(HTMLParser):
	def __init__(self):
		super().__init__()
		self.h1_count = 0
		self.ids = set()
		self.links = []
		self.theme_toggle_count = 0

	def handle_starttag(self, tag, attrs):
		attrs = dict(attrs)
		if tag == "h1":
			self.h1_count += 1
		if "id" in attrs:
			assert attrs["id"] not in self.ids
			self.ids.add(attrs["id"])
		for name in ("href", "src"):
			if name in attrs:
				self.links.append(attrs[name])
		if "data-theme-toggle" in attrs:
				self.theme_toggle_count += 1


def validate_page(path):
	parser = SiteParser()
	parser.feed(path.read_text())
	assert parser.h1_count == 1
	assert parser.theme_toggle_count == 1

	for link in parser.links:
		if link.startswith("#"):
			assert link[1:] in parser.ids
		elif link.startswith(("http://", "https://", "mailto:")):
			continue
		else:
			parsed = urlsplit(link)
			if parsed.path.startswith("assets/"):
				target = path.parent / parsed.path
				if not target.exists():
					target = ROOT / parsed.path
			else:
				target = path.parent / parsed.path
			if target.is_dir():
				target /= "index.html"
			assert target.is_file()
			if parsed.fragment:
				target_parser = SiteParser()
				target_parser.feed(target.read_text())
				assert parsed.fragment in target_parser.ids


def test_documentation_links_and_anchors_exist():
	for page in (ROOT / "docs").glob("*.html"):
		validate_page(page)


def test_sitemap_lists_public_pages():
	tree = ElementTree.parse(ROOT / "docs/sitemap.xml")
	locations = {
		element.text
		for element in tree.getroot().iter("{http://www.sitemaps.org/schemas/sitemap/0.9}loc")
	}
	for page in (ROOT / "docs").glob("*.html"):
		if page.name == "404.html":
			continue
		url = "https://taykhoomdalal.github.io/pype/"
		if page.name != "index.html":
			url += page.name
		assert url in locations


def test_reproduction_command_is_not_indented():
	page = (ROOT / "docs/reproducibility.html").read_text()
	assert (
		"python -m pip install pype-mr openpyxl\n"
		"python reproducibility/reproduce_paper.py \\\n"
		"  --output paper_reproduction"
	) in page


if __name__ == "__main__":
	test_documentation_links_and_anchors_exist()
