import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
import assemble  # noqa: E402
import bundle  # noqa: E402

SKIP = {"tags", "domains"}


class RewriteLinks(unittest.TestCase):
    def rw(self, s):
        return assemble.rewrite_links(s, "io", SKIP)

    def test_markdown_links_and_images(self):
        self.assertEqual(self.rw("[a](/x/y) ![i](/plots/a.png)\n"),
                         "[a](/io/x/y) ![i](/io/plots/a.png)\n")

    def test_root_link(self):
        self.assertEqual(self.rw("[home](/)"), "[home](/io/)")

    def test_already_prefixed_and_global(self):
        self.assertEqual(self.rw("[a](/io/x) [b](/tags/foo/) [c](/domains)"),
                         "[a](/io/x) [b](/tags/foo/) [c](/domains)")

    def test_relative_absolute_protocol_relative_untouched(self):
        s = "[a](x/y) [b](https://e.com/x) [c](//cdn/x) [d](#frag)"
        self.assertEqual(self.rw(s), s)

    def test_html_attributes(self):
        self.assertEqual(self.rw('<a href="/a/b">x</a> <img src=\'/p.png\'>'),
                         '<a href="/io/a/b">x</a> <img src=\'/io/p.png\'>')

    def test_reference_definition(self):
        self.assertEqual(self.rw("[r]: /docs/x\n"), "[r]: /io/docs/x\n")

    def test_fenced_code_and_code_spans_skipped(self):
        s = "```md\n[a](/x)\n```\n`[b](/y)` and [c](/z)\n~~~\n[d](/w)\n~~~\n"
        self.assertEqual(self.rw(s),
                         "```md\n[a](/x)\n```\n`[b](/y)` and [c](/io/z)\n~~~\n[d](/w)\n~~~\n")

    def test_hugo_ref_shortcodes(self):
        self.assertEqual(self.rw('{{< ref "/components" >}} {{< relref "/tags" >}}'),
                         '{{< ref "/io/components" >}} {{< relref "/tags" >}}')

    def test_idempotent(self):
        once = self.rw("[a](/x)")
        self.assertEqual(self.rw(once), once)


def write(p: Path, text: str = "x"):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


def make_repo(root: Path, id="io", api=False):
    write(root / "docs/site.toml",
          f'id = "{id}"\ntitle = "T"\ndescription = "D"\nrepo = "O/R"\n'
          '[[nav]]\ntitle = "Guide"\npath = "guide/"\n')
    write(root / "docs/content/_index.md", "---\ntitle: T\n---\n[g](/guide/)\n")
    write(root / "docs/content/guide/a.md", "![i](/plots/p.png)")
    write(root / "docs/static/plots/p.png", "png")
    if api:
        write(root / "docs/content/generated/api/_index.md", "api")


class Bundle(unittest.TestCase):
    def test_bundle_copies_and_adds_api_nav(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            make_repo(d / "r", api=True)
            meta = bundle.build_bundle(d / "r", d / "out")
            self.assertEqual([n["title"] for n in meta["nav"]], ["Guide", "API"])
            self.assertTrue((d / "out/content/generated/api/_index.md").is_file())
            self.assertTrue((d / "out/static/plots/p.png").is_file())
            self.assertEqual(json.loads((d / "out/bundle.json").read_text())["id"], "io")

    def test_missing_site_toml(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(bundle.BundleError):
                bundle.build_bundle(Path(d), Path(d) / "out")

    def test_bad_id_and_missing_key(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            write(d / "docs/site.toml", 'id = "Bad Id"\ntitle="T"\ndescription="D"\nrepo="O/R"\n')
            with self.assertRaises(bundle.BundleError):
                bundle.load_site_toml(d)
            write(d / "docs/site.toml", 'id = "io"\ntitle="T"\n')
            with self.assertRaises(bundle.BundleError):
                bundle.load_site_toml(d)


class Assemble(unittest.TestCase):
    def test_assemble_local(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            make_repo(d / "r")
            bundle.build_bundle(d / "r", d / "b")
            site = d / "site"
            write(site / "projects.toml",
                  'global_prefixes = ["tags"]\n[[project]]\nid="io"\nrepo="O/R"\ntitle="T"\norder=5\n'
                  '[[project]]\nid="ast"\nrepo="O/A"\ntitle="A"\n')
            rc = assemble.main(["--site-root", str(site), "--local", f"io={d / 'b'}",
                                "--local", f"ast={d / 'missing'}"])
            self.assertEqual(rc, 0)
            self.assertIn("[g](/io/guide/)", (site / "content/io/_index.md").read_text())
            self.assertIn("/io/plots/p.png", (site / "content/io/guide/a.md").read_text())
            self.assertTrue((site / "static/io/plots/p.png").is_file())
            data = json.loads((site / "data/projects/io.json").read_text())
            self.assertEqual(data["order"], 5)
            self.assertFalse((site / "data/projects/ast.json").exists())

    def test_generates_index_when_missing(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            make_repo(d / "r")
            (d / "r/docs/content/_index.md").unlink()
            bundle.build_bundle(d / "r", d / "b")
            site = d / "site"
            write(site / "projects.toml", '[[project]]\nid="io"\nrepo="O/R"\ntitle="T"\n')
            assemble.main(["--site-root", str(site), "--local", f"io={d / 'b'}"])
            self.assertIn('title: "T"', (site / "content/io/_index.md").read_text())


if __name__ == "__main__":
    unittest.main()
