import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
import assemble  # noqa: E402
import bundle  # noqa: E402
import publish  # noqa: E402,F401
import versions  # noqa: E402

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
            # the API reference is versioned, so it leaves content/ for api/dev/
            self.assertFalse((d / "out/content/generated/api").exists())
            self.assertTrue((d / "out/api/dev/_index.md").is_file())
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


VERSION_NAMES = ["dev", "v1.0.0", "v1.2.0", "v1.10.0", "v1.10.1", "v2.0.0-rc.1", "v1.2.1", "junk"]


class Versions(unittest.TestCase):
    def test_order_is_dev_then_semver_descending(self):
        self.assertEqual(versions.order(VERSION_NAMES),
                         ["dev", "v2.0.0-rc.1", "v1.10.1", "v1.10.0", "v1.2.1", "v1.2.0", "v1.0.0"])

    def test_prerelease_ranks_below_its_release(self):
        self.assertEqual(versions.order(["v1.0.0-rc.1", "v1.0.0", "v1.0.0-alpha"]),
                         ["v1.0.0", "v1.0.0-rc.1", "v1.0.0-alpha"])

    def test_latest_skips_prereleases_and_dev(self):
        self.assertEqual(versions.latest(VERSION_NAMES), "v1.10.1")
        self.assertIsNone(versions.latest(["dev", "v1.0.0-rc.1"]))
        self.assertEqual(versions.effective_latest(["dev"]), "dev")
        self.assertEqual(versions.effective_latest(["dev", "v0.1.0"]), "v0.1.0")
        self.assertIsNone(versions.effective_latest([]))

    def test_prune_keeps_newest_patch_of_newest_minors(self):
        self.assertEqual(versions.prune(VERSION_NAMES, 0),
                         ["dev", "v2.0.0-rc.1", "v1.10.1", "v1.2.1", "v1.0.0"])
        self.assertEqual(versions.prune(VERSION_NAMES, 2), ["dev", "v2.0.0-rc.1", "v1.10.1", "v1.2.1"])

    def test_prune_drops_prerelease_older_than_latest(self):
        self.assertEqual(versions.prune(["v1.0.0-rc.1", "v1.0.0"], 0), ["v1.0.0"])

    def test_merge_dev_then_tags(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            existing = d / "branch"

            def build(version, marker):
                make_repo(d / f"r-{version}", api=True)
                write(d / f"r-{version}/docs/content/generated/api/_index.md", marker)
                out = d / f"b-{version}"
                bundle.build_bundle(d / f"r-{version}", out, version)
                return out

            info = versions.merge(existing, build("dev", "dev-1"), "dev")
            self.assertEqual((info["latest"], info["versions"]), (None, ["dev"]))
            versions.merge(existing, build("v1.0.0", "tag-1.0"), "v1.0.0")
            info = versions.merge(existing, build("v1.1.0", "tag-1.1"), "v1.1.0")
            self.assertEqual(info["versions"], ["dev", "v1.1.0", "v1.0.0"])
            self.assertEqual(info["latest"], "v1.1.0")
            # a tag build leaves main's narrative content and other versions alone
            self.assertTrue((existing / "content/guide/a.md").is_file())
            self.assertEqual((existing / "api/v1.0.0/_index.md").read_text(), "tag-1.0")
            # a new dev build replaces dev only
            versions.merge(existing, build("dev", "dev-2"), "dev")
            self.assertEqual((existing / "api/dev/_index.md").read_text(), "dev-2")
            self.assertEqual((existing / "api/v1.0.0/_index.md").read_text(), "tag-1.0")
            self.assertEqual(json.loads((existing / "versions.json").read_text())["latest"], "v1.1.0")

    def test_bundle_rejects_bad_version(self):
        with tempfile.TemporaryDirectory() as d:
            make_repo(Path(d) / "r", api=True)
            with self.assertRaises(bundle.BundleError):
                bundle.build_bundle(Path(d) / "r", Path(d) / "o", "release-1")

    def test_tag_bundle_has_only_the_api_tree(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            make_repo(d / "r", api=True)
            bundle.build_bundle(d / "r", d / "o", "v1.0.0")
            self.assertTrue((d / "o/api/v1.0.0/_index.md").is_file())
            self.assertFalse((d / "o/content").exists())


class ApiVersionLinks(unittest.TestCase):
    def rw(self, s, ver="v1.0.0"):
        return assemble.rewrite_links(s, "io", SKIP, api_version=ver)

    def test_api_links_stay_within_the_version(self):
        self.assertEqual(self.rw("[a](/generated/api/ref/x) [b](/generated/api) [c](/generated/api/)"),
                         "[a](/io/generated/api/v1.0.0/ref/x) [b](/io/generated/api/v1.0.0) "
                         "[c](/io/generated/api/v1.0.0/)")

    def test_html_and_non_api_links(self):
        self.assertEqual(self.rw('<a href="/generated/api/ref/x">x</a> [n](/generated/narrative/d) [o](/guide/)'),
                         '<a href="/io/generated/api/v1.0.0/ref/x">x</a> '
                         '[n](/io/generated/narrative/d) [o](/io/guide/)')

    def test_lookalike_paths_untouched(self):
        self.assertEqual(self.rw("[a](/generated/apiary/x)"), "[a](/io/generated/apiary/x)")

    def test_idempotent(self):
        once = self.rw("[a](/generated/api/ref/x)")
        self.assertEqual(self.rw(once), once)
        self.assertEqual(assemble.rewrite_links(once, "io", SKIP), once)

    def test_code_left_alone(self):
        s = "`[a](/generated/api/x)`\n```\n[a](/generated/api/x)\n```\n"
        self.assertEqual(self.rw(s), s)


class AssembleVersions(unittest.TestCase):
    def build_site(self, d, tags, with_dev=True):
        """A bundle holding the given versions, each with a distinguishable page."""
        branch = d / "branch"
        order_ = (["dev"] if with_dev else []) + tags
        for v in order_:
            make_repo(d / f"r-{v}", api=True)
            write(d / f"r-{v}/docs/content/generated/api/_index.md",
                  f"home {v}\n[ref](/generated/api/ref/x)\n[g](/guide/)\n")
            write(d / f"r-{v}/docs/content/generated/api/ref/x.md", f"ref {v}\n[home](/generated/api/)")
            out = d / f"b-{v}"
            bundle.build_bundle(d / f"r-{v}", out, v)
            versions.merge(branch, out, v)
        site = d / "site"
        write(site / "projects.toml",
              'max_minor_versions = 8\n[[project]]\nid="io"\nrepo="O/R"\ntitle="T"\n')
        self.assertEqual(assemble.main(["--site-root", str(site), "--local", f"io={branch}"]), 0)
        return site

    def test_latest_is_default_and_every_version_has_its_own_tree(self):
        with tempfile.TemporaryDirectory() as d:
            site = self.build_site(Path(d), ["v1.0.0", "v1.1.0"])
            api = site / "content/io/generated/api"
            self.assertIn("home v1.1.0", (api / "_index.md").read_text())
            self.assertIn("[ref](/io/generated/api/ref/x)", (api / "_index.md").read_text())
            for v in ("dev", "v1.0.0", "v1.1.0"):
                idx = (api / v / "_index.md").read_text()
                self.assertIn(f"home {v}", idx)
                self.assertIn(f"[ref](/io/generated/api/{v}/ref/x)", idx)
                self.assertIn("[g](/io/guide/)", idx)
                self.assertIn(f"[home](/io/generated/api/{v}/)", (api / v / "ref/x.md").read_text())
            data = json.loads((site / "data/projects/io.json").read_text())["api"]
            self.assertEqual(data["latest"], "v1.1.0")
            self.assertTrue(data["latest_is_release"])
            self.assertEqual(data["versions"], ["dev", "v1.1.0", "v1.0.0"])

    def test_dev_only_project_serves_dev_as_latest(self):
        with tempfile.TemporaryDirectory() as d:
            site = self.build_site(Path(d), [])
            data = json.loads((site / "data/projects/io.json").read_text())["api"]
            self.assertEqual((data["latest"], data["latest_is_release"], data["versions"]),
                             ("dev", False, ["dev"]))
            self.assertIn("home dev", (site / "content/io/generated/api/_index.md").read_text())
            self.assertTrue((site / "content/io/generated/api/dev/_index.md").is_file())

    def test_site_cap_drops_old_minor_series(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            site = self.build_site(d, ["v1.0.0", "v1.1.0", "v1.2.0"])
            write(site / "projects.toml",
                  'max_minor_versions = 2\n[[project]]\nid="io"\nrepo="O/R"\ntitle="T"\n')
            assemble.main(["--site-root", str(site), "--local", f"io={d / 'branch'}"])
            self.assertFalse((site / "content/io/generated/api/v1.0.0").exists())
            self.assertTrue((site / "content/io/generated/api/v1.1.0").is_dir())

    def test_project_without_api(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            make_repo(d / "r")
            bundle.build_bundle(d / "r", d / "b")
            site = d / "site"
            write(site / "projects.toml", '[[project]]\nid="io"\nrepo="O/R"\ntitle="T"\n')
            assemble.main(["--site-root", str(site), "--local", f"io={d / 'b'}"])
            data = json.loads((site / "data/projects/io.json").read_text())["api"]
            self.assertEqual((data["latest"], data["versions"]), (None, []))
            self.assertFalse((site / "content/io/generated/api").exists())


if __name__ == "__main__":
    unittest.main()
