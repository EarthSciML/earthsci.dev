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
            # a bundle is one version's complete tree, API reference included
            self.assertTrue((d / "out/content/generated/api/_index.md").is_file())
            self.assertEqual(meta["version"], "dev")
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
                write(d / f"r-{version}/docs/content/guide/a.md", marker)
                out = d / f"b-{version}"
                bundle.build_bundle(d / f"r-{version}", out, version)
                return out

            info = versions.merge(existing, build("dev", "dev-1"), "dev")
            self.assertEqual((info["latest"], info["versions"]), (None, ["dev"]))
            versions.merge(existing, build("v1.0.0", "tag-1.0"), "v1.0.0")
            info = versions.merge(existing, build("v1.1.0", "tag-1.1"), "v1.1.0")
            self.assertEqual(info["versions"], ["dev", "v1.1.0", "v1.0.0"])
            self.assertEqual(info["latest"], "v1.1.0")
            # every version is a complete tree
            for v in ("dev", "v1.0.0", "v1.1.0"):
                self.assertTrue((existing / v / "bundle.json").is_file())
                self.assertTrue((existing / v / "content/generated/api/_index.md").is_file())
                self.assertTrue((existing / v / "static/plots/p.png").is_file())
            self.assertEqual((existing / "v1.0.0/content/guide/a.md").read_text(), "tag-1.0")
            # a new dev build replaces dev only
            versions.merge(existing, build("dev", "dev-2"), "dev")
            self.assertEqual((existing / "dev/content/guide/a.md").read_text(), "dev-2")
            self.assertEqual((existing / "v1.0.0/content/guide/a.md").read_text(), "tag-1.0")
            self.assertEqual(json.loads((existing / "versions.json").read_text())["latest"], "v1.1.0")

    def test_merge_migrates_the_legacy_flat_layout(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            existing = d / "branch"
            # the first layout: main's flat tree plus api/<version>/ trees
            write(existing / "bundle.json", json.dumps({"id": "io", "title": "T", "description": "D"}))
            write(existing / "content/guide/a.md", "legacy")
            write(existing / "static/plots/p.png", "png")
            write(existing / "api/dev/_index.md", "legacy api dev")
            write(existing / "api/v0.9.0/_index.md", "legacy api tag")
            write(existing / "versions.json", "{}")
            make_repo(d / "r", api=True)
            bundle.build_bundle(d / "r", d / "b", "v1.0.0")
            info = versions.merge(existing, d / "b", "v1.0.0")
            self.assertEqual(info["versions"], ["dev", "v1.0.0"])
            self.assertEqual((existing / "dev/content/guide/a.md").read_text(), "legacy")
            self.assertEqual((existing / "dev/content/generated/api/_index.md").read_text(),
                             "legacy api dev")
            self.assertTrue((existing / "dev/static/plots/p.png").is_file())
            for gone in ("api", "content", "static", "bundle.json"):
                self.assertFalse((existing / gone).exists(), gone)

    def test_bundle_rejects_bad_version(self):
        with tempfile.TemporaryDirectory() as d:
            make_repo(Path(d) / "r", api=True)
            with self.assertRaises(bundle.BundleError):
                bundle.build_bundle(Path(d) / "r", Path(d) / "o", "release-1")

    def test_tag_bundle_is_a_complete_tree(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            make_repo(d / "r", api=True)
            meta = bundle.build_bundle(d / "r", d / "o", "v1.0.0")
            self.assertEqual(meta["version"], "v1.0.0")
            self.assertTrue((d / "o/content/guide/a.md").is_file())
            self.assertTrue((d / "o/content/generated/api/_index.md").is_file())
            self.assertTrue((d / "o/static/plots/p.png").is_file())


class VersionedLinks(unittest.TestCase):
    def rw(self, s, ver="v1.0.0"):
        return assemble.rewrite_links(s, "io", SKIP, version=ver)

    def test_links_stay_within_the_version(self):
        self.assertEqual(self.rw("[a](/generated/api/ref/x) [b](/guide/) [c](/) ![i](/plots/p.png)"),
                         "[a](/io/v1.0.0/generated/api/ref/x) [b](/io/v1.0.0/guide/) "
                         "[c](/io/v1.0.0/) ![i](/io/v1.0.0/plots/p.png)")

    def test_html_reference_and_shortcode_links(self):
        self.assertEqual(self.rw('<a href="/a/b">x</a>\n[r]: /docs/x\n{{< ref "/c" >}}'),
                         '<a href="/io/v1.0.0/a/b">x</a>\n[r]: /io/v1.0.0/docs/x\n{{< ref "/io/v1.0.0/c" >}}')

    def test_global_prefixes_and_project_links_untouched(self):
        s = "[t](/tags/x/) [d](/domains) [p](/io/guide/)"
        self.assertEqual(self.rw(s), s)

    def test_idempotent(self):
        once = self.rw("[a](/generated/api/ref/x)")
        self.assertEqual(self.rw(once), once)

    def test_code_left_alone(self):
        s = "`[a](/guide/)`\n```\n[a](/guide/)\n```\n"
        self.assertEqual(self.rw(s), s)


class StripFrontMatter(unittest.TestCase):
    KEYS = {"family", "families", "tags", "tag", "url", "aliases"}

    def strip(self, s):
        return assemble.strip_front_matter_keys(s, self.KEYS)

    def test_yaml_scalars_inline_and_block_lists(self):
        s = ("---\ntitle: A\nfamilies:\n  - x\n  - y\ntags: [a, b]\nurl: /foo/\n"
             "aliases:\n- /old/\nweight: 3\n---\nbody [t](/tags/a/)\n")
        self.assertEqual(self.strip(s), "---\ntitle: A\nweight: 3\n---\nbody [t](/tags/a/)\n")

    def test_yaml_keeps_lookalike_keys_and_nested(self):
        s = "---\ntitle: A\nmytags: z\nparams:\n  tags: [k]\n---\nbody\n"
        self.assertEqual(self.strip(s), s)

    def test_toml(self):
        s = ('+++\ntitle = "A"\nfamilies = [\n  "x",\n  "y",\n]\ntags = ["a"]\n'
             'weight = 3\n[params]\ntags = ["k"]\n+++\nbody\n')
        self.assertEqual(self.strip(s), '+++\ntitle = "A"\nweight = 3\n[params]\ntags = ["k"]\n+++\nbody\n')

    def test_no_front_matter_or_body_dashes(self):
        s = "no front matter\n---\ntags: x\n---\n"
        self.assertEqual(self.strip(s), s)


class AssembleVersions(unittest.TestCase):
    def build_site(self, d, tags, with_dev=True, hugo_toml=""):
        """A docs-bundle branch holding the given versions, each with distinguishable pages."""
        branch = d / "branch"
        order_ = (["dev"] if with_dev else []) + tags
        for v in order_:
            make_repo(d / f"r-{v}", api=True)
            write(d / f"r-{v}/docs/content/_index.md", f"---\ntitle: T\n---\nhome {v}\n[g](/guide/)\n")
            write(d / f"r-{v}/docs/content/guide/a.md",
                  f"---\ntitle: A\nfamilies: [f-{v}]\ntags: [t]\nurl: /custom-{v}/\n---\n"
                  f"guide {v} ![i](/plots/p.png) [x](/tags/t/)\n")
            write(d / f"r-{v}/docs/content/generated/api/_index.md",
                  f"api {v}\n[ref](/generated/api/ref/x)\n")
            write(d / f"r-{v}/docs/content/generated/api/ref/x.md",
                  f"ref {v}\n[home](/generated/api/)")
            if v != "dev":  # a nav that differs per version
                write(d / f"r-{v}/docs/site.toml",
                      f'id = "io"\ntitle = "T"\ndescription = "D"\nrepo = "O/R"\n'
                      f'[[nav]]\ntitle = "Guide {v}"\npath = "guide/"\n')
            out = d / f"b-{v}"
            bundle.build_bundle(d / f"r-{v}", out, v)
            versions.merge(branch, out, v)
        site = d / "site"
        write(site / "hugo.toml", hugo_toml or '[taxonomies]\n  family = "families"\n  tag = "tags"\n')
        write(site / "projects.toml",
              'global_prefixes = ["tags", "families"]\nmax_minor_versions = 8\n'
              '[[project]]\nid="io"\nrepo="O/R"\ntitle="T"\n')
        self.assertEqual(assemble.main(["--site-root", str(site), "--local", f"io={branch}"]), 0)
        return site

    def test_latest_is_default_and_every_version_has_its_own_tree(self):
        with tempfile.TemporaryDirectory() as d:
            site = self.build_site(Path(d), ["v1.0.0", "v1.1.0"])
            c = site / "content/io"
            # latest at /io/ (unversioned links), a full tree per version
            self.assertIn("home v1.1.0", (c / "_index.md").read_text())
            self.assertIn("[g](/io/guide/)", (c / "_index.md").read_text())
            self.assertIn("[ref](/io/generated/api/ref/x)", (c / "generated/api/_index.md").read_text())
            for v in ("dev", "v1.0.0", "v1.1.0"):
                self.assertIn(f"home {v}", (c / v / "_index.md").read_text())
                self.assertIn(f"[g](/io/{v}/guide/)", (c / v / "_index.md").read_text())
                self.assertIn(f"[ref](/io/{v}/generated/api/ref/x)",
                              (c / v / "generated/api/_index.md").read_text())
                self.assertIn(f"[home](/io/{v}/generated/api/)",
                              (c / v / "generated/api/ref/x.md").read_text())
                self.assertIn(f"/io/{v}/plots/p.png", (c / v / "guide/a.md").read_text())
                self.assertTrue((site / "static/io" / v / "plots/p.png").is_file())
            self.assertTrue((site / "static/io/plots/p.png").is_file())
            data = json.loads((site / "data/projects/io.json").read_text())
            self.assertEqual(data["versions"]["latest"], "v1.1.0")
            self.assertTrue(data["versions"]["latest_is_release"])
            self.assertEqual(data["versions"]["list"], ["dev", "v1.1.0", "v1.0.0"])

    def test_nav_is_per_version(self):
        with tempfile.TemporaryDirectory() as d:
            site = self.build_site(Path(d), ["v1.0.0"])
            data = json.loads((site / "data/projects/io.json").read_text())
            titles = lambda n: [x["title"] for x in n]
            self.assertEqual(titles(data["nav"]), ["Guide v1.0.0", "API"])  # latest's nav
            self.assertEqual(titles(data["nav_by_version"]["dev"]), ["Guide", "API"])
            self.assertEqual(titles(data["nav_by_version"]["v1.0.0"]), ["Guide v1.0.0", "API"])

    def test_only_latest_feeds_the_taxonomies(self):
        with tempfile.TemporaryDirectory() as d:
            site = self.build_site(Path(d), ["v1.0.0", "v1.1.0"])
            c = site / "content/io"
            latest = (c / "guide/a.md").read_text()
            self.assertIn("families: [f-v1.1.0]", latest)
            self.assertIn("url: /custom-v1.1.0/", latest)
            for v in ("dev", "v1.0.0", "v1.1.0"):
                fm = (c / v / "guide/a.md").read_text().split("---")[1]
                for key in ("families", "tags", "url", "aliases"):
                    self.assertNotIn(key + ":", fm)
                self.assertIn("title: A", fm)

    def test_dev_only_project_serves_dev_as_latest(self):
        with tempfile.TemporaryDirectory() as d:
            site = self.build_site(Path(d), [])
            data = json.loads((site / "data/projects/io.json").read_text())["versions"]
            self.assertEqual((data["latest"], data["latest_is_release"], data["list"]),
                             ("dev", False, ["dev"]))
            self.assertIn("home dev", (site / "content/io/_index.md").read_text())
            self.assertTrue((site / "content/io/dev/_index.md").is_file())

    def test_site_cap_drops_old_minor_series(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            site = self.build_site(d, ["v1.0.0", "v1.1.0", "v1.2.0"])
            write(site / "projects.toml",
                  'max_minor_versions = 2\n[[project]]\nid="io"\nrepo="O/R"\ntitle="T"\n')
            assemble.main(["--site-root", str(site), "--local", f"io={d / 'branch'}"])
            self.assertFalse((site / "content/io/v1.0.0").exists())
            self.assertFalse((site / "static/io/v1.0.0").exists())
            self.assertTrue((site / "content/io/v1.1.0").is_dir())

    def test_single_tree_bundle_is_one_version(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            make_repo(d / "r", api=True)
            bundle.build_bundle(d / "r", d / "b", "v1.0.0")
            site = d / "site"
            write(site / "projects.toml", '[[project]]\nid="io"\nrepo="O/R"\ntitle="T"\n')
            assemble.main(["--site-root", str(site), "--local", f"io={d / 'b'}"])
            data = json.loads((site / "data/projects/io.json").read_text())["versions"]
            self.assertEqual((data["latest"], data["list"]), ("v1.0.0", ["v1.0.0"]))
            self.assertTrue((site / "content/io/v1.0.0/guide/a.md").is_file())

    def test_legacy_flat_bundle_assembles_as_dev(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            old = d / "old"
            write(old / "bundle.json", json.dumps({"id": "io", "title": "T", "description": "D",
                                                   "repo": "O/R", "nav": []}))
            write(old / "content/_index.md", "---\ntitle: T\n---\nlegacy [g](/guide/)\n")
            write(old / "api/dev/_index.md", "legacy api\n[r](/generated/api/ref/x)")
            site = d / "site"
            write(site / "projects.toml", '[[project]]\nid="io"\nrepo="O/R"\ntitle="T"\n')
            assemble.main(["--site-root", str(site), "--local", f"io={old}"])
            data = json.loads((site / "data/projects/io.json").read_text())["versions"]
            self.assertEqual((data["latest"], data["list"]), ("dev", ["dev"]))
            self.assertIn("legacy [g](/io/guide/)", (site / "content/io/_index.md").read_text())
            self.assertIn("[r](/io/dev/generated/api/ref/x)",
                          (site / "content/io/dev/generated/api/_index.md").read_text())

    def test_project_without_api(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            make_repo(d / "r")
            bundle.build_bundle(d / "r", d / "b")
            site = d / "site"
            write(site / "projects.toml", '[[project]]\nid="io"\nrepo="O/R"\ntitle="T"\n')
            assemble.main(["--site-root", str(site), "--local", f"io={d / 'b'}"])
            self.assertFalse((site / "content/io/generated/api").exists())
            self.assertEqual([n["title"] for n in
                              json.loads((site / "data/projects/io.json").read_text())["nav"]], ["Guide"])


if __name__ == "__main__":
    unittest.main()
