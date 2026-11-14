import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools" / "api_docs"))
import bundle  # noqa: E402
import extract_lean as lean  # noqa: E402

SHA = "0123456789abcdef0123456789abcdef01234567"

NOTATION = textwrap.dedent('''\
    import Mathlib.Data.List.Nodup

    /-!
    # Notation tables

    The tables, as data. More text: `Tier` and `Esm.Doc.Math`.

    Second paragraph.
    -/

    namespace Esm.Doc

    /-! ## Operators -/

    /-- How an operator is written. -/
    inductive Tier where
      /-- conventional -/
      | native
      | call
      deriving DecidableEq

    /-- One row. -/
    structure Row where
      op : String
      /-- The tier. -/
      tier : Tier := .native

    open Tier in
    /-- Every row; see `Row`. -/
    def rows : List Row := [
      { op := "+" }, { op := "-- not a comment" }]

    /-- Rows cover. -/
    @[simp] theorem rows_cover (o : String)
        (h : o ∈ ["+"]) : ∃ r ∈ rows, r.op = o := by
      -- sorry in a comment is not an escape
      simp [rows] at *

    theorem undocumented : True := trivial

    /-- A pattern-matching definition. -/
    def tierName : Tier → String
      | .native => "native"
      | .call => "call /- not a comment -/"

    /-- An instance with `where`. -/
    instance : Inhabited Row where
      default := { op := "" }

    end Esm.Doc

    /-- At the root. -/
    def _root_.topLevel (x : Nat) : Nat := x
''')

MATH = textwrap.dedent('''\
    import Esm.Doc.Notation

    /-!
    # Math: a grammar

    Uses `Esm.Doc.Tier`.
    -/

    namespace Esm.Doc.Math

    /-- Not yet proved. -/
    theorem open_problem : 1 = 1 := by
      sorry

    end Esm.Doc.Math
''')


def write(p: Path, text: str):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


class Scan(unittest.TestCase):
    def test_comments_and_strings_are_masked(self):
        src = 'def x := "a -- b /- c" -- tail\n/- outer /- inner -/ still -/ y\n'
        spans, masked = lean.scan(src)
        self.assertEqual(len(masked), len(src))
        self.assertEqual([s.kind for s in spans], ["string", "comment", "comment"])
        self.assertNotIn("inner", masked)
        self.assertNotIn("still", masked)
        self.assertTrue(masked.rstrip().endswith("y"))

    def test_raw_strings_and_chars_and_primes(self):
        src = 'def a\' := r#"x "-- y" z"# ++ \'"\' ++ r"\\q"\n'
        spans, masked = lean.scan(src)
        self.assertEqual([s.kind for s in spans], ["string", "string", "string"])
        self.assertIn("a'", masked)  # a name's prime is not a char literal
        self.assertNotIn("--", masked)

    def test_doc_kinds(self):
        spans, _ = lean.scan("/-! m -/\n/-- d -/\n/- c -/\n")
        self.assertEqual([s.kind for s in spans], ["moddoc", "doc", "comment"])

    def test_comment_body_dedents(self):
        self.assertEqual(lean.comment_body("/-- First.\n\n  * a\n    b -/"), "First.\n\n* a\n  b")


class ParseModule(unittest.TestCase):
    def setUp(self):
        self.mod = lean.parse_module(NOTATION, "Esm.Doc.Notation", "formal/Esm/Doc/Notation.lean")
        self.by_name = {d.full: d for d in self.mod.decls}

    def test_title_summary_and_prose_order(self):
        self.assertEqual(self.mod.title, "Notation tables")
        self.assertEqual(self.mod.summary, "The tables, as data.")
        self.assertEqual(self.mod.imports, ["Mathlib.Data.List.Nodup"])
        kinds = [type(i).__name__ if not isinstance(i, str) else "prose" for i in self.mod.items]
        self.assertEqual(kinds, ["prose", "prose", "Decl", "Decl", "Decl", "Decl", "Decl",
                                 "Decl", "Decl"])
        self.assertTrue(self.mod.items[0].startswith("The tables"))
        self.assertEqual(self.mod.items[1], "## Operators")

    def test_names_are_namespace_qualified(self):
        self.assertEqual(sorted(self.by_name), [
            "", "Esm.Doc.Row", "Esm.Doc.Tier", "Esm.Doc.rows", "Esm.Doc.rows_cover",
            "Esm.Doc.tierName", "topLevel"])
        inst = self.by_name[""]
        self.assertEqual(inst.kind, "instance")
        self.assertEqual(inst.anchor, f"Esm.Doc.Notation.L{inst.line}")

    def test_structures_and_inductives_shown_whole(self):
        tier = self.by_name["Esm.Doc.Tier"]
        self.assertEqual(tier.signature.splitlines()[0], "inductive Tier where")
        self.assertIn("/-- conventional -/", tier.signature)
        self.assertTrue(tier.signature.endswith("deriving DecidableEq"))
        row = self.by_name["Esm.Doc.Row"]
        self.assertTrue(row.signature.endswith("tier : Tier := .native"))

    def test_signature_stops_at_body(self):
        self.assertEqual(self.by_name["Esm.Doc.rows"].signature, "def rows : List Row")
        self.assertEqual(self.by_name["Esm.Doc.rows_cover"].signature,
                         "@[simp] theorem rows_cover (o : String)\n"
                         "    (h : o ∈ [\"+\"]) : ∃ r ∈ rows, r.op = o")
        self.assertEqual(self.by_name["Esm.Doc.tierName"].signature,
                         "def tierName : Tier → String")
        self.assertEqual(self.by_name[""].signature, "instance : Inhabited Row")
        self.assertEqual(self.by_name["topLevel"].signature,
                         "def _root_.topLevel (x : Nat) : Nat")

    def test_lines_and_docs(self):
        cover = self.by_name["Esm.Doc.rows_cover"]
        lines = NOTATION.splitlines()
        self.assertTrue(lines[cover.line - 1].startswith("@[simp] theorem rows_cover"))
        self.assertTrue(lines[cover.end_line - 1].strip().startswith("simp [rows]"))
        self.assertEqual(cover.doc, "Rows cover.")
        self.assertFalse(cover.escape)  # `sorry` only in a comment

    def test_undocumented_theorems_counted(self):
        self.assertEqual(self.mod.undocumented_theorems, 1)

    def test_sorry_in_proof(self):
        mod = lean.parse_module(MATH, "Esm.Doc.Math", "m.lean")
        self.assertTrue(mod.decls[0].escape)


class Modules(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        for name in ("Esm/Doc/Notation.lean", "Esm/Doc/Math/B.lean", "Esm/Doc/Math/A.lean",
                     "Esm/Doc/Math/Deep/C.lean", "Esm/Units/Core.lean"):
            write(self.root / name, "-- x\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_order_globs_and_dedupe(self):
        names, warnings = lean.expand_modules(self.root, [
            "Esm.Doc.Math.B", "Esm.Doc.Math.*", "Esm.Doc.**", "Esm.Missing.**"])
        self.assertEqual(names, ["Esm.Doc.Math.B", "Esm.Doc.Math.A", "Esm.Doc.Math.Deep.C",
                                 "Esm.Doc.Notation"])
        self.assertEqual(len(warnings), 1)
        self.assertIn("Esm.Missing.**", warnings[0])

    def test_missing_module_is_an_error(self):
        with self.assertRaises(lean.LeanError):
            lean.expand_modules(self.root, ["Esm.Doc.Nope"])


class Build(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.repo = Path(self.tmp.name)
        write(self.repo / "docs/site.toml", textwrap.dedent('''\
            id = "ast"
            title = "EarthSciAST"
            description = "d"
            repo = "EarthSciML/EarthSciAST"
        '''))
        write(self.repo / "docs/content/_index.md", "---\ntitle: x\n---\n")
        write(self.repo / "docs/lean/intro.md", "Intro prose.\n")
        write(self.repo / "docs/lean/config.toml", textwrap.dedent('''\
            root = "formal"
            [[chapter]]
            title = "Document formats"
            description = "The formats."
            out = "docs/content/generated/document-formats"
            intro = "docs/lean/intro.md"
            modules = ["Esm.Doc.Notation", "Esm.Doc.**"]
        '''))
        write(self.repo / "formal/Esm/Doc/Notation.lean", NOTATION)
        self.out = self.repo / "docs/content/generated/document-formats"

    def tearDown(self):
        self.tmp.cleanup()

    def build(self):
        chapters, warnings = lean.build(self.repo, sha=SHA)
        self.assertEqual(warnings, [])
        return chapters[0]

    def page(self, slug):
        return (self.out / f"{slug}.md").read_text(encoding="utf-8")

    def test_pages_and_links(self):
        write(self.repo / "formal/Esm/Doc/Math.lean", MATH.replace("  sorry", "  rfl"))
        ch = self.build()
        self.assertEqual([m.name for m in ch.modules], ["Esm.Doc.Notation", "Esm.Doc.Math"])
        self.assertEqual(ch.url, "/generated/document-formats")
        landing = (self.out / "_index.md").read_text(encoding="utf-8")
        self.assertIn('title: "Document formats"\ntype: lean\ndescription: "The formats."', landing)
        self.assertIn("Intro prose.", landing)
        self.assertIn("- [Notation tables](/generated/document-formats/esm-doc-notation) "
                      "(`Esm.Doc.Notation`): The tables, as data.", landing)

        page = self.page("esm-doc-notation")
        self.assertTrue(page.startswith('---\ntitle: "Notation tables"\ntype: lean\n'))
        self.assertIn('<div class="lean-decl lean-theorem" id="Esm.Doc.rows_cover">', page)
        self.assertIn('<span class="lean-badge lean-badge-proved"', page)
        self.assertIn(f"https://github.com/EarthSciML/EarthSciAST/blob/{SHA}/"
                      "formal/Esm/Doc/Notation.lean#L", page)
        self.assertIn("```lean\ndef rows : List Row\n```", page)
        # Code spans that name a declaration or module become links.
        self.assertIn("[`Tier`](#Esm.Doc.Tier)", page)  # from prose outside the namespace
        self.assertIn("[`Row`](#Esm.Doc.Row)", page)
        self.assertIn("[`Esm.Doc.Math`](/generated/document-formats/esm-doc-math)", page)
        self.assertIn("also proves 1 lemma without doc comments", page)
        self.assertNotIn("# Notation tables", page)  # the title is front matter
        math = self.page("esm-doc-math")
        self.assertIn("[`Esm.Doc.Tier`](/generated/document-formats/esm-doc-notation"
                      "#Esm.Doc.Tier)", math)
        self.assertIn('<a href="/generated/document-formats/esm-doc-notation">'
                      "<code>Esm.Doc.Notation</code></a>", math)

    def test_sorry_anywhere_withholds_proved(self):
        write(self.repo / "formal/Esm/Doc/Math.lean", MATH)
        self.build()
        self.assertNotIn("lean-badge-proved", self.page("esm-doc-notation"))
        self.assertIn("lean-badge-sorry", self.page("esm-doc-math"))
        self.assertNotIn("**proved**", (self.out / "_index.md").read_text(encoding="utf-8"))

    def test_rebuild_replaces_output(self):
        write(self.out / "stale.md", "x")
        self.build()
        self.assertFalse((self.out / "stale.md").exists())

    def test_out_outside_content_is_an_error(self):
        cfg = self.repo / "docs/lean/config.toml"
        cfg.write_text(cfg.read_text().replace("docs/content/generated", "elsewhere"))
        with self.assertRaises(lean.LeanError):
            lean.build(self.repo, sha=SHA)

    def test_bundle_warns_when_chapter_not_generated(self):
        self.assertEqual(bundle.missing_lean_chapters(self.repo),
                         ["docs/content/generated/document-formats"])
        self.build()
        self.assertEqual(bundle.missing_lean_chapters(self.repo), [])


if __name__ == "__main__":
    unittest.main()
