# Adding a Supplementary Annotation to fastVEP

A complete, followable guide for adding a new supplementary-annotation (SA) source or field
end-to-end, the same way gnomAD and dbSNP were added.

This guide uses only already-public sources (gnomAD, dbSNP, ClinVar, dbNSFP) as worked
examples. Every file path and symbol below is a real anchor in this repository (verified
against branch `development`).

---

## 1. Overview: the SA system and the three record kinds

fastVEP annotates variants with data from external "supplementary annotation" (SA) sources.
Each source is compiled once into an on-disk index by `fastvep sa-build`, then discovered
automatically at `fastvep annotate` time from a directory you point `--sa-dir` at. There is
no manifest and no per-source annotate flag: discovery is purely by **file extension**.

An SA source produces exactly one of **three record kinds**, and the record kind — not the
source file — determines the on-disk file kind:

| Record kind | On-disk file | Keyed by | Producer example | Defined at |
|---|---|---|---|---|
| `AnnotationRecord` | `.osa` (+ `.osa.idx`) | chrom / pos / ref / alt | `gnomad`, `dbsnp` | `crates/fastvep-sa/src/common.rs:44` |
| `IntervalRecord` | `.osi` | chrom / start / end | only custom BED | `crates/fastvep-sa/src/common.rs:59` |
| `GeneRecord` | `.oga` | gene symbol | `gnomad_genes` | `crates/fastvep-sa/src/common.rs:109` |

Extensions are declared as constants: `OSA_EXT` / `OSI_EXT` / `OGA_EXT` (`common.rs:31`,
`:37`, `:40`).

> A fourth extension, `.osa2`, also exists and is handled by the annotate loader. It is an
> alternate on-disk **encoding** of the allele (`AnnotationRecord`) kind, **not** a fourth
> record kind — you do not choose it when writing a new source; `.osa` is the default.

### Choosing a record kind

- **Per-variant data** (a frequency, a pathogenicity call, a clinical significance, keyed to
  a specific base change) → `AnnotationRecord` → `.osa`. This is gnomAD and dbSNP. Most new
  sources are this kind.
- **Per-position data with no allele** (a conservation score at a genomic position) → still
  `AnnotationRecord` → `.osa`, but built with `is_positional: true` and empty `ref`/`alt`
  (PhyloP/GERP/DANN work this way).
- **Per-gene data** (a constraint metric keyed by symbol, not coordinates) → `GeneRecord` →
  `.oga`. This is gnomAD gene constraint.
- **Per-interval data** (a region/track) → `IntervalRecord` → `.osi`. No built-in source
  under `sources/` emits intervals today; the only `.osi` producer is the generic custom-BED
  path (`crates/fastvep-sa/src/custom.rs:160`). If your data is a BED region set and you
  don't need bespoke parsing, prefer `--source custom_bed` over writing a new parser.

> **The record kind is the single most important decision.** It picks which writer, which
> reader, which CLI dispatch branch, and which projection formatter your source flows
> through. Get it right first.

---

## 2. The pipeline at a glance

A new field (or a new source) crosses these layers in order:

1. **Parser** — `crates/fastvep-sa/src/sources/<name>.rs`. A free function that reads
   VCF/TSV and hand-builds a camelCase JSON string per record. It does **not** open files,
   decompress, declare headers, or write output.
2. **Module registry** — `crates/fastvep-sa/src/sources/mod.rs`. One `pub mod <name>;` line
   so the parser compiles into the crate.
3. **CLI dispatch (registration)** — `crates/fastvep-cli/src/pipeline.rs`. Two match arms:
   an `IndexHeader` arm (declares the `json_key` and lookup flags) and a parser arm (calls
   your parser). Plus the supported-source `bail!` list.
4. **Output / projection seam** — `crates/fastvep-io/src/output.rs`. The `*_FIELDS` table
   (UPPER_CASE label → camelCase subkey) and a `VcfProjectionSpec` that surface the field in
   VCF `FV_*` INFO and the TSV column. JSON output needs no change here — it carries the
   whole object verbatim.
5. **Docs** — `docs/SUPPLEMENTARY_ANNOTATIONS.md`. A verbatim copy of the pipe layout,
   enforced by a unit test that fails CI otherwise.
6. **Tests** — inline `#[cfg(test)]` parser tests (these run in CI) plus, for a new source,
   the test enumeration helpers in `output.rs`.
7. **Build & verify** — `cargo test --workspace --lib`, then `fastvep sa-build …` and
   `fastvep annotate --sa-dir …`.

The rest of the guide is these layers, step by step.

---

## 3. Step-by-step

Two scenarios share most of the steps:

- **Add a FIELD** to an existing source (e.g. a new gnomAD or dbSNP column) → steps
  **a, c, d, e, f** (skip registration; the source already dispatches).
- **Add a whole new SOURCE** → all steps **a–f**.

### 3a. Author or extend the source parser

**File:** `crates/fastvep-sa/src/sources/<name>.rs`

A source parser is a free function. For an allele/interval source it takes the chrom map; for
a **gene** source it omits it:

```rust
// allele/interval source
pub fn parse_<name><R: BufRead>(
    reader: R,
    chrom_map: &HashMap<String, u16>,
) -> anyhow::Result<Vec<AnnotationRecord>>

// gene source (no chrom_map)
pub fn parse_<name><R: BufRead>(reader: R) -> anyhow::Result<Vec<GeneRecord>>
```

**Templates to copy:**
- Simplest full VCF source → `dbsnp.rs` (`parse_dbsnp_vcf` at `dbsnp.rs:11`). Inline INFO
  parse, no column-detection struct.
- Richer VCF source with release-flavor INFO name detection → `gnomad.rs` (`parse_gnomad_vcf`
  at `:179`, `build_gnomad_json` at `:270`).
- TSV with header-driven column detection → copy the detection-struct shape from
  `DbNsfpColumns::from_header` (`dbnsfp.rs:163`) or `GnomadGeneCols` (`gnomad_gene.rs:128`).
- Gene source → `gnomad_gene.rs` (`parse_gnomad_gene_scores` at `:15`).

**The emitted-field / camelCase convention.** There is **no schema file**. A source
"declares" its fields implicitly by hand-formatting JSON key/value fragments and joining
them. Accumulate `parts`, push `"camelCaseKey":value` fragments, then wrap:

```rust
// dbsnp.rs:71 — the whole record JSON, built by hand
let mut parts = vec![format!("\"id\":\"{}\"", rs_id)];
if let Some(f) = freq {
    parts.push(format!("\"globalMaf\":{:.6e}", f));   // camelCase key, scientific float
}
let json = format!("{{{}}}", parts.join(","));
```

```rust
// gnomad.rs:284 / :319 / :322 / :326 — same pattern, ACMG fields
parts.push(format!("\"allAf\":{:.6e}", f));
parts.push(format!("\"grpmaxAf\":{:.6e}", f));
parts.push(format!("\"grpmaxGroup\":\"{}\"", escape_json(g)));  // STRING → escaped AND quoted
parts.push(format!("\"faf95\":{:.6e}", f));
```

Rules for the JSON you emit:
- **Keys are camelCase** (`allAf`, `grpmaxGroup`, `faf95`, `pLI`, `misZ`). This is the
  byte-exact key the projection layer will look up later.
- **Numeric floats use `{:.6e}`** (scientific) so tests can assert exact rendering, e.g.
  `"faf95":1.000000e-4`.
- **String values must be escaped _and_ wrapped in quotes** — note the `\"{}\"` around
  `escape_json(g)` above. `escape_json` (`gnomad.rs:379`) escapes the contents; the literal
  `\"…\"` makes it a JSON string. Emitting `"grpmaxGroup":afr` (unquoted) is invalid JSON.
  For object-shaped values, build via `serde_json::Map` instead (`custom.rs:209`).
- Emitting a **JSON-only** field (present in JSON but deliberately not projected to VCF/TSV)
  is legitimate and common — gnomAD's per-population `nhomalt`, sex-stratified AF, and region
  flags are pushed here and never added to `*_FIELDS`. This is a design choice, documented in
  prose.

**Contig normalization — copy the right one.** Do **not** blindly prepend `chr`. NCBI
releases (dbSNP, dbNSFP) name contigs by RefSeq accession (`NC_000001.11`). Mangling those to
`chrNC_000001.11` yields "0 records parsed." Copy `dbsnp.rs`'s `normalize_chrom` (`:101`),
which special-cases `fastvep_core::looks_like_refseq_accession`. Then look up
`chrom_map.get(chrom)`; records whose chrom is absent are silently dropped (this is the
zero-record failure mode). **When you copy `normalize_chrom`, drop or re-target its
`issue #51` comment** — it references a specific historical fix that won't apply to your
source.

**Multi-allelic handling.** Allele-indexed INFO fields (`AF`/`AC`/`nhomalt`, `Number=A`) are
comma-lists that must be split and indexed per-ALT; `Number=1` fields like `AN` are
single-valued and reused for every ALT. See `gnomad.rs:232-262`. A parser that ignores this
attaches allele-0's frequency to every alt.

**Sort before returning (allele/interval only).** The writer requires sorted input:

```rust
records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));  // gnomad.rs:265, dbsnp.rs:87
```

Gene parsers are **not** sorted — they return in first-seen order.

**Determinism.** If your JSON key order derives from iterating a `HashMap`, output bytes
become non-deterministic. Iterate a fixed slice (as gnomAD does over its population list) or
parse INFO into a `BTreeMap` (see the note at `custom.rs:238`).

### 3b. Register the module (new source only)

**File:** `crates/fastvep-sa/src/sources/mod.rs`

Add one line alongside the existing declarations (`dbsnp` is at `mod.rs:8`):

```rust
pub mod <name>;
```

### 3c. Register the source in the CLI dispatch (new source only)

**File:** `crates/fastvep-cli/src/pipeline.rs`, function `run_sa_build` (`:2200`)

`--source` is a **free `String`** in `main.rs:180` — there is no enum and no clap validation.
The only registration points are the two `match source` blocks plus the `bail!` list. Miss
either arm and you get the "Unknown source" bail or an `unreachable!()` panic at the parser
match.

**Gene sources dispatch separately** through `run_oga_build` (`:2664`). A gene source needs
**three** edits there, not one:
1. add its name to the early guard below so it routes to `run_oga_build` at all,
2. add its `(json_key, name)` header arm in `run_oga_build`'s header match (`:2668`), and
3. add its parser arm in `run_oga_build`'s parser match (`:2691`).
Each match also has its own `bail!` / `unreachable!` fallback (`:2674` / `:2699`), so miss an
arm and you get the same silent-drop / panic failure as the allele path.

```rust
// pipeline.rs:2212 — a new gene source must join this guard or it falls to the .osa path and fails
if matches!(source, "omim" | "gnomad_genes" | "gnomad_gene" | "clinvar_protein") {
    return run_oga_build(source, input, output, assembly);
}
```

**For an allele/positional `.osa` source, add two arms.**

Arm 1 — the header, in the `let header = match source { … }` block (`pipeline.rs:2235`; dbSNP
arm at `:2258`, gnomAD arm at `:2247`):

```rust
"<name>" => IndexHeader {
    schema_version: fastvep_sa::common::SCHEMA_VERSION,
    json_key: "<name>".into(),          // the output key AND downstream lookup key — load-bearing
    name: "<Display>".into(),
    version: "latest".into(),
    description: format!("… for {}", assembly),
    assembly: assembly.into(),
    match_by_allele: true,   // true → matches on ref+alt; false for positional-only sources
    is_array: false,
    is_positional: false,    // true → result wrapped as Positional, ref/alt ignored (scores)
},
```

`match_by_allele` vs `is_positional` decide the lookup semantics — set them wrong and lookups
silently miss or over-match. (`IndexHeader` is defined in `index.rs:27`; positional scores
like phylop use `is_positional: true`, `match_by_allele: false` at `pipeline.rs:2285`.)

Arm 2 — the parser dispatch, in the `let records = match source { … }` block
(`pipeline.rs:2418`, dbSNP at `:2421`):

```rust
"<name>" => fastvep_sa::sources::<name>::parse_<name>(buf_reader, &chrom_map)?,
```

The `_ => unreachable!()` at the end of that block (`:2434`) assumes the source was already
validated upstream — which is what the `bail!` list does.

**Update the supported-sources error.** Add `<name>` to the `bail!` list (`:2378-2379`) so
unknown sources still list the full set:

```rust
_ => anyhow::bail!(
    "Unknown source: {}. Supported: clinvar, gnomad, dbsnp, cosmic, onekg, topmed, mitomap, phylop, gerp, dann, revel, spliceai, primateai, dbnsfp, omim, gnomad_genes, clinvar_protein, custom_vcf, custom_bed, custom",
    …
),
```

> Adding **both** the header arm and the parser arm is mandatory. Miss the header arm and the
> source is written with the wrong/absent `json_key` and silently ignored downstream; miss
> the parser arm and you hit the `unreachable!()`. No change is needed in `main.rs` —
> `--source` is a free string.

### 3d. Wire the output / projection seam

**File:** `crates/fastvep-io/src/output.rs`

**There are three outputs but only two mechanisms.** JSON output (`format_json` at `:1381` /
`format_nirvana_json` at `:1628`) inserts the **entire** per-source object verbatim under its
`json_key` and does **not** consult `*_FIELDS`. So any camelCase key you emitted in step 3a
already appears in JSON with zero change here. Only the VCF `FV_*` INFO field and the TSV
column are positional projections driven by `*_FIELDS` — those are what you edit.

**Adding a FIELD to an existing source.** Append a `(LABEL, jsonKey)` tuple to the source's
`*_FIELDS` constant — **append only, never reorder**:

```rust
// GNOMAD_FIELDS at output.rs:345 — append after ("FILTER", "filter")
    ("NHOMALT", "nhomalt"),   // UPPER_CASE label ; byte-exact camelCase subkey from the parser
```

The smallest example is `DBSNP_FIELDS` at `output.rs:369`:

```rust
const DBSNP_FIELDS: &[(&str, &str)] = &[("ID", "id"), ("GLOBAL_MAF", "globalMaf")];
```

`DBNSFP_FIELDS` (`output.rs:383`) shows a real append — `ALPHAMISSENSE`/`BAYESDEL` were added
after the original entries.

Then extend the matching `VcfProjectionSpec.description` — append `|NHOMALT` to the pipe list
after `Format: `. Keep the description's `Format:` tail and the doc's pipe line (step 3e)
character-identical.

The right-hand element must **byte-match** the producer's camelCase key. The generic
projection loop is `for (_, json_key) in spec.fields { obj.get(json_key) … }` (`:880`,
`:982`, `:1089`); a typo silently yields an empty field via `unwrap_or_default`, not an
error. No emitter code changes — VCF INFO and TSV both read the table positionally.

**Adding a whole NEW source.** In `output.rs`:

1. Define the fields table:
   ```rust
   const NEWSRC_FIELDS: &[(&str, &str)] = &[("FIELD_A", "fieldA"), ("FIELD_B", "fieldB")];
   ```
2. Append a spec to `VCF_PROJECTION_SPECS` (`output.rs:399`):
   ```rust
   VcfProjectionSpec {
       json_key: "<name>",                 // MUST equal the producer's json_key from step 3c
       info_id: "FV_<NAME>",
       description: "fastVEP <Name> annotations. Format: ALLELE|FIELD_A|FIELD_B",
       fields: NEWSRC_FIELDS,
       kind: VcfProjectionKind::AlleleObject,
   },
   ```
3. Pick `kind` (`VcfProjectionKind`, `output.rs:325`) by shape:
   - `AlleleObject` — per-allele JSON object with named subkeys (gnomAD, dbSNP, ClinVar).
   - `AlleleScalar` — per-allele bare number/string; use `SCORE_FIELDS = &[("SCORE","")]`
     (`:381`) with an **empty** `json_key` (PhyloP/GERP/DANN). Don't copy this for object
     sources.
   - `GeneObject` — gene-symbol keyed (OMIM, gnomAD constraint); looked up in `gene_keys`,
     not `sa_keys`.
   - `ClinvarProtein` — the special protein-variants collector.

`json_key` selects the key list: `GeneObject`/`ClinvarProtein` look in `gene_keys`;
`AlleleObject`/`AlleleScalar` look in `sa_keys` (`LoadedSupplementarySpecs::new`, `:603`). If
`json_key` isn't present in the loaded keys, the column/INFO silently never emits.

Note casing is per-registration, not uniformly snake_case: existing keys include `gnomad`,
`dbsnp`, `oneKg`, `primateAI`, `gnomad_genes`. Match whatever the producer literally sets.

### 3e. Document the fields (and satisfy the CI doc-coverage check)

**File:** `docs/SUPPLEMENTARY_ANNOTATIONS.md`

A unit test — `supplementary_annotations_doc_lists_every_projection_spec` at `output.rs:2249`
— `include_str!`s this doc and, for **every** `VCF_PROJECTION_SPECS` entry, asserts both:
- `doc.contains(spec.info_id)`, and
- `doc.contains(<the exact pipe layout after "Format: ">)` — a **verbatim substring** match.

It runs under `cargo test --workspace --lib` (the CI command, `.github/workflows/ci.yml:21`).
It is a unit test, not a `build.rs` step or a YAML grep. A one-character difference fails CI
with, e.g.:

> `docs/SUPPLEMENTARY_ANNOTATIONS.md is missing pipe format for FV_DBSNP: expected ALLELE|ID|GLOBAL_MAF`

**Adding a field** — extend the source's per-source pipe line under `### Allele-level`
(gnomAD at `docs/SUPPLEMENTARY_ANNOTATIONS.md:84`, dbSNP at `:85`) so it matches the spec
description verbatim:

```markdown
- `FV_GNOMAD`: `ALLELE|ALL_AF|…|FILTER|NHOMALT` — …
```

**Adding a source** — add a row to the identifier table (`:17-38`) and a per-source bullet
under `### Allele-level` or `### Gene-level` (allele `:82-96`, gene `:100-104`) with the exact
`Format:` layout.

**JSON-only fields** are invisible to this test (they're not in a spec). Document them in
prose — the existing pattern is the "emitted in **JSON only**" notes on the gnomAD and ClinVar
lines. Nothing enforces this, so don't skip it.

**SpliceAI** is special-cased: not namespaced under `FV_*`, not in `VCF_PROJECTION_SPECS`, and
asserted with a separate hard-coded literal in the same test (`output.rs:2280`). Only relevant
if you touch SpliceAI-shaped output.

### 3f. Add tests

**Inline parser unit tests (these run in CI).** Add a `#[cfg(test)] mod tests` in your new
source file, modeled on `test_parse_dbsnp_vcf` (`dbsnp.rs:117`) or
`test_gnomad_acmg_fields_joint` (`gnomad.rs:449`). Feed a tiny in-memory VCF/TSV string and
assert the JSON contains the expected camelCase keys. Match the exact `{:.6e}` rendering:

```rust
assert!(j.contains("\"faf95\":1.000000e-4"));   // gnomad.rs:472
```

Also cover the RefSeq-accession contig path if your source is NCBI-flavored (see
`test_parse_dbsnp_refseq_accessions`, `dbsnp.rs:141`) — this is exactly the case that produced
the historical zero-record bug.

**Test enumerations (new source only).** For a new source, update the helper lists in
`output.rs` so the column-order test stays complete:
- Add your `json_key` to `all_supplementary_sa_keys` (`output.rs:1922`) or
  `all_supplementary_gene_keys` (`:1941`).
- Add your `FV_<NAME>` to the expected column list inside
  `tab_supplementary_column_names_match_vcf_header_order` (`output.rs:1946`).

**Integration tests do NOT run in CI.** CI runs `--lib` only (`ci.yml:21`), so tests under
`crates/*/tests` (e.g. `crates/fastvep-io/tests/vep_compat.rs`,
`crates/fastvep-sa/tests/gnomad_indel_normalization.rs`) execute only locally with
`cargo test --workspace`. Anything you need CI to enforce must be an inline `#[cfg(test)]`
unit test.

### 3g. Build the SA and verify annotation picks it up

Build tooling: `cargo` lives at `~/.cargo/bin` (not on the default PATH).

```bash
# 1. Recommended locally (hygiene). NOTE: CI itself runs only
#    'cargo check --workspace --lib --bins' and 'cargo test --workspace --lib' —
#    it does NOT run fmt --check or clippy, but keep them clean anyway.
cargo fmt
cargo clippy --workspace
cargo test --workspace --lib      # runs the doc-coverage + column-order tests (the CI gate)

# 2. Build the binary
cargo build --release

# 3. Build the SA index from data
fastvep sa-build \
  --source <name> \
  --input <data.grch38.tsv.gz> \
  --output <dir>/<name> \
  --assembly GRCh38
# → writes <dir>/<name>.osa and <dir>/<name>.osa.idx
#   (or <name>.osi for a BED-derived source; <name>.oga for a gene source)
```

Confirm the `Parsed N records` line is **> 0**. Zero records means the file's contig naming
didn't match the assembly's chrom map (usually a RefSeq-accession normalization miss) — the
annotation column will be silently empty.

```bash
# 4. Annotate — drop the .osa (+ .osa.idx sidecar) in a dir and point --sa-dir at it.
#    Discovery is automatic by file extension; no per-source flag or reader change.
fastvep annotate \
  -i in.vcf \
  -o out.json \
  --gff3 gencode.gff3 \
  --fasta ref.fa \
  --sa-dir <dir>
```

At annotate time, discovery opens every `*.osa` / `*.osi` / `*.osa2` / `*.oga` in `--sa-dir`
by extension (`crates/fastvep-annotate/src/lib.rs:1200-1272`), opens each as an
`AnnotationProvider`, and the per-variant loop calls `annotate_position(chrom, pos, ref, alt)`
(`:646`). The result is attached to each allele under your `json_key`. Verify:

- **JSON output** contains your object under `json_key` with all emitted camelCase keys.
- **VCF output** contains `FV_<NAME>=…` with the projected fields in `*_FIELDS` order.
- **TSV output** has the `FV_<NAME>` column with the same pipe layout.

> A missing `.osa.idx` sidecar, or a wrong extension, makes the file skipped with only a
> `tracing::warn` — not an error. If your column is empty, check the sidecar and the extension
> first, then the `Parsed N records` count.

---

## 4. Conventions & invariants

- **Append-only, no reordering (load-bearing).** VCF `FV_*` INFO values and TSV columns are
  pipe-delimited **positional** strings with no per-field key. Downstream consumers parse by
  ordinal position against the documented `Format:` header. Inserting or reordering a
  `*_FIELDS` tuple silently shifts the meaning of every later field for already-emitted data.
  Always append at the end of the table **and** at the end of the `Format:` description. The
  `Appended (never reordered)` comments in `CLINVAR_FIELDS` (`output.rs:339`), `GNOMAD_FIELDS`
  (`:361`), and `DBNSFP_FIELDS` (`:386`) codify this.
- **Naming.** JSON side is **camelCase** (`allAf`, `grpmaxAf`, `faf95`, `pLI`). VCF/TSV
  positional labels are **UPPER_CASE with underscores** (`ALL_AF`, `GRPMAX_AF`, `FAF95`,
  `PLI`). The tuple's right element must byte-match the producer's camelCase key exactly.
- **JSON is a superset of VCF/TSV.** JSON output carries every producer key verbatim; VCF/TSV
  carry only the subset in `*_FIELDS`. Deliberately-omitted keys are "JSON-only by design"
  and should still be mentioned in the doc prose.
- **Record kind is destiny.** `AnnotationRecord → .osa`, `GeneRecord → .oga`,
  `IntervalRecord → .osi`. It picks the writer, reader, CLI branch, and projection kind. Gene
  sources dispatch through `run_oga_build` and the early guard, not the `.osa` match arms.
- **`json_key` is the single source of truth** tying the parser, the CLI header, the
  projection spec, and the doc together. It also drives special runtime handling (e.g. the
  `gnomad` key triggers allele normalization). Pick it once, use it identically everywhere.
- **Contig normalization** must respect RefSeq accessions (`NC_000001.11`) — never blindly
  prepend `chr`. Copy `dbsnp.rs`'s `normalize_chrom`.
- **Sort allele/interval records** by `(chrom_idx, position)` before returning; gene records
  stay in first-seen order.
- **Determinism.** Emit keys from a fixed slice or a `BTreeMap`, not from `HashMap` iteration
  order, or bytes vary run-to-run.
- **VCF subfield escaping** goes through `json_value_to_vcf` / `escape_vcf_subfield`
  (`output.rs:1101` / `:1125`), which `%`-encode `: ; = % , | &` and whitespace. Don't
  hand-format values that bypass this.
- **CI runs `--lib` only.** Put anything CI must enforce in inline `#[cfg(test)]` unit tests,
  not integration tests.
- **Public-repo hygiene.** This repository is public. Keep the source name, `json_key`, commit
  messages, comments, docs, and branch names free of internal references (ticket IDs,
  downstream-consumer names, roadmap notes). Only add commercially-licensable predictor
  sources.

---

## 5. PR checklist

Copy-paste and tick each item.

```
Record kind
- [ ] Chose the correct record kind (.osa allele/positional / .oga gene / .osi interval)
- [ ] If BED regions with no bespoke parsing: used --source custom_bed instead of a new parser

Parser  (crates/fastvep-sa/src/sources/<name>.rs)
- [ ] parse_<name> signature matches the record kind (chrom_map present for allele/interval, absent for gene)
- [ ] JSON built by hand with camelCase keys; floats {:.6e}; strings escaped AND quoted
- [ ] Contig normalization handles RefSeq accessions (copied from dbsnp.rs), NOT a naive chr-prefix
- [ ] issue #51 comment dropped/re-targeted if normalize_chrom was copied
- [ ] Multi-allelic INFO (Number=A) split and indexed per-ALT; Number=1 reused
- [ ] Records dropped when chrom absent from chrom_map is understood/expected
- [ ] Allele/interval records sorted by (chrom_idx, position); gene records first-seen
- [ ] Key emission is deterministic (fixed slice / BTreeMap, not HashMap iteration)

Registration  (new source only)
- [ ] pub mod <name>; added to crates/fastvep-sa/src/sources/mod.rs
- [ ] IndexHeader arm added in pipeline.rs run_sa_build (json_key, match_by_allele, is_positional, is_array)
      OR, for a gene source: name added to the run_oga_build guard + its (json_key,name) header arm + parser arm
- [ ] Parser arm added in the `let records = match source` block
- [ ] <name> added to the "Unknown source" bail! supported list

Projection seam  (crates/fastvep-io/src/output.rs)
- [ ] Field(s) appended to <SRC>_FIELDS (APPEND-ONLY, no reorder)
- [ ] New source: VcfProjectionSpec added to VCF_PROJECTION_SPECS with correct kind
- [ ] json_key in the spec byte-matches the producer's json_key
- [ ] Right-hand camelCase subkeys byte-match the parser's emitted keys
- [ ] VcfProjectionSpec.description "Format:" tail updated with the appended label(s)

Docs  (docs/SUPPLEMENTARY_ANNOTATIONS.md)
- [ ] Per-source pipe line matches the spec description "Format:" tail VERBATIM
- [ ] New source: identifier-table row added
- [ ] JSON-only fields documented in prose

Tests
- [ ] Inline #[cfg(test)] parser tests assert expected camelCase keys (exact {:.6e} rendering)
- [ ] RefSeq-accession contig path covered if source is NCBI-flavored
- [ ] New source: json_key added to all_supplementary_sa_keys / all_supplementary_gene_keys
- [ ] New source: FV_<NAME> added to tab_supplementary_column_names_match_vcf_header_order

Build & verify
- [ ] cargo fmt && cargo clippy --workspace clean (local hygiene)
- [ ] cargo test --workspace --lib passes (doc-coverage + column-order tests included)
- [ ] fastvep sa-build … reports "Parsed N records" with N > 0
- [ ] .osa (+ .osa.idx sidecar) present in --sa-dir; correct extension
- [ ] fastvep annotate … shows the field in JSON, VCF (FV_<NAME>), and TSV

Hygiene
- [ ] Source name, json_key, commits, comments, docs, and branch name free of internal references
- [ ] Any new predictor source is commercially licensable
```
