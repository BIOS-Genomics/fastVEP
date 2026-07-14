//! gnomAD v4.1 all-sites Allele Number (AN) — per-locus callability.
//!
//! Parses the gnomAD v4.1 `allele_number_all_sites` TSVs (one per callset:
//! exomes and genomes). Each row is a locus (`chr:pos`) with the AN observed
//! across all samples at that site, plus region flags.
//!
//! These are **positional** records (`match_by_allele = false`,
//! `is_positional = true`): AN describes how well a *site* was called,
//! independent of the specific alt allele. That lets a downstream consumer tell
//! a variant at a well-covered locus (real absence) apart from one at a locus
//! no sample could be called at — a distinction the per-variant sites data
//! cannot make, because an unobserved allele simply has no row there.
//!
//! **Distinct from the joint-VCF AN.** This is a *different* quantity than the
//! combined `gnomad.allAn` that the `gnomad` (joint sites VCF) source emits: it
//! is per-callset callability at a *position*, not the AN of an observed
//! variant. It lives under its own source keys (`gnomad_an_exomes` /
//! `gnomad_an_genomes`) and is emitted as `allSitesAn`, so it never overrides
//! `allAn`.
//!
//! **The genome-scale files stream** — the genomes callset has a row per
//! callable base across essentially the whole genome (100M+ loci), so the
//! pipeline drives this through `iter_gnomad_an` + `run_streaming_sa_build`
//! (like gnomAD/dbSNP/TOPMed), never buffering every record. `parse_gnomad_an`
//! (buffer-then-sort) is retained for tests and small inputs. The input must be
//! coordinate-sorted (all gnomAD releases are); the streaming writer bails on
//! out-of-order input rather than sorting.
//!
//! **Rows with `AN == 0` are dropped at build time.** `AN` counts *called*
//! genotypes (reference + alternate), so `AN == 0` means no sample had a
//! confident call at that locus — the site was not covered, which is the same
//! signal as a locus simply absent from the index. Dropping the zeros changes
//! no downstream decision (a missing record already means "not covered") while
//! trimming the index to the genuinely-covered subset — a large saving, since
//! most of the genome sits outside the exome calling intervals (`AN == 0`).
//! Only `AN == 0` is dropped: a low-but-nonzero `AN` is kept so the
//! "well-covered" cutoff stays a *consumer-side* calibration, not a threshold
//! baked into the annotation. The drop is a build-time choice and reversible.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Column indices resolved from the TSV header, so column order isn't hard-coded.
struct AnColumns {
    locus: usize,
    an: usize,
    outside_broad_capture: Option<usize>,
    outside_ukb_capture: Option<usize>,
    outside_broad_calling: Option<usize>,
    outside_ukb_calling: Option<usize>,
}

impl AnColumns {
    fn from_header(header: &str) -> Result<Self> {
        let fields: Vec<&str> = header.split('\t').collect();
        let find = |name: &str| {
            fields
                .iter()
                .position(|f| f.trim().eq_ignore_ascii_case(name))
        };
        Ok(Self {
            locus: find("locus").context("all-sites-AN TSV: missing 'locus' column")?,
            an: find("AN").context("all-sites-AN TSV: missing 'AN' column")?,
            outside_broad_capture: find("outside_broad_capture_region"),
            outside_ukb_capture: find("outside_ukb_capture_region"),
            outside_broad_calling: find("outside_broad_calling_region"),
            outside_ukb_calling: find("outside_ukb_calling_region"),
        })
    }
}

/// A region flag → `Some(true)`/`Some(false)` only for an explicit `true`/`false`
/// value; anything else (blank, `NA`, `.`, a short row) is `None` so we omit the
/// key rather than assert a wrong "inside-region" claim.
fn bool_field(fields: &[&str], idx: Option<usize>) -> Option<bool> {
    let v = idx.and_then(|i| fields.get(i))?.trim();
    if v.eq_ignore_ascii_case("true") {
        Some(true)
    } else if v.eq_ignore_ascii_case("false") {
        Some(false)
    } else {
        None
    }
}

/// Contig normalization: accept a bare `1` as well as `chr1` (every sibling
/// parser does this). The gnomAD loci are already `chr`-prefixed, so this is a
/// no-op for them but keeps the lookup robust.
fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

/// Build the JSON object for one covered locus: `{"allSitesAn": <u64> [, flags]}`.
fn record_json(an: u64, fields: &[&str], cols: &AnColumns) -> String {
    let mut parts = vec![format!("\"allSitesAn\":{}", an)];
    for (key, idx) in [
        ("outsideBroadCaptureRegion", cols.outside_broad_capture),
        ("outsideUkbCaptureRegion", cols.outside_ukb_capture),
        ("outsideBroadCallingRegion", cols.outside_broad_calling),
        ("outsideUkbCallingRegion", cols.outside_ukb_calling),
    ] {
        if let Some(b) = bool_field(fields, idx) {
            parts.push(format!("\"{}\":{}", key, b));
        }
    }
    format!("{{{}}}", parts.join(","))
}

/// Stream a coordinate-sorted gnomAD all-sites-AN TSV as positional
/// `AnnotationRecord`s without buffering the whole file in memory.
///
/// The input must already be sorted by `(chrom, pos)` — all gnomAD releases are;
/// the streaming writer bails on out-of-order input. `AN == 0` rows, unknown
/// contigs, and unparseable rows are skipped; the first non-comment line is the
/// header.
pub fn iter_gnomad_an<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> GnomadAnRecordIter<'_, R> {
    GnomadAnRecordIter {
        lines: reader.lines(),
        chrom_to_idx,
        cols: None,
    }
}

pub struct GnomadAnRecordIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    cols: Option<AnColumns>,
}

impl<R: BufRead> Iterator for GnomadAnRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e).context("Reading all-sites-AN TSV line")),
            };
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            // First non-comment line is the header.
            if self.cols.is_none() {
                match AnColumns::from_header(&line) {
                    Ok(c) => {
                        self.cols = Some(c);
                        continue;
                    }
                    Err(e) => return Some(Err(e)),
                }
            }
            let cols = self.cols.as_ref().expect("header parsed above");

            let fields: Vec<&str> = line.split('\t').collect();
            let locus = match fields.get(cols.locus) {
                Some(l) => *l,
                None => continue,
            };
            let (chrom_str, pos_str) = match locus.split_once(':') {
                Some(cp) => cp,
                None => continue,
            };
            let chrom_idx = match self.chrom_to_idx.get(&normalize_chrom(chrom_str)) {
                Some(&idx) => idx,
                None => continue,
            };
            let pos: u32 = match pos_str.trim().parse() {
                Ok(p) => p,
                Err(_) => continue,
            };
            let an: u64 = match fields.get(cols.an).and_then(|v| v.trim().parse().ok()) {
                Some(a) => a,
                None => continue,
            };
            if an == 0 {
                continue; // uncalled locus — an absent record already means "not covered"
            }

            return Some(Ok(AnnotationRecord {
                chrom_idx,
                position: pos,
                ref_allele: String::new(),
                alt_allele: String::new(),
                json: record_json(an, &fields, cols),
            }));
        }
    }
}

/// Buffered variant — collects every record and sorts in memory.
///
/// Retained for tests and small inputs; the pipeline uses [`iter_gnomad_an`] via
/// `run_streaming_sa_build` for the genome-scale build (see the module docs).
pub fn parse_gnomad_an<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<AnnotationRecord> =
        iter_gnomad_an(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_chrom_map() -> HashMap<String, u16> {
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0);
        m.insert("chr20".into(), 1);
        m
    }

    #[test]
    fn parses_an_and_flags_skips_zero_and_unknown_contig() {
        // Header order matches the real gnomAD exomes all-sites-AN TSV.
        let data = "\
locus\tAN\toutside_broad_capture_region\toutside_ukb_capture_region\toutside_broad_calling_region\toutside_ukb_calling_region
chr1:11719\t0\ttrue\ttrue\tfalse\ttrue
chr20:22584230\t152092\tfalse\tfalse\tfalse\tfalse
chr20:22584231\t730948\ttrue\tfalse\tfalse\tfalse
chr7:100\t500\tfalse\tfalse\tfalse\tfalse
";
        let records = parse_gnomad_an(data.as_bytes(), &test_chrom_map()).unwrap();
        // chr1:11719 dropped (AN==0); chr7 dropped (contig not in map); 2 remain.
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom_idx, 1);
        assert_eq!(records[0].position, 22584230);
        assert_eq!(records[0].ref_allele, "");
        assert_eq!(records[0].alt_allele, "");
        assert_eq!(
            records[0].json,
            "{\"allSitesAn\":152092,\"outsideBroadCaptureRegion\":false,\"outsideUkbCaptureRegion\":false,\"outsideBroadCallingRegion\":false,\"outsideUkbCallingRegion\":false}"
        );
        assert_eq!(records[1].position, 22584231);
        assert!(records[1].json.contains("\"allSitesAn\":730948"));
        assert!(records[1].json.contains("\"outsideBroadCaptureRegion\":true"));
    }

    #[test]
    fn genomes_format_without_capture_flags() {
        // Genomes file may carry only calling-region flags; header detection copes.
        let data = "\
locus\tAN\toutside_broad_calling_region\toutside_ukb_calling_region
chr1:200\t151000\tfalse\tfalse
";
        let records = parse_gnomad_an(data.as_bytes(), &test_chrom_map()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].position, 200);
        assert_eq!(
            records[0].json,
            "{\"allSitesAn\":151000,\"outsideBroadCallingRegion\":false,\"outsideUkbCallingRegion\":false}"
        );
    }

    #[test]
    fn unknown_or_blank_region_flag_is_omitted_not_coerced_false() {
        // A non-true/false flag value (NA / blank) must be omitted, not called false.
        let data = "\
locus\tAN\toutside_broad_capture_region
chr1:300\t900\tNA
";
        let records = parse_gnomad_an(data.as_bytes(), &test_chrom_map()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].json, "{\"allSitesAn\":900}"); // flag omitted, not false
    }

    #[test]
    fn iter_is_lazy_and_matches_buffered() {
        let data = "\
locus\tAN\toutside_broad_calling_region
chr1:10\t5\tfalse
chr20:20\t7\ttrue
";
        let via_iter: Vec<_> = iter_gnomad_an(data.as_bytes(), &test_chrom_map())
            .collect::<Result<_>>()
            .unwrap();
        assert_eq!(via_iter.len(), 2);
        assert_eq!(via_iter[0].position, 10);
        assert_eq!(via_iter[1].chrom_idx, 1);
    }
}
