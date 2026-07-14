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
//! **Rows with `AN == 0` are dropped at build time.** `AN` counts *called*
//! genotypes (reference + alternate), so `AN == 0` means no sample had a
//! confident call at that locus — the site was not covered, which is the same
//! signal as a locus simply absent from the index. Dropping the zeros therefore
//! changes no downstream decision (a missing record already means "not covered")
//! while trimming the index to the genuinely-covered subset — a large saving,
//! since most of the genome sits outside the exome calling intervals (`AN == 0`).
//! Only `AN == 0` is dropped: a low-but-nonzero `AN` is kept so the
//! "well-covered" cutoff stays a *consumer-side* calibration, not a threshold
//! baked into the annotation. The drop is a build-time choice and reversible.
//!
//! JSON per record: `{"an": <u64> [, "outside…Region": <bool>]}`. The AN is the
//! load-bearing field; the region flags are carried for audit (a low AN inside
//! an out-of-capture region is expected, not a coverage concern).

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

fn bool_field(fields: &[&str], idx: Option<usize>) -> Option<bool> {
    idx.and_then(|i| fields.get(i))
        .map(|v| v.trim().eq_ignore_ascii_case("true"))
}

/// Parse a gnomAD all-sites-AN TSV into positional `AnnotationRecord`s.
///
/// Locus format is `chr1:11719` (1-based, `chr`-prefixed). The first non-comment
/// line is the header. Rows with `AN == 0`, an unknown contig, or an unparseable
/// locus/AN are skipped. Records are sorted by `(chrom_idx, position)` as the
/// writer requires.
pub fn parse_gnomad_an<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();
    let mut cols: Option<AnColumns> = None;

    for line in reader.lines() {
        let line = line.context("Reading all-sites-AN TSV line")?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // First non-comment line is the header.
        if cols.is_none() {
            cols = Some(AnColumns::from_header(&line)?);
            continue;
        }
        let cols = cols.as_ref().expect("header parsed above");

        let fields: Vec<&str> = line.split('\t').collect();
        let locus = match fields.get(cols.locus) {
            Some(l) => *l,
            None => continue,
        };
        let (chrom_str, pos_str) = match locus.split_once(':') {
            Some(cp) => cp,
            None => continue,
        };
        let chrom_idx = match chrom_to_idx.get(chrom_str) {
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

        let mut parts = vec![format!("\"an\":{}", an)];
        for (key, idx) in [
            ("outsideBroadCaptureRegion", cols.outside_broad_capture),
            ("outsideUkbCaptureRegion", cols.outside_ukb_capture),
            ("outsideBroadCallingRegion", cols.outside_broad_calling),
            ("outsideUkbCallingRegion", cols.outside_ukb_calling),
        ] {
            if let Some(b) = bool_field(&fields, idx) {
                parts.push(format!("\"{}\":{}", key, b));
            }
        }
        let json = format!("{{{}}}", parts.join(","));

        records.push(AnnotationRecord {
            chrom_idx,
            position: pos,
            ref_allele: String::new(),
            alt_allele: String::new(),
            json,
        });
    }

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
            "{\"an\":152092,\"outsideBroadCaptureRegion\":false,\"outsideUkbCaptureRegion\":false,\"outsideBroadCallingRegion\":false,\"outsideUkbCallingRegion\":false}"
        );
        // Records are sorted by (chrom_idx, position).
        assert_eq!(records[1].position, 22584231);
        assert!(records[1].json.contains("\"an\":730948"));
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
            "{\"an\":151000,\"outsideBroadCallingRegion\":false,\"outsideUkbCallingRegion\":false}"
        );
    }
}
