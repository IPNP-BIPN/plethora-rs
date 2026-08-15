//! The 1000 Genomes `sequence.index`, as the pipeline reads it.
//!
//! Two upstream scripts parse this file, and they disagree about it:
//!
//! ```text
//! # download_sample.pl
//! push @read_counts, $attributes[25];
//!
//! # clean_files.pl
//! $expected_number_of_reads += $data[23];
//! ```
//!
//! Column 23 is `READ_COUNT` and column 25 is `ANALYSIS_GROUP`, so
//! `download_sample.pl` reads a study label where it means a count.
//! The array it fills is never used, which is why nothing ever broke.
//! [`Record::read_count`] is column 23.
//!
//! The columns are named rather than indexed at the call sites, because that
//! is the mistake to make here.

use std::io::BufRead;

/// The header the file carries, used to check the layout is the expected one.
///
/// The 2017 file starts with `FASTQ_FILE`, and the EBI edition names the same
/// column `FASTQ_ENA_PATH`; `preprocessing_1000genomes.R` renames one to the
/// other for exactly that reason.
pub const FIRST_COLUMN: &str = "FASTQ_FILE";

/// Column positions, zero-based, in the order the file lists them.
mod column {
    pub const FASTQ_FILE: usize = 0;
    pub const MD5: usize = 1;
    pub const RUN_ID: usize = 2;
    pub const STUDY_NAME: usize = 4;
    pub const CENTER_NAME: usize = 5;
    pub const SAMPLE_NAME: usize = 9;
    pub const POPULATION: usize = 10;
    pub const INSTRUMENT_PLATFORM: usize = 12;
    pub const INSTRUMENT_MODEL: usize = 13;
    pub const LIBRARY_NAME: usize = 14;
    pub const INSERT_SIZE: usize = 17;
    pub const LIBRARY_LAYOUT: usize = 18;
    pub const WITHDRAWN: usize = 20;
    /// The one the two upstream scripts disagree about.
    pub const READ_COUNT: usize = 23;
    pub const BASE_COUNT: usize = 24;
    pub const ANALYSIS_GROUP: usize = 25;
    /// How many columns a well-formed row has.
    pub const COUNT: usize = 26;
}

/// One row: a FASTQ file and what is known about it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// Where the file can be fetched from.
    pub fastq_file: String,
    /// The checksum `download_sample.pl` verifies against.
    pub md5: String,
    pub run_id: String,
    pub study_name: String,
    pub center_name: String,
    pub sample_name: String,
    pub population: String,
    pub instrument_platform: String,
    pub instrument_model: String,
    pub library_name: String,
    /// Blank or `NA` in the file for some rows, hence the option.
    pub insert_size: Option<u64>,
    pub library_layout: String,
    /// Non-zero for a withdrawn run, which the selection drops.
    pub withdrawn: bool,
    /// Column 23, the real read count.
    pub read_count: Option<u64>,
    pub base_count: Option<u64>,
    pub analysis_group: String,
}

impl Record {
    /// Parses one line, or `None` if it is the header or too short.
    #[must_use]
    pub fn parse(line: &str) -> Option<Self> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < column::COUNT || f[column::FASTQ_FILE] == FIRST_COLUMN {
            return None;
        }

        // "NA" and an empty field both mean absent here.
        let number = |i: usize| -> Option<u64> {
            let v = f[i].trim();
            if v.is_empty() || v == "NA" {
                None
            } else {
                v.parse().ok()
            }
        };
        let text = |i: usize| f[i].to_string();

        Some(Self {
            fastq_file: text(column::FASTQ_FILE),
            md5: text(column::MD5),
            run_id: text(column::RUN_ID),
            study_name: text(column::STUDY_NAME),
            center_name: text(column::CENTER_NAME),
            sample_name: text(column::SAMPLE_NAME),
            population: text(column::POPULATION),
            instrument_platform: text(column::INSTRUMENT_PLATFORM),
            instrument_model: text(column::INSTRUMENT_MODEL),
            library_name: text(column::LIBRARY_NAME),
            insert_size: number(column::INSERT_SIZE),
            library_layout: text(column::LIBRARY_LAYOUT),
            withdrawn: number(column::WITHDRAWN).is_some_and(|v| v != 0),
            read_count: number(column::READ_COUNT),
            base_count: number(column::BASE_COUNT),
            analysis_group: text(column::ANALYSIS_GROUP),
        })
    }

    /// The bare file name, which is what lands in `fastq/<sample>/`.
    ///
    /// `download_sample.pl` derives it with `s/^.+?\/([^\/]+)$/$1/`, so
    /// everything up to the last slash goes.
    #[must_use]
    pub fn file_name(&self) -> &str {
        self.fastq_file
            .rsplit('/')
            .next()
            .unwrap_or(&self.fastq_file)
    }

    /// Which mate of the pair this file holds, if it says.
    ///
    /// The scripts match on `_1(.filt)*.fastq.gz` and `_2(.filt)*.fastq.gz`,
    /// so a file named neither way belongs to no mate and is skipped.
    #[must_use]
    pub fn mate(&self) -> Option<u8> {
        let name = self.file_name();
        // `_1` before `_2`, as the Perl's if/elsif does, and every occurrence
        // is considered: the regex is unanchored, so a name carrying `_1`
        // twice matches on whichever is followed by the extension.
        for (marker, mate) in [("_1", 1_u8), ("_2", 2)] {
            let mut rest = name;
            while let Some(index) = rest.find(marker) {
                let after = &rest[index + marker.len()..];
                if after.starts_with(".fastq") || after.starts_with(".filt.fastq") {
                    return Some(mate);
                }
                rest = after;
            }
        }
        None
    }
}

/// Reads an index, keeping every row.
///
/// # Errors
/// Returns an error if the input cannot be read.
pub fn read<R: BufRead>(input: R) -> std::io::Result<Vec<Record>> {
    let mut out = Vec::new();
    for line in input.lines() {
        if let Some(record) = Record::parse(&line?) {
            out.push(record);
        }
    }
    Ok(out)
}

/// Every row belonging to one sample, in file order.
#[must_use]
pub fn for_sample<'a>(records: &'a [Record], sample: &str) -> Vec<&'a Record> {
    records.iter().filter(|r| r.sample_name == sample).collect()
}

/// The read count `clean_files.pl` compares a sample's files against.
///
/// Summed across every row for the sample, both mates, and **not** halved. The
/// chain compares it against the FASTQ record count over both files and against
/// the BAM record count, which both hold every read; only the BED is per
/// fragment, and [`crate::onekg::clean::plan`] halves it there. Upstream says
/// so where it does the same sum: "note both pairs will be present in alignment
/// file, so we need to count both".
///
/// `None` when the sample is absent or no row carries a count, which upstream
/// treats as "not in the manifest" and refuses to act on.
#[must_use]
pub fn expected_reads(records: &[Record], sample: &str) -> Option<u64> {
    let total: u64 = for_sample(records, sample)
        .iter()
        .filter_map(|r| r.read_count)
        .sum();
    (total > 0).then_some(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The count is over both mates and is not halved: it is checked against
    /// the FASTQ records and the BAM records, which both hold every read.
    ///
    /// Built by editing the real row, so the column positions stay the ones
    /// the file actually uses.
    #[test]
    fn the_expected_count_covers_both_mates() {
        let mate1 = ROW.replace("101519492", "10506");
        let mate2 = mate1.replace("ERR020229_1.fastq.gz", "ERR020229_2.fastq.gz");
        let index = format!("{HEADER}\n{mate1}\n{mate2}\n");

        let records = read(index.as_bytes()).expect("read");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].read_count, Some(10_506));
        assert_eq!(
            expected_reads(&records, "HG00108"),
            Some(21_012),
            "both mates, not halved"
        );
        assert_eq!(expected_reads(&records, "NA19323"), None, "absent sample");
    }

    const HEADER: &str = "FASTQ_FILE\tMD5\tRUN_ID\tSTUDY_ID\tSTUDY_NAME\tCENTER_NAME\tSUBMISSION_ID\tSUBMISSION_DATE\tSAMPLE_ID\tSAMPLE_NAME\tPOPULATION\tEXPERIMENT_ID\tINSTRUMENT_PLATFORM\tINSTRUMENT_MODEL\tLIBRARY_NAME\tRUN_NAME\tRUN_BLOCK_NAME\tINSERT_SIZE\tLIBRARY_LAYOUT\tPAIRED_FASTQ\tWITHDRAWN\tWITHDRAWN_DATE\tCOMMENT\tREAD_COUNT\tBASE_COUNT\tANALYSIS_GROUP";

    /// A row taken from the real file, so the column positions are checked
    /// against the layout rather than against my own header.
    const ROW: &str = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR020/ERR020229/ERR020229_1.fastq.gz\t85dc912bde2389c5341cc102fa0c0bc7\tERR020229\tSRP001294\t1000 Genomes GBR population sequencing\tBGI\tERA014686\t2010-09-26 00:00:00\tSRS006849\tHG00108\tGBR\tERX008348\tILLUMINA\tIllumina HiSeq 2000\tHUMgfdRAADIAAPEI-9\tBGI-A80940ABXX_L1_HUMgfdRAADIAAPEI-9\tNA\t466\tPAIRED\tftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR020/ERR020229/ERR020229_2.fastq.gz\t0\tNA\tNA\t101519492\t18476547544\tlow coverage";

    #[test]
    fn the_header_is_not_a_record() {
        assert!(Record::parse(HEADER).is_none());
    }

    #[test]
    fn a_short_line_is_not_a_record() {
        assert!(Record::parse("chr1\t1\t2").is_none());
    }

    /// The columns that matter, read off a real row.
    #[test]
    fn a_real_row_parses() {
        let r = Record::parse(ROW).unwrap();
        assert_eq!(r.sample_name, "HG00108");
        assert_eq!(r.population, "GBR");
        assert_eq!(r.center_name, "BGI");
        assert_eq!(r.instrument_model, "Illumina HiSeq 2000");
        assert_eq!(r.library_layout, "PAIRED");
        assert_eq!(r.insert_size, Some(466));
        assert!(!r.withdrawn);
        assert_eq!(r.md5, "85dc912bde2389c5341cc102fa0c0bc7");
    }

    /// The mistake the two upstream scripts make between them: column 23 is the
    /// count, column 25 is a label.
    #[test]
    fn the_read_count_is_column_23_not_25() {
        let r = Record::parse(ROW).unwrap();
        assert_eq!(r.read_count, Some(101_519_492), "column 23");
        assert_eq!(r.base_count, Some(18_476_547_544));
        assert_eq!(
            r.analysis_group, "low coverage",
            "column 25, which download_sample.pl reads as a read count"
        );
    }

    #[test]
    fn the_file_name_drops_everything_up_to_the_last_slash() {
        let r = Record::parse(ROW).unwrap();
        assert_eq!(r.file_name(), "ERR020229_1.fastq.gz");
    }

    #[test]
    fn the_mate_comes_from_the_file_name() {
        let r = Record::parse(ROW).unwrap();
        assert_eq!(r.mate(), Some(1));

        let two = ROW.replacen("_1.fastq.gz", "_2.fastq.gz", 1);
        assert_eq!(Record::parse(&two).unwrap().mate(), Some(2));

        // The .filt spelling the scripts also match on.
        let filt = ROW.replacen("_1.fastq.gz", "_1.filt.fastq.gz", 1);
        assert_eq!(Record::parse(&filt).unwrap().mate(), Some(1));

        // Single-end and oddly named files belong to no mate.
        let single = ROW.replacen("ERR020229_1.fastq.gz", "ERR020229.fastq.gz", 1);
        assert_eq!(Record::parse(&single).unwrap().mate(), None);
    }

    /// "NA" and an empty field both mean absent, and the file uses both.
    #[test]
    fn absent_numbers_read_as_none() {
        let na = ROW.replacen("\t466\t", "\tNA\t", 1);
        assert_eq!(Record::parse(&na).unwrap().insert_size, None);

        let empty = ROW.replacen("\t466\t", "\t\t", 1);
        assert_eq!(Record::parse(&empty).unwrap().insert_size, None);
    }

    /// Column 20 carries the withdrawal flag, and the selection drops anything
    /// non-zero.
    #[test]
    fn a_withdrawn_run_is_flagged() {
        assert!(!Record::parse(ROW).unwrap().withdrawn);

        let mut f: Vec<&str> = ROW.split('\t').collect();
        f[20] = "1";
        assert!(Record::parse(&f.join("\t")).unwrap().withdrawn);
    }

    #[test]
    fn reading_skips_the_header_and_groups_by_sample() {
        let text = format!(
            "{HEADER}\n{ROW}\n{}\n",
            ROW.replacen("HG00108", "HG00109", 1)
        );
        let records = read(text.as_bytes()).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(for_sample(&records, "HG00108").len(), 1);
        assert_eq!(for_sample(&records, "HG00109").len(), 1);
        assert_eq!(for_sample(&records, "nobody").len(), 0);
    }
}
