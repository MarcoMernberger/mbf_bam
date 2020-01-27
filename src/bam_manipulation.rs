use crate::BamError;
use bio::io::fastq;
use flate2::read::GzDecoder;
use rust_htslib::bam::{index, record::Aux, Header, Read, Reader, Record, Writer};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufReader;
use std::str;

/// substract all reads from the subtrahend bam file
/// from the minuend bam file,
/// writing to output bam file
/// and indexing it
pub fn py_substract_bam(
    output_filename: &str,
    minuend_filename: &str,
    subtrahend_filename: &str,
) -> Result<(), BamError> {
    let mut minuend = Reader::from_path(minuend_filename)?;
    let mut subtrahend = Reader::from_path(subtrahend_filename)?;
    let header = Header::from_template(minuend.header());

    let mut seen = HashSet::new();

    let mut read: Record = Record::new();
    while let Ok(_) = subtrahend.read(&mut read) {
        if !read.is_unmapped() {
            let q = read.qname().to_owned();
            seen.insert(q);
        }
    }

    {
        let mut output = Writer::from_path(output_filename, &header)?;
        while let Ok(_) = minuend.read(&mut read) {
            if !seen.contains(read.qname()) {
                output.write(&read)?;
            }
        }
    } // output.drop will be called
    index::build(output_filename, None, index::Type::BAI, 4)?; //I see four threads
    Ok(())
}

fn read_gz_or_not(input_filename: &str) -> Result<Box<std::io::Read>, BamError> {
    let file = File::open(input_filename)?;
    if input_filename.ends_with(".gz") {
        return Ok(Box::new(GzDecoder::new(file)));
    } else {
        return Ok(Box::new(file));
    }
}

pub fn py_annotate_barcodes_from_fastq(
    output_filename: &str,
    input_filename: &str,
    fastq2_filenames: Vec<&str>,
    barcodes: Vec<(String, usize, usize)>,
) -> Result<(), BamError> {
    let mut input = Reader::from_path(input_filename)?;
    let header = Header::from_template(input.header());

    let mut qname_to_tags: HashMap<Vec<u8>, Vec<(Vec<u8>, Vec<u8>)>> = HashMap::new();
    for filename in fastq2_filenames {
        let buf = BufReader::new(read_gz_or_not(filename)?);
        let reader = fastq::Reader::new(buf);
        for record in reader.records() {
            let record = record?;
            let mut tags: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
            for (tag, start, end) in barcodes.iter() {
                let barcode = (record.seq()[*start..*end]).to_vec();
                tags.push((tag.as_bytes().to_vec(), barcode));
            }
            qname_to_tags.insert(record.id().as_bytes().to_vec(), tags);
        }
    }
    {
        let mut output = Writer::from_path(output_filename, &header)?;
        let mut read: Record = Record::new();
        while let Ok(_) = input.read(&mut read) {
            let tags = qname_to_tags
                .get(read.qname())
                .ok_or_else(|| BamError::UnknownError {
                    msg: format!(
                        "Read not found in fastq: {}",
                        std::str::from_utf8(read.qname()).unwrap()
                    )
                    .to_string(),
                })?;
            for (tag, value) in tags {
                read.push_aux(tag, &Aux::String(value))?;
            }
            output.write(&read)?;
        }
    } // output.drop will be called
    index::build(output_filename, None, index::Type::BAI, 4)?; //I see four threads
    Ok(())
}
