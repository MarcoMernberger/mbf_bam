use rust_htslib::bam::{Reader, Writer, Header, Read, Record, index};
use crate::{BamError};
use std::collections::HashSet;



/// substract all reads from the subtrahend bam file
/// from the minuend bam file,
/// writing to output bam file
/// and indexing it
pub fn py_substract_bam(
    output_filename: &str,
    minuend_filename: &str,
    subtrahend_filename: &str) -> Result<(), BamError>  {

    let mut minuend = Reader::from_path(minuend_filename)?;
    let mut subtrahend = Reader::from_path(subtrahend_filename)?;
    let header = Header::from_template(minuend.header());


    let mut seen = HashSet::new();

    let mut read: Record = Record::new();
    while let Ok(_) = subtrahend.read(&mut read) {
        let q = read.qname().to_owned();
        seen.insert(q);
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


