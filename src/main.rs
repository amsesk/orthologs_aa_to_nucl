use bio::io::gff;
use clap::{App, Arg};
use csv;
use polars_core::prelude::*;
use polars_io::prelude::*;
use serde::Deserialize;
use std::fs::File;
use std::io;
use std::path::Path;
use std::process;

#[derive(Debug, Deserialize)]
pub struct GffFnaPair {
    gff_path: Box<Path>,
    fna_path: Box<Path>,
}

fn read_ortholog_matrix(path: &str) -> Result<DataFrame> {
    CsvReader::from_path(path)?
        .has_header(true)
        .with_delimiter(b'\t')
        .finish()
}

fn main() {
    let args = App::new("orthologs_aa_to_nucl")
        .version("0.1")
        .author("Kevin Amses")
        .about("Parse a table of protein ortholog names and get their corresponding nucleotide sequences.")
        .arg(
            Arg::new("orthotable")
                .short('o')
                .long("ortholog_table")
                .takes_value(true)
                .required(true),
        ).get_matches();

    let orthotable_path = args.value_of("orthotable").unwrap();

    let mut orthotable =
        read_ortholog_matrix(orthotable_path).expect("Unable to read csv to DataFrame.");

    //println!("{:?}", orthotable);

    let mut csvrdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(io::stdin());
    let mut csvrdr_iter = csvrdr.deserialize::<GffFnaPair>();
    while let Some(record) = csvrdr_iter.next() {
        let pair: GffFnaPair = record.expect("Unable to parse line into GffFnaPair.");

        let mut gff_reader =
            gff::Reader::new(File::open(pair.gff_path).unwrap(), gff::GffType::GFF3);

        /*
        for r in gff_reader.records() {
            println!("{:?}", r)
        }
        */
        /*)
        for i in 1..ncol.1 {
            for row_spl in orthotable
                .select_at_idx(i)
                .expect("Select idx does not exist in DataFrame.")
                .iter()
                .map(|x| x.to_string())
                .map(|x| x.split('|').collect::<Vec<&str>>())
            {
                match row_spl.len() {
                    1 => println!(">>>{:?}", row_spl),
                    2 => println!("{:?}", row_spl),
                    _ => println!("Bah"),
                }
            }
        }
        */
    }
}
