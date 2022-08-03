use calamine::{open_workbook, Xlsx};
use clap::{App, Arg};
use csv;
use serde::Deserialize;
use std::io;
use std::path::Path;
use std::process;

#[derive(Debug, Deserialize)]
pub struct GffFnaPair {
    gff_path: Box<Path>,
    fna_path: Box<Path>,
}

fn read_ortholog_matrix(path: String) {
    let workbook: Xlsx<_> = open_workbook(path).expect("Cannot open workbook");
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

    let mut csvrdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(io::stdin());
    let mut csvrdr_iter = csvrdr.deserialize::<GffFnaPair>();
    while let Some(record) = csvrdr_iter.next() {
        let record: GffFnaPair = record;
    }
    /*
    for result in csvrdr.deserialize() {
        let record: GffFnaPair = match result {
            Ok(r) => r,
            Err(e) => {
                println!("Error reading pairs of gff and fna: {}", e);
                process::exit(1);
            }
        };
        println!("{:?}", record);
    }
    */
}
