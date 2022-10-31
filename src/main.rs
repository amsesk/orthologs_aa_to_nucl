use bio::alphabets::dna::revcomp;
use bio::io::fasta::IndexedReader;
use bio::io::gff;
use bio::utils::Text;
use clap::{App, Arg};
use csv;
use polars::prelude::col;
use polars::prelude::*;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::fs::OpenOptions;
use std::io;
use std::io::prelude::*;
use std::path::Path;

#[derive(Debug, Deserialize)]
pub struct GffFnaPair {
    gff_path: Box<Path>,
    fna_path: Box<Path>,
}

pub struct OutputFileHandles {
    fasta: File,
    metadata: File,
}
impl OutputFileHandles {
    fn new(prefix: &str) -> OutputFileHandles {
        OutputFileHandles {
            fasta: OpenOptions::new()
                .write(true)
                .append(true)
                .create(true)
                .open(format!("{}.fasta", &prefix))
                .expect(&format!(
                    "Unable to create output file at {}.fasta",
                    &prefix
                )),
            metadata: OpenOptions::new()
                .write(true)
                .append(true)
                .create(true)
                .open(format!("{}.metadata.tsv", &prefix))
                .expect(&format!(
                    "Unable to create output file at {}.metadata.tsv",
                    &prefix
                )),
        }
    }
}

fn read_ortholog_matrix(path: &str) -> Result<DataFrame, PolarsError> {
    CsvReader::from_path(path)
        .unwrap()
        .has_header(true)
        .with_delimiter(b'\t')
        .finish()
}

// Create marker list from orthotable DataFrame
// Create a hash map with the marker name as key
// and a file handle as value
fn generate_file_handles(orthotable: DataFrame) -> HashMap<String, OutputFileHandles> {
    let mut per_marker_file_handles: HashMap<String, OutputFileHandles> = HashMap::new();

    orthotable
        .column("marker")
        .unwrap()
        .utf8()
        .unwrap()
        .into_iter()
        .map(|m| m.unwrap().to_owned())
        .for_each(|m| {
            let fhs = OutputFileHandles::new(&m);
            per_marker_file_handles.insert(m, fhs);
        });

    per_marker_file_handles
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

    let orthotable =
        read_ortholog_matrix(orthotable_path).expect("Unable to read csv to DataFrame.");

    let mut per_marker_file_handles = generate_file_handles(orthotable.clone());

    let mut csvrdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(io::stdin());
    let mut csvrdr_iter = csvrdr.deserialize::<GffFnaPair>();
    while let Some(record) = csvrdr_iter.next() {
        let pair: GffFnaPair = record.expect("Unable to parse line into GffFnaPair.");

        let mut gff_reader = gff::Reader::new(
            File::open(&pair.gff_path).expect(&format!(
                "Unable to open input GFF file at {}",
                &pair.gff_path.to_str().unwrap()
            )),
            gff::GffType::GFF3,
        );

        let mut fasta_reader = IndexedReader::from_file(&pair.fna_path).expect(&format!(
            "Unable to open input FASTA file at {}",
            &pair.fna_path.to_str().unwrap()
        ));

        let gff_fname_str = pair
            .gff_path
            .file_name()
            .expect("Unable to pull filename from GFF Path.")
            .to_str()
            .expect("Unable to convert GFF Path to str.");

        let ltp = match &gff_fname_str[gff_fname_str.len() - 4..] {
            ".gff" => Some(gff_fname_str[..&gff_fname_str.len() - 4].to_owned()),
            "gff3" => Some(gff_fname_str[..&gff_fname_str.len() - 5].to_owned()),
            _ => None,
        };

        let ltp = ltp.expect("Unable to convert GFF filename into corresponding Locus Tag Prefix.");

        let orthotable_this_ltp = orthotable
            .select(["marker", &ltp])
            .expect("Unable to find column matching LTP in Ortholog DataFrame.");

        for r in gff_reader.records() {
            match r {
                Ok(r) => {
                    if r.attributes().contains_key("ID") {
                        let prot_id = r.attributes()["ID"].split('|').collect::<Vec<&str>>()[0];
                        let search_pat = format!("^{}[|]{}$", &ltp, &prot_id);
                        let search = orthotable_this_ltp
                            .clone()
                            .lazy()
                            .filter(col(&ltp).str().contains(&search_pat))
                            .collect()
                            .unwrap();
                        match search.shape().0 {
                            0 => (),
                            1 => {
                                // Safe to index 0 and unwrap here because we have just matched the shape
                                // of search to be exactly 1 row
                                let new_match = search.get(0).unwrap();

                                if let (&AnyValue::Utf8(m), &AnyValue::Utf8(p)) =
                                    (&new_match[0], &new_match[1])
                                {
                                    let mut fetched_seq = Text::new();
                                    fasta_reader
                                        .fetch(&r.seqname(), *r.start(), *r.end())
                                        .expect("Unable to find sequence.");
                                    fasta_reader
                                        .read(&mut fetched_seq)
                                        .expect("Unable to read fetched sequence into vector.");

                                    if r.strand().unwrap().strand_symbol() == "-" {
                                        fetched_seq = revcomp(fetched_seq);
                                    }

                                    let fasta_header =
                                        format!("{}|{}::{}::{}-{}", ltp, m, p, r.start(), r.end());

                                    let _wt_fasta =
                                        &per_marker_file_handles.get_mut(m).unwrap().fasta.write(
                                            format!(
                                                ">{}\n{}\n",
                                                fasta_header,
                                                String::from_utf8(fetched_seq).unwrap()
                                            )
                                            .as_bytes(),
                                        );

                                    let _wt_meta = &per_marker_file_handles
                                        .get_mut(m)
                                        .unwrap()
                                        .metadata
                                        .write(
                                            format!(
                                                "{}\t{}\t{}\t{}\t{}\t{}\n",
                                                &m,
                                                &p,
                                                r.start(),
                                                r.end(),
                                                r.strand().unwrap().strand_symbol().to_owned(),
                                                r.frame()
                                            )
                                            .as_bytes(),
                                        );

                                    /*
                                    matches.push((
                                        m.to_owned(),
                                        p.to_owned(),
                                        *r.start(),
                                        *r.end(),
                                        r.strand().unwrap().strand_symbol().to_owned(),
                                    ));
                                    */
                                }

                                //println!("{:?}\t{:?}", new_match[0], new_match[1]);
                                //matches.push(new_match.unwrap());
                                /*
                                println!(
                                    "{}",
                                    search
                                        .column("marker")
                                        .unwrap()
                                        .utf8()
                                        .unwrap()
                                        .into_iter()
                                        .collect::<Vec<Option<&str>>>()[0]
                                        .unwrap()
                                        );
                                        */
                                //match_list.push(search.column("marker").iter().collect()),
                            }
                            _ => println!(
                                "Protein {} corresponds to >1 marker, so skipping:\n {:?}",
                                prot_id, search
                            ),
                        }
                    };
                }
                Err(_) => {
                    // Break because this means we are down at the FASTA portion of the prodigal GFF
                    // i.e., no more gene records
                    break;
                }
            }
        }

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
