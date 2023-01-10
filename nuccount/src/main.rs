use dna::PackedDna;
use std::{process, str::FromStr};
use structopt::StructOpt;

/// Count the number of occurrences of each nucleotide in the provided DNA.
#[derive(Debug, StructOpt)]
struct Opts {
    /// The DNA sequence for which we should retrieve a nucleotide count.
    ///
    /// It is case insensitive but only nucleotides A, C, G and T are supported.
    #[structopt(short = "d", long, required = true)]
    dna: String,
}

fn main() {
    let opts = Opts::from_args();
    let dna = opts.dna;
    println!("Input: {}\n", &dna);
    let packed_dna = PackedDna::from_str(&dna);
    match packed_dna {
        Ok(ref _x) => {
            let nuc_counts = packed_dna.unwrap().get_counts();
            for (nuc, counts) in nuc_counts {
                println!("{} {}", nuc, counts);
            }
        }
        // Doubt: I don't understand why we need to exit the service when someone passes incorrect string.
        Err(e) => {
            println!(
                "Invalid character for Nuclieotide passed in DNA string\nError: {}",
                e
            );
            process::exit(1);
        }
    }
}
