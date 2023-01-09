// TODO: implement a nucleotide counter
//
// command line argument parsing has been provided
// you must use the PackedDna struct you previously implemented
// if there is any functionality you would like to add to PackedDna feel free to do so in the DNA
// crate
//
// If run with `nuccount --dna ACGTTT" it should print the following to stdout:
// ```
// Input: ACGTTT
//
// A: 1
// C: 1
// G: 1
// T: 3
// ```
//
// be sure to exit with informative error messages if the input is invalid

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
        /// Doubt: I don't understand why we need to exit the service when someone passes incorrect string.
        Err(e) => {
            println!("Invalid DNA String Passed \nError: {}", e);
            process::exit(1);
        }
    }
}