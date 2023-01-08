//! A general-purpose genomics crate for dealing with DNA.

#![warn(missing_docs)]

use std::{convert::TryFrom, fmt::Display, iter::FromIterator, str::FromStr};

// TODO: add a packed module with the PackedDna struct
//
// this struct must have the following:
// 1. A representation that is more memory efficient that simply storing a vector of `Nuc`
// 2. A FromStr implementation (should be case insensitive like the `Nuc` impl)
// 3. A `FromIterator` implementation to construct it from an iterator over `Nuc`s
// 4. A `fn get(&self, idx: usize) -> Nuc` getter for a particular nucleotide
//
// Make sure to unit test and document all elements
// Also, the internal representation of the PackedDna struct should be privately scoped

/// A nucleotide
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nuc {
    /// Adenine
    A,
    /// Cytosine
    C,
    /// Guanine
    G,
    /// Thymine
    T,
}

/// An error that can occur when parsing a nucleotide.
#[derive(Debug, thiserror::Error)]
#[error("failed to parse nucleotide from {0}")]
pub struct ParseNucError<T: Display>(T);

impl TryFrom<char> for Nuc {
    type Error = ParseNucError<char>;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value.to_ascii_uppercase() {
            'A' => Ok(Self::A),
            'C' => Ok(Self::C),
            'G' => Ok(Self::G),
            'T' => Ok(Self::T),
            _ => Err(ParseNucError(value)),
        }
    }
}

impl FromStr for Nuc {
    type Err = ParseNucError<String>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let upper = s.to_ascii_uppercase();
        match upper.as_str() {
            "A" => Ok(Self::A),
            "C" => Ok(Self::C),
            "G" => Ok(Self::G),
            "T" => Ok(Self::T),
            _ => Err(ParseNucError(upper)),
        }
    }
}

// Struct for PackedDNA
struct PackedDna {
    packed_dna: Vec<u8>,
}

// Implementation for PackedDNA
impl PackedDna {
    fn new(vec: Vec<u8>) -> PackedDna {
        PackedDna { packed_dna: vec }
    }

    fn print(&self) {
        for x in &self.packed_dna {
            print!("{}", x);
        }
        println!()
    }

    fn get(&self, idx: usize) -> Nuc {
        let nuc = self.packed_dna[idx];
        match nuc {
            0b00 => Nuc::A,
            0b01 => Nuc::C,
            0b10 => Nuc::G,
            _ => Nuc::T,
        }
    }
}

impl FromStr for PackedDna {
    type Err = ParseNucError<String>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut packed_dna = Vec::new();
        let upper = s.to_ascii_uppercase();
        for char in upper.chars() {
            match char {
                'A' => packed_dna.push(0b00),
                'C' => packed_dna.push(0b01),
                'G' => packed_dna.push(0b10),
                'T' => packed_dna.push(0b11),
                _ => return Err(ParseNucError(upper)),
            }
        }
        Ok(PackedDna::new(packed_dna))
    }
}

impl FromIterator<Nuc> for PackedDna {
    fn from_iter<I: IntoIterator<Item = Nuc>>(iter: I) -> Self {
        let mut temp_dna = Vec::new();
        for nuc in iter {
            match nuc {
                Nuc::A => temp_dna.push(0),
                Nuc::C => temp_dna.push(1),
                Nuc::G => temp_dna.push(2),
                Nuc::T => temp_dna.push(3),
            }
        }
        PackedDna::new(temp_dna)
    }
}

#[cfg(test)]
mod tests {
    // TODO: fill in tests

    #[test]
    fn tryfrom_char() {
        assert!(false);
    }

    #[test]
    fn fromstr() {
        assert!(false);
    }
}
