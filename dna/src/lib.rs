//! A general-purpose genomics crate for dealing with DNA.

#![warn(missing_docs)]

use std::{convert::TryFrom, fmt::Display, iter::FromIterator, str::FromStr};

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

/// PackedDNA
#[derive(Debug)]
pub struct PackedDna {
    packed_dna: Vec<u8>,
    last_nuc_set_count: usize,
    a_count: usize,
    c_count: usize,
    g_count: usize,
    t_count: usize,
}

/// Implementation for PackedDNA
impl PackedDna {
    /// Function to get nucleotide at a given index
    pub fn get(&self, idx: usize) -> Result<Nuc, String> {
        let vec_index = (idx - 1) / 4;
        let bit_index = (idx - 1) % 4;
        if (vec_index >= self.packed_dna.len())
            || (vec_index >= self.packed_dna.len() - 1 && bit_index >= self.last_nuc_set_count)
        {
            let error = format!("Index {} is greater than the given DNA Length", idx);
            return Err(error);
        }
        let mut binary_rep = format!("{:08b}", self.packed_dna[vec_index]);
        if (vec_index == self.packed_dna.len() - 1) && (self.last_nuc_set_count != 0) {
            binary_rep = binary_rep
                .chars()
                .rev()
                .take(self.last_nuc_set_count * 2)
                .collect();
            binary_rep = binary_rep.chars().rev().collect();
        }
        let char_vec: Vec<char> = binary_rep.chars().collect();
        let str_rep = format!("{}{}", char_vec[bit_index * 2], char_vec[bit_index * 2 + 1]);
        match str_rep.as_str() {
            "00" => Ok(Nuc::A),
            "01" => Ok(Nuc::C),
            "10" => Ok(Nuc::G),
            "11" => Ok(Nuc::T),
            _ => Err("Invalid String Encountered".to_string()),
        }
    }

    /// Get the counts of individual nucleotides
    pub fn get_counts(&self) -> Vec<(char, usize)> {
        return vec![
            ('A', self.a_count),
            ('C', self.c_count),
            ('G', self.g_count),
            ('T', self.t_count),
        ];
    }
}

/// FromString implementation for PackedDna
/// Takes in the string DNA as the input and stores the DNA in efficient way
impl FromStr for PackedDna {
    type Err = ParseNucError<String>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let string_dna = s.to_ascii_uppercase();
        let extra_nuc = string_dna.len() % 4;
        let mut vec: Vec<u8> = Vec::new();
        let mut curr = 0;
        let (mut a, mut c, mut g, mut t) = (0, 0, 0, 0);
        for (i, char) in string_dna.chars().enumerate() {
            if (i != 0) && (i % 4) == 0 {
                vec.push(curr);
                curr = 0;
            }
            match char {
                'A' => {
                    curr <<= 2;
                    a += 1
                }
                'C' => {
                    curr = curr << 2 | 1;
                    c += 1
                }
                'G' => {
                    curr = curr << 2 | 2;
                    g += 1
                }
                'T' => {
                    curr = curr << 2 | 3;
                    t += 1
                }
                _ => return Err(ParseNucError(string_dna)),
            }
        }
        vec.push(curr);
        Ok(PackedDna {
            packed_dna: vec,
            last_nuc_set_count: extra_nuc,
            a_count: a,
            c_count: c,
            g_count: g,
            t_count: t,
        })
    }
}

/// FromIterator implementation for PackedDna
/// Takes in the vector of Nuc as the input and stores the DNA in efficient way
impl FromIterator<Nuc> for PackedDna {
    fn from_iter<I: IntoIterator<Item = Nuc>>(iter: I) -> Self {
        let mut extra_nuc = 0;
        let mut vec: Vec<u8> = Vec::new();
        let mut curr = 0;
        let (mut a, mut c, mut g, mut t) = (0, 0, 0, 0);
        for (counter, nuc) in iter.into_iter().enumerate() {
            if (counter != 0) && (counter % 4) == 0 {
                vec.push(curr);
                curr = 0;
                extra_nuc = 0;
            }
            match nuc {
                Nuc::A => {
                    curr <<= 2;
                    a += 1
                }
                Nuc::C => {
                    curr = curr << 2 | 1;
                    c += 1
                }
                Nuc::G => {
                    curr = curr << 2 | 2;
                    g += 1
                }
                Nuc::T => {
                    curr = curr << 2 | 3;
                    t += 1
                }
            }
            extra_nuc += 1;
        }
        extra_nuc %= 4;
        vec.push(curr);
        PackedDna {
            packed_dna: vec,
            last_nuc_set_count: extra_nuc,
            a_count: a,
            c_count: c,
            g_count: g,
            t_count: t,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_from_nuc_a_uppercase() {
        let nuc_from_char = Nuc::try_from('A').unwrap();
        assert_eq!(nuc_from_char, Nuc::A);
    }

    #[test]
    fn try_from_nuc_a_lowercase() {
        let nuc_from_char = Nuc::try_from('a').unwrap();
        assert_eq!(nuc_from_char, Nuc::A);
    }

    #[test]
    fn try_from_nuc_c_uppercase() {
        let nuc_from_char = Nuc::try_from('C').unwrap();
        assert_eq!(nuc_from_char, Nuc::C);
    }

    #[test]
    fn try_from_nuc_c_lowercase() {
        let nuc_from_char = Nuc::try_from('c').unwrap();
        assert_eq!(nuc_from_char, Nuc::C);
    }

    #[test]
    fn try_from_nuc_g_uppercase() {
        let nuc_from_char = Nuc::try_from('G').unwrap();
        assert_eq!(nuc_from_char, Nuc::G);
    }

    #[test]
    fn try_from_nuc_g_lowercase() {
        let nuc_from_char = Nuc::try_from('g').unwrap();
        assert_eq!(nuc_from_char, Nuc::G);
    }

    #[test]
    fn try_from_nuc_t_uppercase() {
        let nuc_from_char = Nuc::try_from('T').unwrap();
        assert_eq!(nuc_from_char, Nuc::T);
    }

    #[test]
    fn try_from_nuc_t_lowercase() {
        let nuc_from_char = Nuc::try_from('t').unwrap();
        assert_eq!(nuc_from_char, Nuc::T);
    }

    #[test]
    fn try_from_nuc_negative() {
        let nuc_from_char = Nuc::try_from('X');
        match nuc_from_char {
            Ok(ref _x) => {}
            Err(e) => {
                assert_eq!("failed to parse nucleotide from X", e.to_string());
            }
        }
    }

    #[test]
    fn from_str_nuc_a_uppercase() {
        let nuc_from_string = Nuc::from_str("A").unwrap();
        assert_eq!(nuc_from_string, Nuc::A);
    }

    #[test]
    fn from_str_nuc_a_lowercase() {
        let nuc_from_string = Nuc::from_str("a").unwrap();
        assert_eq!(nuc_from_string, Nuc::A);
    }

    #[test]
    fn from_str_nuc_c_uppercase() {
        let nuc_from_string = Nuc::from_str("C").unwrap();
        assert_eq!(nuc_from_string, Nuc::C);
    }

    #[test]
    fn from_str_nuc_c_lowercase() {
        let nuc_from_string = Nuc::from_str("c").unwrap();
        assert_eq!(nuc_from_string, Nuc::C);
    }

    #[test]
    fn from_str_nuc_g_uppercase() {
        let nuc_from_string = Nuc::from_str("G").unwrap();
        assert_eq!(nuc_from_string, Nuc::G);
    }

    #[test]
    fn from_str_nuc_g_lowercase() {
        let nuc_from_string = Nuc::from_str("g").unwrap();
        assert_eq!(nuc_from_string, Nuc::G);
    }

    #[test]
    fn from_str_nuc_t_uppercase() {
        let nuc_from_string = Nuc::from_str("T").unwrap();
        assert_eq!(nuc_from_string, Nuc::T);
    }

    #[test]
    fn from_str_nuc_t_lowercase() {
        let nuc_from_string = Nuc::from_str("t").unwrap();
        assert_eq!(nuc_from_string, Nuc::T);
    }

    #[test]
    fn from_string_nuc_negative() {
        let nuc_from_string = Nuc::from_str("X");
        match nuc_from_string {
            Ok(ref _x) => {}
            Err(e) => {
                assert_eq!("failed to parse nucleotide from X", e.to_string());
            }
        }
    }

    #[test]
    fn from_string_test_len10_positive() {
        let dna_from_string = PackedDna::from_str("ACGTTGCACT").unwrap();
        let vec = [27, 228, 7];
        assert_eq!(dna_from_string.packed_dna, vec);
        assert_eq!(dna_from_string.last_nuc_set_count, 2);
    }

    #[test]
    fn from_string_test_len8_positive() {
        let dna_from_string = PackedDna::from_str("ACGTTGCA").unwrap();
        let vec = [27, 228];
        assert_eq!(dna_from_string.packed_dna, vec);
        assert_eq!(dna_from_string.last_nuc_set_count, 0);
    }

    #[test]
    fn from_string_test_lowercase_positive() {
        let dna_from_string = PackedDna::from_str("acgttgca").unwrap();
        let vec = [27, 228];
        assert_eq!(dna_from_string.packed_dna, vec);
        assert_eq!(dna_from_string.last_nuc_set_count, 0);
    }

    #[test]
    fn from_string_invalid_nuc() {
        let packed_dna = PackedDna::from_str("acgtx");
        match packed_dna {
            Ok(ref _x) => {}
            Err(e) => {
                assert_eq!("failed to parse nucleotide from ACGTX", e.to_string());
            }
        }
    }

    #[test]
    fn from_iter_test_len8_positive() {
        let dna_from_string = PackedDna::from_iter([
            Nuc::A,
            Nuc::C,
            Nuc::G,
            Nuc::T,
            Nuc::A,
            Nuc::C,
            Nuc::G,
            Nuc::T,
        ]);
        let vec = [27, 27];
        assert_eq!(dna_from_string.packed_dna, vec);
        assert_eq!(dna_from_string.last_nuc_set_count, 0);
    }

    #[test]
    fn from_iter_test_len11_positive() {
        let dna_from_string = PackedDna::from_iter([
            Nuc::A,
            Nuc::C,
            Nuc::G,
            Nuc::T,
            Nuc::A,
            Nuc::C,
            Nuc::G,
            Nuc::T,
            Nuc::T,
            Nuc::T,
            Nuc::T,
        ]);
        let vec = [27, 27, 63];
        assert_eq!(dna_from_string.packed_dna, vec);
        assert_eq!(dna_from_string.last_nuc_set_count, 3);
    }

    #[test]
    fn get_nuc_test_positive_a() {
        let dna_from_string = PackedDna::from_str("ACGTTGCACT").unwrap();
        match dna_from_string.get(1) {
            Ok(x) => assert_eq!(x, Nuc::A),
            Err(e) => println!("Oops! You ran into an error: {e:?}"),
        }
    }

    #[test]
    fn get_nuc_test_positive_c() {
        let dna_from_string = PackedDna::from_str("ACGTTGCACT").unwrap();
        match dna_from_string.get(7) {
            Ok(x) => assert_eq!(x, Nuc::C),
            Err(e) => println!("Oops! You ran into an error: {e:?}"),
        }
    }

    #[test]
    fn get_nuc_test_positive_g() {
        let dna_from_string = PackedDna::from_str("ACGTTGCACT").unwrap();
        match dna_from_string.get(3) {
            Ok(x) => assert_eq!(x, Nuc::G),
            Err(e) => println!("Oops! You ran into an error: {e:?}"),
        }
    }

    #[test]
    fn get_nuc_test_positive_t() {
        let dna_from_string = PackedDna::from_str("ACGTTGCACT").unwrap();
        match dna_from_string.get(10) {
            Ok(x) => assert_eq!(x, Nuc::T),
            Err(e) => println!("Oops! You ran into an error: {e:?}"),
        }
    }

    #[test]
    fn get_nuc_test_negative() {
        let dna_from_string = PackedDna::from_str("ACGTTGCACT").unwrap();
        let get11 = dna_from_string.get(11);
        match get11 {
            Ok(_x) => {}
            Err(e) => assert_eq!("Index 11 is greater than the given DNA Length", e),
        }
    }

    #[test]
    fn get_nuc_count_test() {
        let dna_from_string = PackedDna::from_str("ACGTTGCACT").unwrap();
        let counts = dna_from_string.get_counts();
        assert_eq!(counts, vec![('A', 2), ('C', 3), ('G', 2), ('T', 3)]);
    }
}
