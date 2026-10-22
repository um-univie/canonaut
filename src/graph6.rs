//! Parser for the graph6 encoding (Brendan McKay, nauty).
//!
//! graph6 encodes undirected simple graphs. Encoding:
//!   - n ≤ 62:       1 byte  = n + 63
//!   - n ≤ 258047:   byte 126, then 3 bytes each = 6-bit chunk + 63
//!   - larger n:     bytes 126,126, then 6 bytes each = 6-bit chunk + 63
//!
//! Adjacency: upper triangle (j=1..n-1, i=0..j-1), LSB-right packed 6 bits/byte,
//! each byte stored as value + 63.

use crate::structs::graph::{add_one_edge, DenseGraph};
use crate::utilities::words_needed;

#[derive(Debug, PartialEq)]
pub enum Graph6Error {
    EmptyInput,
    NotGraph6Format,
    InvalidByte(u8),
    UnexpectedEnd,
    TrailingData,
    InvalidPadding,
    UnsupportedOrder,
}

impl std::fmt::Display for Graph6Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "empty input"),
            Self::NotGraph6Format => write!(f, "not graph6 (sparse6/digraph6 not supported)"),
            Self::InvalidByte(b) => write!(f, "invalid graph6 byte: {b}"),
            Self::UnexpectedEnd => write!(f, "unexpected end of data"),
            Self::TrailingData => write!(f, "extra data after graph6 adjacency"),
            Self::InvalidPadding => write!(f, "nonzero graph6 padding bits"),
            Self::UnsupportedOrder => write!(f, "graph order exceeds the supported limit"),
        }
    }
}

impl std::error::Error for Graph6Error {}

/// Parse one graph6 string (bytes, no trailing newline required) into a `DenseGraph`.
pub fn parse_graph6(input: &[u8]) -> Result<DenseGraph, Graph6Error> {
    // Strip optional `>>graph6<<` header
    let input = input.strip_prefix(b">>graph6<<").unwrap_or(input);

    // Strip trailing CR/LF
    let input = match input.iter().rposition(|&byte| byte != b'\n' && byte != b'\r') {
        Some(last) => &input[..=last],
        None => return Err(Graph6Error::EmptyInput),
    };

    if input.is_empty() {
        return Err(Graph6Error::EmptyInput);
    }

    // Reject sparse6 (':') and digraph6 ('&')
    match input[0] {
        b':' | b'&' => return Err(Graph6Error::NotGraph6Format),
        _ => {}
    }

    let (size_start, size_end) = if input[0] != 126 {
        (0, 1)
    } else if input.get(1) == Some(&126) {
        (2, 8)
    } else {
        (1, 4)
    };
    let size_bytes = input.get(size_start..size_end).ok_or(Graph6Error::UnexpectedEnd)?;
    let mut vertex_count = 0usize;
    for &byte in size_bytes {
        if !(63..=126).contains(&byte) {
            return Err(Graph6Error::InvalidByte(byte));
        }
        vertex_count = (vertex_count << 6) | (byte - 63) as usize;
    }
    if vertex_count > (crate::consts::ALGORITHM_INFINITY - 2) as usize {
        return Err(Graph6Error::UnsupportedOrder);
    }
    let adjacency_bits = vertex_count.checked_mul(vertex_count.saturating_sub(1))
        .ok_or(Graph6Error::UnsupportedOrder)? / 2;
    let required_bytes = adjacency_bits.div_ceil(6);
    let adjacency_data = &input[size_end..];
    if adjacency_data.len() < required_bytes {
        return Err(Graph6Error::UnexpectedEnd);
    }
    if adjacency_data.len() > required_bytes {
        return Err(Graph6Error::TrailingData);
    }
    // Validate each byte once before allocating, rather than up to six times
    // in the adjacency loop. Invalid/truncated input cannot trigger large allocations.
    if let Some(&byte) = adjacency_data.iter().find(|&&byte| !(63..=126).contains(&byte)) {
        return Err(Graph6Error::InvalidByte(byte));
    }
    let padding = (6 - adjacency_bits % 6) % 6;
    if padding != 0 && (adjacency_data[required_bytes - 1] - 63) & ((1 << padding) - 1) != 0 {
        return Err(Graph6Error::InvalidPadding);
    }
    let words_per_vertex = words_needed(vertex_count) as u32;
    let graph_words = vertex_count.checked_mul(words_per_vertex as usize)
        .filter(|&words| words <= isize::MAX as usize / size_of::<usize>())
        .ok_or(Graph6Error::UnsupportedOrder)?;
    let mut adjacency = vec![0usize; graph_words];
    let mut bit_index = 0usize;

    // Upper triangle: column from 1..vertex_count, row from 0..column
    for column in 1..vertex_count {
        for row in 0..column {
            let data_byte_index = bit_index / 6;
            let bit_offset_in_byte = 5 - (bit_index % 6); // MSB-first within each 6-bit group

            if ((adjacency_data[data_byte_index] - 63) >> bit_offset_in_byte) & 1 == 1 {
                add_one_edge(&mut adjacency, row, column, words_per_vertex);
            }

            bit_index += 1;
        }
    }

    Ok(DenseGraph {
        adjacency,
        number_of_vertices: vertex_count as u32,
        words_per_vertex,
        directed: false,
        colors: None,
    })
}

/// Parse all graphs from a graph6 file (one graph per line).
pub fn parse_graph6_file(contents: &[u8]) -> Vec<Result<DenseGraph, Graph6Error>> {
    contents
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.is_empty() && *line != b"\r")
        .map(parse_graph6)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structs::graph::get_graph_row;
    use crate::structs::BitSetOpsMSB0;

    fn edge(graph: &DenseGraph, row: usize, column: usize) -> bool {
        get_graph_row(&graph.adjacency, row, graph.words_per_vertex).is_bit_set(column)
    }

    fn edge_count(graph: &DenseGraph) -> usize {
        let vertex_count = graph.number_of_vertices as usize;
        (0..vertex_count)
            .flat_map(|row| (row + 1..vertex_count).map(move |column| (row, column)))
            .filter(|&(row, column)| edge(graph, row, column))
            .count()
    }

    // --- n encoding ---

    #[test]
    fn empty_graph_n0() {
        // '?' = 63 = n + 63 → n=0
        let graph = parse_graph6(b"?").unwrap();
        assert_eq!(graph.number_of_vertices, 0);
        assert!(graph.adjacency.is_empty());
    }

    #[test]
    fn single_vertex_n1() {
        // '@' = 64 → n=1, no adjacency bits
        let graph = parse_graph6(b"@").unwrap();
        assert_eq!(graph.number_of_vertices, 1);
        assert_eq!(edge_count(&graph), 0);
    }

    // --- small well-known graphs ---

    #[test]
    fn k2_complete() {
        // 'A' → n=2; '_' = 95, 95-63=32=0b100000, bit5=1 → edge(0,1)
        let graph = parse_graph6(b"A_").unwrap();
        assert_eq!(graph.number_of_vertices, 2);
        assert!(edge(&graph, 0, 1));
        assert!(edge(&graph, 1, 0));
        assert_eq!(edge_count(&graph), 1);
    }

    #[test]
    fn empty_n2() {
        // 'A' → n=2; '?' = 63, bits=0 → no edges
        let graph = parse_graph6(b"A?").unwrap();
        assert_eq!(graph.number_of_vertices, 2);
        assert!(!edge(&graph, 0, 1));
        assert_eq!(edge_count(&graph), 0);
    }

    #[test]
    fn k3_complete() {
        // 'B' → n=3; 'w' = 119, 119-63=56=0b111000 → all 3 edges set
        let graph = parse_graph6(b"Bw").unwrap();
        assert_eq!(graph.number_of_vertices, 3);
        assert!(edge(&graph, 0, 1));
        assert!(edge(&graph, 0, 2));
        assert!(edge(&graph, 1, 2));
        assert_eq!(edge_count(&graph), 3);
    }

    #[test]
    fn path_p3() {
        // 'B' → n=3; 'g' = 103, 103-63=40=0b101000 → edges (0,1) and (1,2)
        let graph = parse_graph6(b"Bg").unwrap();
        assert_eq!(graph.number_of_vertices, 3);
        assert!(edge(&graph, 0, 1));
        assert!(!edge(&graph, 0, 2));
        assert!(edge(&graph, 1, 2));
        assert_eq!(edge_count(&graph), 2);
    }

    #[test]
    fn k4_complete() {
        // 'C' → n=4; '~' = 126, 126-63=63=0b111111 → all 6 edges
        let graph = parse_graph6(b"C~").unwrap();
        assert_eq!(graph.number_of_vertices, 4);
        for row in 0..4 {
            for column in 0..4 {
                assert_eq!(edge(&graph, row, column), row != column, "K4 edge({row},{column})");
            }
        }
        assert_eq!(edge_count(&graph), 6);
    }

    #[test]
    fn cycle_c4() {
        // 'C' → n=4; 'l' = 108, 108-63=45=0b101101
        // bits: (0,1)=1, (0,2)=0, (1,2)=1, (0,3)=1, (1,3)=0, (2,3)=1
        // edges: 0-1, 1-2, 0-3, 2-3  → 4-cycle 0-1-2-3-0
        let graph = parse_graph6(b"Cl").unwrap();
        assert_eq!(graph.number_of_vertices, 4);
        assert!(edge(&graph, 0, 1));
        assert!(!edge(&graph, 0, 2));
        assert!(edge(&graph, 1, 2));
        assert!(edge(&graph, 0, 3));
        assert!(!edge(&graph, 1, 3));
        assert!(edge(&graph, 2, 3));
        assert_eq!(edge_count(&graph), 4);
    }

    #[test]
    fn undirected_symmetry() {
        // Every parsed graph must have a symmetric adjacency matrix
        for g6_string in [b"A_".as_ref(), b"Bw", b"Bg", b"C~", b"Cl"] {
            let graph = parse_graph6(g6_string).unwrap();
            let vertex_count = graph.number_of_vertices as usize;
            for row in 0..vertex_count {
                for column in 0..vertex_count {
                    assert_eq!(
                        edge(&graph, row, column),
                        edge(&graph, column, row),
                        "symmetry failed for {:?} at ({row},{column})",
                        std::str::from_utf8(g6_string)
                    );
                }
            }
        }
    }

    #[test]
    fn header_stripped() {
        let with_header = b">>graph6<<A_";
        let without = b"A_";
        let parsed_with_header = parse_graph6(with_header).unwrap();
        let parsed_without_header = parse_graph6(without).unwrap();
        assert_eq!(parsed_with_header.number_of_vertices, parsed_without_header.number_of_vertices);
        assert_eq!(parsed_with_header.adjacency, parsed_without_header.adjacency);
    }

    #[test]
    fn trailing_newline_stripped() {
        let with_newline = parse_graph6(b"Bw\n").unwrap();
        let without_newline = parse_graph6(b"Bw").unwrap();
        assert_eq!(with_newline.adjacency, without_newline.adjacency);
    }

    #[test]
    fn multiline_file() {
        // Three graphs: n=2 empty, K2, K3
        let file_contents = b"A?\nA_\nBw\n";
        let results = parse_graph6_file(file_contents);
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].as_ref().unwrap().number_of_vertices, 2);
        assert_eq!(edge_count(results[0].as_ref().unwrap()), 0);
        assert_eq!(edge_count(results[1].as_ref().unwrap()), 1);
        assert_eq!(edge_count(results[2].as_ref().unwrap()), 3);
    }

    // --- error cases ---

    #[test]
    fn err_empty() {
        assert_eq!(parse_graph6(b""), Err(Graph6Error::EmptyInput));
        assert_eq!(parse_graph6(b"\n"), Err(Graph6Error::EmptyInput));
    }

    #[test]
    fn err_sparse6_prefix() {
        assert_eq!(
            parse_graph6(b":Bw"),
            Err(Graph6Error::NotGraph6Format)
        );
    }

    #[test]
    fn err_digraph6_prefix() {
        assert_eq!(
            parse_graph6(b"&Bw"),
            Err(Graph6Error::NotGraph6Format)
        );
    }

    #[test]
    fn err_invalid_byte() {
        // 'B' → n=3; space (0x20=32) is < 63 → invalid adjacency byte
        assert_eq!(
            parse_graph6(b"B\x20"),
            Err(Graph6Error::InvalidByte(0x20))
        );
    }

    #[test]
    fn reject_truncated_or_extra_adjacency() {
        for input in [b"A".as_slice(), b"D~", b"~??~"] {
            assert_eq!(parse_graph6(input), Err(Graph6Error::UnexpectedEnd));
        }
        for input in [b"A?trailing".as_slice(), b"??", b"@?"] {
            assert_eq!(parse_graph6(input), Err(Graph6Error::TrailingData));
        }
        assert_eq!(parse_graph6(b"A@"), Err(Graph6Error::InvalidPadding));
    }

    #[test]
    fn reject_invalid_order_bytes_without_arithmetic_overflow() {
        for byte in (0..=255u8).filter(|b| !(63..=126).contains(b)) {
            if byte != b'\n' && byte != b'\r' && byte != b':' && byte != b'&' {
                assert_eq!(parse_graph6(&[byte]), Err(Graph6Error::InvalidByte(byte)));
            }
            // A payload byte after the size keeps CR/LF size bytes from being
            // stripped as line terminators.
            for size in [vec![126, 63, 63, 126], vec![126, 126, 63, 63, 63, 63, 63, 126]] {
                let start = if size.len() == 4 { 1 } else { 2 };
                for i in start..size.len() {
                    // A '~' second byte switches to the long encoding, but all
                    // invalid bytes in this test are outside the '~' range.
                    let mut input = size.clone();
                    input[i] = byte;
                    input.push(63);
                    assert_eq!(parse_graph6(&input), Err(Graph6Error::InvalidByte(byte)));
                }
            }
        }
        assert_eq!(parse_graph6(b"~?? "), Err(Graph6Error::InvalidByte(32)));
        assert_eq!(parse_graph6(b"~~????? "), Err(Graph6Error::InvalidByte(32)));
        assert_eq!(parse_graph6(b"~~~~~~~~"), Err(Graph6Error::UnsupportedOrder));
    }

    #[test]
    fn valid_extended_order_and_zero_padding() {
        let mut input = b"~??~".to_vec(); // n=63
        input.resize(4 + (63usize * 62 / 2).div_ceil(6), b'?');
        let graph = parse_graph6(&input).unwrap();
        assert_eq!(graph.number_of_vertices(), 63);
        assert!(graph.adjacency.iter().all(|&word| word == 0));
        // A valid long order with no payload must fail before allocating.
        assert_eq!(parse_graph6(b"~~???~??"), Err(Graph6Error::UnexpectedEnd));
    }
}
