// Reads one uncompressed MRT RIB dump on stdin and writes the pair counts
// that the former bgpdump | get_neigh.py / get_pairs.py | awk pipelines in
// parsebview.sh produced, in one pass and several times faster:
//   <as_as out>:     "a b count" for consecutive ASes on each AS path
//   <prefix_as out>: "prefix asn count" for each prefix and its origin AS(es)
// Lines are unsorted; parsebview.sh sums them over all dumps and sorts them.
//
// The output matches the old pipeline byte for byte, except that IPv6
// prefixes are written in RFC 5952 form. bgpdump wrote a single zero group
// as "::" (e.g. 2001:668::3:ffff:1:0:2270 instead of 2001:668:0:3:ffff:1:0:2270).
//
// Usage: zcat bview.gz | bviewpairs <as_as out> <prefix_as out>
use bgpkit_parser::BgpkitParser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

fn clean(s: &str) -> String {
    // sed -e 's/[{},]/ /g'
    s.chars().map(|c| if c == '{' || c == '}' || c == ',' { ' ' } else { c }).collect()
}

fn write_counts(path: &str, counts: &HashMap<String, u64>) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    for (k, v) in counts {
        writeln!(w, "{} {}", k, v)?;
    }
    w.flush()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: zcat bview.gz | {} <as_as out> <prefix_as out>", args[0]);
        std::process::exit(2);
    }
    let stdin = std::io::stdin();
    let parser = BgpkitParser::from_reader(std::io::BufReader::with_capacity(1 << 20, stdin.lock()));

    let mut as_as: HashMap<String, u64> = HashMap::new();
    let mut prefix_as: HashMap<String, u64> = HashMap::new();

    for elem in parser {
        let prefix = elem.prefix.to_string();
        let path = elem.as_path.as_ref().map(|p| p.to_string()).unwrap_or_default();

        // get_neigh.py: pairs of consecutive tokens, smaller (as strings) first
        let cleaned = clean(&path);
        let toks: Vec<&str> = cleaned.split_whitespace().collect();
        for w in toks.windows(2) {
            if w[0] == w[1] {
                continue;
            }
            let (a, b) = if w[0] < w[1] { (w[0], w[1]) } else { (w[1], w[0]) };
            *as_as.entry(format!("{} {}", a, b)).or_insert(0) += 1;
        }

        // awk: first and last token of "prefix path", then get_pairs.py
        let line = format!("{} {}", prefix, path);
        let raw: Vec<&str> = line.split_whitespace().collect();
        if raw.is_empty() {
            continue;
        }
        let s = clean(&format!("{} {}", raw[0], raw[raw.len() - 1]));
        let toks: Vec<&str> = s.split_whitespace().collect();
        if toks.is_empty() {
            continue;
        }
        for asn in &toks[1..] {
            if *asn == toks[0] {
                continue;
            }
            *prefix_as.entry(format!("{} {}", toks[0], asn)).or_insert(0) += 1;
        }
    }

    write_counts(&args[1], &as_as).expect("writing as_as");
    write_counts(&args[2], &prefix_as).expect("writing prefix_as");
}
