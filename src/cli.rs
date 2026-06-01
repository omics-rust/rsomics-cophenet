use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, HelpSpec, Origin};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-cophenet", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    pub linkage: PathBuf,

    /// Original condensed (or --square) distance matrix; enables the cophenetic correlation coefficient.
    #[arg(long = "dist")]
    dist: Option<PathBuf>,

    /// Treat --dist as a square matrix rather than a condensed vector.
    #[arg(long = "dist-square")]
    dist_square: bool,

    /// Emit a square cophenetic matrix instead of the condensed vector.
    #[arg(long)]
    square: bool,

    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        rsomics_cophenet::run(
            &self.linkage,
            self.dist.as_deref(),
            self.dist_square,
            self.square,
            &mut out,
        )
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Cophenetic distances from a linkage matrix; cophenetic correlation when given the original distances.",
    origin: Some(Origin {
        upstream: "scipy.cluster.hierarchy.cophenet",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &[
        "<linkage.tsv> [--dist orig.tsv [--dist-square]] [--square] [-o cophenetic.tsv]",
    ],
    sections: &[],
    examples: &[
        Example {
            description: "Condensed cophenetic distance vector from a UPGMA linkage matrix",
            command: "rsomics-cophenet Z.tsv -o cophenetic.tsv",
        },
        Example {
            description: "Also report the cophenetic correlation coefficient against the original distances",
            command: "rsomics-cophenet Z.tsv --dist original.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
