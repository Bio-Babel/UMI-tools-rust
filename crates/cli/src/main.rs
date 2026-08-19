//! The `umi_tools` binary (Python source: the `umi_tools = "umi_tools.umi_tools:main"`
//! console-script entry point declared in `pyproject.toml`).

use std::process::ExitCode;

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().collect();
    match umi_tools_cli::run(&argv) {
        Ok(code) => ExitCode::from(u8::try_from(code).unwrap_or(1)),
        Err(e) => ExitCode::from(u8::try_from(umi_tools_cli::report_fatal(&e)).unwrap_or(1)),
    }
}
