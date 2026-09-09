fn main() {
    println!("cargo:rerun-if-changed=benches/c_nauty_helper.c");

    // The C reference implementation is only needed for the comparison
    // benchmarks; a plain `cargo build` must not require libnauty.
    if std::env::var_os("CARGO_FEATURE_C_NAUTY_BENCH").is_none() {
        return;
    }

    println!("cargo:rerun-if-env-changed=NAUTY_INCLUDE_DIR");
    println!("cargo:rerun-if-env-changed=NAUTY_LIB_DIR");
    let include =
        std::env::var("NAUTY_INCLUDE_DIR").unwrap_or_else(|_| "/usr/local/include".into());
    let library = std::env::var("NAUTY_LIB_DIR").unwrap_or_else(|_| "/usr/local/lib".into());
    cc::Build::new()
        .file("benches/c_nauty_helper.c")
        .include(include)
        .opt_level(3)
        .compile("c_nauty_helper");

    println!("cargo:rustc-link-search=native={library}");
    println!("cargo:rustc-link-lib=dylib=nauty");
}
