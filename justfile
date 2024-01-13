view-cov: coverage
  firefox target/debug/coverage/index.html

coverage:
  rm -r target/debug/coverage/
  CARGO_INCREMENTAL=0 RUSTFLAGS='-Cinstrument-coverage' LLVM_PROFILE_FILE='cargo-test-%p-%m.profraw' cargo test
  grcov . -s . --binary-path ./target/debug/ -t html --branch --ignore-not-existing -o ./target/debug/coverage/
  rm *.profraw

ebnf:
  ebnf2railroad grammar/peptidoglycan.ebnf -t PGLang
  firefox grammar/peptidoglycan.html
