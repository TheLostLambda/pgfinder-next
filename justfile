watch:
  watchexec -e rs,toml,kdl just test

test:
  cargo test --workspace

# FIXME: The --workspace flag appears to be broken at the moment: https://github.com/mitsuhiko/insta/issues/396
review:
  cargo insta test --review
  sh -c 'for crate in `ls crates`; do (cd "crates/$crate" && cargo insta test --review); done'

bench:
  cargo bench --workspace

lint:
  cargo fmt --check
  cargo clippy --workspace --tests

# FIXME: Get rid of these -A flags
annoy:
  cargo clippy --workspace --tests -- -W clippy::nursery -W clippy::pedantic -W clippy::cargo -A clippy::missing_errors_doc -A clippy::cargo_common_metadata -A clippy::multiple_crate_versions

cov:
  cargo llvm-cov --workspace --branch --open

ci-cov:
  cargo llvm-cov --workspace --codecov --output-path codecov.json

check-wasm:
  sh -c 'for crate in `ls crates`; do (cd "crates/$crate" && cargo check --target wasm32-unknown-unknown); done'

ebnf:
  ebnf2railroad grammar/peptidoglycan.ebnf -t PGLang --write-style
  firefox grammar/peptidoglycan.html
