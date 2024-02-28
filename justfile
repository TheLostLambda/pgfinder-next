watch:
  watchexec -e rs,toml,kdl just test

test:
  cargo test --workspace

review:
  (cd polychem; cargo insta test --review)

bench:
  cargo bench --workspace

lint:
  cargo clippy --workspace --tests

# FIXME: Get rid of these -A flags
annoy:
  cargo clippy --workspace --tests -- -W clippy::nursery -W clippy::pedantic -W clippy::cargo -A clippy::missing_errors_doc -A clippy::cargo_common_metadata -A clippy::multiple_crate_versions

cov:
  cargo llvm-cov --workspace --open

ebnf:
  ebnf2railroad grammar/peptidoglycan.ebnf -t PGLang --write-style
  firefox grammar/peptidoglycan.html
