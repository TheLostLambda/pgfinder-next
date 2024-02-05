watch:
  watchexec -e rs,toml just test

test:
  cargo test --workspace

review:
  cargo insta test --workspace --review

lint:
  cargo clippy --workspace --tests

annoy:
  cargo clippy --workspace --tests -- -W clippy::nursery -W clippy::pedantic -W clippy::cargo

cov:
  cargo llvm-cov --workspace --open

ebnf:
  ebnf2railroad grammar/peptidoglycan.ebnf -t PGLang --write-style
  firefox grammar/peptidoglycan.html
