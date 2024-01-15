bully:
  cargo clippy --workspace -- -W clippy::pedantic -W clippy::nursery -W clippy::cargo

cov:
  cargo llvm-cov --open

ebnf:
  ebnf2railroad grammar/peptidoglycan.ebnf -t PGLang
  firefox grammar/peptidoglycan.html
