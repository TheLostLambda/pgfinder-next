lint:
  cargo clippy --tests

cov:
  cargo llvm-cov --open

ebnf:
  ebnf2railroad grammar/peptidoglycan.ebnf -t PGLang --write-style
  firefox grammar/peptidoglycan.html
