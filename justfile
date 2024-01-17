lint:
  cargo clippy

cov:
  cargo llvm-cov --open

ebnf:
  ebnf2railroad grammar/peptidoglycan.ebnf -t PGLang
  firefox grammar/peptidoglycan.html
