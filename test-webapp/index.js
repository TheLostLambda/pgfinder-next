import * as wasm from "../pkg/smithereens";

let pg = new wasm.Peptidoglycan("gm-AEJA");
console.log(`Monoisotopic Mass: ${pg.monoisotopic_mass()}`);
console.log(`Fragments:\n${wasm.pg_to_fragments(pg)}`);
