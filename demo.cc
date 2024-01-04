#include "include/smithereens.rs.h"
#include <iostream>

// I'm so sorry you have to see this...
int main() {
  auto pg = build_pg("gm-AEJA");
  auto mass = std::string(pg->monoisotopic_mass());
  std::cout << "Monoisotopic Mass: " << mass << std::endl;
  return 0;
}
