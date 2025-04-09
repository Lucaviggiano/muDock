#include <filesystem>
#include <iomanip>
#include <iostream>
#include <mudock/format.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
using namespace OpenBabel;
using namespace std;

int main(int argc, char** argv) {
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <input.pdb>" << endl;
    return 1;
  }

  std::filesystem::path pdb{argv[1]};

  const mudock::ob_mol_wrapper ob_mol = mudock::parser(pdb);

  for (auto atom_it = ob_mol->BeginAtoms(); atom_it < ob_mol->EndAtoms(); ++atom_it) {
    const auto atom = *atom_it;

    // Get basic info
    int idx          = atom->GetIdx();
    string elSymbol  = OBElements::GetSymbol(atom->GetAtomicNum());
    string atomType  = atom->GetType();
    double vdwRadius = OBElements::GetVdwRad(atom->GetAtomicNum());

    // Check donor/acceptor properties
    bool isDonor    = atom->IsHbondDonor();
    bool isAcceptor = atom->IsHbondAcceptor();

    // Get residue info if available
    string resInfo     = "N/A";
    OBResidue* residue = atom->GetResidue();
    if (residue) {
      resInfo = residue->GetName() + " " + to_string(residue->GetNum());
    }

    // Print atom info
    cout << left << setw(6) << idx << setw(5) << atomType << setw(5) << elSymbol << fixed << setprecision(2)
         << setw(8) << vdwRadius << setw(10) << (isDonor ? "Yes" : "No") << setw(10)
         << (isAcceptor ? "Yes" : "No") << setw(15) << resInfo << endl;
  }

  return 0;
}
