/******************************************************************************
 OSRA: Optical Structure Recognition

 This is a U.S. Government work (2007-2010) and is therefore not subject to
 copyright. However, portions of this work were obtained from a GPL or
 GPL-compatible source.
 Created by Igor Filippov, 2007-2010 (igorf@helix.nih.gov)

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
 St, Fifth Floor, Boston, MA 02110-1301, USA
 *****************************************************************************/

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h>
#include <openbabel/alias.h>

#include <sstream> // std:ostringstream
#include <iostream> // std::cerr

#include "osra_common.h"
#include "osra_openbabel.h"
#include "osra.h"
#include "mcdlutil.h"

using namespace OpenBabel;

#define HYDROGEN_ATOMIC_NUM     1
#define CARBON_ATOMIC_NUM       6
#define OXYGEN_ATOMIC_NUM       8
#define FLUORINE_ATOMIC_NUM     9
#define CHLORINE_ATOMIC_NUM     17
#define ARGON_ATOMIC_NUM        18
#define BROMINE_ATOMIC_NUM      35
#define IODINE_ATOMIC_NUM       53

// Function: create_atom()
//
// For the atom represented by its OCR'ed label create a new atom in the given molecule.
//
// Parameters:
//      mol - the current molecule
//      atom - atom to add to molecule (will be updated)
//      scale - scale factor / coordinates multiplier
//      superatom - superatom dictionary, that maps atom labels to corresponding SMILES
//      verbose - print debug information
//
// Returns:
//      0 in case when no bonds have been added to molecule or the index of the first added bond
//      
int create_atom(OBMol &mol, atom_t &atom, double scale, const map<string, string> &superatom, bool verbose)
{
  if (atom.label.empty() || atom.label == " ")
    {
      atom.anum = CARBON_ATOMIC_NUM;
    }
  else
    {
      // Lookup in superatom dictionary:
      map<string, string>::const_iterator it = superatom.find(atom.label);

      if (it != superatom.end())
        {
          // "superatom" case (e.g. "COOH")
          const string &smiles_superatom = it->second;

          OBConversion conv;
          OBMol superatom_mol;

          conv.SetInFormat("SMI");
          conv.ReadString(&superatom_mol, smiles_superatom);

          if (verbose)
            cout << "Considering superatom " << atom.label << "->" << smiles_superatom <<
                    " vector: " << atom.x * scale << "x" << -atom.y * scale << endl;

          // This is the index of first atom in superatom in molecule:
          atom.n = mol.NumAtoms() + 1;

          OBAtomIterator atom_iter;

          // Transfer all atoms from "superatom" molecule to current molecule.
          for (OBAtom *atom = superatom_mol.BeginAtom(atom_iter); atom; atom = superatom_mol.NextAtom(atom_iter))
            {
              if (verbose)
                cout << "Adding atom #" << mol.NumAtoms() + 1 << ", anum: " << atom->GetAtomicNum() << endl;

              OBAtom *a = mol.NewAtom();

              a->SetAtomicNum(atom->GetAtomicNum());
              a->SetFormalCharge(atom->GetFormalCharge());
            }

          // Correct first atom meta-info:
          OBAtom *first_superatom = mol.GetAtom(atom.n);

          first_superatom->SetVector(atom.x * scale, -atom.y * scale, 0);

          atom.anum = first_superatom->GetAtomicNum();

          if (atom.anum == 0)
            {
              // Unknown atom?
              AliasData* ad = new AliasData();
              ad->SetAlias(atom.label);
              ad->SetOrigin(external);
              first_superatom->SetData(ad);
            }

          int first_bond_index = mol.NumBonds();

          OBBondIterator bond_iter;

          // Transfer all bonds from "superatom" molecule to current molecule:
          for (OBBond *bond = superatom_mol.BeginBond(bond_iter); bond; bond = superatom_mol.NextBond(bond_iter))
            {
              if (verbose)
                cout << "Adding bond #" << mol.NumBonds() << " " << bond->GetBeginAtomIdx() + atom.n - 1 << "->"
                     << bond->GetEndAtomIdx() + atom.n - 1 << ", flags: " << bond->GetFlags() << endl;

              mol.AddBond(bond->GetBeginAtomIdx() + atom.n - 1, bond->GetEndAtomIdx() + atom.n - 1,
                          bond->GetBondOrder(), bond->GetFlags());
            }

          // Return the first bond index if at least one bond was added:
          return first_bond_index == mol.NumBonds() ? 0 : first_bond_index;
        }

      // If not found, lookup the atom number in periodic table of elements:
      atom.anum = etab.GetAtomicNum(atom.label.c_str());
    }

  atom.n = mol.NumAtoms() + 1;

  if (verbose)
    cout << "Creating atom #" << atom.n << " \"" << atom.label << "\", anum: " << atom.anum << endl;

  OBAtom *a = mol.NewAtom();

  a->SetAtomicNum(atom.anum);
  a->SetVector(atom.x * scale, -atom.y * scale, 0);

  if (atom.charge != 0)
    a->SetFormalCharge(atom.charge);

  if (atom.anum == 0)
    {
      // Unknown atom?
      AliasData* ad = new AliasData();
      ad->SetAlias(atom.label);
      ad->SetOrigin(external);
      a->SetData(ad);
    }

  return 0;
}

// Function: confidence_function()
//
// Calculates confidence estimate based on molecular counts provided by <create_molecule()>
//
// Parameters:
//      C_Count, N_Count, O_Count, F_Count, S_Count, Cl_Count, Br_Count - number of carbon, nitrogen, oxygen, fluorine, sulfur, chlorine, and bromine atoms
//      R_Count - number of recognized Markush atomic labels, such as R1, R2....
//      Xx_Count - number of unrecognized atomic labels from <osra.cpp::remove_small_terminal_bonds()>
//      num_rings - number of rings
//      num_aromatic - number of aromatic rings
//      num_fragments - number of fragments
//      Num_Rings - vector of counts for number of 3,4,5,6,7-member rings
//
// Returns:
//      confidence estimate
double confidence_function(int C_Count, int N_Count, int O_Count, int F_Count, int S_Count, int Cl_Count, int Br_Count,
                           int R_Count, int Xx_Count, int num_rings, int num_aromatic, int num_fragments, const vector<int> &Num_Rings)
{
  double confidence = 0.316030 //
                      - 0.016315 * C_Count //
                      + 0.034336 * N_Count //
                      + 0.066810 * O_Count //
                      + 0.035674 * F_Count //
                      + 0.065504 * S_Count //
                      + 0.04 * Cl_Count //
                      + 0.066811 * Br_Count //
                      + 0.01 * R_Count //
                      - 0.02 * Xx_Count //
                      - 0.212739 * num_rings //
                      + 0.071300 * num_aromatic //
                      + 0.329922 * Num_Rings[5] //
                      + 0.342865 * Num_Rings[6] //
                      - 0.037796 * num_fragments;

  return (confidence);
}

// Function: create_molecule()
//
// Converts vectors of atoms and bonds into a molecular object and calculates the molecule statistics.
// Note: this function changes the atoms!
//
// Parameters:
//      atom - vector of <atom_s> atoms
//      bond - vector of <bond_s> bonds
//      n_bond - total number of bonds
//      avg_bond_length - average bond length as measured from the image (to be included into output if provided)
//      molecule_statistics - the molecule statistics (returned to the caller)
//      generate_2D_coordinates - generate 2D coordinates for chemical groups
//      confidence - confidence score (returned to the caller if provided)
//      superatom - dictionary of superatom labels mapped to SMILES
//      verbose - print debug info
void create_molecule(OBMol &mol, vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond, double avg_bond_length, molecule_statistics_t &molecule_statistics,
                     bool generate_2D_coordinates, double * const confidence, const map<string, string> &superatom, bool verbose)
{
  string str;
  double scale = CC_BOND_LENGTH / avg_bond_length;
  vector<int> atomN, bondN;
  int anum;

  mol.SetDimension(2);
  mol.BeginModify();
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && i < MAX_ATOMS - 1 && bond[i].a < MAX_ATOMS - 1 && bond[i].b < MAX_ATOMS - 1)
      {
        atom_t* bond_atoms[] = { &atom[bond[i].a], &atom[bond[i].b] };

        for (int j = 0; j < 2; j++)
          if (bond_atoms[j]->n == 0)
            {
              if (create_atom(mol, *bond_atoms[j], scale, superatom, verbose) > 0)
                {
                  atomN.push_back(bond_atoms[j]->n);
                  bondN.push_back(mol.NumBonds());
                }
            }

        if (bond[i].hash)
          {
            if (verbose)
              cout << "Creating bond #" << mol.NumBonds() << " " << atom[bond[i].a].n << "<->"
                   << atom[bond[i].b].n << ", flags: " << OB_HASH_BOND << endl;

            if (atom[bond[i].a].anum == OXYGEN_ATOMIC_NUM || atom[bond[i].a].anum == HYDROGEN_ATOMIC_NUM || atom[bond[i].a].anum == FLUORINE_ATOMIC_NUM
                || atom[bond[i].a].anum == IODINE_ATOMIC_NUM || atom[bond[i].a].anum == CHLORINE_ATOMIC_NUM || atom[bond[i].a].anum == BROMINE_ATOMIC_NUM
                || atom[bond[i].a].anum == ARGON_ATOMIC_NUM || atom[bond[i].a].terminal)
              mol.AddBond(atom[bond[i].b].n, atom[bond[i].a].n, bond[i].type, OB_HASH_BOND);
            else
              mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_HASH_BOND);
          }
        else
          {
            int bond_flags = 0;

            if (bond[i].arom)
              bond_flags = OB_AROMATIC_BOND;
            else if (bond[i].wedge)
              bond_flags = OB_WEDGE_BOND;
            else if (bond[i].up)
              bond_flags = OB_TORUP_BOND;
            else if (bond[i].down)
              bond_flags = OB_TORDOWN_BOND;

            if (verbose)
              cout << "Creating bond #" << mol.NumBonds() << " " << atom[bond[i].a].n << "->"
                   << atom[bond[i].b].n << ", flags: " << bond_flags << endl;

            mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, bond_flags);
          }
      }
  mol.EndModify();

  mol.FindRingAtomsAndBonds();

  // Clear the counters of created OBAtom objects:
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists)
      {
        atom[bond[i].a].n = 0;
        atom[bond[i].b].n = 0;
      }

  // The logic below calculates the information both for molecule statistics and for confidence function:

  OBBondIterator bond_iter;

  // This block modifies the molecule:
  for (OBBond *b = mol.BeginBond(bond_iter); b; b = mol.NextBond(bond_iter))
    {
      if (b->IsInRing())
        {
          // Clear any indication of "/" and "\" double bond stereochemistry:
          b->UnsetUp();
          b->UnsetDown();
        }
      else
        // Clear all aromaticity information for the bond (affects "num_aromatic" variable below):
        b->UnsetAromatic();
    }

  vector<int> Num_Rings(8, 0); // number of rings of the given size (e.g. "Num_Rings[2]" = number of rings of size 2)
  int num_rings = 0, // total number of rings
      num_aromatic = 0; // total number of aromatic rings

  // Get the Smallest Set of Smallest Rings:
  vector<OBRing*> vr = mol.GetSSSR();

  for (vector<OBRing*>::iterator iter = vr.begin(); iter != vr.end(); iter++)
    {
      num_rings++;
      if ((*iter)->IsAromatic())
        num_aromatic++;
      if ((*iter)->Size() < 8)
        Num_Rings[(*iter)->Size()]++;
    }

  // Get a list of contiguous fragments sorted by size from largest to smallest:
  std::vector<std::vector<int> > cfl;
  mol.ContigFragList(cfl);

  molecule_statistics.rotors = mol.NumRotors();
  molecule_statistics.fragments = cfl.size();
  molecule_statistics.rings56 = Num_Rings[5] + Num_Rings[6];
  molecule_statistics.num_atoms = mol.NumAtoms();
  molecule_statistics.num_bonds = mol.NumBonds();

  if (confidence)
    {
      int C_Count = 0;
      int N_Count = 0;
      int O_Count = 0;
      int F_Count = 0;
      int S_Count = 0;
      int Cl_Count = 0;
      int Br_Count = 0;
      int R_Count = 0;
      int Xx_Count = 0;

      OBAtomIterator atom_iter;

      for (OBAtom *a = mol.BeginAtom(atom_iter); a; a = mol.NextAtom(atom_iter))
        {
          if (a->IsCarbon())
            C_Count++;
          else if (a->IsNitrogen())
            N_Count++;
          else if (a->IsOxygen())
            O_Count++;
          else if (a->IsSulfur())
            S_Count++;
          else if (a->GetAtomicNum() == FLUORINE_ATOMIC_NUM)
            F_Count++;
          else if (a->GetAtomicNum() == CHLORINE_ATOMIC_NUM)
            Cl_Count++;
          else if (a->GetAtomicNum() == BROMINE_ATOMIC_NUM)
            Br_Count++;
          else if (a->GetAtomicNum() == 0)
            {
              AliasData *ad;
              ad = (AliasData *) a->GetData(OBGenericDataType::SetData);
              if (ad != NULL && ad->GetAlias() != "Xx")
                R_Count++;
              else
                Xx_Count++;
            }
        }

      *confidence = confidence_function(C_Count, N_Count, O_Count, F_Count, S_Count, Cl_Count, Br_Count, R_Count,
                                        Xx_Count, num_rings, num_aromatic, molecule_statistics.fragments, Num_Rings);
    }

  if (generate_2D_coordinates)
    {
      if (verbose && !atomN.empty())
        cout << "Generating 2D coordinates for atoms " << atomN << " and bonds " << bondN << endl;

      for (unsigned int i = 0; i < atomN.size(); i++)
        {
          groupRedraw(&mol, bondN[i], atomN[i], true);
        }
    }
}

molecule_statistics_t caclulate_molecule_statistics(vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond, double avg_bond_length, const map<string, string> &superatom, bool verbose)
{
  molecule_statistics_t molecule_statistics;

  #pragma omp critical
  {
    OBMol mol;
    create_molecule(mol, atom, bond, n_bond, avg_bond_length, molecule_statistics, false, NULL, superatom, false);
    mol.Clear();
  }

  if (verbose)
    cout << "Molecule fragments: " << molecule_statistics.fragments << '.' << endl;

  return molecule_statistics;
}

const string get_formatted_structure(vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond, const string &format, const string &embedded_format, molecule_statistics_t &molecule_statistics,
                                     double &confidence, bool show_confidence, double avg_bond_length, double scaled_avg_bond_length, bool show_avg_bond_length, const int * const resolution,
                                     const int * const page, const box_t * const surrounding_box,
                                     const map<string, string> &superatom, bool verbose)
{
  ostringstream strstr;
  #pragma omp critical
  {
    OBMol mol;
    create_molecule(mol, atom, bond, n_bond, avg_bond_length, molecule_statistics, format == "sdf", &confidence, superatom, verbose);

    // Add hydrogens to the entire molecule to fill out implicit valence spots:
    mol.AddHydrogens(true, false); // polarOnly, correctForPh
    // Find all chiral atom centers:
    mol.FindChiralCenters();

    // Clear any indication of 2D "wedge" notation:
    OBBondIterator bond_iter;

    for (OBBond *b = mol.BeginBond(bond_iter); b; b = mol.NextBond(bond_iter))
      {
        if (!b->GetBeginAtom()->IsChiral() && !b->GetEndAtom()->IsChiral())
          {
            //b->UnsetHash();
            b->UnsetWedge();
          }
      }

    // Add single bonds based on atom proximity:
    mol.ConnectTheDots();
    // Copies each disconnected fragment as a separate OBMol:
    //mol.Separate();
    // Deletes all atoms except for the largest contiguous fragment:
    mol.StripSalts(MIN_A_COUNT);

    if (show_confidence)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Confidence_estimate");
        ostringstream cs;
        cs << confidence;
        label->SetValue(cs.str());
        mol.SetData(label);
      }

    if (show_avg_bond_length)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Average_bond_length");
        ostringstream cs;
        cs << scaled_avg_bond_length;
        label->SetValue(cs.str());
        mol.SetData(label);
      }

    if (resolution)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Resolution");
        ostringstream cs;
        cs << *resolution;
        label->SetValue(cs.str());
        mol.SetData(label);
      }

    if (page)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Page");
        ostringstream cs;
        cs << *page;
        label->SetValue(cs.str());
        mol.SetData(label);
      }

    if (surrounding_box)
      {
        OBPairData *label = new OBPairData;
        label->SetAttribute("Surrounding_box");
        ostringstream cs;
        cs << surrounding_box->x1 << 'x' << surrounding_box->y1 << '-' << surrounding_box->x2 << 'x' << surrounding_box->y2;
        label->SetValue(cs.str());
        mol.SetData(label);
      }

    if (!embedded_format.empty())
      {
        string value;

        OBConversion conv;

        conv.SetOutFormat(embedded_format.c_str());

        value = conv.WriteString(&mol);
        trim(value);

        if (!value.empty())
          {
            OBPairData *label = new OBPairData;
            label->SetAttribute(embedded_format);
            label->SetValue(value.c_str());
            mol.SetData(label);
            if (embedded_format == "inchi")
              {
                conv.SetOptions("K", conv.OUTOPTIONS);

                value = conv.WriteString(&mol);
                trim(value);

                label = new OBPairData;
                label->SetAttribute("InChI_key");
                label->SetValue(value.c_str());
                mol.SetData(label);
              }
          }
      }

    OBConversion conv;

    conv.SetOutFormat(format.c_str());
    conv.Read(&mol);

    strstr << conv.WriteString(&mol, true);

    if (format == "smi" || format == "can")
      {
        if (show_avg_bond_length)
          strstr << " " << scaled_avg_bond_length;
        if (resolution)
          strstr << " " << *resolution;
        if (show_confidence)
          strstr << " " << confidence;
        if (page)
          strstr << " " << *page;
        if (surrounding_box)
          strstr << " "<< surrounding_box->x1 << 'x' << surrounding_box->y1 << '-' << surrounding_box->x2 << 'x' << surrounding_box->y2;
      }

    strstr << endl;

    mol.Clear();
  }

  if (verbose)
    cout << "Structure length: " << strstr.str().length() << ", molecule fragments: " << molecule_statistics.fragments << '.' << endl;

  return (strstr.str());
}
