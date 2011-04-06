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

#include "osra_openbabel.h"
#include "osra.h"
#include "mcdlutil.h"

using namespace OpenBabel;

int abbreviation_to_mol(OBMol &mol, int &n, int &bondn, const string &smiles_superatom)
{
  OBMol mol1;
  OBConversion conv;
  OBAtom *atom, *a1;
  OBBond *bond;

  conv.SetInFormat("SMI");
  conv.ReadString(&mol1, smiles_superatom);

  a1 = mol1.GetFirstAtom();

  if (a1 == NULL)
    {
      cerr << "Unable to parse the SMILES " << smiles_superatom << '.' << endl << "That means:" << endl
           << "(a) You have no /usr/lib/openbabel/x.x.x/smilesformat.so library installed. Install the format libraries / check http://openbabel.org/docs/dev/Installation/install.html#environment-variables" << endl
           << "(b) The format libraries are installed, but do not correspond to /usr/lib/libopenbabel.so.y.y.y. Check they correspond to the same OpenBabel version." << endl
           << "OSRA will produce unpredictable/wrong results." << endl;

      return 0;
    }

  unsigned int anum = a1->GetAtomicNum();
  int firstatom = a1->GetIdx();
  int prevatms = mol.NumAtoms();
  int numatms = mol1.NumAtoms();

  for (unsigned int i = mol1.NumAtoms(); i >= 1; i--)
    {
      atom = mol1.GetAtom(i);
      if (atom != NULL)
        {
          OBAtom *a = mol.CreateAtom();
          a->SetAtomicNum(atom->GetAtomicNum());
          a->SetFormalCharge(atom->GetFormalCharge());
          a->SetIdx(mol.NumAtoms() + 1);
          // This operation copies the given atom into new one before addition,
          // so the caller is still responsible to free the memory:
          mol.AddAtom(*a);
          delete a;
          n++;
        }
    }

  for (unsigned int j = 0; j <= mol1.NumBonds(); j++)
    {
      bond = mol1.GetBond(j);
      if (bond != NULL)
        {
          int b1 = (numatms - bond->GetBeginAtomIdx() + 1) - firstatom + 1 + prevatms;
          int b2 = (numatms - bond->GetEndAtomIdx() + 1) - firstatom + 1 + prevatms;
          mol.AddBond(b1, b2, bond->GetBO(), bond->GetFlags());
          bondn++;
        }
    }

  return (anum);
}

/*
void mol_to_abbr() {
	OBMol mol1, mol;
	OBConversion conv;
	OBAtom *atom, *a, *a1;
	OBBond *bond;
	int n = 1, bondn = 0;

	conv.SetOutFormat("SMI");

	addTHPO(&mol1, &n, &bondn);
	a = mol1.CreateAtom();
	a->SetAtomicNum(8);

	// This operation copies the given atom into new one before addition,
	// so the caller is still responsible to free the memory:
	mol1.AddAtom(*a);
	delete a;

	a1 = mol1.GetFirstAtom();

	int firstatom = a1->GetIdx();
	int prevatms = mol.NumAtoms();
	int numatms = mol1.NumAtoms();

	for (unsigned int i = mol1.NumAtoms(); i >= 1; i--) {
		atom = mol1.GetAtom(i);
		if (atom != NULL) {
			a = mol.CreateAtom();
			a->SetAtomicNum(atom->GetAtomicNum());
			a->SetFormalCharge(atom->GetFormalCharge());
			a->SetIdx(mol.NumAtoms() + 1);
			// This operation copies the given atom into new one before addition,
			// so the caller is still responsible to free the memory:
			mol.AddAtom(*a);
			delete a;
			n++;
		}
	}

	for (unsigned int j = 0; j <= mol1.NumBonds(); j++) {
		bond = mol1.GetBond(j);
		if (bond != NULL) {
			int b1 = (numatms - bond->GetBeginAtomIdx() + 1) - firstatom + 1 + prevatms;
			int b2 = (numatms - bond->GetEndAtomIdx() + 1) - firstatom + 1 + prevatms;
			mol.AddBond(b1, b2, bond->GetBO(), bond->GetFlags());
			bondn++;
		}
	}
	cout << conv.WriteString(&mol, true) << endl;
	exit(0);
}
*/

int get_atomic_num(const string &s, OBMol &mol, int &n, int &bondn, const map<string, string> &superatom)
{

  if (s.empty() || s == " ") return 6;
  map<string, string>::const_iterator it = superatom.find(s);

  if (it != superatom.end())
    {
      return (abbreviation_to_mol(mol, n, bondn, it->second));
    }

  int isotope;
  int anum = etab.GetAtomicNum(s.c_str(), isotope);
//  if (anum == 0) return(6);
  return (anum);
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
//
// Returns:
//      calculated molecule statistics
void create_molecule(OBMol &mol, vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond, double avg_bond_length, molecule_statistics_t &molecule_statistics,
                     bool generate_2D_coordinates, double * const confidence, const map<string, string> &superatom)
{
  string str;
  int n = 1;
  double scale = CC_BOND_LENGTH / avg_bond_length;
  vector<int> atomN, bondN;
  int bondn = 0;
  int anum;

  mol.SetDimension(2);

  mol.BeginModify();
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && i < MAX_ATOMS - 1 && bond[i].a < MAX_ATOMS - 1 && bond[i].b < MAX_ATOMS - 1)
      {
        if (atom[bond[i].a].n == 0)
          {
            int oldn = n;
            anum = get_atomic_num(atom[bond[i].a].label, mol, n, bondn, superatom);
            if (oldn != n)
              {
                if (n != oldn+1)
                  {
                    atomN.push_back(n - 1);
                    bondN.push_back(bondn);
                  }
                atom[bond[i].a].n = n - 1;
                OBAtom *a = mol.GetAtom(n - 1);
                a->SetVector(atom[bond[i].a].x * scale, -atom[bond[i].a].y * scale, 0);
                if (anum == 0)
                  {
                    AliasData* ad = new AliasData();
                    ad->SetAlias(atom[bond[i].a].label);
                    ad->SetOrigin(external);
                    a->SetData(ad);
//                      ad->Expand(mol, anum);
                  }
              }
            else
              {
                OBAtom *a = mol.CreateAtom();
                a->SetAtomicNum(anum);
                if (atom[bond[i].a].charge != 0)
                  a->SetFormalCharge(atom[bond[i].a].charge);
                a->SetVector(atom[bond[i].a].x * scale, -atom[bond[i].a].y * scale, 0);
                if (anum == 0)
                  {
                    AliasData* ad = new AliasData();
                    ad->SetAlias(atom[bond[i].a].label);
                    ad->SetOrigin(external);
                    a->SetData(ad);
//                      ad->Expand(mol, anum); //Make chemically meaningful, if possible.
                  }
                // This operation copies the given atom into new one before addition,
                // so the caller is still responsible to free the memory:
                mol.AddAtom(*a);
                delete a;
                atom[bond[i].a].n = n;
                n++;
              }
            atom[bond[i].a].anum = anum;
          }
        if (atom[bond[i].b].n == 0)
          {
            int oldn = n;
            anum = get_atomic_num(atom[bond[i].b].label, mol, n, bondn, superatom);
            if (oldn != n)
              {
                if (n != oldn+1)
                  {
                    atomN.push_back(n - 1);
                    bondN.push_back(bondn);
                  }
                atom[bond[i].b].n = n - 1;
                OBAtom *b = mol.GetAtom(n - 1);
                b->SetVector(atom[bond[i].b].x * scale, -atom[bond[i].b].y * scale, 0);
                if (anum == 0)
                  {
                    AliasData* ad = new AliasData();
                    ad->SetAlias(atom[bond[i].b].label);
                    ad->SetOrigin(external);
                    b->SetData(ad);
//                      ad->Expand(mol, anum);
                  }
              }
            else
              {
                OBAtom *b = mol.CreateAtom();
                b->SetAtomicNum(anum);
                if (atom[bond[i].b].charge != 0)
                  b->SetFormalCharge(atom[bond[i].b].charge);
                b->SetVector(atom[bond[i].b].x * scale, -atom[bond[i].b].y * scale, 0);
                if (anum == 0)
                  {
                    AliasData* ad = new AliasData();
                    ad->SetAlias(atom[bond[i].b].label);
                    ad->SetOrigin(external);
                    b->SetData(ad);
//                      ad->Expand(mol, anum); // Make chemically meaningful, if possible.
                  }
                // This operation copies the given atom into new one before addition,
                // so the caller is still responsible to free the memory:
                mol.AddAtom(*b);
                delete b;
                atom[bond[i].b].n = n;
                n++;
              }
            atom[bond[i].b].anum = anum;
          }

        if (bond[i].arom)
          {
            mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, 5);
            bondn++;
          }
        else if (bond[i].hash)
          {
            if (atom[bond[i].a].anum == 8 || atom[bond[i].a].anum == 1 || atom[bond[i].a].anum == 9
                || atom[bond[i].a].anum == 53 || atom[bond[i].a].anum == 17 || atom[bond[i].a].anum == 35
                || atom[bond[i].a].anum == 18 || atom[bond[i].a].terminal)
              mol.AddBond(atom[bond[i].b].n, atom[bond[i].a].n, bond[i].type, OB_HASH_BOND);
            else
              mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_HASH_BOND);
            bondn++;
          }
        else if (bond[i].wedge)
          {
            mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_WEDGE_BOND);
            bondn++;
          }
        else if (bond[i].up)
          {
            mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_TORUP_BOND);
            bondn++;
          }
        else if (bond[i].down)
          {
            mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_TORDOWN_BOND);
            bondn++;
          }
        else
          {
            mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type);
            bondn++;
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

  // This block modifies the molecule:
  for (unsigned int j = 0; j <= mol.NumBonds(); j++)
    {
      OBBond *b = mol.GetBond(j);
      if (b != NULL)
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

      for (unsigned int i = 1; i <= mol.NumAtoms(); i++)
        {
          OBAtom *a = mol.GetAtom(i);
          if (a->IsCarbon())
            C_Count++;
          else if (a->IsNitrogen())
            N_Count++;
          else if (a->IsOxygen())
            O_Count++;
          else if (a->IsSulfur())
            S_Count++;
          else if (a->GetAtomicNum() == 9)
            F_Count++;
          else if (a->GetAtomicNum() == 17)
            Cl_Count++;
          else if (a->GetAtomicNum() == 35)
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
      for (unsigned int i = 0; i < atomN.size(); i++)
        {
          groupRedraw(&mol, bondN[i], atomN[i], true);
        }
    }
}

molecule_statistics_t caclulate_molecule_statistics(vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond, double avg_bond_length, const map<string, string> &superatom)
{
  molecule_statistics_t molecule_statistics;
  #pragma omp critical
  {
    OBMol mol;
    create_molecule(mol, atom, bond, n_bond, avg_bond_length, molecule_statistics, false, NULL, superatom);
    mol.Clear();
  }
  return molecule_statistics;
}

const string get_formatted_structure(vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond, const string &format, const string &embedded_format, molecule_statistics_t &molecule_statistics,
                                     double &confidence, bool show_confidence, double avg_bond_length, double scaled_avg_bond_length, bool show_avg_bond_length, const int * const resolution,
                                     const int * const page, const box_t * const surrounding_box,
                                     const map<string, string> &superatom)
{
  ostringstream strstr;
  #pragma omp critical
  {
    OBMol mol;
    create_molecule(mol, atom, bond, n_bond, avg_bond_length, molecule_statistics, format == "sdf", &confidence, superatom);

    // Add hydrogens to the entire molecule to fill out implicit valence spots:
    mol.AddHydrogens(true, false); // polarOnly, correctForPh
    // Find all chiral atom centers:
    mol.FindChiralCenters();

    // Clear any indication of 2D "wedge" notation:
    for (unsigned int j = 0; j <= mol.NumBonds(); j++)
      {
        OBBond *b = mol.GetBond(j);
        if (b != NULL && !b->GetBeginAtom()->IsChiral() && !b->GetEndAtom()->IsChiral())
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
        OBConversion conv;
        string value;

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
          strstr << " "<<surrounding_box->x1 << 'x' << surrounding_box->y1 << '-' << surrounding_box->x2 << 'x' << surrounding_box->y2;
      }

    strstr << endl;

    mol.Clear();
  }
  return (strstr.str());
}
