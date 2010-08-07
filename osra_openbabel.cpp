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

#include <map>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/builder.h> 
#include <openbabel/alias.h>

#include "osra.h"
#include "mcdlutil.h"

using namespace OpenBabel;

int abbreviation_to_mol(OBMol &mol, int &n, int &bondn, const string &smiles_superatom) {
	OBMol mol1;
	OBConversion conv;
	OBAtom *atom, *a1;
	OBBond *bond;

	conv.SetInFormat("SMI");
	conv.ReadString(&mol1, smiles_superatom);
	a1 = mol1.GetFirstAtom();
	unsigned int anum = a1->GetAtomicNum();


	int firstatom = a1->GetIdx();
	int prevatms = mol.NumAtoms();
	int numatms = mol1.NumAtoms();

	for (unsigned int i = mol1.NumAtoms(); i > 1; i--) {
		atom = mol1.GetAtom(i);
		if (atom != NULL) {
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

	for (unsigned int j = 0; j <= mol1.NumBonds(); j++) {
		bond = mol1.GetBond(j);
		if (bond != NULL) {
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

int getAnum(const string &s, OBMol &mol, int &n, int &bondn, const map<string, string> &superatom) {
	//mol_to_abbr();
  /*	if (s == "Xx" || s == "X" || s == "R" || s == "Y" || s == "Z" || s == "R1" || s == "R2" || s == "R3" || s == "R4"
			|| s == "R5" || s == "R6" || s == "R7" || s == "R8" || s == "R9" || s == "R10" || s == "Y2")
		return (0);
  */
	map<string, string>::const_iterator it = superatom.find(s);

	if (it != superatom.end()) {
		return (abbreviation_to_mol(mol, n, bondn, it->second));
	}

	int isotope;
	int anum = etab.GetAtomicNum(s.c_str(), isotope);

	if (anum != 0)
		return (anum);
	return (6);
}

const string get_smiles(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, int &rotors, double &confidence,
		int &num_fragments, int &r56, double avg, const string &format, int resolution, bool conf, bool guess,
		bool showpage, int page, const map<string, string> &superatom, bool showbond) {
	stringstream strstr;

#pragma omp critical
	{
		OBMol mol;
		string str;
		OBConversion conv;
		int n = 1;
		double scale = CC_BOND_LENGTH / avg;
		vector<int> atomB, atomN, bondB, bondN;
		int bondn = 0;
		int anum;

		conv.SetOutFormat(format.c_str());
		conv.Read(&mol);
		mol.SetDimension(2);

		mol.BeginModify();
		for (int i = 0; i < n_bond; i++)
			if (bond[i].exists && i < MAX_ATOMS - 1 && bond[i].a < MAX_ATOMS - 1 && bond[i].b < MAX_ATOMS - 1) {
				//cout << i << endl;
				if (atom[bond[i].a].n == 0) {
					int oldn = n;
					int oldbond = bondn;
					//cout << i << " " << bond[i].a << " " << atom[bond[i].a].label << endl;
					anum = getAnum(atom[bond[i].a].label, mol, n, bondn, superatom);
					if (oldn != n) {
						atomB.push_back(oldn);
						atomN.push_back(n - 1);
						bondB.push_back(oldbond);
						bondN.push_back(bondn);
						atom[bond[i].a].n = n - 1;
						OBAtom *a = mol.GetAtom(n - 1);
						a->SetVector(atom[bond[i].a].x * scale, -atom[bond[i].a].y * scale, 0);
					} else {
						OBAtom *a = mol.CreateAtom();
						a->SetAtomicNum(anum);
						if (atom[bond[i].a].charge != 0)
							a->SetFormalCharge(atom[bond[i].a].charge);
						a->SetVector(atom[bond[i].a].x * scale, -atom[bond[i].a].y * scale, 0);
						if (anum == 0) {
							AliasData* ad = new AliasData();
							ad->SetAlias(atom[bond[i].a].label);
							ad->SetOrigin(external);
							a->SetData(ad);
							ad->Expand(mol, anum); //Make chemically meaningful, if possible.
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
				if (atom[bond[i].b].n == 0) {
					int oldn = n;
					int oldbond = bondn;
					anum = getAnum(atom[bond[i].b].label, mol, n, bondn, superatom);
					if (oldn != n) {
						atomB.push_back(oldn);
						atomN.push_back(n - 1);
						bondB.push_back(oldbond);
						bondN.push_back(bondn);
						atom[bond[i].b].n = n - 1;
						OBAtom *b = mol.GetAtom(n - 1);
						b->SetVector(atom[bond[i].b].x * scale, -atom[bond[i].b].y * scale, 0);
					} else {
						OBAtom *b = mol.CreateAtom();
						b->SetAtomicNum(anum);
						if (atom[bond[i].b].charge != 0)
							b->SetFormalCharge(atom[bond[i].b].charge);
						b->SetVector(atom[bond[i].b].x * scale, -atom[bond[i].b].y * scale, 0);
						if (anum == 0) {
							AliasData* ad = new AliasData();
							ad->SetAlias(atom[bond[i].b].label);
							ad->SetOrigin(external);
							b->SetData(ad);
							ad->Expand(mol, anum); // Make chemically meaningful, if possible.
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

				if (bond[i].arom) {
					mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, 5);
					bondn++;
				} else if (bond[i].hash) {
					if (atom[bond[i].a].anum == 8 || atom[bond[i].a].anum == 1 || atom[bond[i].a].anum == 9
							|| atom[bond[i].a].anum == 53 || atom[bond[i].a].anum == 17 || atom[bond[i].a].anum == 35
							|| atom[bond[i].a].anum == 18 || atom[bond[i].a].terminal)
						mol.AddBond(atom[bond[i].b].n, atom[bond[i].a].n, bond[i].type, OB_HASH_BOND);
					else
						mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_HASH_BOND);
					bondn++;
				} else if (bond[i].wedge) {
					mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_WEDGE_BOND);
					bondn++;
				} else if (bond[i].up) {
					mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_TORUP_BOND);
					bondn++;
				} else if (bond[i].down) {
					mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type, OB_TORDOWN_BOND);
					bondn++;
				} else {
					mol.AddBond(atom[bond[i].a].n, atom[bond[i].b].n, bond[i].type);
					bondn++;
				}
			}
		mol.EndModify();

		mol.FindRingAtomsAndBonds();
		int num_double = 0;
		int num_triple = 0;
		for (unsigned int j = 0; j <= mol.NumBonds(); j++) {
			OBBond *b = mol.GetBond(j);
			if (b != NULL) {
				if (b->IsInRing()) {
					b->UnsetUp();
					b->UnsetDown();
				} else
					b->UnsetAromatic();
				if (b->IsDouble())
					num_double++;
				if (b->IsTriple())
					num_triple++;
			}
		}
		int C_Count = 0;
		int N_Count = 0;
		int O_Count = 0;
		int F_Count = 0;
		int S_Count = 0;
		int Cl_Count = 0;
		int Br_Count = 0;
		int R_Count = 0;
		int Xx_Count = 0;
		for (unsigned int i = 1; i <= mol.NumAtoms(); i++) {
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
			else if (a->GetAtomicNum() == 0) {
				AliasData *ad;
				ad = (AliasData *) a->GetData(OBGenericDataType::SetData);
				if (ad != NULL && ad->GetAlias() != "Xx")
					R_Count++;
				else
					Xx_Count++;
			}
		}

		vector<OBRing*> vr = mol.GetSSSR();
		vector<OBRing*>::iterator iter;
		vector<int> Num_Rings(8, 0);
		int num_rings = 0, num_aromatic = 0;
		for (iter = vr.begin(); iter != vr.end(); iter++) {
			num_rings++;
			if ((*iter)->IsAromatic())
				num_aromatic++;
			if ((*iter)->Size() < 8)
				Num_Rings[(*iter)->Size()]++;
		}
		rotors = mol.NumRotors();

		std::vector<std::vector<int> > cfl;
		mol.ContigFragList(cfl);
		num_fragments = cfl.size();

		confidence = confidence_function(C_Count, N_Count, O_Count, F_Count, S_Count, Cl_Count, Br_Count, R_Count,
				Xx_Count, num_rings, num_aromatic, num_fragments, Num_Rings, num_double, num_triple);

		r56 = Num_Rings[5] + Num_Rings[6];

		if (conf) {
			OBPairData *label = new OBPairData;
			label->SetAttribute("Confidence_Estimate");
			stringstream cs;
			cs << confidence;
			label->SetValue(cs.str());
			//label->SetOrigin(userInput); // set by user, not by Open Babel
			mol.SetData(label);
		}
		if (guess) {
			OBPairData *label = new OBPairData;
			label->SetAttribute("Resolution");
			stringstream cs;
			cs << resolution;
			label->SetValue(cs.str());
			//label->SetOrigin(userInput); // set by user, not by Open Babel
			mol.SetData(label);
		}

		if (showpage) {
			OBPairData *label = new OBPairData;
			label->SetAttribute("Page");
			stringstream cs;
			cs << page;
			label->SetValue(cs.str());
			//label->SetOrigin(userInput); // set by user, not by Open Babel
			mol.SetData(label);
		}

		if (showbond) {
			OBPairData *label = new OBPairData;
			label->SetAttribute("Average_bond_length");
			stringstream cs;
			cs << avg;
			label->SetValue(cs.str());
			//label->SetOrigin(userInput); // set by user, not by Open Babel
			mol.SetData(label);
		}

		if (format == "sdf") {
			for (unsigned int i = 0; i < atomN.size(); i++) {
				//groupRedrawBeginEnd(&mol, atomB[i], atomN[i], bondB[i], bondN[i]);
				groupRedraw(&mol, bondN[i], atomN[i], true);
			}
		}

		mol.AddHydrogens(true, false); // polarOnly, correctForPh
		mol.FindChiralCenters();

		for (unsigned int j = 0; j <= mol.NumBonds(); j++) {
			OBBond *b = mol.GetBond(j);
			if (b != NULL) {
				if (!b->GetBeginAtom()->IsChiral()) {
					b->UnsetHash();
					b->UnsetWedge();
				}
			}
		}

		mol.ConnectTheDots();
		//mol.Separate();
		mol.StripSalts(MIN_A_COUNT);

		if (format != "empty")
			str = conv.WriteString(&mol, true);
		else
			str = "";

		for (int i = 0; i < n_bond; i++)
			if (bond[i].exists) {
				atom[bond[i].a].n = 0;
				atom[bond[i].b].n = 0;
			}

		strstr << str;
		if (format == "smi" || format == "can") {
			if (guess)
				strstr << " " << resolution;
			if (conf)
				strstr << " " << confidence;
			if (showpage)
				strstr << " " << page;
			if (showbond)
				strstr << " " << avg;
		}
		strstr << endl;

		mol.Clear();
	}

	return (strstr.str());
}
