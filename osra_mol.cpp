/*********************************************************************
  OSRA: Optical Structure Recognition
  
  This is a U.S. Government work (year) and is therefore not subject to copyright.  
  However, portions of this work were obtained from a GPL or GPL-compatiple source.   
  Created by Igor Filippov, 2007-2008 (igorf@helix.nih.gov)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
  USA

*********************************************************************/

#include "osra.h"

#include "openbabel/mol.h"
#include "openbabel/obconversion.h" 
using namespace OpenBabel;


void addMeX(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
}


void addCF(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(9);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
}

void addCF3(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(9);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(9);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(9);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-3,(*n),1);
  mol->AddBond((*n)-2,(*n),1);
  mol->AddBond((*n)-1,(*n),1);
}

void addF3CN(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(9);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(9);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(9);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-3,(*n),1);
  mol->AddBond((*n)-2,(*n),1);
  mol->AddBond((*n)-1,(*n),1);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
}

void addPh(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-5,(*n)-4,5);
  mol->AddBond((*n)-4,(*n)-3,5);
  mol->AddBond((*n)-3,(*n)-2,5);
  mol->AddBond((*n)-2,(*n)-1,5);
  mol->AddBond((*n)-1,(*n),5);
  mol->AddBond((*n)-5,(*n),5);
}

void addNO2(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-2,(*n),2);
  mol->AddBond((*n)-1,(*n),2);
}

void addSO3H(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-3,(*n),1);
  mol->AddBond((*n)-2,(*n),2);
  mol->AddBond((*n)-1,(*n),2);
}


void addNC(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(7);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),3);
}

void addnBu(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
}

void addiPr(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
}

void addEtO(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
}

void addOiBu(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
 a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n)-2,1);
  mol->AddBond((*n)-2,(*n),1);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
}

void addtBu(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-3,(*n),1);
  mol->AddBond((*n)-2,(*n),1);
  mol->AddBond((*n)-1,(*n),1);
}

/*void addtBu(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n)-3,1);
  mol->AddBond((*n)-3,(*n),1);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  }*/

void addCOOH(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-2,(*n),2);
  mol->AddBond((*n)-1,(*n),1);
}

void addAc(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),2);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
}

void addAcO(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),2);
  mol->AddBond((*n)-2,(*n),1);
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
}

int getAnum(string s, OBMol *mol,int *n)
{
  if (s=="C") return(6);
  if (s=="N") return(7);
  if (s=="H") return(1);
  if (s=="O") return(8);
  if (s=="F") return(9);
  if (s=="P") return(15);
  if (s=="S") return(16);
  if (s=="I") return(53);
  if (s=="Cl") return(17);
  if (s=="Br") return(35);
  if (s=="X") return(0);
  if (s=="Ar") return(18);
  if (s=="MeO") 
    {
      addMeX(mol,n);
      return(8);
    }
  if (s=="CF")
    {
      addCF(mol,n);
      return(6);
    }
  if (s=="CF3")
    {
      addCF3(mol,n);
      return(6);
    }
  if (s=="F3CN")
    {
      addF3CN(mol,n);
      return(7);
    }
  if (s=="CN")
    {
      addNC(mol,n);
      return(6);
    }
  if (s=="nBu")
    {
      addnBu(mol,n);
      return(6);
    }
  if (s=="EtO")
    {
      addEtO(mol,n);
      return(8);
    }
  if (s=="OiBu")
    {
      addOiBu(mol,n);
      return(8);
    }
  if (s=="iPr")
    {
      addiPr(mol,n);
      return(6);
    }
  if (s=="tBu")
    {
      addtBu(mol,n);
      return(6);
    }
  if (s=="COOH")
    {
      addCOOH(mol,n);
      return(6);
    }
  if (s=="Ac")
    {
      addAc(mol,n);
      return(6);
    }
  if (s=="AcO")
    {
      addAcO(mol,n);
      return(6);
    }
  if (s=="NO2")
    {
      addNO2(mol,n);
      return(7);
    }
  if (s=="Ph")
    {
      addPh(mol,n);
      return(6);
    }
  if (s=="MeS") 
    {
      addMeX(mol,n);
      return(16);
    }
  if (s=="MeN") 
    {
      addMeX(mol,n);
      return(7);
    }
  if (s=="SO3H") 
    {
      addSO3H(mol,n);
      return(16);
    }
  return(6);
}

int getValency(string s)
{
  if (s=="C") return(4);
  if (s=="N") return(5);
  if (s=="H") return(1);
  if (s=="O") return(2);
  if (s=="F") return(1);
  if (s=="P") return(5);
  if (s=="S") return(6);
  if (s=="I") return(1);
  if (s=="Cl") return(1);
  if (s=="Br") return(1);
  if (s=="Ar") return(1);
  return(4);
}

int count_fragments(string input)
{
  int r=1;
  for(string::size_type i = input.find(".", 0); i != string::npos; i = input.find(".", i))
    {
      r++;
      i++;
    }
  return(r);
}

string get_smiles(atom_t *atom, bond_t *bond, int n_bond, int &rotors, 
		  double &confidence, int &num_fragments, int &r56)
{
 OBMol mol;
 OBAtom *a,*b;
 string str;
 //stringstream ss;
 //OBConversion conv(NULL,&ss);
 OBConversion conv;
 int n=1;
 int anum;

 conv.SetOutFormat("can");
 conv.Read(&mol);
 mol.SetDimension(2);
 for (int i=0;i<n_bond;i++)
   if (bond[i].exists) 
     {
       if (atom[bond[i].a].n==0)
	 {
	   a=mol.CreateAtom();
	   anum=getAnum(atom[bond[i].a].label,&mol,&n);
	   a->SetAtomicNum(anum);
	   if (atom[bond[i].a].charge!=0)
	     a->SetFormalCharge(atom[bond[i].a].charge);
	   //a->SetVector(atom[bond[i].a].x,atom[bond[i].a].y,0);
	   mol.AddAtom(*a);
	   atom[bond[i].a].n=n++;
	 }
       if (atom[bond[i].b].n==0)
	 {
	   b=mol.CreateAtom();
	   anum=getAnum(atom[bond[i].b].label,&mol,&n);
	   b->SetAtomicNum(anum);
	   if (atom[bond[i].b].charge!=0)
	     b->SetFormalCharge(atom[bond[i].b].charge);
	   //b->SetVector(atom[bond[i].b].x,atom[bond[i].b].y,0);
	   mol.AddAtom(*b);
	   atom[bond[i].b].n=n++;
	 }
       if (bond[i].arom)
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,5);
	 }
       else if (bond[i].hash)
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,bond[i].type,OB_HASH_BOND);
	 }
       else if (bond[i].wedge)
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,bond[i].type,OB_WEDGE_BOND);
	 }
       else if (bond[i].up)
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,bond[i].type,OB_TORUP_BOND);
	 }
       else if (bond[i].down)
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,bond[i].type,OB_TORDOWN_BOND);
	 }
       else
	 mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,bond[i].type);
     }
 mol.FindRingAtomsAndBonds();
 int j=0;
 for (int i=0;i<n_bond;i++)
   if (bond[i].exists)
     {
       OBBond *b=mol.GetBond(j);
       if (b!=NULL && b->IsInRing())
	 {
	   //b->UnsetHash();
	   //b->UnsetWedge();
	   b->UnsetUp();
	   b->UnsetDown();
	 }
       else if (b!=NULL && !b->IsInRing())
	 b->UnsetAromatic();
       j++;
     }
 int C_Count=0;
 int N_Count=0;
 int O_Count=0;
 int F_Count=0;
 int S_Count=0;
 int Cl_Count=0;
 int R_Count=0;
 for (unsigned int i=1;i<=mol.NumAtoms();i++)
   {
     OBAtom *a=mol.GetAtom(i);
     if (a->IsCarbon()) C_Count++;
     else if (a->IsNitrogen()) N_Count++;
     else if (a->IsOxygen()) O_Count++;
     else if (a->IsSulfur()) S_Count++;
     else if (a->GetAtomicNum()==9) F_Count++;
     else if (a->GetAtomicNum()==17) Cl_Count++;
     else if (a->GetAtomicNum()==0) R_Count++;
   }

 vector<OBRing*> vr=mol.GetSSSR();
 vector<OBRing*>::iterator iter;
 vector<int> Num_Rings(8,0);
 int num_rings=0,num_aromatic=0;
 for (iter = vr.begin();iter!=vr.end();iter++)
   {
     num_rings++;
     if ((*iter)->IsAromatic())
       num_aromatic++;
     if ((*iter)->Size()<8)
       Num_Rings[(*iter)->Size()]++;
   }
 rotors=mol.NumRotors();
 str=conv.WriteString(&mol,true);
 num_fragments=count_fragments(str);
 confidence=0.316030
   -0.016315*C_Count
   +0.034336*N_Count
   +0.066810*O_Count
   +0.035674*F_Count
   +0.065504*S_Count
   +0.198795*Cl_Count
   +0.1*R_Count
   -0.212739*num_rings
   +0.071300*num_aromatic
   +0.339289*Num_Rings[3]
   +0.422291*Num_Rings[4]
   +0.329922*Num_Rings[5]
   +0.342865*Num_Rings[6]
   +0.350747*Num_Rings[7]
   -0.037796*num_fragments;

 r56=Num_Rings[5]+Num_Rings[6];

 for (int i=0;i<n_bond;i++)
   if (bond[i].exists) 
     {
       atom[bond[i].a].n=0;
       atom[bond[i].b].n=0;
     }
 return(str);
}
