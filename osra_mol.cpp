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


string fix_atom_name(string s,int n)
{
  string r=s;
  //cout<<s<<endl;
  if (s.length()==1) r=toupper(s.at(0));
  if (s=="Ci" || s=="Cf" || s=="Cll") r="Cl";
  if (s=="H" && n>1) r="N";
  if (s=="HN" || s=="NH" || s=="M" || s=="Hm" || s=="MN" || s=="N2"
      || s=="NM" || s=="NH2" || s=="H2N" || s=="NHZ" || s=="HZN" || s=="NH3")   r="N";
  if (s=="OH" || s=="oH" || s=="Ho" || s=="HO" || s=="ol"
      || s=="On" || s=="on" || s=="no" || s=="nO") r="O";
  if (s=="Meo" || s=="oMe" || s=="oMg" || s=="omg" || s=="Mgo"
      || s=="leo" || s=="ohle" || s=="lleo" || s=="olllle" 
      || s=="OMe" || s=="OM8" || s=="OMo" || s=="OMB")   r="MeO";
  if (s=="FC")  r="CF";
  if (s=="NC")  r="CN";
  if ((s=="nBU") || (s=="neU") ||(s=="ngU")) r="nBu";
  if ((s=="Eto") || (s=="oEt") || (s=="Elo") || (s=="oEl")
      || s=="ElO") r="EtO";
  if ((s=="olgU") || (s=="oleU") || s=="OlBU") r="OiBu";
  if ((s=="npr") || (s=="llpll") || (s=="lpl") || (s=="npl")
      || s=="lPl") r="iPr";
  if ((s=="tBU") || (s=="BU") || (s=="llBU") || (s=="lBU")) r="tBu";
  if (s=="CooH" || s=="HooC" || s=="Co2H" || s=="CO2H") r="COOH";
  if (s=="AC") r="Ac";
  if (s=="ACo" || s=="opC" || s=="pcO") r="AcO";
  if (s=="Bl" || s=="el") r="Br";
  if (s=="CH3" || s=="H3C") r="C";
  if (s=="R" || s=="Rl" || s=="Rlo" || s=="R2" || s=="R3" || s=="Rg"
      || s=="R4" || s=="R5" || s=="R6" || s=="R7" || s=="R8" || s=="Z" 
      || s=="Y" || s=="2") r="X";
  if (s=="pl" || s=="nl") r="Ar";
  if (s=="oX") r="Ox";
  if (s=="NoZ" || s=="o2N" || s=="No2" || s=="No" || s=="O2N") r="NO2";
  if (s=="ph") r="Ph";
  if (s=="F3C") r="CF3";
  if (s=="F3Co") r="F3CN";
  if (s=="S3" || s=="Se") r="S";
  return(r);
}



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

string get_smiles(atom_t *atom, bond_t *bond, int n_bond, int &rotors, double &confidence)
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
   if (bond[i].exists) //&& atom[bond[i].a].label!="H" && atom[bond[i].b].label!="H")
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
	   b->UnsetHash();
	   b->UnsetWedge();
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
 for (unsigned int i=1;i<=mol.NumAtoms();i++)
   {
     OBAtom *a=mol.GetAtom(i);
     if (a->IsCarbon()) C_Count++;
     else if (a->IsNitrogen()) N_Count++;
     else if (a->IsOxygen()) O_Count++;
     else if (a->IsSulfur()) S_Count++;
     else if (a->GetAtomicNum()==9) F_Count++;
     else if (a->GetAtomicNum()==17) Cl_Count++;
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
 confidence=0.252867-0.018512*C_Count+0.032893*N_Count+0.065282*O_Count+0.034093*F_Count+0.069519*S_Count+0.202079*Cl_Count-0.183504*num_rings+0.083272*num_aromatic+0.324456*Num_Rings[3]+0.397955*Num_Rings[4]+0.306363*Num_Rings[5]+0.327634*Num_Rings[6]+0.325800*Num_Rings[7];

 str=conv.WriteString(&mol,true);
// conv.Write(&mol);
 for (int i=0;i<n_bond;i++)
   if (bond[i].exists) 
     {
       atom[bond[i].a].n=0;
       atom[bond[i].b].n=0;
     }
 rotors=mol.NumRotors();
 return(str);
}
