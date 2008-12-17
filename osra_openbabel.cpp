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
#include <openbabel/builder.h> 
#include "mcdlutil.h"
using namespace OpenBabel;
//#include <openbabel/generic.h>

void addMeX(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addOR(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(0);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addCF(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(9);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addCF3(OBMol *mol,int *n,int *bondn)
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
  (*bondn)++;
  mol->AddBond((*n)-2,(*n),1);
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addF3CN(OBMol *mol,int *n,int *bondn)
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
  (*bondn)++;
  mol->AddBond((*n)-2,(*n),1);
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addPh(OBMol *mol,int *n,int *bondn)
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
  (*bondn)++;
  mol->AddBond((*n)-4,(*n)-3,5);
  (*bondn)++;
  mol->AddBond((*n)-3,(*n)-2,5);
  (*bondn)++;
  mol->AddBond((*n)-2,(*n)-1,5);
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),5);
  (*bondn)++;
  mol->AddBond((*n)-5,(*n),5);
  (*bondn)++;
}

void addBzO(OBMol *mol,int *n,int *bondn)
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
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  mol->AddBond((*n)-5,(*n)-4,5);
  (*bondn)++;
  mol->AddBond((*n)-4,(*n)-3,5);
  (*bondn)++;
  mol->AddBond((*n)-3,(*n)-2,5);
  (*bondn)++;
  mol->AddBond((*n)-2,(*n)-1,5);
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),5);
  (*bondn)++;
  mol->AddBond((*n)-5,(*n),5);
  (*bondn)++;
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addTHPO(OBMol *mol,int *n,int *bondn)
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
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  mol->AddBond((*n)-5,(*n)-4,1);
  (*bondn)++;
  mol->AddBond((*n)-4,(*n)-3,1);
  (*bondn)++;
  mol->AddBond((*n)-3,(*n)-2,1);
  (*bondn)++;
  mol->AddBond((*n)-2,(*n)-1,1);
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
  mol->AddBond((*n)-5,(*n),1);
  (*bondn)++;
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addNO2(OBMol *mol,int *n,int *bondn)
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
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),2);
  (*bondn)++;
}

void addNO(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),2);
  (*bondn)++;
}


void addNOHCH3(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-2,(*n),1);
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addSO3H(OBMol *mol,int *n,int *bondn)
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
  (*bondn)++;
  mol->AddBond((*n)-2,(*n),2);
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),2);
  (*bondn)++;
}


void addNC(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(7);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),3);
  (*bondn)++;
}

void addnBu(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addiPr(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addEtO(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addOiBu(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-2,(*n),1);
  (*bondn)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addtBu(OBMol *mol,int *n,int *bondn)
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
  (*bondn)++;
  mol->AddBond((*n)-2,(*n),1);
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}


void addCOOH(OBMol *mol,int *n,int *bondn)
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
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addAc(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-2,(*n),2);
  (*bondn)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

void addAcO(OBMol *mol,int *n,int *bondn)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(8);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),2);
  (*bondn)++;
  mol->AddBond((*n)-2,(*n),1);
  (*bondn)++;
  a=mol->CreateAtom();
  a->SetAtomicNum(6);
  mol->AddAtom(*a);
  (*n)++;
  mol->AddBond((*n)-1,(*n),1);
  (*bondn)++;
}

int getAnum(string s, OBMol *mol,int *n, int *bondn)
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
      addMeX(mol,n,bondn);
      return(8);
    }
  if (s=="CF")
    {
      addCF(mol,n,bondn);
      return(6);
    }
  if (s=="CF3")
    {
      addCF3(mol,n,bondn);
      return(6);
    }
  if (s=="F3CN")
    {
      addF3CN(mol,n,bondn);
      return(7);
    }
  if (s=="CN")
    {
      addNC(mol,n,bondn);
      return(6);
    }
  if (s=="nBu")
    {
      addnBu(mol,n,bondn);
      return(6);
    }
  if (s=="EtO")
    {
      addEtO(mol,n,bondn);
      return(8);
    }
  if (s=="OiBu")
    {
      addOiBu(mol,n,bondn);
      return(8);
    }
  if (s=="iPr")
    {
      addiPr(mol,n,bondn);
      return(6);
    }
  if (s=="tBu")
    {
      addtBu(mol,n,bondn);
      return(6);
    }
  if (s=="COOH")
    {
      addCOOH(mol,n,bondn);
      return(6);
    }
  if (s=="Ac")
    {
      addAc(mol,n,bondn);
      return(6);
    }
  if (s=="AcO")
    {
      addAcO(mol,n,bondn);
      return(8);
    }
  if (s=="NO2")
    {
      addNO2(mol,n,bondn);
      return(7);
    }
  if (s=="NO")
    {
      addNO(mol,n,bondn);
      return(7);
    }
  if (s=="Ph")
    {
      addPh(mol,n,bondn);
      return(6);
    }
  if (s=="MeS") 
    {
      addMeX(mol,n,bondn);
      return(16);
    }
  if (s=="MeN") 
    {
      addMeX(mol,n,bondn);
      return(7);
    }
  if (s=="SO3H") 
    {
      addSO3H(mol,n,bondn);
      return(16);
    }
  if (s=="OR") 
    {
      addOR(mol,n,bondn);
      return(8);
    }
  if (s=="BzO")
    {
      addBzO(mol,n,bondn);
      return(8);
    }
  if (s=="N(OH)CH3")
    {
      addNOHCH3(mol,n,bondn);
      return(7);
    }
  if (s=="THPO")
    {
      addTHPO(mol,n,bondn);
      return(8);
    }
  return(6);
}



string get_smiles(atom_t *atom, bond_t *bond, int n_bond, int &rotors, 
		  double &confidence, int &num_fragments, int &r56, double avg,
		  string format,int resolution,bool conf, bool guess)
{
 OBMol mol;
 OBAtom *a,*b;
 string str;
 OBConversion conv;
 int n=1;
 int anum;
 double scale=CC_BOND_LENGTH/avg;
 vector<int> atomB,atomN,bondB,bondN;
 int bondn=0;

 conv.SetOutFormat(format.c_str());
 conv.Read(&mol);
 mol.SetDimension(2);
 mol.BeginModify();
 for (int i=0;i<n_bond;i++)
   if (bond[i].exists) 
     {
       if (atom[bond[i].a].n==0)
	 {
	   a=mol.CreateAtom();
	   int oldn=n;
	   int oldbond=bondn;
	   anum=getAnum(atom[bond[i].a].label,&mol,&n,&bondn);
	   if (oldn!=n) 
	     {
	       atomB.push_back(oldn);
	       atomN.push_back(n);
	       bondB.push_back(oldbond);
	       bondN.push_back(bondn);
	     }
	   a->SetAtomicNum(anum);
	   if (atom[bond[i].a].charge!=0)
	     a->SetFormalCharge(atom[bond[i].a].charge);
	   a->SetVector(atom[bond[i].a].x*scale,-atom[bond[i].a].y*scale,0);
	   mol.AddAtom(*a);
	   atom[bond[i].a].n=n++;
	 }
       if (atom[bond[i].b].n==0)
	 {
	   b=mol.CreateAtom();
	   int oldn=n;
	   int oldbond=bondn;
	   anum=getAnum(atom[bond[i].b].label,&mol,&n,&bondn);
	   if (oldn!=n) 
	     {
	       atomB.push_back(oldn);
	       atomN.push_back(n);
	       bondB.push_back(oldbond);
	       bondN.push_back(bondn);
	     }
	   b->SetAtomicNum(anum);
	   if (atom[bond[i].b].charge!=0)
	     b->SetFormalCharge(atom[bond[i].b].charge);
	   b->SetVector(atom[bond[i].b].x*scale,-atom[bond[i].b].y*scale,0);
	   mol.AddAtom(*b);
	   atom[bond[i].b].n=n++;
	 }

       if (bond[i].arom)
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,5);
	   bondn++;
	 }
       else if (bond[i].hash)
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,bond[i].type,OB_HASH_BOND);
	   bondn++;
	 }
       else if (bond[i].wedge)
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,bond[i].type,OB_WEDGE_BOND);
	   bondn++;
	 }
       else if (bond[i].up)
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,bond[i].type,OB_TORUP_BOND);
	   bondn++;
	 }
       else if (bond[i].down)
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,bond[i].type,OB_TORDOWN_BOND);
	   bondn++;
         }
       else
	 {
	   mol.AddBond(atom[bond[i].a].n,atom[bond[i].b].n,bond[i].type);
	   bondn++;
	 }
     }
 mol.EndModify();
 mol.FindRingAtomsAndBonds();
 for (unsigned int j=1;j<=mol.NumBonds();j++)
     {
       OBBond *b=mol.GetBond(j);
       if (b!=NULL)
       {
        if (b->IsInRing())
	 {
	   b->UnsetUp();
	   b->UnsetDown();
	 }
	 else
	   b->UnsetAromatic();
       }
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
 
 std::vector< std::vector< int > > cfl;
 mol.ContigFragList(cfl);
 num_fragments=cfl.size();

 confidence=confidence_function(C_Count,N_Count,O_Count,F_Count,S_Count,Cl_Count,
				num_rings,num_aromatic,num_fragments,&Num_Rings);

 r56=Num_Rings[5]+Num_Rings[6];

 if (conf)
   {
     OBPairData *label = new OBPairData;
     label->SetAttribute("Confidence_Estimate");
     stringstream cs;
     cs<<confidence;
     label->SetValue(cs.str());
     //label->SetOrigin(userInput); // set by user, not by Open Babel
     mol.SetData(label);
   }
 if (guess)
   {
     OBPairData *label = new OBPairData;
     label->SetAttribute("Resolution");
     stringstream cs;
     cs<<resolution;
     label->SetValue(cs.str());
     //label->SetOrigin(userInput); // set by user, not by Open Babel
     mol.SetData(label);
   }

if (format=="sdf")
   {
     for (unsigned int i=0;i<atomN.size();i++)
       {
	 //	 groupRedrawBeginEnd(&mol,atomB[i],atomN[i],bondB[i],bondN[i]);
         groupRedraw(&mol,bondN[i],atomN[i],true);
       }
   }

 mol.AddHydrogens(true,false); // polarOnly, correctForPh

 mol.FindChiralCenters();
 for (unsigned int j=1;j<=mol.NumBonds();j++)
     {
       OBBond *b=mol.GetBond(j);
       if (b!=NULL)
       {
	 if (!b->GetBeginAtom()->IsChiral())
	   {
	     b->UnsetHash(); 
	     b->UnsetWedge();
	   }
       }
     }

 mol.ConnectTheDots();
 // mol.Separate();

 str=conv.WriteString(&mol,true);

 for (int i=0;i<n_bond;i++)
   if (bond[i].exists) 
     {
       atom[bond[i].a].n=0;
       atom[bond[i].b].n=0;
     }
 stringstream strstr;
 strstr<<str;
 if (format=="smi" || format=="can")
   {
     if (guess)
       strstr<<" "<<resolution;
     if (conf)
       strstr<<" "<<confidence;
   }
 strstr<<endl;
 return(strstr.str());
}
