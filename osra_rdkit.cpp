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

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>
using namespace RDKit;



void addMeX(RWMol *mol,int *n)
{
  Atom *a=new Atom(6);
  mol->addAtom(a);
  (*n)++;
  mol->addBond((*n)-1,(*n),Bond::SINGLE);
}
/*
void addOR(OBMol *mol,int *n)
{
  OBAtom *a;
  a=mol->CreateAtom();
  a->SetAtomicNum(0);
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
  if (s=="OR") 
    {
      addOR(mol,n);
      return(8);
    }
  return(6);
}
*/

int getAnum(string s,RWMol *mol,int *n)
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
  return(6);
}

string get_smiles(atom_t *atom, int real_atoms,bond_t *bond, int n_bond, int &rotors, 
		  double &confidence, int &num_fragments, int &r56)
{
  RWMol *mol=new RWMol();
  int n=0,bondid=0;
  int anum;
  Conformer *conf = new Conformer(real_atoms);	
  for (int i=0;i<n_bond;i++)
    if (bond[i].exists) 
      {
	if (atom[bond[i].a].n==0)
	  {
  	    RDGeom::Point3D pos;
  	    pos.x=atom[bond[i].a].x;
  	    pos.y=atom[bond[i].a].y;
  	    pos.z=0;
	    anum=getAnum(atom[bond[i].a].label,mol,&n);
	    Atom *a=new Atom(anum);
	    if (atom[bond[i].a].charge!=0)
	      a->setFormalCharge(atom[bond[i].a].charge);
	    unsigned int aid=mol->addAtom(a);
	    conf->setAtomPos(aid, pos);
	    atom[bond[i].a].n=n++;
	  }
	if (atom[bond[i].b].n==0)
	  {
	    RDGeom::Point3D pos;
	    pos.x=atom[bond[i].b].x;
	    pos.y=atom[bond[i].b].y;
	    pos.z=0;
	    anum=getAnum(atom[bond[i].b].label,mol,&n);
	    Atom *b=new Atom(anum);
	    if (atom[bond[i].b].charge!=0)
	      b->setFormalCharge(atom[bond[i].b].charge);
	    unsigned int aid=mol->addAtom(b);
	    conf->setAtomPos(aid, pos);
	    atom[bond[i].b].n=n++;
	  }
	if (bond[i].arom)
	  bondid=mol->addBond(atom[bond[i].a].n,atom[bond[i].b].n,Bond::AROMATIC)-1;
	else if (bond[i].type==2)
	  bondid=mol->addBond(atom[bond[i].a].n,atom[bond[i].b].n,Bond::DOUBLE)-1;
	else if (bond[i].type==3)
	  bondid=mol->addBond(atom[bond[i].a].n,atom[bond[i].b].n,Bond::TRIPLE)-1;
	else
	  bondid=mol->addBond(atom[bond[i].a].n,atom[bond[i].b].n,Bond::SINGLE)-1;
	
	if(bond[i].up)
	  mol->getBondWithIdx(bondid)->setBondDir(Bond::ENDUPRIGHT);
	if(bond[i].down)
	  mol->getBondWithIdx(bondid)->setBondDir(Bond::ENDDOWNRIGHT);
        if(bond[i].hash)
          mol->getBondWithIdx(bondid)->setBondDir(Bond::BEGINDASH);
        if(bond[i].wedge)
          mol->getBondWithIdx(bondid)->setBondDir(Bond::BEGINWEDGE);

      }
    mol->addConformer(conf, true);
    for(RWMol::AtomIterator atomIt=mol->beginAtoms();atomIt!=mol->endAtoms();atomIt++) 
     (*atomIt)->calcExplicitValence();
    RDKit::MolOps::cleanUp(*mol);
    const Conformer &conf2 = mol->getConformer();
    DetectAtomStereoChemistry(*mol, &conf2);
  
                              
  RDKit::MolOps::sanitizeMol(*mol);
  RDKit::MolOps::assignBondStereoCodes(*mol);
             
      
  RingInfo *ringInfo = mol->getRingInfo();

  for (unsigned int i=0;i<mol->getNumBonds();i++)
      {
	Bond *b=mol->getBondWithIdx(i);
	if (b!=NULL && ringInfo->numBondRings(i)!=0 &&
	    (b->getBondDir()==Bond::ENDUPRIGHT || b->getBondDir()==Bond::ENDDOWNRIGHT))
	  b->setBondDir(Bond::NONE);
	else if (b!=NULL && ringInfo->numBondRings(i)==0 && b->getIsAromatic())
	  b->setIsAromatic(false);
     }

 int C_Count=0;
 int N_Count=0;
 int O_Count=0;
 int F_Count=0;
 int S_Count=0;
 int Cl_Count=0;
 int R_Count=0;
 for (unsigned int i=0;i<mol->getNumAtoms();i++)
   {
     int anum=mol->getAtomWithIdx(i)->getAtomicNum();
     if (anum==6) C_Count++;
     if (anum==7) N_Count++;
     if (anum==8) O_Count++;
     if (anum==9) F_Count++;
     if (anum==16) S_Count++;
     if (anum==17) Cl_Count++;
     if (anum==0) R_Count++;
   }

 int num_rings=ringInfo->numRings();
 ROMol *pattern_rotors=SmartsToMol("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]");
 
 std::vector<MatchVectType> matches;
 rotors=SubstructMatch(*mol,*pattern_rotors,matches);

  vector<int> Num_Rings(8,0);
  VECT_INT_VECT atomRings; // VECT_INT_VECT is vector< vector<int> >
  atomRings=ringInfo->atomRings();
  for(VECT_INT_VECT_CI ringIt=atomRings.begin();ringIt!=atomRings.end();++ringIt)
    if(ringIt->size()<8) Num_Rings[ringIt->size()]++;
  VECT_INT_VECT bondRings; // VECT_INT_VECT is vector< vector<int> >
  bondRings=ringInfo->bondRings();
  unsigned int num_aromatic=0;
  for(VECT_INT_VECT_CI ringIt=bondRings.begin();ringIt!=bondRings.end();++ringIt)
   {
     bool isAromatic=true;
     for(INT_VECT_CI bondIt=ringIt->begin();bondIt!=ringIt->end();++bondIt)
      if(!mol->getBondWithIdx(*bondIt)->getIsAromatic())
         {
          isAromatic=false;
          break;
         }
     if(isAromatic) num_aromatic++;
   }
  std::string smiles;
  smiles = MolToSmiles(*(static_cast<ROMol *>(mol)),true,false); 
  num_fragments=count_fragments(smiles);

  confidence=confidence_function(C_Count,N_Count,O_Count,F_Count,S_Count,Cl_Count,
				 num_rings,num_aromatic,num_fragments,&Num_Rings);

  r56=Num_Rings[5]+Num_Rings[6];


  for (int i=0;i<n_bond;i++)
    if (bond[i].exists) 
      {
	atom[bond[i].a].n=0;
	atom[bond[i].b].n=0;
      }
  return(smiles);
}

