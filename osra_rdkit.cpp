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
#include <GraphMol/RDKitQueries.h>
#include <vector>
#include <algorithm>
using namespace RDKit;



void addMeX(RWMol *mol,unsigned int aid)
{
  Atom *a=new Atom(6);
  unsigned int aid1=mol->addAtom(a);
  mol->addBond(aid,aid1,Bond::SINGLE);
}

void addOR(RWMol *mol,unsigned int aid)
{
  QueryAtom *query=new QueryAtom(0);
  query->setQuery(makeAtomNullQuery());
  query->setProp("dummyLabel",std::string("?"));
  unsigned int aid1=mol->addAtom(query);
  mol->addBond(aid,aid1,Bond::SINGLE);
}


void addCF(RWMol *mol,unsigned int aid)
{
  Atom *a=new Atom(9);
  unsigned int aid1=mol->addAtom(a);
  mol->addBond(aid,aid1,Bond::SINGLE);
}

void addCF3(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(9);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(9);
  unsigned int aid2=mol->addAtom(a2);
  Atom *a3=new Atom(9);
  unsigned int aid3=mol->addAtom(a3);
  mol->addBond(aid1,aid,Bond::SINGLE);
  mol->addBond(aid2,aid,Bond::SINGLE);
  mol->addBond(aid3,aid,Bond::SINGLE);
}

void addF3CN(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(9);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(9);
  unsigned int aid2=mol->addAtom(a2);
  Atom *a3=new Atom(9);
  unsigned int aid3=mol->addAtom(a3);
  Atom *a4=new Atom(6);
  unsigned int aid4=mol->addAtom(a4);
  mol->addBond(aid1,aid4,Bond::SINGLE);
  mol->addBond(aid2,aid4,Bond::SINGLE);
  mol->addBond(aid3,aid4,Bond::SINGLE);
  mol->addBond(aid4,aid,Bond::SINGLE);
}

void addPh(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(6);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(6);
  unsigned int aid2=mol->addAtom(a2);
  Atom *a3=new Atom(6);
  unsigned int aid3=mol->addAtom(a3);
  Atom *a4=new Atom(6);
  unsigned int aid4=mol->addAtom(a4);
  Atom *a5=new Atom(6);
  unsigned int aid5=mol->addAtom(a5);
  mol->addBond(aid1,aid2,Bond::AROMATIC);
  mol->addBond(aid2,aid3,Bond::AROMATIC);
  mol->addBond(aid3,aid4,Bond::AROMATIC);
  mol->addBond(aid4,aid5,Bond::AROMATIC);
  mol->addBond(aid5,aid,Bond::AROMATIC);
  mol->addBond(aid1,aid,Bond::AROMATIC);
}

void addBzO(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(6);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(6);
  unsigned int aid2=mol->addAtom(a2);
  Atom *a3=new Atom(6);
  unsigned int aid3=mol->addAtom(a3);
  Atom *a4=new Atom(6);
  unsigned int aid4=mol->addAtom(a4);
  Atom *a5=new Atom(6);
  unsigned int aid5=mol->addAtom(a5);
  Atom *a6=new Atom(6);
  unsigned int aid6=mol->addAtom(a6);
  mol->addBond(aid1,aid2,Bond::AROMATIC);
  mol->addBond(aid2,aid3,Bond::AROMATIC);
  mol->addBond(aid3,aid4,Bond::AROMATIC);
  mol->addBond(aid4,aid5,Bond::AROMATIC);
  mol->addBond(aid5,aid6,Bond::AROMATIC);
  mol->addBond(aid1,aid6,Bond::AROMATIC);
  Atom *a7=new Atom(6);
  unsigned int aid7=mol->addAtom(a7);
  mol->addBond(aid1,aid7,Bond::SINGLE);
  mol->addBond(aid,aid7,Bond::SINGLE);
}

void addTHPO(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(6);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(8);
  unsigned int aid2=mol->addAtom(a2);
  Atom *a3=new Atom(6);
  unsigned int aid3=mol->addAtom(a3);
  Atom *a4=new Atom(6);
  unsigned int aid4=mol->addAtom(a4);
  Atom *a5=new Atom(6);
  unsigned int aid5=mol->addAtom(a5);
  Atom *a6=new Atom(6);
  unsigned int aid6=mol->addAtom(a6);
  mol->addBond(aid1,aid2,Bond::SINGLE);
  mol->addBond(aid2,aid3,Bond::SINGLE);
  mol->addBond(aid3,aid4,Bond::SINGLE);
  mol->addBond(aid4,aid5,Bond::SINGLE);
  mol->addBond(aid5,aid6,Bond::SINGLE);
  mol->addBond(aid1,aid6,Bond::SINGLE);
  mol->addBond(aid,aid1,Bond::SINGLE);
}

void addNO2(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(8);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(8);
  unsigned int aid2=mol->addAtom(a2);
  mol->addBond(aid1,aid,Bond::DOUBLE);
  mol->addBond(aid2,aid,Bond::DOUBLE);
}

void addNOHCH3(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(6);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(8);
  unsigned int aid2=mol->addAtom(a2);
  mol->addBond(aid1,aid,Bond::SINGLE);
  mol->addBond(aid2,aid,Bond::SINGLE);
}

void addSO3H(RWMol *mol,unsigned int aid)
{

  Atom *a1=new Atom(8);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(8);
  unsigned int aid2=mol->addAtom(a2);
  Atom *a3=new Atom(8);
  unsigned int aid3=mol->addAtom(a3);
  mol->addBond(aid1,aid,Bond::SINGLE);
  mol->addBond(aid2,aid,Bond::DOUBLE);
  mol->addBond(aid3,aid,Bond::DOUBLE);
}


void addNC(RWMol *mol,unsigned int aid)
{
  Atom *a=new Atom(7);
  unsigned int aid1=mol->addAtom(a);
  mol->addBond(aid1,aid,Bond::TRIPLE);
}

void addnBu(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(6);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(6);
  unsigned int aid2=mol->addAtom(a2);
  Atom *a3=new Atom(6);
  unsigned int aid3=mol->addAtom(a3);
  mol->addBond(aid2,aid3,Bond::SINGLE);
  mol->addBond(aid1,aid2,Bond::SINGLE);
  mol->addBond(aid3,aid,Bond::SINGLE);
}

void addiPr(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(6);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(6);
  unsigned int aid2=mol->addAtom(a2);
  mol->addBond(aid1,aid2,Bond::SINGLE);
  mol->addBond(aid2,aid,Bond::SINGLE);
}

void addEtO(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(6);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(6);
  unsigned int aid2=mol->addAtom(a2);
  mol->addBond(aid1,aid2,Bond::SINGLE);
  mol->addBond(aid2,aid,Bond::SINGLE);
}

void addOiBu(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(6);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(6);
  unsigned int aid2=mol->addAtom(a2);
  Atom *a3=new Atom(6);
  unsigned int aid3=mol->addAtom(a3);
  Atom *a4=new Atom(6);
  unsigned int aid4=mol->addAtom(a4);
  mol->addBond(aid,aid1,Bond::SINGLE);
  mol->addBond(aid1,aid2,Bond::SINGLE);
  mol->addBond(aid2,aid3,Bond::SINGLE);
  mol->addBond(aid2,aid4,Bond::SINGLE);
}

void addtBu(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(6);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(6);
  unsigned int aid2=mol->addAtom(a2);
  Atom *a3=new Atom(6);
  unsigned int aid3=mol->addAtom(a3);
  mol->addBond(aid1,aid,Bond::SINGLE);
  mol->addBond(aid2,aid,Bond::SINGLE);
  mol->addBond(aid3,aid,Bond::SINGLE);
}


void addCOOH(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(8);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(8);
  unsigned int aid2=mol->addAtom(a2);
  mol->addBond(aid1,aid,Bond::DOUBLE);
  mol->addBond(aid2,aid,Bond::SINGLE);
}

void addAc(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(8);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(6);
  unsigned int aid2=mol->addAtom(a2);
  mol->addBond(aid1,aid2,Bond::DOUBLE);
  mol->addBond(aid2,aid,Bond::SINGLE);
}

void addAcO(RWMol *mol,unsigned int aid)
{
  Atom *a1=new Atom(8);
  unsigned int aid1=mol->addAtom(a1);
  Atom *a2=new Atom(8);
  unsigned int aid2=mol->addAtom(a2);
  Atom *a3=new Atom(6);
  unsigned int aid3=mol->addAtom(a3);
  mol->addBond(aid2,aid3,Bond::DOUBLE);
  mol->addBond(aid1,aid3,Bond::SINGLE);
  mol->addBond(aid3,aid,Bond::SINGLE);
}

int getAnum(string s)
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
  if (s=="MeO") return(8);
  if (s=="CF") return(6);
  if (s=="CF3") return(6);
  if (s=="F3CN") return(7);
  if (s=="CN")  return(6);
  if (s=="nBu") return(6);
  if (s=="EtO") return(8);
  if (s=="OiBu") return(8);
  if (s=="iPr")  return(6);
  if (s=="tBu") return(6);
  if (s=="COOH") return(6);
  if (s=="Ac") return(6);
  if (s=="AcO") return(8);
  if (s=="NO2") return(7);
  if (s=="Ph") return(6);
  if (s=="MeS") return(16);
  if (s=="MeN") return(7);
  if (s=="SO3H") return(16);
  if (s=="OR")  return(8);
  if (s=="BzO") return(8);
  if (s=="N(OH)CH3") return(7);
  if (s=="THPO")     return(8);
  return(6);
}

void superatom(string s,RWMol *mol,unsigned int n)
{
  if (s=="MeO") addMeX(mol,n);
  if (s=="CF")  addCF(mol,n);
  if (s=="CF3") addCF3(mol,n);
  if (s=="F3CN") addF3CN(mol,n);
  if (s=="CN")   addNC(mol,n);
  if (s=="nBu")  addnBu(mol,n);
  if (s=="EtO")  addEtO(mol,n);
  if (s=="OiBu") addOiBu(mol,n);
  if (s=="iPr")  addiPr(mol,n);
  if (s=="tBu")  addtBu(mol,n);
  if (s=="COOH") addCOOH(mol,n);
  if (s=="Ac")   addAc(mol,n);
  if (s=="AcO")  addAcO(mol,n);
  if (s=="NO2")  addNO2(mol,n);
  if (s=="Ph")   addPh(mol,n);
  if (s=="MeS")  addMeX(mol,n);
  if (s=="MeN")  addMeX(mol,n);
  if (s=="SO3H") addSO3H(mol,n);
  if (s=="OR")   addOR(mol,n);
  if (s=="BzO")  addBzO(mol,n);
  if (s=="N(OH)CH3")  addNOHCH3(mol,n);
  if (s=="THPO")      addTHPO(mol,n);
}

string get_smiles(atom_t *atom, int real_atoms,bond_t *bond, int n_bond, int &rotors, 
		  double &confidence, int &num_fragments, int &r56, double avg)
{
  RWMol *mol=new RWMol();
  int bondid=0;
  int anum;
  double scale=CC_BOND_LENGTH/avg;

  Conformer *conf = new Conformer(real_atoms);	
  std::string smiles="";
  rotors=0;
  confidence=-1000;
  num_fragments=0;
  r56=0;
  vector<int> bondid_to_i(MAX_ATOMS,-1);

  for (int i=0;i<n_bond;i++)
    if (bond[i].exists) 
      {
	atom[bond[i].a].n=-1;
	atom[bond[i].b].n=-1;
      }
  for (int i=0;i<n_bond;i++)
    if (bond[i].exists) 
      {
	if (atom[bond[i].a].n<0)
	  {
  	    RDGeom::Point3D pos;
  	    pos.x=atom[bond[i].a].x*scale;
  	    pos.y=atom[bond[i].a].y*scale;
  	    pos.z=0;
	    anum=getAnum(atom[bond[i].a].label);
	    Atom *a=new Atom(anum);
	    if (atom[bond[i].a].charge!=0)
	       a->setFormalCharge(atom[bond[i].a].charge);
	    if (anum==0)
	    {
	      QueryAtom *query=new QueryAtom(0);
	      query->setQuery(makeAtomNullQuery());
	      delete a;
	      a=query;
	      a->setProp("dummyLabel",std::string("?"));
             }  
	    unsigned int aid=mol->addAtom(a);                         
	    superatom(atom[bond[i].a].label,mol,aid);
	    conf->setAtomPos(aid, pos);
	    atom[bond[i].a].n=aid;
	  }
	if (atom[bond[i].b].n<0)
	  {
	    RDGeom::Point3D pos;
	    pos.x=atom[bond[i].b].x*scale;
	    pos.y=atom[bond[i].b].y*scale;
	    pos.z=0;
	    anum=getAnum(atom[bond[i].b].label);
	    Atom *b=new Atom(anum);
	    if (atom[bond[i].b].charge!=0)
	      b->setFormalCharge(atom[bond[i].b].charge);
            if (anum==0)
               {
                 QueryAtom *query=new QueryAtom(0);
                 query->setQuery(makeAtomNullQuery());
                 delete b;
                 b=query;
                 b->setProp("dummyLabel",std::string("?"));
	       }
	    unsigned int aid=mol->addAtom(b);
	    superatom(atom[bond[i].b].label,mol,aid);
	    conf->setAtomPos(aid, pos);
	    atom[bond[i].b].n=aid;
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
	bondid_to_i[bondid]=i;
      }
  MolOps::findSSSR(*mol);
  for(ROMol::BondIterator bondIt=mol->beginBonds();bondIt!=mol->endBonds();++bondIt)
    if( ((*bondIt)->getIsAromatic() || (*bondIt)->getBondType()==Bond::AROMATIC)
        && !mol->getRingInfo()->numBondRings((*bondIt)->getIdx()) )
      {
	(*bondIt)->setIsAromatic(false);
	(*bondIt)->setBondType(Bond::SINGLE);
	int i=bondid_to_i[(*bondIt)->getIdx()];
	if (i>=0)
	  if (bond[i].type==2)
	    (*bondIt)->setBondType(Bond::DOUBLE);
	  else if (bond[i].type==3)
	    (*bondIt)->setBondType(Bond::TRIPLE);
      }
   bool doStereo=true;
   try {
    mol->addConformer(conf, true);
   }
   catch (...)
     {
       doStereo=false;
     }
  
    for(RWMol::AtomIterator atomIt=mol->beginAtoms();atomIt!=mol->endAtoms();atomIt++) 
    {
     try {
       (*atomIt)->calcExplicitValence();
     }
     catch (...)
     {
       string symbol=(*atomIt)->getSymbol();
       int id=(*atomIt)->getIdx();
       QueryAtom *query=new QueryAtom(0);
       query->setQuery(makeAtomNullQuery());
       query->setProp("dummyLabel",symbol);
       mol->replaceAtom(id,query);
     }
     try {
       (*atomIt)->calcImplicitValence();
     }
     catch (...)
     {
      (*atomIt)->setNoImplicit(true);
     }
    }
    if (doStereo)
    {
     RDKit::MolOps::cleanUp(*mol);
     const Conformer &conf2 = mol->getConformer();
     DetectAtomStereoChemistry(*mol, &conf2);
    }
    
    try {                            
      RDKit::MolOps::sanitizeMol(*mol);
    }
    catch (...)
      {
	delete mol;
	return(smiles);
      }

    try {
      RDKit::MolOps::assignBondStereoCodes(*mol);
    } catch (...)
      {
	// do nothing?
      }
             
    
  RingInfo *ringInfo = mol->getRingInfo();

    for (unsigned int i=0;i<mol->getNumBonds();i++)
      {
	Bond *b=mol->getBondWithIdx(i);
	if (b!=NULL && ringInfo->numBondRings(i)!=0 &&
	    (b->getBondDir()==Bond::ENDUPRIGHT || b->getBondDir()==Bond::ENDDOWNRIGHT))
	  b->setBondDir(Bond::NONE);
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
 try 
   {
     rotors=SubstructMatch(*mol,*pattern_rotors,matches,
       true,false,false,false);
   }
 catch (...)
   {
     rotors=-1;
   }

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
  
  smiles = MolToSmiles(*(static_cast<ROMol *>(mol)),true,false); 
  num_fragments=count_fragments(smiles);

  confidence=confidence_function(C_Count,N_Count,O_Count,F_Count,S_Count,Cl_Count,
				 num_rings,num_aromatic,num_fragments,&Num_Rings);

  r56=Num_Rings[5]+Num_Rings[6];


  return(smiles);
}

