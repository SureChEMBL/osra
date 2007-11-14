/*********************************************************************
  OSRA: Optical Structure Recognition
  
  Created by Igor Filippov, 2007 (igorf@helix.nih.gov)

  This program is free software; the part of the software that was written 
  at the National Cancer Institute is in the public domain.  This does not
  preclude, however, that components such as specific libraries used in the
  software may be covered by specific licenses, including but not limited
  to the GNU General Public License as published by the Free Software Foundation; 
  either version 2 of the License, or (at your option) any later version; 
  which may impose specific terms for redistribution or modification.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
  USA

*********************************************************************/


#include <iostream>



using namespace std;
extern "C" {
#include "potracelib.h"
}

struct atom_s {
  double x,y;
  string label;
  int n;
  potrace_path_t *curve;
  bool exists,corner;
  int charge;
};
typedef struct atom_s atom_t;

struct bond_s {
  int a,b,type;
  potrace_path_t *curve;
  bool exists;
  bool hash;
  bool wedge;
  bool up;
  bool down;
  bool small;
};
typedef struct bond_s bond_t;

struct letters_s {
  double x,y,r;
  char a;
  bool free;
};
typedef struct letters_s letters_t;

struct label_s {
  double x1,y1,r1,x2,y2,r2;
  string a;
};
typedef struct label_s label_t;


struct box_s {
  int x1,y1,x2,y2;
};
typedef struct box_s box_t;

struct lbond_s{
  int a,b;
  double x;
  bool exists;
};
typedef struct lbond_s lbond_t;

struct dash_s {
  double x,y;
  bool free;
  potrace_path_t *curve;
};
typedef struct dash_s dash_t;

string fix_atom_name(string s,int n);
int getValency(string s);
string get_smiles(atom_t *atom, bond_t *bond, int n_bond, int &rotors);
