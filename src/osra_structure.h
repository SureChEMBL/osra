#include <math.h> // fabs(double)
#include "osra_labels.h"
#include "osra_common.h"
#include "osra_ocr.h"
#include "osra_openbabel.h"
#include "osra.h"

#define PI 3.14159265358979323846


//struct: dash_s
// used to identify dashed bonds in <osra.cpp::find_dashed_bonds()> and small bonds in <osra.cpp::find_small_bonds()>
struct dash_s
{
  //double: x,y
  //coordinates
  double x, y;
  //bool: free
  // is this dash available for a perspective dashed bond?
  bool free;
  //pointer: curve
  // pointer to the curve found by Potrace
  const potrace_path_t *curve;
  //int: area
  // area occupied by the dash
  int area;
};
//typedef: dash_t
//defines dash_t type based on dash_s struct
typedef struct dash_s dash_t;


void remove_disconnected_atoms(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond);
void remove_zero_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom);
void collapse_doubleup_bonds(vector<bond_t> &bond, int n_bond);
double skeletize(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, const Image &image, double threshold,const ColorGray &bgColor, double dist, double avg);
double dist_double_bonds(const vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, double avg);
int double_triple_bonds(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, double avg, int &n_atom,double max_dist_double_bond);
void extend_terminal_bond_to_label(vector<atom_t> &atom, const vector<letters_t> &letters, int n_letters, const vector<bond_t> &bond, int n_bond, const vector<label_t> &label, int n_label, double avg, double maxh,
                                   double max_dist_double_bond);
void extend_terminal_bond_to_bonds(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, double avg, double maxh,double max_dist_double_bond);
void assign_charge(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond, const map<string, string> &fix,const map<string, string> &superatom, bool debug);
int find_bonds(vector<atom_t> &atom, vector<bond_t> &bond, int b_atom, int n_atom, int n_bond, const potrace_path_t * const p);
int find_atoms(const potrace_path_t *p, vector<atom_t> &atom, vector<bond_t> &bond, int *n_bond);
int find_dashed_bonds(const potrace_path_t *p, vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int *n_bond,int max, double avg, const Image &img,
                      const ColorGray &bg, double THRESHOLD, bool thick, double dist);
int find_small_bonds(const potrace_path_t *p, vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int *n_bond,
                     double max_area, double Small, double thickness);
int resolve_bridge_bonds(vector<atom_t> &atom, int n_atom, vector<bond_t> &bond, int n_bond, double thickness,
                         double avg_bond_length, const map<string, string> &superatom);
void collapse_atoms(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond, double dist);
void collapse_bonds(vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond, double dist);
int fix_one_sided_bonds(vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom, double thickness, double avg);
double find_wedge_bonds(const Image &image, vector<atom_t> &atom, int n_atom, vector<bond_t> &bond, int n_bond,
                        const ColorGray &bgColor, double THRESHOLD_BOND, double max_dist_double_bond, double avg, int limit, int dist = 0);
void collapse_double_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double dist);
void find_up_down_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double thickness);
void find_old_aromatic_bonds(const potrace_path_t *p, vector<bond_t> &bond, int n_bond, vector<atom_t> &atom,int n_atom, double avg);
void flatten_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double maxh);
void mark_terminal_atoms(const vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, int n_atom);
void find_limits_on_avg_bond(double &min_bond, double &max_bond, const vector<vector<double> > &pages_of_avg_bonds,
                             const vector<vector<double> > &pages_of_ind_conf);


