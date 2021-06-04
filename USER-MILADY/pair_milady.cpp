/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Thomas Swinburne (CNRS/CINaM)
                         Mihai-Cosmin Marinica (CEA Saclay)
------------------------------------------------------------------------- */


#include "pair_milady.h"

// as for pair_eam
#include <cmath>
#include <cstring>

#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "update.h"

#include "tokenizer.h"
#include "potential_file_reader.h"

// extra - do I need them?
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "neigh_request.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

#define NPERATOMVEC 0

// # entries in peratom_property_vector

/* ---------------------------------------------------------------------- */

PairMILADY::PairMILADY(LAMMPS *lmp) : Pair(lmp)
{

  one_coeff = 1;
  manybody_flag = 1;

  nmax = 0;

  maxneigh = 0;
  allocated = 0;

  MiladyParams = NULL;

  // set comm size needed by this Pair
  comm_forward = 0;
  comm_reverse = 0;

}

/* ----------------------------------------------------------------------
   free all arrays
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairMILADY::~PairMILADY()
{
  int i;
  if(MiladyParams) {
    for (i = 0; i < MiladyParams->nelements; i++) delete [] MiladyParams->elements[i];
    delete [] MiladyParams->elements;
    delete [] MiladyParams->mass;
    delete [] MiladyParams->details;
    memory->destroy(MiladyParams->weights);
    delete MiladyParams;
  }
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] fmap;
    map = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void PairMILADY::compute(int eflag, int vflag)
{

  int i,j,ii,n,inum_half,inum_full,errorflag;
  int *ilist_half,*numneigh_half,**firstneigh_half;
  int *ilist_full,*numneigh_full,**firstneigh_full;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  // grow local arrays if necessary

  if (atom->nmax > nmax) nmax = atom->nmax;

  // neighbor list info

  inum_half = listhalf->inum;
  ilist_half = listhalf->ilist;
  numneigh_half = listhalf->numneigh;
  firstneigh_half = listhalf->firstneigh;

  inum_full = listfull->inum;
  ilist_full = listfull->ilist;
  numneigh_full = listfull->numneigh;
  firstneigh_full = listfull->firstneigh;

  // strip neighbor lists of special bond flags before use with MILADY
  // necessary before doing neigh_f2c and neigh_c2f conversions each step

  if (neighbor->ago == 0) {
    neigh_strip(inum_half,ilist_half,numneigh_half,firstneigh_half);
    neigh_strip(inum_full,ilist_full,numneigh_full,firstneigh_full);
  }

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  // monitor max neighbor - not used yet by this routine...
  n = 0;
  for (ii = 0; ii < inum_half; ii++) n += numneigh_half[ilist_half[ii]];
  if (n > maxneigh) maxneigh = n;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int ntype = atom->ntypes;

  // change neighbor list indices to Fortran indexing

  neigh_c2f(inum_half,ilist_half,numneigh_half,firstneigh_half);
  neigh_c2f(inum_full,ilist_full,numneigh_full,firstneigh_full);


  int ifort;
  int offset = 0;
  errorflag = 0;



  // vptr is first value in vatom if it will be used by milady_force_()
  // else vatom may not exist, so pass dummy ptr
  double *vptr;
  if (vflag_atom) vptr = &vatom[0][0];
  else vptr = &cutmax; // odd choice of dummy....
  for (ii = 0; ii < inum_half; ii++) {
    i = ilist_half[ii];
    ifort = i+1;
    // MILADY computation of descriptors and thus forces
    milady_force_(&ifort,&nmax,&eflag_either,&eflag_global,&eflag_atom,
                &vflag_atom,&eng_vdwl,eatom,&ntype,type,fmap,MiladyParams->mass,&x[0][0],
                &numneigh_half[i],firstneigh_half[i],
                &numneigh_full[i],firstneigh_full[i],
                &f[0][0],vptr,&errorflag,&escale,&fscale);
    if (errorflag) {
      char str[128];
      sprintf(str,"MILADY library error %d",errorflag);
      error->one(FLERR,str);
    }
    offset += numneigh_half[i];
  }

  // change neighbor list indices back to C indexing

  neigh_f2c(inum_half,ilist_half,numneigh_half,firstneigh_half);
  neigh_f2c(inum_full,ilist_full,numneigh_full,firstneigh_full);


  if (vflag_fdotr) virial_fdotr_compute();
}

/* --------------------------------------------------------------------- */

void PairMILADY::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  map = new int[n+1];
  fmap = new int[n];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMILADY::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMILADY::coeff(int narg, char **arg)
{
  int i,j,m,n;

  if (!allocated) allocate();

  /*
  We use same coeff syntax principle as MEAM:
  pair_coeff	* * file.milady ele_1 ele_2 ... map_ele_1 map_ele_2...
  [ele_i] has nelement entries, where nelement = # elements in MILADY file
  [m_ele_i] has atom->ntypes entries, mapping MILADY elements to LAMMPS

  So monatomic example would be
  pair_coeff	* * file.milady ele

  */

  if (narg < 4) error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *
  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (MiladyParams) {
    for (i = 0; i < MiladyParams->nelements; i++)
      delete [] MiladyParams->elements[i];
    delete [] MiladyParams->elements;
    memory->destroy(MiladyParams->mass);
    memory->destroy(MiladyParams->details);
    memory->destroy(MiladyParams->weights);
    memory->destroy(MiladyParams->parameters);
    delete MiladyParams;
  }
  MiladyParams = new MiladyStruct();
  read_file(arg[2]);

  milady_setup_(&MiladyParams->nelements,&MiladyParams->ninteractions,
                &MiladyParams->nparameters,&MiladyParams->parameters[0],
                &MiladyParams->rcutvec[0],&MiladyParams->ndetails,
                &MiladyParams->details[0],&MiladyParams->weights[0],
                &MiladyParams->descriptor_type);
  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < MiladyParams->nelements; j++)
      if (strcmp(arg[i],MiladyParams->elements[j]) == 0) break;
    if (j < MiladyParams->nelements) map[i-2] = j;
    else error->all(FLERR,"No matching element in MILADY potential file");
  }

  // clear setflag since coeff() called once with I,J = * *
  n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j
  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(FLERR,i,MiladyParams->mass[map[i]]);
        count++;
      }
    }
  }
  escale = 1.0;
  fscale = 1.0;
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMILADY::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style MILADY requires newton pair on");

  // need full and half neighbor list

  int irequest_full = neighbor->request(this,instance_me);
  neighbor->requests[irequest_full]->id = 1;
  neighbor->requests[irequest_full]->half = 0;
  neighbor->requests[irequest_full]->full = 1;
  int irequest_half = neighbor->request(this,instance_me);
  neighbor->requests[irequest_half]->id = 2;

  // setup Fortran-style mapping array needed by MILADY package
  // fmap is indexed from 1:ntypes by Fortran and stores a Fortran index
  // if type I is not a MILADY atom, fmap stores a 0
  for (int i = 1; i <= atom->ntypes; i++) fmap[i-1] = map[i] + 1;


}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   half or full
------------------------------------------------------------------------- */

void PairMILADY::init_list(int id, NeighList *ptr)
{
  if (id == 1) listfull = ptr;
  else if (id == 2) listhalf = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMILADY::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  int ii = std::min(map[i]+1,map[j]+1);
  int jj = std::max(map[i]+1,map[j]+1);

  // e.g. nelements=3 : 11 22 33 12 13 23
  int ij=ii-1;
  // ii<jj: ij = sum(k=1,ii)(nelements-(k-1)) + jj-(ii+1)
  if (ii!=jj) ij=ii*MiladyParams->nelements-ii*(ii+1)/2+jj-1;
  //cutghost[i][j] = rcutvec[ij];
  //cutghost[j][i] = cutghost[i][j];
  //cutsq[i][j] = rcutvec[ij]*rcutvec[ij];
  //cutsq[j][i] = cutsq[i][j];
  return MiladyParams->rcutvec[ij];
}

/* ---------------------------------------------------------------------- */

void PairMILADY::read_file(char *filename)
{
  MiladyStruct *file = MiladyParams;

  // read potential file
  if(comm->me == 0) {
    PotentialFileReader reader(lmp, filename, "milady", unit_convert_flag);

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY,
                                                            unit_convert);
    try {
      // line 1-6: comment lines
      for(int i=0; i<6; i++) reader.skip_line();


      // line 7: nelements element_1 weight_1 element_2 weight_2 ...
      ValueTokenizer values = reader.next_values(1);
      file->nelements = values.next_int();
      file->ninteractions = ( file->nelements * (file->nelements+1) ) / 2;


      //values = reader.next_values(2*file->nelements);

      // extract element names from nelements line

      if ((int)values.count() != 2*file->nelements+1)
        error->one(FLERR,"Incorrect number of element names,weights in MILADY potential file");

      // Following convention in pair_eam.cpp - use of new, not memory->create
      file->elements = new char*[file->nelements];
      memory->create(file->weights,file->nelements,"pair:weights");

      for (int i = 0; i < file->nelements; i++) {
        const std::string word = values.next_string();
        const int n = word.length() + 1;
        file->elements[i] = new char[n];
        strcpy(file->elements[i], word.c_str());
        const double w = values.next_double();
        file->weights[i] = w;
      }


      // line 8: descriptor type, nparameters, ndetails
      values = reader.next_values(3);
      file->descriptor_type = values.next_int();
      file->nparameters = values.next_int();
      file->ndetails = values.next_int();

      memory->create(file->mass,file->nelements,"pair:mass");
      memory->create(file->details,file->ndetails,"pair:details");
      memory->create(file->rcutvec,file->ninteractions,"pair:rcutvec");
      memory->create(file->parameters,file->ninteractions*file->nparameters+1,"pair:parameters");

      // line 9: mass_1,...mass_nelements ..... check nelements?
      values = reader.next_values(file->nelements);

      if ((int)values.count() != file->nelements)
        error->one(FLERR,"Incorrect element masses in MILADY potential file");

      for (int i = 0; i < file->nelements; i++)
        file->mass[i] = values.next_double();

      // line 10: detail_1,...detail_ndetails (e.g.  # radials #angulars)
      values = reader.next_values(file->ndetails);
      if ((int)values.count() != file->ndetails)
        error->one(FLERR,"Incorrect details array in MILADY potential file");
      for (int i = 0; i < file->ndetails; i++)
        file->details[i] = values.next_double();

      // line 11: r_cut_ij for i=1,N j=(i,N) where N= number of elements;
      values = reader.next_values(file->ninteractions);
      if ((int)values.count() != file->ninteractions)
        error->one(FLERR,"Incorrect rcut array in MILADY potential file");
      for (int i = 0; i < file->ninteractions; i++)
        file->rcutvec[i] = values.next_double();

      // remaining lines w_ij  for i=1,N j=(i,N) ( i.e. w_11 w_12 w_22 per line for N=2)
      double *dline = new double[file->ninteractions];
      for (int i = 0; i < file->nparameters; i++) {
        reader.next_dvector(dline,file->ninteractions);
        for (int j = 0; j < file->ninteractions; j++)
          file->parameters[j*file->nparameters+i] = dline[j];
      }

      if (unit_convert)
        for (int i = 1; i <= file->ninteractions * file->nparameters; ++i)
          file->parameters[i] *= conversion_factor;

    } catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }

  // broadcast potential information
  MPI_Bcast(&file->nelements, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->ninteractions, 1, MPI_INT, 0, world);
  MPI_Bcast(&file->descriptor_type,1,MPI_INT,0,world);
  MPI_Bcast(&file->nparameters,1,MPI_INT,0,world);
  MPI_Bcast(&file->ndetails,1,MPI_INT,0,world);

  if (comm->me != 0) {
    file->elements = new char*[file->nelements];
    for (int i = 0; i < file->nelements; i++) file->elements[i] = nullptr;
    memory->create(file->mass,file->nelements,"pair:mass");
    memory->create(file->details,file->ndetails,"pair:details");
    memory->create(file->rcutvec,file->ninteractions,"pair:rcutvec");
    memory->create(file->weights,file->nelements,"pair:weights");
    memory->create(file->parameters,file->ninteractions*file->nparameters+1,"pair:parameters");
  }

  // broadcast file->elements string array
  for (int i = 0; i < file->nelements; i++) {
    int n;
    if (comm->me == 0) n = strlen(file->elements[i]) + 1;
    MPI_Bcast(&n, 1, MPI_INT, 0, world);
    if (comm->me != 0) file->elements[i] = new char[n];
    MPI_Bcast(file->elements[i], n, MPI_CHAR, 0, world);
  }
  // broadcast arrays
  MPI_Bcast(&file->mass[0],file->nelements,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->weights[0],file->nelements,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->details[0],file->ndetails,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->rcutvec[0],file->ninteractions,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->parameters[0],file->ninteractions*file->nparameters+1,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

int PairMILADY::pack_forward_comm(int n, int *list, double *buf,
                                int pbc_flag, int *pbc)
{
  int i,j,k,m;
  m = 0;
  return m;
}

/* ---------------------------------------------------------------------- */

void PairMILADY::unpack_forward_comm(int n, int first, double *buf)
{
  int i,k,m,last;
  m = 0;
  last = first + n;
}

/* ---------------------------------------------------------------------- */

int PairMILADY::pack_reverse_comm(int n, int first, double *buf)
{
  int i,k,m,last;

  m = 0;
  last = first + n;
  return m;
}

/* ---------------------------------------------------------------------- */

void PairMILADY::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,k,m;
  m = 0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairMILADY::memory_usage()
{
  double bytes =  ( maxeatom + maxvatom * 6 ) * nmax * sizeof(double);
  //bytes += 3 * maxneigh * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   strip special bond flags from neighbor list entries
   are not used with MILADY
   need to do here so Fortran lib doesn't see them
   done once per reneighbor so that neigh_f2c and neigh_c2f don't see them
------------------------------------------------------------------------- */

void PairMILADY::neigh_strip(int inum, int *ilist,
                           int *numneigh, int **firstneigh)
{
  int i,j,ii,jnum;
  int *jlist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (j = 0; j < jnum; j++) jlist[j] &= NEIGHMASK;
  }
}

/* ----------------------------------------------------------------------
   toggle neighbor list indices between zero- and one-based values
   needed for access by MILADY Fortran library
------------------------------------------------------------------------- */

void PairMILADY::neigh_f2c(int inum, int *ilist, int *numneigh, int **firstneigh)
{
  int i,j,ii,jnum;
  int *jlist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (j = 0; j < jnum; j++) jlist[j]--;
  }
}

void PairMILADY::neigh_c2f(int inum, int *ilist, int *numneigh, int **firstneigh)
{
  int i,j,ii,jnum;
  int *jlist;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (j = 0; j < jnum; j++) jlist[j]++;
  }
}


void *PairMILADY::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"escale") == 0) return (void *) &escale;
  if (strcmp(str,"fscale") == 0) return (void *) &fscale;
  return nullptr;
}
