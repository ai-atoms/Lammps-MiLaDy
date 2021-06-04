/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef PAIR_CLASS

PairStyle(milady,PairMILADY)

#else

#ifndef LMP_PAIR_MILADY_H
#define LMP_PAIR_MILADY_H

extern "C" {
  /*milady_force(i, nmax,   &
            eflag_either, eflag_global, eflag_atom, vflag_atom, &
            eng_vdwl, eatom, ntype, type, fmap, map_mass, x, &
            numneigh, firstneigh, numneigh_full, firstneigh_full, &
            f, vatom, errorflag, escale, fscale)*/
  void milady_force_(int *, int *, int *, int *, int *, int *,
                   double *, double *, int *, int *, int *, double *,
                   double *, int *, int *, int *, int *,
                   double *, double *, int *, double *, double *);
  void milady_setup_(int *, int *, int *, double *, double *, int *, double *, double *, int *);
}


#include "pair.h"
#include <iostream>

namespace LAMMPS_NS {

class PairMILADY : public Pair {
 public:
  PairMILADY(class LAMMPS *);
  ~PairMILADY();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void *extract(const char *, int &);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  double cutmax;      // max cutoff for all elements
  double fscale, escale; // for thermodynamic integration


  struct MiladyStruct {
    int nelements; // number of elements
    int ninteractions; // = nelements * (nelements+1)/2
    double *rcutvec;    // nvector of nelement cutoff radii

    char **elements; // vector of nelement chemical symbols
    double *mass; // vector of nelement masses

    int descriptor_type; // MILADY descriptor type
    double *weights; // vector of nelement weights

    int nparameters; // number of parameters per interaction
    double *parameters; // (flattened) nparameters*ninteractions matrix

    int ndetails; // # of MILADY details (e.g. number of radial channels)
    double *details;// array of MILADY details

  };

  MiladyStruct *MiladyParams;

  int *map;           // mapping from atom types to elements
  int *fmap;          // Fortran version of map array for MILADY lib

  int maxneigh;
  int nmax;


  double **peratom_property_vector;    // per atom property vector

  void allocate();
  void grab(FILE *, int, double *);
  void read_file(char *);
  void neigh_strip(int, int *, int *, int **);
  void neigh_f2c(int, int *, int *, int **);
  void neigh_c2f(int, int *, int *, int **);
  void data_c2f(int, int *, int *, int **);
};

}

#endif
#endif
