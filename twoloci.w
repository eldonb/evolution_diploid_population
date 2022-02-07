\pdfoutput=1
\documentclass[a4paper,12pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[lf]{Baskervaldx}
\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
\usepackage[round,numbers,super]{natbib}
\usepackage{color}
\usepackage{a4wide,fullpage}
\usepackage{setspace}
\usepackage{hyperref}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage[right]{lineno}
\usepackage{verbatim}
\usepackage{tabto}
\usepackage{lipsum}
\usepackage{orcidlink}
\setstretch{1.5}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\be}{\begin{equation}}%
\newcommand{\ee}{\end{equation}}%
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\G}{\ensuremath{\mathbb{G}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushright}
  \textbf{\LARGE \@@title}

  \@@author
\end{flushright}\egroup
}
\makeatother
\title{Fixation at two loci}
\author{Bjarki Eldon\footnote{MfN Berlin, Germany} \footnote{Supported by 
  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
%% „funded by the Deutsche Forschungsgemein-schaft (DFG, German Research Foundation) –Projektnummer(n)“.
%% 
  through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17-2 to Wolfgang Stephan; acknowledge  funding by the Icelandic Centre of Research through an
Icelandic Research Fund Grant of Excellence no.\
185151-051 to  Einar \'Arnason, Katr\'in Halld\'orsd\'ottir, Alison M.\ Etheridge,   WS, and BE. BE also acknowledges Start-up module grants through SPP 1819  with Jere Koskela and Maite Wilke-Berenguer, and  with Iulia Dahmer. \\ \today} \orcidlink{https://orcid.org/0000-0001-9354-2391} }

\begin{document}
\maketitle

\rule{\textwidth}{.8pt}


\begin{abstract}
 This C++ code generates excursions of the evolution of a diploid
      population partitioned into two genetic types at two loci, with viability
      weight determined by $W = e^{-s(z_0 - z)^2}$, where $z$ is the
      trait value of a  given individual, and $z_0$ is the optimal
      trait value, and $s > 0$ is the strength of selection.  The population
      evolves according to a model of random sweepstakes and viability
      selection and randomly occurring bottlenecks.  We estimate the probability of fixation of the type
      conferring advantage, and the expected time to fixation
      conditional on fixation of the type conferring  selective advantage at the two loci.  
\end{abstract}

\tableofcontents


@* {\bf Copyright}. 

Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline

This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.


@* {\bf Compilation,  output and execution}. 
\label{compile}

 This CWEB
      \citep{knuth1994cweb} document (the {\tt .w} file) can be
      compiled with {\tt cweave} to generate a {\tt .tex} file, and
      with {\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c}
      file.

One can use {\tt cweave} to generate a {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file. To compile the C++ code (the {\tt
.c} file), one needs the GNU Scientific Library.   
Using a Makefile can be helpful, calling this file {\tt iguana.w}


 {\tt
iguana.pdf : iguana.tex \\
\tab\quad\quad\quad\quad cweave iguana.w \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        bibtex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        ctangle iguana \\
\tab\quad\quad\quad\quad        c++ -Wall -Wextra -pedantic -O3 -march=native -m64 iguana.c -lm -lgsl -lgslcblas \\
        
clean :  \\
\tab\quad\quad\quad\quad        rm -vf iguana.c iguana.tex \\
}


Use {\tt valgrind} to check for memory leaks:

{\tt valgrind -v --leak-check=full --show-leak-kinds=all <program call>}




@* {\bf introduction}. 
\label{intro}


We consider a diploid population of maximum size $2N$ diploid individuals.  Let $X^N, X_1^N, \ldots, X_N^N$ be i.i.d.\ discrete random variables taking values in $\{2, \ldots, \Psi_N\}$; the $X_1^N, \ldots, X_N^N$ denote the random number of diploid  juveniles independently   produced  in a given generation according to
\be
\label{PXN}
   \prb{X^N = k} = \frac{ (\Psi_N +1)^\alpha }{ (\Psi_N + 1)^\alpha -
 (1/2)^{\alpha} } \left( \frac{1}{k^\alpha} - \frac{1}{(k+1)^{\alpha}} \right),
\quad 1 \leq k \leq \Psi_N.  \ee   The mass in Eq \eqref{PXN} is normalised so that $\prb{ 2 \leq
X^N \leq \Psi_N} =1 $, and $\prb{X^N = k} \ge \prb{X^N = k+1}$. Given
a pool of at least $N$ juveniles, we sample $N$ juveniles for the next
generation.  Leaving out an atom at zero and one gives $X_1^N + \cdots + X_N^N
\ge 2N$ almost surely, guaranteeing that we always have at least $2N$
juveniles to choose from in each generation. 



Write $X_{1} \sim L(\alpha,\Psi_{N})$ if $X_{1}$ is distributed
according to Eq~\eqref{PXN} for given values of $\alpha$ and
$\Psi_{N}$. Let $0 < \alpha_{1} < 2$ and $\alpha_{2} > 2$ be fixed and
consider the mixture distribution\cite{dahmer_coales}
\begin{equation}
\label{eq:3}
X_{1}, \ldots , X_{N} \sim
\begin{cases}
L(\alpha_{1}, \Psi_{N}) & \text{with probability $\varepsilon_{N}$,} \\
L(\alpha_{2}, \Psi_{N}) & \text{with probability $1 - \varepsilon_{N}$.} \\
\end{cases}
\end{equation}
Similarly, by identifying the appropriate scaling of $\varepsilon_{N}$
one can keep
$\alpha$ fixed and varied $\Psi_{N}$\cite{chetwyn-diggle_beta}.  

Each diploid  juvenile inherits two alleles, one from each parent, and is assigned a
viability weight $z$ accoring to the two-locus type; the wild type is assigned
the weight $z=e^{-sf(g)}$ for some fixed $s > 0$ and $f(g)$ is a
function for how the two-locus type affects  the weight.  If the total
number of juveniles at any given time  exceeds  $2N$   we sample an exponential with rate
the given viability weight,  and  $2N$ juveniles with the smallest
exponential replace the parents.  In any given generation a bottleneck
of a fixed size $N_{b}$
occurs with a fixed probability. If a bottleneck occurs we sample
$N_{b}$ individuals independently and uniformly at random without
replacement. The surviving individuals then produce juveniles, and if
the total number of juveniles is less than  the capacity $2N$ all the
juveniles survive, otherwise we assign weights and sample $2N$
juveniles according to the weights as just described.



At both loci there are two types ($0$, $1$),  so there are three
genotypes at each locus, and nine two-locus genotypes.  
Let $Y_{t} \equiv \{Y_{t} : t \ge 0 \}$ denote the frequency of the
two-locus  type configuration in the population, i.e.\  $Y_{t}$ takes
values in $[0,2N]^{9}$,     and write  $T_{k}(y) :=  \min\{ t
\ge 0 :  Y_{t} = k, Y_{0} = y \}$.    We are interested in the
quantitites
\begin{equation}
\label{eq:4}
\begin{split}
p_{N}(y_{0}) & :=   \prb{ T_{N}(y_0) < T_{0}(y_0) } \\
\tau_{N}(y_{0}) & := \EE{ T_{N}(y_0) :  T_{N}(y_0) < T_{0}(y_0) } \\
\end{split}
\end{equation}
where $y_{0}$ is the starting configuration, i.e.\ the  number of
diploid individuals of each two-locus genotype at time 0. The quantity
$p_{N}(y_{0})$ is the  probability of all $2N$ diploid  individuals  reaching the two-locus type
conferring maximum  advantage when starting from configuration
$y_{0}$, and  $\tau_{N}(y_{0})$ is the expected time to do so
conditional on the population  reaching the configuration when
starting from $y_{0}$.  For example, if type $1$ confers advantage at
both loci, then  the two-locus type in question  would  be  the 
type $1/1-1/1$  where individuals are homozygous for type 1 at both
loci. 



@*{\bf Code}. 
\label{SEC:code}

We collect the key containers and constants  into a struct
\S~\ref{SEC:structM},  we use the GSL random  number generator
\S~\ref{SEC:rng}, in \S~\ref{SEC:cdf} we compute the cumulative
density function for sampling random numbers of juveniles according to
the inverse CDF method,  in \S~\ref{SEC:random_number_juveniles} we
sample a random number of juveniles, in \S~\ref{SEC:comp} we define a
comparison function for sorting the exponentials in \S~\ref{SEC:nth},
in \S~\ref{SEC:samplepool} we sample a pool of juveniles and assign
weight to them in \S~\ref{SEC:assignweight},  in
\S~\ref{sec:survivehypergeom} we sample the number of individuals  of the advantageous
type surviving  a bottleneck,  in \S~\ref{sec:surviveweight} we count
the number of advantageous type surviving selection according to their
weight,  in \S~\ref{sec:onestep} we step through one generation by
checking if a bottleneck occurs and then produce juveniles if neither
fixation nor loss of the advantageous type occurs,  in
\S~\ref{sec:trajectory} we generate one excursion until fixation or
loss of the advantageous type starting with one copy of the
advantageous type, the main module \S~\ref{SEC:main} generates a given
number of trajectories, \S~\ref{sec:examples} holds examples of
trajectories to fixation of the advantageous type. 

@*1 {Includes}. 
\label{SEC:includes}

The included libraries.

@<includes@>=@#
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
#include <cstddef>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <list>
#include <string>
#include <fstream>
#include <chrono>
#include <forward_list>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "diploid_excursions_random_bottlenecks.hpp"




@*1 {\bf the random number generator}. 
\label{SEC:rng}


@<gslrng@>=@#
gsl_rng * rngtype ;
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}


@*1 {\bf the number of diploid individuals by index}.
\label{sec:numberbyindex}

return the number of diploid individuals by index. Let $0$ denote the
homozygous type $0/0$, $1$ the heterozygous type $0/1$, and 2 the
homozygous type $1/1$ at each locus.    The population is
an array with indexes zero to nine  with the configuration
\begin{table}[htp]
\centering
\caption{index}
\label{tab:popconfiguration}
\begin{tabular}[htp]{lllllllllll}
\hline
index &  0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8  \\
two-locus type & $(0,0)$ & $(0,1)$ & $(0,2)$ & $(1,0)$ & $(1,1)$ & $(1,2)$ & $(2,0)$ & $(2,1)$ & $(2,2)$ \\  
\hline
\end{tabular}
\end{table}


@<numberbyindex@>=@#
static  unsigned int number_diploid_individuals_by_index(const std::vector<unsigned>& population,    const size_t c_index)
  {
    return population[ c_index ] ;
  }


@*1 {\bf check if lost type}.
\label{sec:checklosttype}

Check if lost  type 1  conferring  advantage  at either locus; see
Table~\ref{tab:popconfiguration}.

@<checklosttype@>=@#
static  size_t check_if_lost_type( const std::vector<unsigned>& population )
  {
    /* return 1 if lost type at either locus, otherwise 0 \newline  */
    return( (population[1] + population[2] + population[4] + population[5] + population[7] + population[8] < 1) || (population[3] + population[4] + population[5] + population[6] + population[7] + population[8] < 1) ? 1 : 0 ); 
  }


@*1 {\bf clear the container for juveniles }.
\label{sec:clearjuveniles}

clear the container containing the juveniles; each juvenile is stored
as a pair of genotype index (see Table~\ref{tab:popconfiguration}),
and viability weight. 

@<clearjuveniles@>=@#
static  void removealljuveniles( std::vector< std::pair< size_t, double>>& juveniles )
  { juveniles.clear() ;
    juveniles.shrink_to_fit();
    assert( juveniles.size() < 1 );
  }


@*1 {\bf total number of juveniles}.
\label{sec:totalnumberjuvs}

return the total number of juveniles

@<totalnumberjuvs@>=@#
static  size_t totalnumberjuveniles(  const std::vector< std::pair< size_t, double>>& juveniles )
  { return juveniles.size() ; }
 

@*1 {\bf add a juvenile with a given type index and weight}. 
\label{sec:addjuv}

add a juvenile with a given viability weight and two-locus genotype index
@<addjuv@>=@#
static  void add_juvenile(  std::vector< std::pair< size_t, double>>& juveniles,  const size_t g, const double weight)
  { juveniles.push_back( std::make_pair(g, weight) ) ;
  }


@*1 {\bf number of individuals with given types}.
\label{sec:numberwithgiventype}

see Table~\ref{tab:popconfiguration}. 

@<numberwithgiventype@>=@#
unsigned int number_type( const std::vector<unsigned>& population,    const size_t c_lone, const size_t c_ltwo)
  { return population[ (3*c_lone) + c_ltwo ]; }


@*1 {\bf set population to zero}.
\label{sec:setpopzero}

set number of individuals in the  population to zero for all types

@<setpopzero@>=@#
static void set_population_zero( std::vector<unsigned>& population)
  {
    std::fill( std::begin(population), std::end(population), 0);
    assert( std::accumulate( std::begin(population), std::end(population), 0) == 0); 
    /*
    for( |size_t| i = 0; i < 9 ; ++i){
      population[i] = 0; }
    */
  }


@*1 {\bf current number of diploid individuals}.
\label{sec:currentnumberind}

return the current number of individuals in the population

@<currentind@>=@#
static  unsigned int current_number_individuals(const std::vector<unsigned>& population)
  {
    return std::accumulate( std::begin(population), std::end(population), 0); 
    /*
    |size_t| s = 0 ;
    for( |size_t| i = 0 ; i < 9 ; ++i){
      s += population[i] ; }
    return (s); 
    */
  }


@*1 {\bf update the number of individuals of a given type }. 
\label{sec:updatecount}

add or subtract by one  the number of diploid individuals with a given
type Table~\ref{tab:popconfiguration}.

@<updatecount@>=@#
static void update_count_type( std::vector<unsigned>& population,  const size_t c_type, const size_t add_subtract )
  {
  assert( population[c_type] > 0 );
  population[c_type ] +=  (add_subtract < 1 ? 1 : -1) ; }


@*1 {\bf return the most numerous type  }. 
\label{sec:mosttype}

@<mosttype@>=@#
static   size_t type_most_copies( const std::vector<unsigned>& population)
  { return std::distance( population.begin(), std::max_element(population.begin(), population.end()) );  }



@*1 {\bf all juveniles survive }. 
\label{sec:allsurvive}

all juveniles survive if the total number of juveniles produced at any
given time does not exceed the carrying capacity.

@<allsurvive@>=@#
static  void update_population_all_juveniles(const std::vector< std::pair< size_t, double>>& juveniles, std::vector<unsigned>& population)
  {
   /* set the population to zero \S~\ref{sec:setpopzero} \newline */
    set_population_zero(population) ;
    assert( current_number_individuals(population) < 1 );
    /* add a juvenile to the population by updating the corresponding
    number of individuals with the type of the juvenile \newline */
    for( const auto &j : juveniles){
      population[ std::get<0>( j ) ] += 1;}
  }


@*1 {\bf initialise the containers}.
\label{sec:initconts}

initialise the containers for the population and  the juveniles and
the cumulative density functions; the initial population configuration
$y_{0}$ is with  $2N-2$ diploid individuals as double homozygous for
the wild type allele, one individual  $(0,1)$ and the other as
$(1,0)$. 

@<initconts@>=@#
static void init_containers( std::vector< unsigned>& population, std::vector<double>& cdf_one, std::vector<double>& cdf_two   )
  {

    population.clear() ;
    population.assign( 9, 0); 
    assert( current_number_individuals(population) < 1);
    assert( population.size() == 9 ); 
    population[0] = GLOBAL_CONST_II - 2 ;
    population[1] = 1;
    population[3] = 1;
    assert( current_number_individuals( population) == GLOBAL_CONST_II );

    cdf_one.clear() ;
    cdf_two.clear() ;
    cdf_one.reserve( GLOBAL_CONST_CUTOFF_ONE + 2) ;
    cdf_two.reserve( GLOBAL_CONST_CUTOFF_TWO + 2) ;
    /* set $\prb{X^{N} < 2} = 0$ Eq~\eqref{PXN} \newline */
    cdf_one.push_back( 0.);
    cdf_one.push_back( 0.);
    cdf_two.push_back(0.);
    cdf_two.push_back(0.);
    assert( cdf_one.size() == 2 );
    assert( cdf_two.size() == 2 );
  }



@*1 {\bf initialise for a trajectory}.
\label{sec:inittraj}

initialise the population  for a new trajectory 

@<inittraj@>=@#
static  void init_for_trajectory( std::vector< unsigned>& population)
  {
    set_population_zero(population) ;
    assert( current_number_individuals(population) < 1);
    population[0] = GLOBAL_CONST_II -2 ;
    population[1] = 1 ;
    population[3] = 1 ;
    
    assert( current_number_individuals(population) == GLOBAL_CONST_II);
  }


@*1 {\bf comparison function for sorting juveniles }. 
\label{sec:comp}

comparison function for sorting juveniles according to viability weight

@<fcom@>=@#
static bool comp( const std::pair<size_t, double> a, const std::pair<size_t, double> b) 
  { 
    return ( std::get<1>(a) < std::get<1>(b) );
  }


@*1 {\bf sort the juveniles }. 
\label{sec:nth}

partially sort the juveniles and return the $2N$th  sorted viability  weight
sorted in ascending order


@<nth@>=@#
static  double nthelm( std::vector< std::pair< size_t, double>>& juveniles )
  {
  /* partially sort the weights using \S~\ref{sec:comp} \newline */
 std::nth_element( juveniles.begin(), juveniles.begin() + (GLOBAL_CONST_II - 1), juveniles.end(), comp);
    
    return( std::get<1>(juveniles[ GLOBAL_CONST_II - 1]) );
  }



@*1 {\bf sample juveniles according to weight }. 
\label{sec:samplejuvsweight}

sample juveniles surviving selection by sampling according to
viability weight

@<samplejuvsweight@>=@#
static  void sample_juveniles_according_to_weight(std::vector<unsigned>& population, const std::vector< std::pair< size_t, double>>& juveniles,   const double c_nth)
  {
    assert( c_nth > 0. );
    set_population_zero(population);
    /* check that the population is correctly initialised
    \S~\ref{sec:setpopzero} \newline */
    assert( current_number_individuals(population) < 1); 
    
    size_t j = 0 ;
    while( j < GLOBAL_CONST_II ){
      assert( j < GLOBAL_CONST_II) ;
      population[ std::get<0>(juveniles[j]) ] += std::get<1>(juveniles[j]) <= c_nth ? 1 : 0 ;
      ++j ;
    }
    /* check that we have sampled correct number of juveniles
    \S~\ref{sec:currentnumberind} \newline */
    assert( current_number_individuals(population) == GLOBAL_CONST_II );
  }


@*1 {\bf genotype from an index }.
\label{sec:secondg}

return the second genotype 0,1, or 2 from a given genotype index as in Table~\ref{tab:popconfiguration}

@<secondg@>=@#
static unsigned int second_locus_genotype_from_index(const int c_i)
{
  assert( c_i > -1); 
  unsigned int x {};
  switch( c_i){
  case 0 : {
    x = 0 ;
    break ;}
  case 1 : {
    x = 1 ;
    break ; }
  case 2 : {
    x = 2 ;
    break ; }
  case 3 : {
    x = 0 ;
    break ; }
  case 4 : {
    x = 1 ;
    break ; }
  case 5 : {
    x = 2;
    break ; }
  case 6 : {
    x = 0 ;
    break ; }
  case 7 : {
    x = 1 ;
    break ; }
  case 8 : {
    x = 2 ;
    break ; }
  default : break ; 
  }
  return x ;
}




/* sample genotype of one parent */
/* return the index of the genotype sampled */

@*1 {\bf sample a parent}.
\label{sec:sampleparent}

return the genotype index (Table~\ref{tab:popconfiguration})  of a
sampled parent

@<sampleparent@>=@#
static  int sample_genotype_parent( std::vector<unsigned>& p,   gsl_rng *r )
{
  int i  = 0 ;
  unsigned int nothers = current_number_individuals(p) - number_type(p,  0,0) ;
  unsigned int x =  gsl_ran_hypergeometric( r, number_diploid_individuals_by_index(p, 0), nothers, 1);

  while( (x < 1) && (i < 7) ){
    ++i ;
    /* update the number of remaining individuals
    \S~\ref{sec:numberbyindex} \newline */
    nothers -= number_diploid_individuals_by_index(p, i);
    x =  gsl_ran_hypergeometric( r, number_diploid_individuals_by_index(p, i), nothers,  1);
  }
  i += (x < 1 ? 1 : 0); 
  /* update  the number of remaining parents \S~\ref{sec:updatecount}  \newline */
  update_count_type( p, i, 1);
  /* return the index  of the genotype of the parent \newline  */
  /* index is between 0 and 8 \newline  */
  return i ;
}


@*1 {\bf assign a genotype to juvenile }.
\label{sec:genotypej}


 assign a single locus  genotype to juvenile given single locus
 genotypes in parents by sampling one allele from each parent
 independently and uniformly at random 
 @<genotypej@>=@#
static int assign_type_juvenile( const int gone, const int gtwo, gsl_rng *r)
{
  int g {} ;
  const double u = gsl_rng_uniform(r) ;
  switch(gone){
  case 0 : {
    g = (gtwo < 1 ? 0 : (gtwo < 2 ? (u < 0.5 ? 0 : 1) : 1) );
    break ;}
  case 1 : {
    g = (gtwo < 1 ? (u < .5 ? 0 : 1) : (gtwo < 2 ? (u < 0.25 ? 0 : ( u < 0.75 ? 1 : 2)) : (u < 0.5 ? 1 : 2) ) ) ;
    break ; }
  case 2 : {
    g = (gtwo < 1 ? 1 : (gtwo < 2 ? (u < .5 ? 1 : 2) : 2) ) ;
    break ; }
  default : break ; } 

  return g ;
}


@*1 {\bf sample random number juveniles}. 
\label{sec:samplerandomjuvs}

sample a random number of juveniles using the inverse CDF method,
i.e.\ return
\begin{equation}
\label{eq:1}
\min \{j \ge 2 :  F(j) \ge u \}
\end{equation}
where $F$ is the CDF
and $u$ is a random uniform on the unit interval

@<samplerandomjuvs@>=@#
static size_t sample_random_number_juveniles( const size_t c_twoone, const std::vector<double>& cdfone, const std::vector<double>& cdftwo,   gsl_rng *r)
{
  const double u = gsl_rng_uniform(r);
  size_t j = 2 ;
  if( c_twoone < 2 ){
    while( u > cdfone[j] ){ ++j ;}}
  else{
    while( u > cdftwo[j]){ ++j ;}}

  assert( j > 1) ;
  return j ;
}


@*1 {\bf compute the viability weight }. 
\label{sec:computeweight}

compute  the viability weight  as  $\exp(-sf(g))$ where $s \ge 0$ is
the strength of selection and  $f(g)$ is the genotype to phenotype
map, i.e.\ we interpret $f(g)$ as a trait value, a phenotype, for the
given two-locus genotype  $g$.  We take $f(g) = \left( h_{1}(g_{1}) +
h_{2}(g_{2}) \right)/2 $ where  $h_{1}$ and $h_{2}$ determine the
contribution of the genotypes at the two loci resp.   
@<computew@>=@#
static double computeweight(const int g_one, const int g_two,   gsl_rng *r)
{
  /* | (g_one < 1 ? 2 : (g_one < 2 ? 1 : 0) ); | \newline  */
  /* complete dominance, ie the wild type is recessive \newline  */
  /*   $h_{1}(g_{1}) =  2\one{g_{1} < 1} $ \newline   */
  const double ggone =  (g_one < 1 ? 2 : 0 ) ;
   /*   $h_{2}(g_{2}) =  2\one{g_{2} = 0} +  \one{g_{2} = 1} $ \newline   */
  const double ggtwo =  (g_two < 1 ? 2 : (g_two < 2 ? 1 : 0 ) ) ;


 /* return a random exponential with rate $\exp(-s f(g))$ used for
 sorting the juveniles according to viability weight  \newline */
 
  return ( gsl_ran_exponential( r,  1./exp( (-GLOBAL_CONST_SELECTION)*pow( (ggone + ggtwo)/2. , 2.) )) ) ;
}


@*1 {\bf sample a litter of juveniles }.
\label{sec:litter}

sample a litter of juveniles, i.e.\ a random number of juveniles with
alleles, from a pair of parents with given genotypes


@<litter@>=@#
static void add_juveniles_for_given_parent_pair( const std::vector<double>& cdfone, const std::vector<double>& cdftwo, std::vector< std::pair<size_t, double>>& jvs,   const int gone, const int gtwo, const size_t conetwo,  gsl_rng *r)
{
/* |gone| and |gtwo|  are the two-locus genotype indexes for the two
parents Table~\ref{tab:popconfiguration}  \newline  */

 /* first sample the number of juveniles produced by the parent pair
 \S~\ref{sec:samplerandomjuvs} \newline */
  const size_t numberj = sample_random_number_juveniles(  conetwo, cdfone, cdftwo, r);
  assert( numberj > 1 ) ;
  int g_locus_one {} ;
  int g_locus_two {} ;
  for( size_t j = 0; j < numberj ; ++j){
/* for each juvenile in the litter sample the two alleles \S~\ref{sec:genotypej} \newline */
    g_locus_one = assign_type_juvenile( (gone < 3 ? 0 : (gone < 6 ? 1 : 2)), (gtwo < 3 ? 0 : (gtwo < 6 ? 1 : 2)),  r) ;
    g_locus_two = assign_type_juvenile(  second_locus_genotype_from_index( gone), second_locus_genotype_from_index( gtwo), r);
    assert( g_locus_one == 0 || g_locus_one == 1 || g_locus_one == 2) ;
    assert( g_locus_two == 0 || g_locus_two == 1 || g_locus_two == 2)
    ;
    /* given alleles add juvenile with viability weight
    \S~\ref{sec:computeweight} \newline */
    add_juvenile( jvs, (3*g_locus_one) + g_locus_two, computeweight( g_locus_one, g_locus_two, r) ) ; }
}


@*1 {\bf pool of juveniles }.
\label{sec:pool}

generate a pool of juveniles for all parent pairs

@<pool@>=@#
static void generate_pool_juveniles( std::vector< std::pair<size_t, double>>& jvs, std::vector<unsigned>& p,  const std::vector<double>& cdfone, const std::vector<double>& cdftwo,   gsl_rng *r)
{
/* clear the container of juveniles \S~\ref{sec:clearjuveniles}
\newline */
  removealljuveniles(jvs) ;
  int gone {} ;
  int gtwo {} ;
  /* sample distribution of number of juveniles Eq~\eqref{eq:3} \newline  */
  const size_t conetwo = (gsl_rng_uniform(r) < GLOBAL_CONST_EPSILON ? 1 : 2) ; 
  /* |i| runs over number of pairs that can be formed from the current
  number of individuals; if $n$ the current number of individuals \S~\ref{sec:currentnumberind}  can
 produce   $\lfloor n/2 \rfloor$ pairs of two-locus genotypes  \newline */
  assert( current_number_individuals(p) < GLOBAL_CONST_I + 1) ;
  const double currenti = current_number_individuals(p) ;
  for ( double i = 0 ; i < floor( currenti / 2. ) ; ++i){
  /* sample a parent genotype \S~\ref{sec:sampleparent} \newline */
    gone = sample_genotype_parent(p, r);
     /* sample another  parent genotype \S~\ref{sec:sampleparent} \newline */
    gtwo = sample_genotype_parent(p, r);
    assert(gone > -1);
    assert( gtwo > -1);
    /* given the parent genotypes add a litter \S~\ref{sec:litter}
    \newline */
    add_juveniles_for_given_parent_pair(cdfone, cdftwo, jvs,   gone, gtwo, conetwo,  r) ;}

  assert( totalnumberjuveniles(jvs) >= static_cast<size_t>( currenti ) ) ;
}


@*1 {\bf bottleneck}. 
\label{sec:bottleneck}


 sample diploid individuals surviving a bottleneck; we sample
 uniformly at random without replacement $N_{b}$ diploid individuals
 by sampling the number of each two-locus genotype surviving a
 bottleneck; we therefore sample a hypergeometric by  updating the
 relevent numbers each time 

@<bottle@>=@#
static void sample_surviving_bottleneck( std::vector<unsigned>& p,   gsl_rng * r)
{
  size_t i = 0 ;

 /* |p| is the population indexed as in
 Table~\ref{tab:popconfiguration} \newline */
 /* |nothers| is the number of individuals in the pot  of the colour
 not being sampled; |p[i]| is the number of the colour being sampled  \newline */
  unsigned int nothers = current_number_individuals(p) - p[i] ;
  unsigned newn = gsl_ran_hypergeometric( r,  p[i],  nothers,  GLOBAL_CONST_BOTTLENECK);
  unsigned int remaining = GLOBAL_CONST_BOTTLENECK - newn ;
  
  /* update count of individuals of type index |i| surviving
  bottleneck \newline  */
  p[i] = newn ;
  while( (i < 7) ){
    ++i ;
    nothers -= p[i] ;
    newn = (remaining > 0 ? gsl_ran_hypergeometric( r, p[i], nothers, remaining) : 0) ;
    p[i] = newn ;
    remaining -= newn ;
  }
  /* update for index 8 with the remaining to sample \newline  */
  p[8] =  (remaining < GLOBAL_CONST_BOTTLENECK ? remaining : GLOBAL_CONST_II );

  assert( current_number_individuals(p) >= GLOBAL_CONST_BOTTLENECK ); 
}



@*1 {\bf take one step}. 
\label{sec:onestep}

step through one generation by first checking if a bottleneck occurs,
and then sample juveniles if the type neither lost nor fixed 

@<onestep@>=@#
static void onestep( std::vector<unsigned>& p,  const std::vector<double>& cdfone, const std::vector<double>& cdftwo, std::vector< std::pair<size_t, double>>& jvs,  gsl_rng *r )
{
  double nth {} ;
  /* check if bottleneck */
  if( gsl_rng_uniform( r) < GLOBAL_CONST_PROBABILITY_BOTTLENECK ){
    /* bottleneck occurs ; sample surviving types \S~\ref{sec:bottleneck}  and update
    population |p| */
    
   sample_surviving_bottleneck(p, r) ; }
  /* first check if lost type at either loci
  \S~\ref{sec:checklosttype} \newline */
  if ( check_if_lost_type(p) < 1 ){
    /* not lost type; check if fixed at both
    \S~\ref{sec:numberbyindex} \newline   */
    if( number_diploid_individuals_by_index(p, 8) < GLOBAL_CONST_II ){
	/* not all  individuals of type 2, so sample juveniles \S~\ref{sec:pool}
        \newline  */
      generate_pool_juveniles(jvs, p, cdfone, cdftwo, r) ;
	if( totalnumberjuveniles(jvs) <= GLOBAL_CONST_II )
	  {
	    /* total number of juveniles not over capacity so all
            survive \S~\ref{sec:allsurvive} \newline */
	    update_population_all_juveniles(jvs, p) ;
	    assert( current_number_individuals(p) >= GLOBAL_CONST_BOTTLENECK ); 
	  }
	else{
	  /* need to sort juveniles \S~\ref{sec:nth} and sample
          according to weight \S~\ref{sec:samplejuvsweight} \newline */
	  nth = nthelm(jvs) ;
	  sample_juveniles_according_to_weight(p, jvs, nth) ;
	  assert( current_number_individuals(p) >= GLOBAL_CONST_II ); 
	}
      }
      /* mutation has fixed at both loci */
    }
    /* mutation has been lost */ 
}


@*1 {\bf trajectory }.
\label{sec:trajectory}

generate one trajectory by stepping through the generations one step 
at a time \S~\ref{sec:onestep} until either lost a type or fixed at
both loci, the current optimal configuration 


@<trajectory@>=@#
static void trajectory( std::vector<unsigned>& p,  const std::vector<double>& cdfone, const std::vector<double>& cdftwo, std::vector< std::pair<size_t, double>>& jvs,   const int numer,   gsl_rng *r)
{
/* initialise for a trajectory \S~\ref{sec:inittraj} \newline */
  init_for_trajectory(p) ;

  std::vector< double > excursion_to_fixation {} ;
  int timi = 0 ;
  while( ( check_if_lost_type(p) < 1 ) && ( current_number_individuals(p) - p[8] > 0) ){
    /* record the number of diploid individuals  homozygous  |1/1|  at both loci (sites) over current number of diploid individuals  */
    
    assert( current_number_individuals(p) >= GLOBAL_CONST_BOTTLENECK) ;
    excursion_to_fixation.push_back( static_cast< double>( number_diploid_individuals_by_index(p, 8) ) / static_cast<double>( current_number_individuals(p) ) ) ;
    ++ timi ;
    onestep(p, cdfone, cdftwo, jvs, r) ; }

  const std::string eskra = "twolocitrajectory" +  std::to_string( numer) + ".txt" ;

/* check if lost the type after an excursion
\S~\ref{sec:checklosttype} \newline */
  if( check_if_lost_type(p) < 1 ){
    /* did not lose the type; check that fixed at both loci for the
    type conferring selective advantage \S~\ref{sec:mosttype} \newline */
    assert( type_most_copies(p) == 8 );
    /* fixation occurs so  print excursion to file \newline  */
    std::cout << eskra << '\n';
    std::ofstream f(eskra, std::ofstream::app);
    assert( f.is_open() ); 
    for( const auto& y: excursion_to_fixation){ f << y << ' ' ;  }
    f << '\n' ;
    f.close() ;
  }
  std::cout << ( check_if_lost_type(p) > 0 ? 0 : 1) << ' ' << timi << '\n' ;
}


@*1 {\bf the mass function Eq~\eqref{PXN} }. 
\label{sec:pxn}

The mass function in Eq~\eqref{PXN} used for sampling random number of
juveniles; see \S~\ref{sec:initcdf}

@<pxn@>=@#
static double px(const double k, const double calpha, const double ccutoff)
{
  return ( (pow( 1./k, calpha) - pow( 1./(k + 1.), calpha) )/( pow( .5, calpha) - pow( 1./(ccutoff + 1.), calpha) ) ) ;
}



@*1 {\bf initialise the CDF from Eq~\eqref{PXN}}. 
\label{sec:initcdf}

Initialise the CDF for the distribution of random number of juveniles
Eq~\eqref{PXN} \S~\ref{sec:pxn}; the population evolves according to
Eq~\eqref{eq:3} so need two versions of the CDF for different values
of $\alpha$  or the cutoff $\Psi_{N}$ 


@<initcdf@>=@#
static void initialise_cdf( std::vector<double>& cdfo, std::vector<double>& cdft )
{
  
  for( double i = 2; i <=  GLOBAL_CONST_PSI_ONE ; ++i){
    cdfo.push_back( cdfo.back() + px( i, GLOBAL_CONST_ALPHA_ONE, GLOBAL_CONST_PSI_ONE) ) ;}

  for( double j = 2; j <= GLOBAL_CONST_PSI_TWO; ++j){
    cdft.push_back( cdft.back() + px( j, GLOBAL_CONST_ALPHA_TWO, GLOBAL_CONST_PSI_TWO)) ; }
}


@*1 {\bf generate a given number of trajectories}.
\label{sec:runsims}


generate a given number of trajectories 

@<runsims@>=@#
static void runsims( const int cnumer,   gsl_rng *r)
{
  
  std::vector< unsigned int> population (9, 0) ;
  std::vector< std::pair< size_t, double>> juveniles {} ;

  std::vector<double> cdf_one {} ;
  std::vector<double> cdf_two {} ;

  /* initialise the main  containers for  the objects we  need
  \S~\ref{sec:initconts} \newline */
  init_containers(population, cdf_one, cdf_two) ;
  /* initialise both CDFs \S~\ref{sec:initcdf} \newline */ 
  initialise_cdf(cdf_one, cdf_two) ;
    
  int z = GLOBAL_CONST_NUMBER_EXPERIMENTS + 1;
  while( --z > 0){
  /* sample a trajectory \S~\ref{sec:trajectory} \newline */
    trajectory(population, cdf_one, cdf_two, juveniles, cnumer, r);
  }
}


@*1 {\bf the main module}. 
\label{SEC:main}


The |main| function

@C

@<includes@>@#
@<gslrng@>@#
@<numberbyindex@>@#
@<checklosttype@>@#
@<clearjuveniles@>@#
@<totalnumberjuvs@>@#
@<addjuv@>@#
@<numberwithgiventype@>@#
@<setpopzero@>@#
@<currentind@>@#
@<updatecount@>@#
@<mosttype@>@#
@<allsurvive@>@#
@<initconts@>@#
@<inittraj@>@#
@<fcom@>@#
@<nth@>@#
@<samplejuvsweight@>@#
@<secondg@>@#
@<sampleparent@>@#
@<genotypej@>@#
@<samplerandomjuvs@>@#
@<computew@>@#
@<litter@>@#
@<pool@>@#
@<bottle@>@#
@<onestep@>@#
@<trajectory@>@#
@<pxn@>@#
@<initcdf@>@#
@<runsims@>@#

int main(int argc, char *argv[])
{
  
  setup_rng(  static_cast<unsigned long int>(atoi(argv[1])) );

  /* run a given number of trajectories \S~\ref{sec:runsims} \newline
  */
  runsims( atoi(argv[1]), rngtype ); 

/* free the random number generator in \S~\ref{SEC:rng}  \newline */
  gsl_rng_free( rngtype ); 
  return 0 ;
}


@* {\bf examples}. 
\label{sec:examples}


\begin{figure}[htp]
\includegraphics[scale=1]{twolocigraphs-crop}
\caption{$2N = 10^{6}$, $\Psi_{N} = 2N$,
$\alpha_{2} = 3$, $\varepsilon_{N} = 0$ (a,b) resp.\ $0.1$;   $s=0.5$, bottleneck size
$10^{2}$ (a,c) resp.\  $10^4$, probability of a bottleneck in any given generation $0.01$;
$h_{i}(g) =  2\one{g < 1}$ ie  the wild type is recessive, and
viability weight  $\exp\left( -s\left(  \left( h_{1}(g_{1} ) +
h_{2}(g_{2}) \right)/2 \right)^{2} \right) $
%% A :  14.0000000 754.1428571 337.5909808   0.4476486
%% B : bottleneck 1e4
%% B :    28.0000000 5722.6071429 2892.8004545    0.50550
%% bottleneck 1e2, p 0.01, e 0.1 
%% E : 8.0000000 360.6250000  76.4851760   0.212090 from 1e2 expr
from $10^{2}$ replicates. The excursions are shown as the frequency at time $t$ of
diploid individual homozygous for the beneficial type at both loci 
relative to the number of diploid individuals in the population at time $t$. }
\label{fig:twoloci1}
\end{figure}


\clearpage
\pagebreak
\newpage


\begin{figure}[htp]
\centering
\includegraphics[scale=1]{twolocigraphsCHD-crop}
\caption{$2N = 10^{6}$, $\Psi_{N} = 2N$,
$\alpha_{2} = 3$, $\varepsilon_{N} = 0$ (a,b) resp.\ $0.1$;   $s=0.5$, bottleneck size
$10^{2}$ (a,c) resp.\  $10^4$, probability of a bottleneck in any given generation $0.01$;
$h_{1}(g) =  2\one{g < 1}$  ie  the wild type is recessive,
$h_{2}(g) =  2\one{g = 0} + \one{g=1} $  and
viability weight  $\exp\left( -s\left(  \left( h_{1}(g_{1} ) +
h_{2}(g_{2}) \right)/2 \right)^{2} \right) $
%%
%% A :  bottleneck 1e2, p 0.01, e 0
%% A :  21.0000000 513.8095238 199.1059063   0.387509 from 1e2 expr
%% B : bottleneck 1e4, p 0.01, e 0, s 0.5
%% B : 21.0000000 3545.4761905 1435.6746365    0.404931 from 1e2 expr
%% bottleneck 1e2, p 0.01, e 0.1, alpha1 0.75,  s 0.5
%% E :   8.0000000 329.5000000 112.3336866   0.34092 from 1e2 expr
%% bottleneck 1e4, p 0.01, e 0.1, alpha 0.75, s 0.5
%% F : 2.0000000 385.5000000 154.8563851   0.4017 from 1e2 expr
%%
from $10^{2}$ replicates. The excursions are shown as the frequency at time $t$ of
diploid individual homozygous for the beneficial type at both loci 
relative to the number of diploid individuals in the population at time $t$.
}
\label{fig:chd}
\end{figure}



@* {\bf conclusion}. 
\label{sec:concl}





@* {\bf references}.
\label{sec:refs}

\bibliographystyle{plain}
\bibliography{/home/bjarki/verk/master_bibfile/refs}



@
\end{document}