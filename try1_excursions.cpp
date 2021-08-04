/* *****************************
Copyright (C) 2021 Bjarki Eldon


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

******************* */

/* *********
the optimum is defined as all genotypes are 0/1; 
i.e. fixed for the heterozygote genotype 0/1 at all chromosomes 
so if 2N diploid individuals then 4N gene copies at each chromosome,
and  we define Y = \sum_chromosomes #{type 1 at chromosome};
but need to track at all L loci, so is multidimensional (Y_1, ..., Y_L), 
since if lose type at any locus then stop
 * *********** */

/* ****************   compile 
c++ -Wall -m64 -O3 -DNDEBUG -mtune=corei7 -march=native try2_excursions.cpp -lm -lgsl -lgslcblas
* ************** */

#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <list>
#include <string>
#include <fstream>
#include <forward_list>
#include <assert.h>
#include <cmath>
/* **********
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/beta_distribution.hpp>
************ */
#include <gsl/gsl_matrix.h>

using namespace std;

 /* obtain a seed out of thin air for the random number engine */
  std::random_device randomseed;
  /* Standard mersenne twister  random number engine seeded with rng() */
  std::mt19937_64 rng(randomseed());
 /* construct a uniform on the unit interval; call with uni01( rng ) */
std::uniform_real_distribution<double> uni01(0.0, 1.0);
/* construct a random exponential */
/*  timi = timi + ( -log(1. - uni01(rng) )/(ct) ); */

/* the main containers are organised in |struct G|. */
struct G {
private : 
  /* popsize is N; there are 2N diploid individuals in the population */
  static const size_t popsize = 500;
  static const size_t cutoff =  2*popsize;
  const double __dcutoff = static_cast<double>( cutoff); 
  static const size_t chroms = 1;
  static constexpr double dchroms = static_cast<double>(chroms); 
  static constexpr double s = 10.0;
  /* should have  alphaone <= alphatwo */ 
  static constexpr double alphaone = 0.01;
  static constexpr double alphatwo = 3.0;
  static constexpr double epsilon = 1.0/static_cast<double>(popsize);
  const string excursions_file_name =  "excursions_N" + std::to_string(popsize) + "_s" + std::to_string(s) + "_alphaone" + std::to_string(alphaone) + "_alphatwo" + std::to_string(alphatwo) + "_resout.txt" ;
  size_t goodtypes {} ;
  size_t losttypes {} ;
  size_t alphaindex {} ;
  /* count individuals with all loci of best genotype */
  std::vector<size_t> counts {};
  /* count types at each locus to check if lost type */
  std::vector<size_t> counttypes {};
  /* matrix for the genotypes of the population; rows are the chromosomes, and columns the diploid individuals; each genotype stored as $0 = 0/0$, $1 = 0/1$, $2 = 1/1$ */
  gsl_matrix_uint * pop = gsl_matrix_uint_calloc( chroms, 2*popsize);
  // std::list< std::vector<size_t> >::iterator it;
  std::vector<size_t> indexes {};
  /* list of juveniles with genotypes */
  std::list<std::vector<size_t>> juveniles ;
  std::vector< double> z {};
  std::vector< size_t > ellindex {};

  /* excursion records excursion to fixation of heterozygous genotype of a single locus */
  std::vector<size_t> excursion {} ;

   // = (double *)calloc( cutoff+1, sizeof(double));
  std::vector<double> cdfxialphaone {};
  std::vector<double> cdfxialphatwo {}; 
  
/* the |public|  functions for  working with the containers */
public :

     void init(){

       excursion.clear();
       excursion.resize(0);
       excursion.shrink_to_fit();
       std::vector<size_t>().swap( excursion);
       
       cdfxialphatwo.clear();
       cdfxialphatwo.resize(0);
	 cdfxialphatwo.shrink_to_fit();
       std::vector<double>().swap( cdfxialphatwo);
       cdfxialphatwo.reserve( cutoff + 1);
       cdfxialphatwo.assign( cutoff+1, 0.0);

       cdfxialphaone.clear();
       cdfxialphaone.resize(0);
       cdfxialphaone.shrink_to_fit();
       std::vector<double>().swap( cdfxialphaone);
       cdfxialphaone.reserve(cutoff+1);
       cdfxialphaone.assign(cutoff+1, 0.0);
    counttypes.resize(0);
    counttypes.shrink_to_fit();
    counttypes.reserve(chroms);
    counttypes.assign( chroms, 0);
    counts.resize(0);
    counts.shrink_to_fit();
    counts.resize(chroms); 
    counts.assign( chroms , 0);
    indexes.resize(0);
    indexes.shrink_to_fit();
    indexes.reserve(2*popsize); 
    indexes.assign( 2*popsize, 0);
    ellindex.resize(0);
    ellindex.shrink_to_fit();
    ellindex.reserve(chroms); 
    ellindex.assign( chroms, 0);
    std::iota( std::begin(indexes), std::end(indexes), 0);
    std::iota( std::begin(ellindex), std::end(ellindex), 0);
    /* |for( const unsigned& i: ellindex){}|  */
    /* the number of juveniles varies from generation to generation */ }

/* initialise the population */
 void initG(  )
{
  gsl_matrix_uint_set_zero( pop);
  /* set only one individual as  heteroz at all loci; all others homozygos 0/0 */
  /* was N individuals homozy 1/1 at all loci; other N homoz 0/0 */
  for( unsigned j = 0 ; j < 1; ++j){
  for( const size_t& ell: ellindex){
    /* ell ranges from chroms - 1, ... , 0 */
    gsl_matrix_uint_set( pop, ell, j + ell + 1, 1); }}
}


/* sample a genotype given the parent genotypes */
size_t  sample_genotype( const  unsigned gone, const  unsigned gtwo )
{
  /* assign genotype at a locus */
  /* |gone| is genotype of one parent */
  /* |gtwo| is genotype of other parent */
  
  size_t  genotype = 0; 
  const double u = uni01(rng); 
  switch( gone) {
  case 0 : {
    /* 0/0 :  */
    genotype = (gtwo < 1 ? 0 : ( gtwo < 2 ? (u < 0.5 ? 0 : 1) : 1 )) ;
    break ; }
  case 1 : {
    /* 0/1 : */
    genotype = ( gtwo < 1 ? (u < 0.5 ? 0 : 1) : (gtwo < 2 ? ( u < 0.5 ? 1 : (u < 0.75 ? 0 : 2)) : ( u < 0.5 ? 1 : 2) ) ) ;
    break ;}
  case 2 : {
    /* 1/1 : */
    genotype = ( gtwo < 1 ? 1 : ( gtwo < 2 ? (uni01(rng) < .5 ? 1 : 2) : 2) ); 
    break ; }
  default :
    break ;}
  assert( genotype < 3);
  return ( genotype );
}

static  double rexp( const double lambda )
{
  assert( lambda > 0);
  /*  \ref{SEC:included}   */
  return( -log( 1. - uni01(rng))/lambda );
}

/* loop over loci and assign genotypes to juveniles  */
 void assigngenome( const unsigned pone, const  unsigned ptwo )
{
  /* |pone| and |ptwo| are indexes of parents */
  std::vector<size_t> g {};
  g.reserve(chroms); 
  g.assign( chroms, 0);
  size_t sumg = 0 ;
  /* loop over the chromosomes and sample genotypes */
  for( const size_t& ell: ellindex){
      g[ell] =  sample_genotype( gsl_matrix_uint_get( pop, ell, pone), gsl_matrix_uint_get( pop, ell, ptwo) );
      /* compute the genotype value */
      // sumg += ( g[ell] < 2 ? 0.0 : 2.0 );
      sumg +=  g[ell]; }
  
  /* add a juvenile with genome g */
  juveniles.push_back(g);
  /* record the $z$ score $exp( -s(znull - z)^2 )$ for juvenile */
  /* now optimum trait value is $z0 = 1$; heteroz $0/1$ at all loci */
  z.push_back( rexp( exp( ( -s )*pow( 1.0  -  ( static_cast<double>(sumg) / dchroms ), 2. ) ) ) ) ;
}


/* sample a random number of juveniles */
size_t sample_random_number_juveniles(  )
{
  size_t j = 2;
  const double __rU = uni01( rng ) ;
  while( __rU  > ( alphaindex < 1 ? cdfxialphaone[j] : cdfxialphatwo[j]) ){
    ++j ; }
    assert( j > 1);
    assert( j <= cutoff);
  return( j );
}

/* compute the 2Nth largest trait value among  SN >= 2N trait values using partial sorting */
  double nthelmnt( )
{
  assert( z.size()  >=  2*popsize );
  ( std::nth_element( z.begin(), z.begin() + ((2*popsize) -1 ), z.end()) );

  return( z[ ( (2*popsize) - 1) ] );
}


/* probability of a given number of juveniles; probability is monotonically decreasing;  */
  double probxi( const double k,  const double a)
{
 /* $\prb{X = k} =  C\left( \frac{1}{k^\alpha} -  \frac{1}{(1+k)^\alpha} \right), \quad 2 \le k \le u_N$ where $u_N$ is the cutoff */
 /* and $\prb{X = k} = 0$ for $k \notin \{2, 3, \ldots, u_N\}$ */
  return( (pow(  1./k, a) - pow( 1./(k+1.), a))/( pow( .5, a)  -  pow( 1./(__dcutoff + 1.), a) ) );
}


/* initializing  probability array for sampling juveniles */
 void generate_cdf_for_sampling_juveniles( )
{
  size_t i = 1;
  while( ++i < cutoff ){
    assert( i > 1);
    assert( i <= cutoff);
    assert( cdfxialphatwo[i] == 0.0);
    cdfxialphatwo[ i ] = cdfxialphatwo[ i-1 ]  +  probxi( static_cast<double>(i), alphatwo);
    assert( cdfxialphatwo[i] > 0.0 );
    assert( cdfxialphatwo[i] <= 1.0);
    assert( cdfxialphaone[i] == 0.0);
    cdfxialphaone[ i ] = cdfxialphaone[ i-1 ]  +  probxi( static_cast<double>(i), alphaone);
    assert( cdfxialphaone[i] > 0.0);
    assert( cdfxialphaone[i] <= 1.0);
  }
  cdfxialphaone[cutoff] = 1.0;
  cdfxialphatwo[cutoff] = 1.0;
}

/* check if lost type at any locus  from population */
/* \label{func:lossoftype} */
  bool lossoftype( )
{
  /* self->popsize is N; number of diploid individuals is 2N; max number of type 1 is 4N at any chromosome (locus) */
  const size_t  k = 4*(popsize);
  /* if 1/1 best then only applies if lost type 1 */
  /*  */
  /* if lost type 1 at a locus from the population then all are 0/0 at that locus */
  /* now heteroz are optimal so check if lost type 1 or 0 at  any chromosome */
  /* counttypes is count of type 1 alleles at each chromosome */
  /* optimum is homoz 1/1 at all loci so only need to check if lost type 1 at any locus */
  return( std::any_of( counttypes.begin(), counttypes.end(), [k]( size_t i){ return ( (i < 1) || (i == k) ); } ) ) ;
}


/* count types among surviving juveniles, i.e. juveniles with trait value less than reference trait value */
size_t count_types( const double nth )
{

  size_t ell {};
  size_t  i {};
  size_t j = 0;
  size_t hlocus {};
  
  i = 0;
  assert( nth > 0.0);

  /* counts counts favored genotypes */
  std::fill( counts.begin(), counts.end(), 0);
  std::fill( counttypes.begin(), counttypes.end(), 0);
  /* count total number of sampled juveniles  with  all good genotypes */
  /* optimal genotype is heterozygote 0/1 */
  goodtypes = 0; 
  for( const std::vector<size_t>& v: juveniles){
    /* check if juvenile is sampled */
    if ( z[j] <= nth ){
      ell = 0;
      hlocus = 0;
      for( const size_t& g: v){
	/* g is genotype at juvenile v */
	/* when heterozygous  is better */
	/* count number of 0/1 loci; stop when all 2NChroms loci are 0/1  */
	/* stop when all chroms are genotype 1 equiv  heteroz 0/1 */
	hlocus += ( g == 1 ? 1 : 0);
	/* count type at locus to check if lost types; counting among surviving juveniles; if lost genotype among juveniles then automatically have lost it in the population */
	/* counts number type 1 {0, 1, or 2 }; adds up number of type 1 among surviving juveniles */
	counttypes[ell] += g; 
	/* count 0/1  genotypes at locus */
	/* counts[ell] += (g == 1 ? 1 : 0) ; */
	gsl_matrix_uint_set( pop, ell, i, g);
	++ell;}
      /* counting number of individuals with all genotypes 0/1 */
      goodtypes += (hlocus <  chroms ? 0 : 1);
      ++i;}
    if( j > popsize*2 ){ break; }
    ++j; }
  // assert( i == (popsize)*2);
  return goodtypes ;
}

/* update population - sample juveniles, assign genotypes,  sort the juveniles, and record the surviving juveniles  */
/* returns number of individuals heteroz */
 size_t  updatepop( )
{

  /* shuffle the parent indexes 0 to 2N */
  std::shuffle( indexes.begin(), indexes.end(), rng);
  int i = popsize ;
  size_t Xi {};

  size_t number_individuals_heteroz {} ;

  z.clear();
  /* 
  z.resize(0);
  z.shrink_to_fit();
  z.reserve( 10*2*(popsize));
  */
  for( auto &j: juveniles){
    j.clear();
    j.resize(0);
    j.shrink_to_fit();
    std::vector<size_t>().swap(j);}
  juveniles.resize(0);
  std::list< std::vector<size_t> >().swap(juveniles);

  /* toss a coin for  choosing between alphas */
  /* if $alphaindex = 0$ then use smaller alphaone; otherwise alphatwo */
  alphaindex = ( uni01(rng) < epsilon ? 0 : 1); 
  while( --i >= 0){
    assert( i >= 0); 
    /* pair ( i, i + N) produce juveniles */
    Xi =  sample_random_number_juveniles(  ); 
    /* assign chromosomes to Xi juveniles */
    assert( Xi > 1);
    while( Xi-- > 0 ){
      assigngenome( indexes[i], indexes[i + popsize ] );}
    }
  /* partial sort for  the 2N-th element of the z scor and return the 2Nth largest element */
  // nth =  nthelmnt( );
  /* update pop with 2N new juveniles */
  /* init the counts of types and heteroz loci */
  // std::fill( self->counts.begin(), self->counts.end(), 0);
  number_individuals_heteroz =   count_types( nthelmnt() );

  return number_individuals_heteroz ;
}


/* run the simulations using function within struct */
void run_within_struct()
{

  init();
  generate_cdf_for_sampling_juveniles( );

  size_t number_heteroz {} ;

  size_t timi = 0;
  int r = 1e3 + 1;
  while( --r > 0){
    timi = 0; 
  initG();
  excursion.clear();
  /* starting with one heteroz individual */
  do{
    excursion.push_back( number_heteroz );
    ++timi; 
    number_heteroz = updatepop();}
  while(  ( !lossoftype() ) &&  ( (goodtypes) <  (popsize)*2 ) );
  std::cout << (goodtypes < (popsize) * 2 ? 0 : 1) << ' ' << timi << '\n';
  /* print excursion to file if fixation */
  if( goodtypes <( (popsize) * 2 ) ){}
  else{
    /* all individuals with optimal heterozygous genotype */
    std::ofstream excursion_file( excursions_file_name, std::ofstream::out | std::ofstream::app );
    for ( size_t k = 0 ; k < timi ; ++k){
      excursion_file << excursion[k] << ' ' ; }
  }
  
  }
  freememory();
}

/* \newline clear the memory of the containers */
    void freememory()
    {
      excursion.clear();
      excursion.resize(0);
      excursion.shrink_to_fit();
      std::vector<size_t>().swap( excursion);

      
      gsl_matrix_uint_free( pop);
      indexes.resize(0);
      indexes.shrink_to_fit();
      std::vector<size_t>().swap(indexes);
      for( auto&j: juveniles){
	j.clear();
	j.resize(0);
	j.shrink_to_fit();
	std::vector<size_t>().swap(j);}
      juveniles.resize(0);
      std::list< std::vector< size_t > >().swap(juveniles);
      z.resize(0);
      z.shrink_to_fit();
      std::vector<double>().swap(z);
      ellindex.resize(0);
      std::vector<size_t >().swap(ellindex);
      counts.resize(0);
      std::vector<size_t>().swap(counts);
      counttypes.resize(0);
      std::vector<size_t>().swap(counttypes);

    
      cdfxialphatwo.clear();
      cdfxialphatwo.resize(0);
	 cdfxialphatwo.shrink_to_fit();
       std::vector<double>().swap( cdfxialphatwo);

       cdfxialphaone.clear();
       cdfxialphaone.resize(0);
       cdfxialphaone.shrink_to_fit();
       std::vector<double>().swap( cdfxialphaone);
    
    }
} ;


static void run()
{
  std::unique_ptr< G> __uptr = std::make_unique<G>();

  /* initialise containers */
  // __uptr->init();
  // __uptr->generate_cdf_for_sampling_juveniles( );
  __uptr->run_within_struct(); 

/* clear |__uptr| */
// __uptr->freememory() ;
 std::unique_ptr<G>().swap(__uptr);
 __uptr.reset();
 assert( !__uptr);
}

int main()
{

run();

return GSL_SUCCESS; 
}

