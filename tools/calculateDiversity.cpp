/****************************************************************************
 * calculateDiversity.cpp                                                   *
 * Written by Patrick Reilly                                                *
 * Version 1.0 written 2020/01/01                                           *
 * Version 1.1 written 2020/01/07 Output sample sizes, remove D_a, add      *
 *                                shared_poly                               *
 * Version 1.2 written 2020/01/09 Bugfixes, including usable fraction no    *
 *                                longer output if site not usable          *
 * Version 1.3 written 2024/05/19 Bugfix for invalid pop files              *
 *                                                                          *
 * Description:                                                             *
 * This script takes in pseudoreference FASTAs and a TSV describing which   *
 *  pseudoreferences are from which population, and calculates various      *
 *  measures of diversity (H_exp, pi, Dxy, H_exp_total, pi_total, SP)       *
 * When indicated in the third column of the populations TSV, each          *
 *  pseudoreference may be individually treated as an inbred line, so that  *
 *  sample will have a single allele chosen at random for each site.        *
 * Note: The estimates are all meaninglessly 0 if at least one population   *
 *  has 0 alleles sampled.                                                  *
 *                                                                          *
 * Syntax: calculateDiversity [options]                                     *
 ****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <cctype>
#include <vector>
#include <map>
#include <sstream>
#include <array>
#include <set>
#include <unordered_map>
#include <algorithm>

//Define constants for getopt:
#define no_argument 0
#define required_argument 1
#define optional_argument 2

//Version:
#define VERSION "1.3"

//Define number of bases, shouldn't ever change:
#define NUM_BASES 4

//Usage/help:
#define USAGE "calculateDiversity\nUsage:\n calculateDiversity [options]\nOptions:\n -h,--help\tPrint this help\n -v,--version\tPrint the version of this program\n -d,--debug\tEnable debugging output to STDERR\n -p,--popfile\tTSV file of FASTA name, and population number\n\t\t(optional 3rd column is effective ploidy, 1 for inbred, 2 for outbred)\n -i,--inbred\tParse the third column of the pop TSV file\n -r,--prng_seed\tSet PRNG seed for random allele selection in inbred lines\n\t\tDefault: 42\n --usable_fraction,-u:\tFourth column represents fraction of unmasked bases\n"

using namespace std;

vector<string> splitString(string line_to_split, char delimiter) {
   vector<string> line_vector;
   string element;
   istringstream line_to_split_stream(line_to_split);
   while (getline(line_to_split_stream, element, delimiter)) {
      line_vector.push_back(element);
   }
   return line_vector;
}

bool openFASTAs(vector<ifstream*> &input_FASTAs, vector<string> input_FASTA_paths) {
   for (auto path_iterator = input_FASTA_paths.begin(); path_iterator != input_FASTA_paths.end(); ++path_iterator) {
      ifstream *input_FASTA = new ifstream(path_iterator->c_str(), ios::in);
      if (!input_FASTA->is_open()) { //Check to make sure the input file was validly opened
         cerr << "Error opening input FASTA: " << *path_iterator << "." << endl;
         return 0;
      }
      input_FASTAs.push_back(input_FASTA);
   }
   return 1;
}

void closeFASTAs(vector<ifstream*> &input_FASTAs) {
   for (auto FASTA_iterator = input_FASTAs.begin(); FASTA_iterator != input_FASTAs.end(); ++FASTA_iterator) {
      (*FASTA_iterator)->close();
      delete *FASTA_iterator;
   }
}

bool readFASTAs(vector<ifstream*> &input_FASTAs, vector<string> &FASTA_lines) {
   unsigned long which_input_FASTA = 0;
   bool ifstream_notfail;
   for (auto FASTA_iterator = input_FASTAs.begin(); FASTA_iterator != input_FASTAs.end(); ++FASTA_iterator) {
      if ((*FASTA_iterator)->fail()) {
         cerr << "Ifstream " << which_input_FASTA+1 << " failed." << endl;
         return 0;
      }
      string FASTA_line;
      ifstream_notfail = (bool)getline(**FASTA_iterator, FASTA_line);
      if (!ifstream_notfail) {
         return ifstream_notfail;
      }
      if (FASTA_lines.size() < input_FASTAs.size()) {
         FASTA_lines.push_back(FASTA_line);
      } else {
         FASTA_lines[which_input_FASTA++] = FASTA_line;
      }
   }
   return ifstream_notfail;
}

void countAlleles(vector<string> &FASTA_sequences, unsigned long i, map<unsigned long, unsigned long> &which_pop, vector<array<unsigned long, NUM_BASES+2>> &allele_counts, map<unsigned long, unsigned long> &inbred) {
   unsigned long num_sequences = FASTA_sequences.size();
   for (unsigned long j = 0; j < num_sequences; j++) {
      switch (FASTA_sequences[j][i]) {
         case 'A':
         case 'a':
            allele_counts[which_pop[j]-1][0] += inbred[j] == 1 ? 1 : 2; //Add 2 A alleles
            allele_counts[which_pop[j]-1][5] += inbred[j] == 1 ? 1 : 2; //Add 2 non-N alleles
            break;
         case 'C':
         case 'c':
            allele_counts[which_pop[j]-1][1] += inbred[j] == 1 ? 1 : 2; //Add 2 C alleles
            allele_counts[which_pop[j]-1][5] += inbred[j] == 1 ? 1 : 2; //Add 2 non-N alleles
            break;
         case 'G':
         case 'g':
            allele_counts[which_pop[j]-1][2] += inbred[j] == 1 ? 1 : 2; //Add 2 G alleles
            allele_counts[which_pop[j]-1][5] += inbred[j] == 1 ? 1 : 2; //Add 2 non-N alleles
            break;
         case 'K': //G/T het site
         case 'k':
            if (inbred[j] == 1) { //Randomly choose one of the alleles
               allele_counts[which_pop[j]-1][rand() <= (RAND_MAX-1)/2 ? 2 : 3]++;
            } else {
               allele_counts[which_pop[j]-1][2]++; //Add 1 G allele
               allele_counts[which_pop[j]-1][3]++; //Add 1 T allele
            }
            allele_counts[which_pop[j]-1][5] += inbred[j] == 1 ? 1 : 2; //Add 2 non-N alleles
            break;
         case 'M': //A/C het site
         case 'm':
            if (inbred[j] == 1) { //Randomly choose one of the alleles
               allele_counts[which_pop[j]-1][rand() <= (RAND_MAX-1)/2 ? 0 : 1]++;
            } else {
               allele_counts[which_pop[j]-1][0]++; //Add 1 A allele
               allele_counts[which_pop[j]-1][1]++; //Add 1 C allele
            }
            allele_counts[which_pop[j]-1][5] += inbred[j] == 1 ? 1 : 2; //Add 2 non-N alleles
            break;
         case 'R': //A/G het site
         case 'r':
            if (inbred[j] == 1) { //Randomly choose one of the alleles
               allele_counts[which_pop[j]-1][rand() <= (RAND_MAX-1)/2 ? 0 : 2]++;
            } else {
               allele_counts[which_pop[j]-1][0]++; //Add 1 A allele
               allele_counts[which_pop[j]-1][2]++; //Add 1 G allele
            }
            allele_counts[which_pop[j]-1][5] += inbred[j] == 1 ? 1 : 2; //Add 2 non-N alleles
            break;
         case 'S': //C/G het site
         case 's':
            if (inbred[j] == 1) { //Randomly choose one of the alleles
               allele_counts[which_pop[j]-1][rand() <= (RAND_MAX-1)/2 ? 1 : 2]++;
            } else {
               allele_counts[which_pop[j]-1][1]++; //Add 1 C allele
               allele_counts[which_pop[j]-1][2]++; //Add 1 G allele
            }
            allele_counts[which_pop[j]-1][5] += inbred[j] == 1 ? 1 : 2; //Add 2 non-N alleles
            break;
         case 'T':
         case 't':
            allele_counts[which_pop[j]-1][3] += inbred[j] == 1 ? 1 : 2; //Add 2 T alleles
            allele_counts[which_pop[j]-1][5] += inbred[j] == 1 ? 1 : 2; //Add 2 non-N alleles
            break;
         case 'W': //A/T het site
         case 'w':
            if (inbred[j] == 1) { //Randomly choose one of the alleles
               allele_counts[which_pop[j]-1][rand() <= (RAND_MAX-1)/2 ? 0 : 3]++;
            } else {
               allele_counts[which_pop[j]-1][0]++; //Add 1 A allele
               allele_counts[which_pop[j]-1][3]++; //Add 1 T allele
            }
            allele_counts[which_pop[j]-1][5] += inbred[j] == 1 ? 1 : 2; //Add 2 non-N alleles
            break;
         case 'Y': //C/T het site
         case 'y':
            if (inbred[j] == 1) { //Randomly choose one of the alleles
               allele_counts[which_pop[j]-1][rand() <= (RAND_MAX-1)/2 ? 1 : 3]++;
            } else {
               allele_counts[which_pop[j]-1][1]++; //Add 1 C allele
               allele_counts[which_pop[j]-1][3]++; //Add 1 T allele
            }
            allele_counts[which_pop[j]-1][5] += inbred[j] == 1 ? 1 : 2; //Add 2 non-N alleles
            break;
         default: //Assume that any case not handled here is an N
            allele_counts[which_pop[j]-1][4] += inbred[j] == 1 ? 1 : 2; //Add 2 N alleles
            break;
      }
   }
   return;
}

//A silly little optimization for memoization, since heterozygosity, pi,
// and Dxy will all evaluate the same for permutations of base identities,
// so we implement a sorting network as per:
//https://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-int-array
static inline void sort4_ulong(unsigned long * d) {
#define min(x, y) (x<y?x:y)
#define max(x, y) (x<y?y:x)
#define SWAP(x, y) { const unsigned long a = min(d[x], d[y]); \
                     const unsigned long b = max(d[x], d[y]); \
                     d[x] = a; d[y] = b; }
   SWAP(0, 1);
   SWAP(2, 3);
   SWAP(0, 2);
   SWAP(1, 3);
   SWAP(1, 2);
#undef SWAP
#undef min
#undef max
}

//Simple function just for debugging allele frequencies:
string join_freqs(array<double, NUM_BASES> &allele_freqs) {
   string output = "";
   if (allele_freqs.size() == NUM_BASES) {
      if (NUM_BASES > 0) {
         output += to_string(allele_freqs[0]);
         for (unsigned long i = 1; i < NUM_BASES; i++) {
            output += "," + to_string(allele_freqs[i]);
         }
      }
   }
   return output;
}

string join_counts(unsigned long *counts, const unsigned int n_counts) {
   string output = "";
   if (n_counts > 0) {
      output += to_string(counts[0]);
      for (unsigned int i = 1; i < n_counts; i++) {
         output += "," + to_string(counts[i]);
      }
   }
   return output;
}

string hetKey(array<unsigned long, NUM_BASES+2> &allele_counts, bool do_sort) {
   unsigned long counts[NUM_BASES] = {allele_counts[0], allele_counts[1], allele_counts[2], allele_counts[3]};
   if (do_sort) {
      sort4_ulong(counts);
   }
   string hetkey = join_counts(counts, NUM_BASES);
   return hetkey;
}

//A slight optimization here would be to sort the two populations,
// since order of populations doesn't matter to Dxy.
//We explicitly do not sort within populations, because order matters
string dxyKey(array<unsigned long, NUM_BASES+2> &pop1_allele_counts, array<unsigned long, NUM_BASES+2> &pop2_allele_counts) {
   string dxykey = hetKey(pop1_allele_counts, 0) + "," + hetKey(pop2_allele_counts, 0);
   return dxykey;
}

bool is_shared_poly(array<unsigned long, NUM_BASES+2> &pop1_allele_counts, array<unsigned long, NUM_BASES+2> &pop2_allele_counts) {
   unsigned int shared_poly = 0;
   //If any two alleles have non-zero product of frequencies across the two populations, the site is a shared polymorphism
   for (unsigned int i = 0; i < NUM_BASES; i++) {
      if (pop1_allele_counts[i] * pop2_allele_counts[i] > 0) {
         shared_poly++;
      }
   }
   return shared_poly > 1;
}

void calculateFreqs(array<unsigned long, NUM_BASES+2> &counts, array<double, NUM_BASES> &freqs) {
   //Zero out the initial frequencies:
   for (unsigned long i = 0; i < NUM_BASES; i++) {
      if (counts[NUM_BASES+1] > 0) {
         freqs[i] = (double)counts[i]/(double)counts[NUM_BASES+1];
      } else {
         freqs[i] = 0.0;
      }
   }
   return;
}

double calculateHet(array<double, NUM_BASES> &freqs) {
   double H_hat = 0.0;
   for (unsigned long i = 0; i < NUM_BASES-1; i++) {
      for (unsigned long j = i+1; j < NUM_BASES; j++) {
         H_hat += 2*freqs[i]*freqs[j];
      }
   }
   return H_hat;
}

double calculatePi(double H_hat, unsigned long nonN_count) {
   double pi_hat = H_hat;
   //Do the \frac{n}{n-1} correction, which now makes this Nei (1987) Eqn. 10.5
   if (nonN_count >= 2) { //Avoid divide-by-zero
      pi_hat *= (double)nonN_count/(double)(nonN_count-1);
   }
   return pi_hat;
}

double calculateDxy(array<double, NUM_BASES> &pop1_freqs, array<double, NUM_BASES> &pop2_freqs) {
   double dxy_hat = 0.0;
   for (unsigned long i = 0; i < NUM_BASES; i++) {
      for (unsigned long j = 0; j < NUM_BASES; j++) {
         if (i != j) { //d_{ij} = 1 for i != j, see Nei (1987) Eqn. 10.20
            //Temporary debugging:
            //cerr << "calculateDxy pre" << '\t' << i << ',' << j << '\t' << to_string(pop1_freqs[i]) << '\t' << to_string(pop2_freqs[j]) << '\t' << to_string(dxy_hat) << endl;
            dxy_hat += pop1_freqs[i]*pop2_freqs[j]; //\hat{x}_{i}\hat{y}_{j} in Nei (1987) Eqn. 10.20
            //Temporary debugging:
            //cerr << "calculateDxy post" << '\t' << i << ',' << j << '\t' << to_string(pop1_freqs[i]) << '\t' << to_string(pop2_freqs[j]) << '\t' << to_string(dxy_hat) << endl;
         }
      }
   }
   return dxy_hat;
}

//Nei's G_{ST} involves:
//H_{pooled} (H_{exp} when allele counts are combined between pops)
//H_{within} (average H_{exp} of the two pops)
//Hudson's K_{ST} involves:
//pi_{pooled} (Tajima's pi when allele counts are combined between pops)
//pi_{within} (weighted average of Tajima's pi of the two pops)
//The optimal choice of weight is w=n1/(n1+n2), and the weighting formula is:
//pi_w = w*pi_1 + (1-w)*pi_2 = (n1*pi_1 + n2*pi_2)/(n1+n2)
//Make sure to output n1 and n2 for use, although it's unclear
// how to decide the weight when doing ratio of averages.
//Still not sure why Hudson, Boos, and Kaplan didn't call their K pi...
//Hudson, Slatkin, and Maddison's F_{ST} involves:
//pi_{within} (average Tajima's pi of the two pops)
//pi_{between} (equivalent to D_{xy} between the populations)
//HSM call pi_{within} \hat{H}_{w}, and pi_{between} \hat{H}_{b}
// but the text is clear that they are pi and D_{xy}

void processScaffold(vector<string> &FASTA_headers, vector<string> &FASTA_sequences, map<unsigned long, unsigned long> &population_map, unsigned long num_pops, unordered_map<string, double> &memoized_het, unordered_map<string, double> &memoized_dxy, map<unsigned long, unsigned long> &inbred, unsigned int debug, bool usable) {
   cerr << "Processing scaffold " << FASTA_headers[0].substr(1) << " of length " << FASTA_sequences[0].length() << endl;
   //Do all the processing for this scaffold:
   //Polymorphism estimator: Given base frequencies at site:
   //\hat{\pi} = \(\frac{n}{n-1}\)\sum_{i=1}^{3}\sum_{j=i+1}^{4} 2\hat{p_{i}}\hat{p_{j}}
   //For biallelic sites, this reduces to the standard estimator: \(\frac{n}{n-1}\)2\hat{p}\hat{q}, since
   //p_{3} = p_{4} = 0
   //Dxy estimator is from Nei (1987) Eqn. 10.20 (\hat{d}_{XY} = \Sum_{i,j} \hat{x}_{i} \hat{y}_{j} d_{i,j}
   //D_{a} is a waste to calculate per-site, variances are too high
   unsigned long scaffold_length = FASTA_sequences[0].length();
   
   unsigned long num_pop_pairs = num_pops*(num_pops-1)/2;
   //Containers for various site statistics:
   array<unsigned long, NUM_BASES+2> init_base_frequency = { {0, 0, 0, 0, 0, 0} }; //Store the count of A, C, G, T, N, nonN for each site
   vector<array<unsigned long, NUM_BASES+2>> pop_site_counts; //Store population-specific allele counts
   pop_site_counts.reserve(num_pops);
   array<double, NUM_BASES> init_p_hats = { {0.0, 0.0, 0.0, 0.0} }; //Store the estimated allele frequencies for each site
   vector<array<double, NUM_BASES>> pop_p_hats; //Store population-specific estimated allele frequencies
   pop_p_hats.reserve(num_pops);
   vector<unsigned long> pop_ns; //Store population-specific haploid sample sizes
   pop_ns.reserve(num_pops);
   vector<double> pop_H_hats; //Store population-specific heterozygosity estimates
   pop_H_hats.reserve(num_pops);
   vector<double> pop_pi_hats; //Store population-specific polymorphism estimates
   pop_pi_hats.reserve(num_pops);
   vector<double> D_xys; //Store population pair-specific D_{xy} estimates (absolute, not net divergence)
   D_xys.reserve(num_pop_pairs);
   vector<bool> SPs; //Store shared polymorphisms
   SPs.reserve(num_pop_pairs);
   vector<double> H_totals; //Store population pair-specific H_{total} estimates
   H_totals.reserve(num_pop_pairs);
   vector<double> pi_totals; //Store population pair-specific \pi_{total} estimates for F_{ST}
   pi_totals.reserve(num_pop_pairs);
   
   for (unsigned long i = 0; i < scaffold_length; i++) {
      //Use site if all populations have at least 2 alleles:
      bool use_site = 1;
   
      //Initialize some of these containers:
      pop_site_counts.clear();
      pop_p_hats.clear();
      pop_ns.clear();
      pop_H_hats.clear();
      pop_pi_hats.clear();
      D_xys.clear();
      SPs.clear();
      H_totals.clear();
      pi_totals.clear();
      //Count up the alleles for each site:
      for (unsigned long j = 0; j < num_pops; j++) {
         pop_site_counts.push_back(init_base_frequency);
         pop_p_hats.push_back(init_p_hats);
         pop_ns.push_back(0);
         pop_H_hats.push_back(0.0);
         pop_pi_hats.push_back(0.0);
      }
      if (debug) {
         cerr << "Counting alleles for site " << i+1 << "." << endl;
      }
      countAlleles(FASTA_sequences, i, population_map, pop_site_counts, inbred);
      
      //For verbose debugging, output the allele counts:
      if (debug > 1) {
         for (unsigned long j = 0; j < num_pops; j++) {
            //Output comma-separated allele counts, and then
            // comma-separated N and non-N counts:
            cerr << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << j+1 << '\t' << pop_site_counts[j][0] << ',' << pop_site_counts[j][1] << ',' << pop_site_counts[j][2] << ',' << pop_site_counts[j][3] << '\t' << pop_site_counts[j][4] << ',' << pop_site_counts[j][5] << endl;
         }
      }
      //Calculate the total and population-specific allele frequencies:
      if (debug) {
         cerr << "Estimating allele frequencies for site " << i+1 << "." << endl;
      }
      unsigned long nonN_bases = 0;
      unsigned long total_bases = 0;
      for (unsigned long pop_index = 0; pop_index < num_pops; pop_index++) {
         //Only generate results if there is at least non-N alleles sampled
         // from each population:
         use_site = use_site && pop_site_counts[pop_index][NUM_BASES+1] >= 1;
         //Keep track of the number of non-N alleles vs. total alleles
         // to calculate the usable fraction of sequences at this site:
         nonN_bases += pop_site_counts[pop_index][NUM_BASES+1];
         total_bases += pop_site_counts[pop_index][NUM_BASES] + pop_site_counts[pop_index][NUM_BASES+1];
         //Keep track of the haploid sample size for this population:
         pop_ns[pop_index] = pop_site_counts[pop_index][NUM_BASES+1];
         //Determine the key for searching the memo:
         string het_key = hetKey(pop_site_counts[pop_index], 1);
         //Verbose debugging:
         if (debug > 1) {
            cerr << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << pop_index << '\t' << het_key << endl;
         }
         //Calculate allele frequencies for each population:
         calculateFreqs(pop_site_counts[pop_index], pop_p_hats[pop_index]);
         //Short-circuit the H and pi calculations if there aren't any
         // non-N alleles:
         if (pop_site_counts[pop_index][NUM_BASES+1] > 0) {
            //Calculate \pi_{i} for each population (or find in memo if seen before):
            if (memoized_het.find(het_key) != memoized_het.end()) {
               pop_H_hats[pop_index] = memoized_het[het_key];
               pop_pi_hats[pop_index] = calculatePi(memoized_het[het_key], pop_site_counts[pop_index][5]);
               //Verbose debugging/profiling:
               if (debug > 1) {
                  cerr << "H memo hit: " << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << het_key << '\t' << to_string(memoized_het[het_key]) << endl;
               }
            } else {
               //Calculate \pi_{i} for each population:
               if (debug) {
                  cerr << "Estimating pi for each population at site " << i+1 << "." << endl;
               }
               double H_hat = calculateHet(pop_p_hats[pop_index]);
               pop_H_hats[pop_index] = H_hat;
               pop_pi_hats[pop_index] = calculatePi(H_hat, pop_site_counts[pop_index][5]);
               //Store the calculated pi for that allele count distribution in the memo:
               memoized_het[het_key] = H_hat;
               //Verbose debugging
               if (debug > 1) {
                  cerr << "H memo miss: " << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << het_key << '\t' << to_string(pop_H_hats[pop_index]) << '\t' << to_string(pop_pi_hats[pop_index]) << '\t' << join_freqs(pop_p_hats[pop_index]) << endl;
               }
            }
         } else {
            pop_H_hats[pop_index] = 0.0;
            pop_pi_hats[pop_index] = 0.0;
            //Verbose debugging
            if (debug > 1) {
               cerr << "H pi no nonN: " << FASTA_headers[0].substr(1) << '\t' << i+1 << endl;
            }
         }
      }
      
      //Calculate D_{xy}, H_T, and pi_T for each pair of populations:
      if (debug) {
         cerr << "Estimating D_xy, H_T, and pi_T for site " << i+1 << "." << endl;
      }
      for (unsigned long pop1_index = 0; pop1_index < num_pops; pop1_index++) {
         for (unsigned long pop2_index = pop1_index+1; pop2_index < num_pops; pop2_index++) {
            //Short-circuit if either population has 0 non-N alleles
            if (pop_site_counts[pop1_index][NUM_BASES+1] > 0 && pop_site_counts[pop2_index][NUM_BASES+1] > 0) {
               //Search in memo for this pair of allele count distributions:
               string dxy_key = dxyKey(pop_site_counts[pop1_index], pop_site_counts[pop2_index]);
               if (memoized_dxy.find(dxy_key) != memoized_dxy.end()) {
                  D_xys.push_back(memoized_dxy[dxy_key]);
                  //Verbose debugging/profiling:
                  if (debug > 1) {
                     cerr << "Dxy memo hit: " << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << dxy_key << '\t' << to_string(memoized_dxy[dxy_key]) << endl;
                  }
               } else {
                  //If not found in memo, calculate it and store it in memo too:
                  double d_xy = calculateDxy(pop_p_hats[pop1_index], pop_p_hats[pop2_index]);
                  D_xys.push_back(d_xy);
                  memoized_dxy[dxy_key] = d_xy;
                  //Verbose debugging
                  if (debug > 1) {
                     cerr << "Dxy memo miss: " << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << dxy_key << '\t' << to_string(D_xys.size()) << '\t' << to_string(d_xy) << '\t' << join_freqs(pop_p_hats[pop1_index]) << '\t' << join_freqs(pop_p_hats[pop2_index]) << endl;
                  }
               }
               //For now, we won't memoize shared polymorphism calculations:
               SPs.push_back(is_shared_poly(pop_site_counts[pop1_index], pop_site_counts[pop2_index]));
               //Combine the allele counts for use in calculating \pi_{total}:
               array<unsigned long, NUM_BASES+2> total_base_counts = { {0, 0, 0, 0, 0, 0} }; //Store the count of A, C, G, T, N, nonN for each site totalled for the two populations
               transform(pop_site_counts[pop1_index].begin(), pop_site_counts[pop1_index].end(), pop_site_counts[pop2_index].begin(), total_base_counts.begin(), plus<unsigned long>());
               //Check the memo for this allele count distribution, else
               // calculate pi from the frequencies:
               string het_t_key = hetKey(total_base_counts, 1);
               double het_total, pi_total;
               if (memoized_het.find(het_t_key) != memoized_het.end()) {
                  het_total = memoized_het[het_t_key];
                  pi_total = calculatePi(het_total, total_base_counts[NUM_BASES+1]);
                  //Verbose debugging/profiling:
                  if (debug > 1) {
                     cerr << "H_T memo hit: " << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << het_t_key << '\t' << to_string(memoized_het[het_t_key]) << endl;
                  }
               } else {
                  //Calculate \pi_{total} for each population pair:
                  if (debug) {
                     cerr << "Estimating pi_total for each population pair at site " << i+1 << "." << endl;
                  }
                  array<double, NUM_BASES> total_freqs;
                  calculateFreqs(total_base_counts, total_freqs);
                  het_total = calculateHet(total_freqs);
                  pi_total = calculatePi(het_total, total_base_counts[NUM_BASES+1]);
                  memoized_het[het_t_key] = het_total;
                  //Verbose debugging
                  if (debug > 1) {
                     cerr << "H_T memo miss: " << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << het_t_key << '\t' << to_string(het_total) << '\t' << to_string(pi_total) << '\t' << join_freqs(total_freqs) << endl;
                  }
               }
               H_totals.push_back(het_total);
               pi_totals.push_back(pi_total);
            } else {
               D_xys.push_back(0.0);
               SPs.push_back(0.0);
               H_totals.push_back(0.0);
               pi_totals.push_back(0.0);
               //Verbose debugging
               if (debug > 1) {
                  cerr << "Dxy SP H_T pi_T no nonN: " << FASTA_headers[0].substr(1) << '\t' << i+1 << endl;
               }
            }
         }
      }
      //Calculate the usable fraction of alleles for this site:
      double usable_fraction = (double)nonN_bases/(double)total_bases;
      
      if (use_site) {
         //Output elements: Scaffold, position, pi_{1}, omit site, H_{i}, pi_{i}, n_{i}, D_{ij}, SP_{ij}, H_{T}, and pi_{T} values
         cout << FASTA_headers[0].substr(1) << '\t' << i+1;
         //Output pi_{1} and omit_site:
         if (!usable) { //If we don't want to output the usable fraction, just output 0
            cout << '\t' << pop_pi_hats[0] << '\t' << "0";
         } else {
            cout << '\t' << pop_pi_hats[0] << '\t' << usable_fraction;
         }
         //Output H_{i}, \pi_{i}, and n_{i} values:
         for (unsigned long j = 0; j < num_pops; j++) {
            cout << '\t' << pop_H_hats[j];
            cout << '\t' << pop_pi_hats[j];
            cout << '\t' << pop_ns[j];
         }
         //Output D_{XY}, SP, H_T, and pi_T values:
         unsigned long num_pairs = D_xys.size();
         for (unsigned long pair_index = 0; pair_index < num_pairs; pair_index++) {
            cout << '\t' << D_xys[pair_index];
            cout << '\t' << SPs[pair_index];
            cout << '\t' << H_totals[pair_index];
            cout << '\t' << pi_totals[pair_index];
         }
         cout << endl;
      } else { //Do not output any estimators where any pop has no alleles
         if (!usable) { //If we don't want to output the usable fraction, just output 1
            cout << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << "0" << '\t' << 1;
         } else { //We output 0 for the usable fraction to not give false info
            //This basically is acting as a flag when 0
            cout << FASTA_headers[0].substr(1) << '\t' << i+1 << '\t' << "0" << '\t' << "0";
         }
         for (unsigned long j = 1; j <= num_pops; j++) {
            cout << '\t' << "0"; //Output 0 (NA)s for H_{i}
            cout << '\t' << "0"; //Output 0 (NA)s for \pi_{i}
            cout << '\t' << "0"; //Output 0 (NA)s for n_{i}
         }
         for (unsigned long j = 1; j <= num_pops; j++) {
            for (unsigned long k = j+1; k <= num_pops; k++) {
               cout << '\t' << "0"; //Output 0 (NA)s for D_{XY}
               cout << '\t' << "0"; //Output 0 (NA)s for SP
               cout << '\t' << "0"; //Output 0 (NA)s for H_total
               cout << '\t' << "0"; //Output 0 (NA)s for pi_total
            }
         }
         cout << endl;
      }
   }
}

int main(int argc, char **argv) {
   //Variables for processing the FASTAs:
   vector<string> input_FASTA_paths;
   vector<ifstream*> input_FASTAs;
   
   //Path to file describing which indivs are in which populations:
   string popfile_path;
   
   //Handle inbred pseudoreferences as haploids:
   bool parse_inbred = 0;
   unsigned int prng_seed = 42;

   //Option to output debugging info on STDERR
   unsigned int debug = 0;

   //Option to output fraction of usable sites:
   bool usable = 0;
   
   //Variables for getopt_long:
   int optchar;
   int structindex = 0;
   extern int optind;
   //Create the struct used for getopt:
   const struct option longoptions[] {
      {"popfile", required_argument, 0, 'p'},
      {"inbred", no_argument, 0, 'i'},
      {"prng_seed", required_argument, 0, 'r'},
      {"usable_fraction", no_argument, 0, 'u'},
      {"debug", no_argument, 0, 'd'},
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'}
   };
   //Read in the options:
   while ((optchar = getopt_long(argc, argv, "p:ir:udvh", longoptions, &structindex)) > -1) {
      switch(optchar) {
         case 'p':
            cerr << "Using population TSV file " << optarg << endl;
            popfile_path = optarg;
            break;
         case 'i':
            cerr << "Reading in inbreeding status from third column of population TSV." << endl;
            parse_inbred = 1;
            break;
         case 'r':
            cerr << "Setting PRNG seed to " << optarg << endl;
            prng_seed = atoi(optarg);
            break;
         case 'u':
            cerr << "Outputting fraction of usable sites rather than omit column" << endl;
            usable = 1;
            break;
         case 'd':
            cerr << "Outputting debug information." << endl;
            debug += 1;
            break;
         case 'v':
            cerr << "calculateDiversity version " << VERSION << endl;
            return 0;
            break;
         case 'h':
            cerr << USAGE;
            return 0;
            break;
         default:
            cerr << "Unknown option " << (unsigned char)optchar << " supplied." << endl;
            cerr << USAGE;
            return 1;
            break;
      }
   }
   
   //Set the seed of the PRNG:
   srand(prng_seed);
   
   //Open the population TSV file:
   ifstream pop_file;
   pop_file.open(popfile_path);
   if (!pop_file) {
      cerr << "Error opening population TSV file " << popfile_path << endl;
      return 9;
   }
   
   //Read the population members into a map from FASTA to population:
   map<unsigned long, unsigned long> population_map;
   set<unsigned long> populations;
   map<unsigned long, unsigned long> inbred;
   string popline;
   if (debug) {
      cerr << "Reading population TSV file." << endl;
   }
   unsigned long fasta_index = 0;
   while (getline(pop_file, popline)) {
      vector<string> line_vector;
      line_vector = splitString(popline, '\t');
      if ((parse_inbred && line_vector.size() < 3) || (!parse_inbred && line_vector.size() < 2)) {
         cerr << "Invalid population TSV file, either not enough columns or the delimiter is not a tab." << endl;
         pop_file.close();
         return 11;
      }
      try {
         population_map[fasta_index] = stoul(line_vector.at(1));
      } catch (const invalid_argument& e) {
         cerr << "Invalid population ID in second column of population TSV." << endl;
         cerr << "Must be a positive integer." << endl;
         pop_file.close();
         return 12;
      }
      if (parse_inbred) {
         try {
            inbred[fasta_index] = stoul(line_vector.at(2));
         } catch (const invalid_argument& e) {
            cerr << "Invalid inbred boolean in third column of population TSV." << endl;
            cerr << "Must be a non-negative integer." << endl;
            pop_file.close();
            return 10;
         }
      } else { //Default to treating as normal diploid
         inbred[fasta_index] = 2;
      }
      fasta_index++;
      input_FASTA_paths.push_back(line_vector.at(0));
      populations.insert(stoul(line_vector.at(1)));
   }
   pop_file.close();
   unsigned long num_pops = populations.size();
   if (debug) {
      cerr << "Read in " << num_pops << " populations." << endl;
   }
   //Output the header line:
   if (!usable) {
      cout << "Scaffold" << '\t' << "Position" << '\t' << "pi_1" << '\t' << "omit_position";
   } else {
      cout << "Scaffold" << '\t' << "Position" << '\t' << "pi_1" << '\t' << "site_weight";
   }
   for (unsigned long i = 1; i <= num_pops; i++) {
      cout << '\t' << "H_" << i; //H_w can be calculated post-hoc from these
      cout << '\t' << "pi_" << i; //pi_w can be calculated post-hoc from these
      cout << '\t' << "n_" << i; //Can be used for calculating weight for K_{ST}
   }
   for (unsigned long i = 1; i <= num_pops; i++) {
      for (unsigned long j = i+1; j <= num_pops; j++) {
         cout << '\t' << "D_" << i << ',' << j; //aka pi_{between}
         cout << '\t' << "SP_" << i << ',' << j; //Shared polymorphisms
         cout << '\t' << "H_T_" << i << "," << j; //aka H_{pooled}
         cout << '\t' << "pi_T_" << i << "," << j; //aka pi_{pooled}
      }
   }
   cout << endl;
   
   //Open the input FASTAs:
   bool successfully_opened = openFASTAs(input_FASTAs, input_FASTA_paths);
   if (!successfully_opened) {
      closeFASTAs(input_FASTAs);
      cerr << "Unable to open at least one of the FASTAs provided." << endl;
      return 2;
   }
   cerr << "Opened " << input_FASTAs.size() << " input FASTA files out of " << input_FASTA_paths.size() << " paths provided." << endl;
   
   //Set up the vector to contain each line from the n FASTA files:
   vector<string> FASTA_lines;
   FASTA_lines.reserve(input_FASTA_paths.size());
   
   //Set up maps to memoize pi and Dxy:
   unordered_map<string, double> memoized_het;
   unordered_map<string, double> memoized_dxy;

   //Iterate over all of the FASTAs synchronously:
   vector<string> FASTA_headers;
   FASTA_headers.reserve(input_FASTA_paths.size());
   vector<string> FASTA_sequences;
   FASTA_sequences.reserve(input_FASTA_paths.size());
   while (readFASTAs(input_FASTAs, FASTA_lines)) {
      //Check if we're on a header line:
      bool all_header_lines = 1;
      bool any_header_lines = 0;
      for (auto line_iterator = FASTA_lines.begin(); line_iterator != FASTA_lines.end(); ++line_iterator) {
         all_header_lines = all_header_lines && ((*line_iterator)[0] == '>');
         any_header_lines = any_header_lines || ((*line_iterator)[0] == '>');
      }
      if (all_header_lines) { //Process the previous scaffold if at a header
         if (!FASTA_sequences.empty()) {
            processScaffold(FASTA_headers, FASTA_sequences, population_map, num_pops, memoized_het, memoized_dxy, inbred, debug, usable);
            FASTA_sequences.clear();
         }
         FASTA_headers = FASTA_lines;
         for (auto header_iterator = FASTA_headers.begin()+1; header_iterator != FASTA_headers.end(); ++header_iterator) {
            if (*header_iterator != FASTA_headers[0]) {
               cerr << "Error: FASTAs are not synchronized, headers differ." << endl;
               cerr << FASTA_headers[0] << endl;
               cerr << *header_iterator << endl;
               closeFASTAs(input_FASTAs);
               return 3;
            }
         }
         //For simplicity, we will trim FASTA headers to ignore anything
         // after a space:
         for (auto header_iterator = FASTA_headers.begin(); header_iterator != FASTA_headers.end(); ++header_iterator) {
            *header_iterator = header_iterator->substr(0, header_iterator->find(" "));
         }
      } else if (any_header_lines) {
         cerr << "Error: FASTAs are not synchronized, or not wrapped at the same length." << endl;
         closeFASTAs(input_FASTAs);
         return 4;
      } else { //Read in the sequences if no files are at header lines
         unsigned long sequence_index = 0;
         for (auto line_iterator = FASTA_lines.begin(); line_iterator != FASTA_lines.end(); ++line_iterator) {
            if (FASTA_sequences.size() < FASTA_lines.size()) { //If vector is empty, push_back
               FASTA_sequences.push_back(*line_iterator);
            } else { //If vector isn't empty, append to existing sequence
               FASTA_sequences[sequence_index++] += *line_iterator;
            }
         }
      }
   }
   
   //Catch any IO errors that kicked us out of the while loop:
   unsigned long which_input_FASTA = 0;
   for (auto FASTA_iterator = input_FASTAs.begin(); FASTA_iterator != input_FASTAs.end(); ++FASTA_iterator) {
      if ((*FASTA_iterator)->bad()) {
         cerr << "Error reading input FASTA: " << input_FASTA_paths[which_input_FASTA] << endl;
         cerr << "Fail bit: " << ((*FASTA_iterator)->rdstate() & ifstream::failbit) << " Bad bit: " << (*FASTA_iterator)->bad() << " EOF bit: " << (*FASTA_iterator)->eof() << endl;
         closeFASTAs(input_FASTAs);
         return 5;
      }
      which_input_FASTA++;
   }
   //If no errors kicked us out of the while loop, process the last scaffold:
   processScaffold(FASTA_headers, FASTA_sequences, population_map, num_pops, memoized_het, memoized_dxy, inbred, debug, usable);
   
   //Close the input FASTAs:
   closeFASTAs(input_FASTAs);
   
   return 0;
}
