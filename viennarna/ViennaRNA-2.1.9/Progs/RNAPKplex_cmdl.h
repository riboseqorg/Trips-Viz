/** @file RNAPKplex_cmdl.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.5
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef RNAPKPLEX_CMDL_H
#define RNAPKPLEX_CMDL_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef PKPLEX_CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define PKPLEX_CMDLINE_PARSER_PACKAGE "RNAPKplex"
#endif

#ifndef PKPLEX_CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define PKPLEX_CMDLINE_PARSER_PACKAGE_NAME "RNAPKplex"
#endif

#ifndef PKPLEX_CMDLINE_PARSER_VERSION
/** @brief the program version */
#define PKPLEX_CMDLINE_PARSER_VERSION VERSION
#endif

/** @brief Where the command line options are stored */
struct PKplex_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *detailed_help_help; /**< @brief Print help, including all details and hidden options, and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  float cutoff_arg;	/**< @brief Report only base pairs with an average probability > cutoff in the dot plot
  
 (default='0.01').  */
  char * cutoff_orig;	/**< @brief Report only base pairs with an average probability > cutoff in the dot plot
  
 original value given at command line.  */
  const char *cutoff_help; /**< @brief Report only base pairs with an average probability > cutoff in the dot plot
  
 help description.  */
  double temp_arg;	/**< @brief Rescale energy parameters to a temperature of temp C. Default is 37C.
  
.  */
  char * temp_orig;	/**< @brief Rescale energy parameters to a temperature of temp C. Default is 37C.
  
 original value given at command line.  */
  const char *temp_help; /**< @brief Rescale energy parameters to a temperature of temp C. Default is 37C.
  
 help description.  */
  int noTetra_flag;	/**< @brief Do not include special stabilizing energies for certain tetra-loops. Mostly for testing.
  
 (default=off).  */
  const char *noTetra_help; /**< @brief Do not include special stabilizing energies for certain tetra-loops. Mostly for testing.
  
 help description.  */
  int noLP_flag;	/**< @brief Produce structures without lonely pairs (helices of length 1).
 (default=off).  */
  const char *noLP_help; /**< @brief Produce structures without lonely pairs (helices of length 1).
 help description.  */
  int noGU_flag;	/**< @brief Do not allow GU pairs
  
 (default=off).  */
  const char *noGU_help; /**< @brief Do not allow GU pairs
  
 help description.  */
  int noClosingGU_flag;	/**< @brief Do not allow GU pairs at the end of helices
  
 (default=off).  */
  const char *noClosingGU_help; /**< @brief Do not allow GU pairs at the end of helices
  
 help description.  */
  int noconv_flag;	/**< @brief Do not automatically substitude nucleotide \"T\" with \"U\"
  
 (default=off).  */
  const char *noconv_help; /**< @brief Do not automatically substitude nucleotide \"T\" with \"U\"
  
 help description.  */
  char * nsp_arg;	/**< @brief Allow other pairs in addition to the usual AU,GC,and GU pairs.
 (default='empty').  */
  char * nsp_orig;	/**< @brief Allow other pairs in addition to the usual AU,GC,and GU pairs.
 original value given at command line.  */
  const char *nsp_help; /**< @brief Allow other pairs in addition to the usual AU,GC,and GU pairs.
 help description.  */
  double energyCutoff_arg;	/**< @brief Energy cutoff or pseudoknot initiation cost. Minimum energy gain of a pseudoknot interaction for it to be returned. Pseudoknots with smaller energy gains are rejected.
  
 (default='-8.10').  */
  char * energyCutoff_orig;	/**< @brief Energy cutoff or pseudoknot initiation cost. Minimum energy gain of a pseudoknot interaction for it to be returned. Pseudoknots with smaller energy gains are rejected.
  
 original value given at command line.  */
  const char *energyCutoff_help; /**< @brief Energy cutoff or pseudoknot initiation cost. Minimum energy gain of a pseudoknot interaction for it to be returned. Pseudoknots with smaller energy gains are rejected.
  
 help description.  */
  char * paramFile_arg;	/**< @brief Read energy parameters from paramfile, instead of using the default parameter set.
.  */
  char * paramFile_orig;	/**< @brief Read energy parameters from paramfile, instead of using the default parameter set.
 original value given at command line.  */
  const char *paramFile_help; /**< @brief Read energy parameters from paramfile, instead of using the default parameter set.
 help description.  */
  int verbose_flag;	/**< @brief print verbose output
 (default=off).  */
  const char *verbose_help; /**< @brief print verbose output
 help description.  */
  double subopts_arg;	/**< @brief print suboptimal structures whose energy difference of the pseudoknot to the optimum pseudoknot is smaller than the given value.
 (default='0.0').  */
  char * subopts_orig;	/**< @brief print suboptimal structures whose energy difference of the pseudoknot to the optimum pseudoknot is smaller than the given value.
 original value given at command line.  */
  const char *subopts_help; /**< @brief print suboptimal structures whose energy difference of the pseudoknot to the optimum pseudoknot is smaller than the given value.
 help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int detailed_help_given ;	/**< @brief Whether detailed-help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int cutoff_given ;	/**< @brief Whether cutoff was given.  */
  unsigned int temp_given ;	/**< @brief Whether temp was given.  */
  unsigned int noTetra_given ;	/**< @brief Whether noTetra was given.  */
  unsigned int noLP_given ;	/**< @brief Whether noLP was given.  */
  unsigned int noGU_given ;	/**< @brief Whether noGU was given.  */
  unsigned int noClosingGU_given ;	/**< @brief Whether noClosingGU was given.  */
  unsigned int noconv_given ;	/**< @brief Whether noconv was given.  */
  unsigned int nsp_given ;	/**< @brief Whether nsp was given.  */
  unsigned int energyCutoff_given ;	/**< @brief Whether energyCutoff was given.  */
  unsigned int paramFile_given ;	/**< @brief Whether paramFile was given.  */
  unsigned int verbose_given ;	/**< @brief Whether verbose was given.  */
  unsigned int subopts_given ;	/**< @brief Whether subopts was given.  */

} ;

/** @brief The additional parameters to pass to parser functions */
struct PKplex_cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure PKplex_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure PKplex_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *PKplex_args_info_purpose;
/** @brief the usage string of the program */
extern const char *PKplex_args_info_usage;
/** @brief all the lines making the help output */
extern const char *PKplex_args_info_help[];
/** @brief all the lines making the detailed help output (including hidden options and details) */
extern const char *PKplex_args_info_detailed_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int PKplex_cmdline_parser (int argc, char **argv,
  struct PKplex_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use PKplex_cmdline_parser_ext() instead
 */
int PKplex_cmdline_parser2 (int argc, char **argv,
  struct PKplex_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int PKplex_cmdline_parser_ext (int argc, char **argv,
  struct PKplex_args_info *args_info,
  struct PKplex_cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int PKplex_cmdline_parser_dump(FILE *outfile,
  struct PKplex_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int PKplex_cmdline_parser_file_save(const char *filename,
  struct PKplex_args_info *args_info);

/**
 * Print the help
 */
void PKplex_cmdline_parser_print_help(void);
/**
 * Print the detailed help (including hidden options and details)
 */
void PKplex_cmdline_parser_print_detailed_help(void);
/**
 * Print the version
 */
void PKplex_cmdline_parser_print_version(void);

/**
 * Initializes all the fields a PKplex_cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void PKplex_cmdline_parser_params_init(struct PKplex_cmdline_parser_params *params);

/**
 * Allocates dynamically a PKplex_cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized PKplex_cmdline_parser_params structure
 */
struct PKplex_cmdline_parser_params *PKplex_cmdline_parser_params_create(void);

/**
 * Initializes the passed PKplex_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void PKplex_cmdline_parser_init (struct PKplex_args_info *args_info);
/**
 * Deallocates the string fields of the PKplex_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void PKplex_cmdline_parser_free (struct PKplex_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int PKplex_cmdline_parser_required (struct PKplex_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* RNAPKPLEX_CMDL_H */
