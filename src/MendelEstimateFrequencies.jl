"""
This module organizes estimation of allele frequencies from pedigree data.
"""
module MendelEstimateFrequencies
#
# Other OpenMendel modules.
#
using MendelBase
# using DataStructures
# using ModelConstruction
# using ElstonStewartPreparation
# using ElstonStewartEvaluation
# using Optimize
#
# External modules.
#
using DataFrames    # From package DataFrames.
using Distributions # From package Distributions.

export EstimateFrequencies

"""
This is the wrapper function for the AIM Selection analysis option.
"""
function EstimateFrequencies(control_file = ""; args...)

  const ESTIMATE_FREQUENCIES_VERSION :: VersionNumber = v"0.1.0"
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println("  Estimate Frequencies analysis option")
  println("        version ", ESTIMATE_FREQUENCIES_VERSION)
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{ASCIIString, Any}())
  #
  # Keywords unique to this analysis may be defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = value
  #
  # Process the run-time user-specified keywords that will control the analysis.
  # This will also initialize the random number generator.
  #
  process_keywords!(keyword, control_file, args)
  #
  # Check that the correct analysis option was specified.
  #
  lc_analysis_option = lowercase(keyword["analysis_option"])
  if (lc_analysis_option != "" &&
      lc_analysis_option != "estimatefrequencies")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
  end
  keyword["analysis_option"] = "EstimateFrequencies"
  #
  # Read the genetic data from the external files named in the keywords.
  #
  (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Execute the specifed analysis.
  #
  println(" \nAnalyzing the data.\n")
  execution_error = false
  skipped_loci = estimate_frequencies_option(pedigree, person, nuclear_family,
      locus, locus_frame, phenotype_frame, pedigree_frame, keyword)
  if execution_error
    println(" \n \nERROR: Mendel terminated prematurely!\n")
  else
    println(" \n \nMendel's analysis is finished.\n")
  end
  #
  # Finish up by closing and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing

end # function EstimateFrequencies

function estimate_frequencies_option(pedigree::Pedigree, person::Person,
  nuclear_family::NuclearFamily, locus::Locus, locus_frame::DataFrame, 
  phenotype_frame::DataFrame, pedigree_frame::DataFrame,
  keyword::Dict{ASCIIString, Any})
  io = keyword["output_unit"]
  skipped_loci = 0
  #
  # Eliminate genotypes but do not lump alleles at each locus.
  #
  keyword["eliminate_genotypes"] = true
  #
  # Prepare to estimate allele frequencies at each locus.
  #
  (model_loci, locus.model_loci) = (locus.model_loci, 1)
  model_locus = similar(locus.model_locus)
  copy!(model_locus, locus.model_locus)
  locus.model_locus = zeros(Int, 1)
  #
  # Loop over all loci.
  #
  for loc = 1:locus.loci
    locus.model_locus[1] = loc
    #
    # Construct the parameter data structure.
    #
    keyword["constraints"] = 1
    keyword["parameters"] = locus.alleles[loc]
    keyword["title"] = "Estimate Frequencies analyis for " * locus.name[loc]
    parameter = set_parameter_defaults(keyword)
    parameter = initialize_optimization(locus, parameter, keyword)
    #
    # Fetch the instructions for conducting the Elston-Stewart algorithm.
    #
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree, 
      person, nuclear_family, locus, keyword)
    if any(elston_stewart_count .>  keyword["complexity_threshold"])
      println("One or more pedigrees exceeds the complexity threshold.")
      println("$locus.name[loc] is being skipped.")
      skipped_loci = skipped_loci + 1
      continue
    end
    #
    # Pass the variables to optimize for maximum likelihood estimation.
    #
    function fun(par)
      copy!(parameter.par, par)
      f = elston_stewart_loglikelihood(pedigree, person, locus, parameter, 
        instruction, keyword)
      return (f, nothing, nothing)
    end # function fun
    (best_par, best_value) = optimize(fun, parameter)
  end
  return skipped_loci
end # function estimate_frequencies_option 

end # module MendelEstimateFrequencies
