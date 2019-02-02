__precompile__()

"""
This module organizes estimation of allele frequencies from pedigree data.
"""
module MendelEstimateFrequencies
#
# Required OpenMendel packages and modules.
#
using MendelBase
# namely: DataStructures, ModelConstructions,
# ElstonStewartPreparations, ElstonStewartEvaluations
using MendelSearch
#
# Required external modules.
#
using DataFrames
using Distributions

export EstimateFrequencies
"""
This is the wrapper function for the Estimate Allele Frequency analysis option.
"""
function EstimateFrequencies(control_file = ""; args...)

  ESTIMATE_FREQUENCIES_VERSION :: VersionNumber = v"0.5.0"
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
  keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
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
  # Execute the specified analysis.
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
  # Finish up by closing, and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing

end # function EstimateFrequencies

function estimate_frequencies_option(pedigree::Pedigree, person::Person,
  nuclear_family::NuclearFamily, locus::Locus, locus_frame::DataFrame,
  phenotype_frame::DataFrame, pedigree_frame::DataFrame,
  keyword::Dict{AbstractString, Any})

  io = keyword["output_unit"]
  skipped_loci = 0
  #
  # Eliminate genotypes but do not lump alleles at each locus.
  #
  keyword["eliminate_genotypes"] = true
  #
  # Add a frequency field to the locus file frame
  # and an index n for the current position in the locus frame.
  #
  locus_frame[:PedFrequency] = zeros(size(locus_frame, 1))
  n = 0
  #
  # Prepare to estimate allele frequencies at each locus.
  #
  (model_loci, locus.model_loci) = (locus.model_loci, 1)
  model_locus = similar(locus.model_locus)
  copyto!(model_locus, locus.model_locus)
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
    keyword["goal"] = "maximize"
    keyword["parameters"] = locus.alleles[loc]
    keyword["title"] = "Estimate Frequencies analyis for " * locus.name[loc]
    parameter = set_parameter_defaults(keyword)
    parameter =
      initialize_optimization_estimate_frequencies!(locus, parameter, keyword)
    #
    # Fetch the instructions for conducting the Elston-Stewart algorithm.
    #
    (instruction, elston_stewart_count) = orchestrate_likelihood(pedigree,
      person, nuclear_family, locus, keyword)
    if any(elston_stewart_count .>  keyword["complexity_threshold"])
      println("One or more pedigrees exceeds the complexity threshold.")
      println("$locus.name[loc] is being skipped.")
      skipped_loci = skipped_loci + 1
      n = n + locus.alleles[loc]
      continue
    end
    #
    # Pass the variables to search for maximum likelihood estimation.
    #
    function fun(par)
      copyto!(parameter.par, par)
      f = elston_stewart_loglikelihood(penetrance_estimate_frequencies,
        prior_estimate_frequencies, transmission_estimate_frequencies,
        pedigree, person, locus, parameter, instruction, keyword)
      return (f, nothing, nothing)
    end # function fun
    (best_par, best_value) = mendel_search(fun, parameter)
    for i = 1:locus.alleles[loc]
      n = n + 1
      locus_frame[n, :PedFrequency] = best_par[i]
    end
  end
  show(locus_frame)
  return skipped_loci
end # function estimate_frequencies_option

"""
Supply a penetrance for individual i.
"""
function penetrance_estimate_frequencies(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int)

  pen = 1.0
  return pen
end # function penetrance_estimate_frequencies

"""
Supply a prior probability for founder i.
"""
function prior_estimate_frequencies(person::Person, locus::Locus,
  multi_genotype::Matrix{Int}, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int)

  prior_prob = 1.0
  for l = start:finish
    loc = locus.model_locus[l]
    allele = multi_genotype[1, l]
    frequency = par[allele]
    prior_prob = prior_prob * frequency
    if !locus.xlinked[loc] || !person.male[i]
      allele = multi_genotype[2, l]
      frequency = par[allele]
      prior_prob = prior_prob * frequency
    end
  end
  return prior_prob
end # function prior_estimate_frequencies

"""
Supply the transmission probability that a parent i with a particular
genotype transmits a particular gamete to his or her child j.
"""
function transmission_estimate_frequencies(person::Person, locus::Locus,
  gamete::Vector{Int}, multi_genotype::Matrix{Int}, par::Vector{Float64},
  keyword::Dict{AbstractString, Any}, start::Int, finish::Int, i::Int, j::Int)
  #
  # For male to male inheritance at an x-linked locus,
  # set the transmission probability equal to 1.
  #
  loc = locus.model_locus[start]
  xlinked = locus.xlinked[loc]
  if xlinked && person.male[i] && person.male[j]
    return 1.0
  end
  #
  # Apply Mendel's segregation rules.
  #
  k = multi_genotype[1, start]
  l = multi_genotype[2, start]
  m = gamete[start]
  trans = 0.0
  if k == m; trans = trans + 0.5; end
  if l == m; trans = trans + 0.5; end
  return trans
end # function transmission_estimate_frequencies

"""
Initialize the optimization problem.
"""
function initialize_optimization_estimate_frequencies!(locus::Locus,
  parameter::Parameter, keyword::Dict{AbstractString, Any})
  #
  # Initialize, bound, and name the parameters.
  #
  loc = locus.model_locus[1]
  for i = 1:parameter.parameters
    parameter.par[i] = locus.frequency[loc][1, i]
    parameter.min[i] = 1e-5
    parameter.name[i] = "freq"*" $i"
    parameter.name[i] = rpad(parameter.name[i], 8, ' ')
  end
  #
  # Force the parameter frequencies to sum to 1.0.
  #
  parameter.constraint[1, :] .= 1.0
  parameter.constraint_level[1] = 1.0
  return parameter
end # function initialize_optimization_estimate_frequencies!
#
# Method to obtain path to this package's data files
# so they can be used in the documentation and testing routines.
# For example, datadir("Control file.txt") will return
# "/path/to/package/data/Control file.txt"
#
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module MendelEstimateFrequencies
