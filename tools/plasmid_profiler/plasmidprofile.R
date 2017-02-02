#' RScript capable
#' Rscript plasmidprofile.R -b blast_runG.tsv -s srst2_runG.tsv -u 0.75 -l 10000 -t "This is a test" -a


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Plasmidprofiler"))
options(bitmapType='cairo')

cl_arguments <- function(){
  # CL arguments ####
  option_list = list(
    make_option(c("-b", "--blastfile"), type="character", default=NULL,
                help="BLAST TSV file name", metavar="character"),
    make_option(c("-s", "--srst2file"), type="character", default=NULL,
                help="SRST2 TSV file name", metavar="character"),
    make_option(c("-u", "--sureness"), type="numeric", default=NA,
                help="Sureness cut off [default = %default]", metavar="numeric"),
    make_option(c("-c", "--coverage"), type="numeric", default=NA,
                help="Percent coverage cut off", metavar="numeric"),
    make_option(c("-l", "--length"), type="numeric", default=NA,
                help="Plasmid length cut off", metavar="numeric"),
    make_option(c("-a", "--anonymize"), action="store_true", default=NA,
                help="Anonymize plasmid and sample names"),
    make_option(c("-o", "--outfile"), type="character", default="P2Run_",
                help="Output filename prefix [default=  %default]", metavar="character"),
    make_option(c("-t", "--title"), type="character", default="Plasmid Profiles",
                help="Title of image [default = %default]", metavar="character"),
    make_option(c("-C", "--combineincs"), action="store_true", default=NA,
                help="Combine very closely related incompatibility groups. eg. ")
    # make_option(c("-T", "--Test"), action="store_true", default=NA,
    #             help="Test filecache")

  );

  opt_parser <- OptionParser(option_list=option_list);
  opt <- parse_args(opt_parser);

  if (is.null(opt$blastfile) | is.null(opt$srst2file)){
    print_help(opt_parser)
    stop("SRST2 and BLAST files must be supplied.", call.=FALSE)
  }
  opt
}


opt <- cl_arguments()

filecache <<- new.env(parent = .GlobalEnv)
assign("name", opt$outfile, envir = filecache)
assign("mods", "Subsampling applied: ", envir = filecache)

main(blast.file = opt$blastfile,
     srst2.file = opt$srst2file,
     coverage.filter = opt$coverage,
     sureness.filter = opt$sureness,
     length.filter = opt$length,
     anonymize = opt$anonymize,
     combine.inc = opt$combineincs,
     main.title = opt$title)
