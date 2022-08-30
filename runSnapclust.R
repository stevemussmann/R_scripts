#!/usr/bin/Rscript

library("optparse")
library("adegenet")

option_list = list(
	make_option(
		c("-f", "--file"), 
		type="character", 
		default=NULL, 
		help="File in single-line-per-individual Structure format.", 
		metavar="character"
	),
	make_option(
		c("-o", "--out"), 
		type="character", 
		default="Output", 
		help="File prefix for outputs. Default = 'Output'", 
		metavar="character"
	),
	make_option(
		c("-b", "--bx"),
		action="store_true",
		default=FALSE,
		help="Turn on testing for backcrosses.",
	)
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
	print_help(opt_parser)
	stop("Input file must be provided.", call.=FALSE)
}

#make output file names
tableout <- paste(opt$out, "csv", sep=".")
pdfout <- paste(opt$out, "pdf", sep=".")

# data file input requirements:
# no header line of loci
# one line per individual
data = read.table(opt$file)
ninds=nrow(data) #number of rows = number of individuals
nloci=(ncol(data)-1)/2 #number of columns - line labels / 2 = number of loci

genodata <- read.structure(
	opt$file,
	n.ind=ninds,
	n.loc=nloci,
	onerowperind=TRUE,
	col.lab=1,
	#col.lab=2,
	NA.char="-9",
	ask=FALSE
)

if(opt$bx == TRUE){
	res.hyb<-snapclust(genodata, k=2, hybrids=TRUE, hybrid.coef = c(0.25, 0.5))
}else{
	res.hyb<-snapclust(genodata, k=2, hybrids=TRUE)
}

#write classifications output to table
write.table(res.hyb, 
	tableout, 
	append=FALSE, 
	sep=",", 
	row.names = TRUE, 
	col.names = TRUE
)

#write compoplot to file
pdf(pdfout)
compoplot(res.hyb, col.pal = hybridpal(), n.col = 2)
dev.off()

quit()
