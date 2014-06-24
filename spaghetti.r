#!/usr/bin/Rscript
options(warn=-1)
Sys.setenv(LANG="EN")
suppressMessages(library(data.table))
suppressMessages(library(gdata))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(optparse))



option_list <- list(
    make_option(c("-b", "--PdfWidth"),type="double", default=12,
                help="Width of .pdf output file, [default %default]"),
    make_option(c("-c", "--PdfHeight"), type="double",default=7,
                help="Height of .pdf output file,, [default %default]"),
    make_option(c("-d", "--EpsWidth"), type="double", default=6,
                help="Width of .eps output file, [default %default]"),
    make_option(c("-e", "--EpsHeight"), type="double", default=4,
                help="Height of .eps output file, [default %default]"),
    make_option(c("-a", "--Angles"), action="store_true", default=FALSE,
                help="Use angles instead of bonds")			
    )

parser    <- OptionParser(usage = "spaghetti.r [options] directory", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt       <- arguments$options

if (length(arguments$args) != 1)
{
    print("Incorrect number of required positional arguments")
    print_help(parser)
    stop()
}





directory <-   arguments$args
if (opt$Angles) {pattern <- "Angle_Restrain(_Breakable)?"
                 name.su <- "_angles"
             }else{
                 pattern = "Distance_Restrain_(Breakable|Morse)"
                 name.su <- "_bonds"}

na.to.0 <- function(x) {if (is.na(x)) 0 else x}

extract.data <- function (filename) ###Индия!
    {
        xs <- readLines(filename)
#        xs <- gsub("'.*","",xs) #kill comments
        k1 <- unlist(strsplit(grep("penalties_weighting_K1",xs, value=TRUE), "[[:space:]]+"))
        k1 <- as.numeric(k1[length(k1)])
        rwp <-  unlist(strsplit(grep("r_wp", xs, value=TRUE)[1],  "[[:space:]]+"))
        rwp <- as.numeric(rwp[which(rwp=="r_wp") +1])
        ys <- grep(pattern,xs, value=TRUE)
        ys <- ys[!grepl("H\\d",ys, perl=TRUE)]
        ys <- strsplit(ys, "([ \t]*,[][ \t]*)|[()][ \t]*")
        l.ys  <- length(ys)
        Delta <- numeric(l.ys)
        Bond  <- numeric(l.ys)
        Error <- numeric(l.ys)
        for (i in 1:l.ys){
            y <- ys[[i]]
            num <- unlist(strsplit(y[4], "(`_|`)"))
            Delta[i] <- as.numeric(num[1]) - as.numeric(y[3])
            Error[i] <- na.to.0(as.numeric(num[2]))
            Bond[i] <- y[2]
        }
        qq <- quantile(Delta, c(0.25,0.75), names=FALSE)
        Low.Q=qq[1]-1.5*diff(qq)
        High.Q=qq[2]+1.5*diff(qq)
        data.frame(Bond=Bond, Delta=Delta,
                   Error=Error, Av.Error=mean(Error),
                   Low.Q=Low.Q, High.Q=High.Q, 
                   K1=k1, Rwp=rwp, Outlier = (Delta > (High.Q + 2 * Error)) |
                                             (Delta < (Low.Q - 2 * Error)),
                   Bad.Bond=FALSE)            
    }

cat("Reading data..\n")

topas.outs <- list.files(directory, full.names=TRUE, pattern="*.(out|OUT)")
errors.table <- rbindlist(Map(extract.data, topas.outs))

cat("Some calculations...\n")


errors.table <- rbindlist(by(errors.table, errors.table$Bond, function(x) {x[["Bad.Bond"]] <-  any(x[["Outlier"]]); return(x)}))

errors.summary <- rbindlist(by(errors.table, errors.table$K1, 
	function(x){ data.frame(K1=x[["K1"]][1], Rwp=x[["Rwp"]][1],
	RMS=round(sqrt(sum((x[["Delta"]])^2/length(x[["Delta"]]))), digits=4),
                                Outliers=sum(x[["Outlier"]]))}))



## Plotting settings

PlotBreaks <- function(x){
    breaks <- seq(0, 100, by=10)
    names(breaks) <- attr(breaks,"labels")
    breaks
}
 
p <- ggplot(data=errors.table, aes(x=K1, y=Delta)) + theme_bw()   # +  scale_x_continuous (breaks=PlotBreaks)
if (sum(errors.table$Av.Error) >  0) {
	p <- p +
		geom_ribbon(aes(ymin=High.Q-2*Av.Error, ymax=High.Q+2*Av.Error, y=NULL), fill="#999999") +
		geom_ribbon(aes(ymin=Low.Q-2*Av.Error, ymax=Low.Q+2*Av.Error, y=NULL), fill="#999999") 
	}else {
	cat("WARNING: You'd better use do_errors keyword in your *.inp files! \n")
	p <- p +  geom_line(aes(y = High.Q), linetype="dashed") + geom_line(aes(y = Low.Q), linetype="dashed") }
	
 if (sum(errors.table$Outlier) > 0){
     p <- p +
         geom_line(data=errors.table[errors.table$Bad.Bond==FALSE,], aes(group=Bond), color="black") +
             geom_line(data=errors.table[errors.table$Bad.Bond==TRUE,], aes(group=Bond, color = Bond)) +
                 geom_point(data=errors.table[errors.table$Outlier==TRUE,], aes(group=Bond, color = Bond), size=2) +
				 scale_colour_discrete(name="Outlying bonds")
 } else {
	p <- p + geom_line(aes(group=Bond) , colour="black")  }

## Output

cat("Plotting and tables writing...\n")

write.table(errors.table[order(errors.table$K1,errors.table$Bond),], paste(directory, name.su, "_rtable.dat", sep=""), row.names=FALSE, sep="\t")	
write.fwf(errors.summary[order(errors.summary$K1)], paste(directory, name.su,"_summary.txt", sep=""))

ggsave(p, file= paste(directory, name.su, "_plot.pdf", sep = ""), width=opt$PdfWidth, height=opt$PdfHeight)
ggsave(p, file= paste(directory, name.su, "_plot.eps", sep = ""), width=opt$EpsWidth, height=opt$EpsHeight)
