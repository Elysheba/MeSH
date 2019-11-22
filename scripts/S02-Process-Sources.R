library(here)
library(XML)
library(parallel)
# library(rdflib)
library(dplyr)
library(tidyr)
library(tibble)

setwd(file.path(here(),"scripts/"))
source(here("..", "00-Utils/writeLastUpdate.R"))

##
mc.cores <- 55
sdir <- "../sources"
ddir <- "../data"

###############################################################################@
## Source information ----
###############################################################################@

sfi <- read.table(
   file.path(sdir, "ARCHIVES/ARCHIVES.txt"),
   sep="\t",
   header=T,
   stringsAsFactors=FALSE
)
Mesh_sourceFiles <- sfi[which(sfi$inUse), c("url", "current", "file")]

# devtools::install_github("hrbrmstr/whatamesh")
library(whatamesh)
# list_mesh_files()
ascii <- read_mesh_file(here("sources", Mesh_sourceFiles$file), wide = T)
Mesh_sourceFiles <- select(Mesh_sourceFiles, url, current)

###############################################################
## Parent/child
## Only keep diseases --> MN IDs starting with C (https://www.nlm.nih.gov/mesh/2019/download/2019New_Mesh_Tree_Hierarchy.txt)
hier <- ascii[,grep(paste("UI","MN",sep = "|"),colnames(ascii))]
hier <- do.call(rbind,sapply(grep("UI",invert = T, colnames(hier), value = T),
                          function(df){
                            d <- hier[,df]
                            names(d) <- "meshid"
                            toRet <- data.frame("id" = hier$UI,
                                                "meshid" = d,
                                                stringsAsFactors = F)
                            toRet <- toRet[!is.na(toRet$meshid),]
                            toRet$meshid <- gsub("[:|:].*","",toRet$meshid)
                            return(toRet)
                          },
                          simplify = FALSE,
                          USE.NAMES = FALSE))
hier <- hier[grep("^C",hier$meshid),]
hier$parent <- ifelse(!grepl("[.]",hier$meshid),
                  NA,
                  substr(hier$meshid,1,nchar(hier$meshid) - 4))
hier$parentid <- hier$id[match(hier$parent,hier$meshid)]
##
parentId <- hier[!is.na(hier$parentid),c("id","parentid")]
names(parentId) <- c("id","parent")
parentId$DB <- parentId$pDB <- "MeSH"
parentId$origin <- "MeSH"

## Add levels
getAncestors <- function(id){
  direct <- termParents[[id]]
  parents <- direct
  level <- 0
  dLev <- c()
  for(d in direct){
    dPar <- getAncestors(d)
    dLev <- c(dLev, dPar$level)
    parents <- c(parents, dPar$parents)
  }
  if(length(dLev)>0){
    level <- max(dLev)+1
  }
  return(list(parents=unique(parents), level=level))
}


parentList <- unstack(parentId, parent~id)
termParents <- parentList
library(BiocParallel)
bpparam <- MulticoreParam(workers = 30)

termAncestors <- bplapply(
  parentId$id,
  getAncestors,
  BPPARAM = bpparam
)
names(termAncestors) <- parentId$id

##################################################
## syn/labels
syn <- ascii[ascii$UI %in% hier$id,grep(paste("UI","ENTRY",sep = "|"),colnames(ascii))]
syn <- do.call(rbind,sapply(grep("UI",invert = T, colnames(syn), value = T),
            function(df){
              d <- syn[,df]
              names(d) <- "syn"
              toRet <- data.frame("id" = syn$UI,
                                   "syn" = d,
                                   stringsAsFactors = F)
              toRet <- toRet[!is.na(toRet$syn),]
              toRet$syn <- gsub("[:|:].*","",toRet$syn)
              return(toRet)
            },
            simplify = FALSE,
            USE.NAMES = FALSE)) 
syn$canonical <- FALSE
##
label <- ascii[ascii$UI %in% hier$id,grep(paste("UI","\\bMH\\b",sep = "|"),colnames(ascii))]
names(label) <- c("syn","id")
label$canonical <- TRUE
##
idNames <- syn %>%
  as_tibble() %>%
  bind_rows(label) %>%
  mutate(DB = "MeSH") %>%
  arrange(canonical) %>%
  distinct() 
# idNames <- idNames[order(idNames$canonical,decreasing = T),]
# idNames <- unique(idNames)
dim(idNames)

## Check characters for \t, \n, \r and put to ASCII
idNames$syn <- iconv(x = idNames$syn,to="ASCII//TRANSLIT")
table(unlist(sapply(idNames$syn, strsplit, split = "")))
idNames$syn <- gsub(paste("\n","\t","\r", sep = "|")," ",idNames$syn)
table(unlist(sapply(idNames$syn, strsplit, split = "")))
idNames$syn <- gsub("\"","'",idNames$syn)
table(unlist(sapply(idNames$syn, strsplit, split = "")))

##################################################################
## entryId
def <- ascii[ascii$UI %in% hier$id,grep(paste("UI","MS",sep = "|"),colnames(ascii))]
entryId <- data.frame(DB = "MeSH",
                      id = unique(hier$id),
                      def = def$MS[match(unique(hier$id),def$UI)],
                      stringsAsFactors = F)

## Empty definition to NA
nc <- nchar(entryId$def)
head(table(nc), n = 20)
# entryId[which(nc < 4),]
# entryId[which(nc < 4),"def"] <- NA
## Check characters for \t, \n, \r and put to ASCII
entryId$def <- iconv(x = entryId$def,to="ASCII//TRANSLIT")
table(unlist(sapply(entryId$def, strsplit, split = "")))
entryId$def <- gsub(paste("\n","\t","\r", sep = "|")," ",entryId$def)
table(unlist(sapply(entryId$def, strsplit, split = "")))
entryId$def <- gsub("\"","'",entryId$def)
entryId$def <- gsub("\\\\","",entryId$def)
table(unlist(sapply(entryId$def, strsplit, split = "")))

## add level
entryId <- entryId %>%
  mutate(
    level=unlist(lapply(termAncestors, function(x) x$level))[entryId$id]
  ) %>%
  mutate(level = case_when(is.na(level) ~ 0,
                           TRUE ~ level))



############################
Mesh_parentId <- parentId[,c("DB","id","pDB","parent","origin")]
Mesh_entryId <- entryId[,c("DB","id","def","level")]
Mesh_idNames <- idNames[,c("DB","id","syn","canonical")]

###############################################################################@
## Writing tables ----
###############################################################################@
message("Writing tables...")
message(Sys.time())
toSave <- grep("^Mesh[_]", ls(), value=T)
for(f in toSave){
  message(paste("Saving", f))
  ## Ensure unicity
  assign(f, get(f))
  if(length(names(f))==0){
    f <- unique(f)
  }
  ##
  write.table(
    get(f),
    file=file.path(ddir, paste(f, ".txt", sep="")),
    sep="\t",
    row.names=FALSE, col.names=TRUE,
    quote=TRUE,
    qmethod = "double"
  )
}
message(Sys.time())
message("... Done\n")
writeLastUpdate()

##############################################################
## Check model
source("../../00-Utils/autoCheckModel.R")
