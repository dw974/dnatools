#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

extract_region <- function(ref=NULL,qry_fld=NULL,temp_dir=NULL,len_thresh=NULL) {
  if (is.null(ref)) stop("Please provide a fasta file containing reference sequence data")
  if (is.null(temp_dir)) temp_dir=tempdir()

  ref=normalizePath(ref)
  temp_dir=normalizePath(temp_dir)

  print("This function extracts DNA regions that best match the provided reference sequence from all sequences in a specified folder, using BLASTN.")
  print("The folder should contain files with the extensions .fa .fas .fna .fasta")
  print(paste0("Temporary files will be written to this location: ",temp_dir))

  lst=list.files(path=qry_fld,pattern=".fa$|.fas$|.fasta|.fna",full.names = T)
  nm=length(lst)
  print(paste0("The data folder contains ", nm," files to be queried."))

  print("Reading reference fasta file.")
  tryCatch(expr={fas=Biostrings::readDNAStringSet(ref)},error=function(e){stop("ERROR: Please provide a valid fasta-format reference file.")})
  print(paste0("The reference file contains ",length(fas)," sequence(s)."))

  print("Generating BLAST database")
  pt=normalizePath(paste0(temp_dir,"/sequences"),winslash = "\\",mustWork = F)
  print(pt)
  pt2=normalizePath(paste0(temp_dir,"/results"),winslash = "\\",mustWork = F)
  print(pt2)
    system(paste0("makeblastdb -dbtype nucl -input_type fasta -in ",ref," -out ",pt))
    print("Done.")

  print("BLASTing all query sequences, and extracting the best hits")

  res=list()
  for (x in 1:nm){
    print(paste0("BLASTing sequence ",x))
      print(paste0("blastn -num_threads 4 -query ",normalizePath(lst[x])," -db ",pt," -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -out ",pt2))

      system(paste0("blastn -num_threads 4 -query ",normalizePath(lst[x])," -db ",pt," -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq\" -out ",pt2))

#
    if(file.info(pt2)$size>0){
      tab=read.table(pt2,header=F,stringsAsFactors = F)
      colnames(tab)=c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq")
      res[[x]]= tab %>% dplyr::group_by(qaccver) %>% dplyr::summarise(best_hit=saccver[order(evalue,-rank(bitscore))[1]],qseq=qseq[order(evalue,-rank(bitscore))[1]],score=evalue[order(evalue,-rank(bitscore))[1]],length=length[order(evalue,-rank(bitscore))[1]],file=basename(lst[x]))
    }else{
      print(paste0("Sequence ",x," (",basename(lst[x]),") has no BLAST hits"))
    }
  }
  res=do.call(rbind,res)
  if(!is.null(len_thresh)){
    res=res[res$length>len_thresh,]
  }
  return(res)
}


align_df=function(df=NULL,temp_dir=NULL){
  if (is.null(temp_dir)) temp_dir=tempdir()

  print("This function aligns the sequences in the output from function extract_region using system calls to mafft.")
  print(paste0("Temporary files will be written to this location: ",temp_dir))
  print("Writing sequence data to disk")
  seqinr::write.fasta(sequences = as.list(gsub("-| ","",df$qseq)),names = as.list(df$file),file.out = paste0(temp_dir,"/aln.fasta"),open="w")
  print("Calling mafft")
  setwd(temp_dir)
  print(paste0("mafft --adjustdirection --thread 4 aln.fasta > aln2.fasta"))
  system(paste0("mafft --adjustdirection --thread 4 aln.fasta > aln2.fasta"),intern=T)
  dna=Biostrings::readDNAStringSet("aln2.fasta")
  df=data.frame(file=names(dna),id=df$qaccver,aln=paste(dna))
  un=unique(df$aln)
  df$allele=match(df$aln,un)
  return(df)
}

get_source_data=function(ids=NULL){

  md=lapply(1:length(ids),function(y){
    print(paste0(y," of ",length(ids)))
    URL=paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=",ids[y],"&rettype=gb&retmode=text")
    data <- tryCatch(scan(file = URL, what = "", sep = "\n", quiet = TRUE),error=function(cond){NA},warning=function(cond){NA})
    if(!is.na(data)){
    pos1=which(unlist(sapply(1:length(data),function(x) length(grep("     source",data[x]))==1)))+1
    pos2=which(unlist(sapply(1:length(data),function(x) length(grep("                     ",data[x]))==0)) & sapply(1:length(data),function(x) x>pos1))[1]-1
    tmp=lapply(pos1:pos2,function(x){
      tmp=str_split(str_match(data[x], "^                     /(.*?)\"$")[[2]],"=\"",simplify = T)
      dftmp=data.frame(tmp=tmp[2])
      colnames(dftmp)=tmp[1]
      return(dftmp)

    })
    return(cbind(data.frame(id=ids[y],stringsAsFactors = F),bind_cols(tmp)))
  }else{
    return(data.frame(id=ids[y],stringsAsFactors = F))
  }
  })
  md=bind_rows(md)
  return(md[,!grepl("V[0-9]",colnames(md))])

}

write_alleles=function(df=NULL,file.out=NULL){
  un=unique(df$aln)
  cor=data.frame(id=df$id,allele=paste0("allele_",match(df$aln,un)),stringsAsFactors = F)
  write.csv(cor,paste0(dirname(file.out),"allele_correspondence.csv"))
  write.fasta(as.list(un),lapply(1:length(un),function(x) paste0("allele_",x)),file.out)
}

make_contigs = function(tab=NULL,outfld=NULL){

  tab=read.table(tab,header=T)
  res=sapply(1:dim(tab)[1],function(x){
      print(paste0("Attempting to generate contig for sequence ",tab$name[x]," using method M2"))
      sancon=sangeranalyseR::SangerContig(parentDirectory = tab$folder[x],
                                          suffixForwardRegExp = tab$forward[x],
                                          suffixReverseRegExp = tab$reverse[x],
                                          TrimmingMethod = "M2",
                                          processorsNum=2,
                                          M1TrimmingCutoff = NULL,
                                          M2CutoffQualityScore = 40,
                                          M2SlidingWindowSize = 30,
                                          contigName = tab$name[x])
    if(length(sancon@contigSeq)==0){
      print("Method M2 was NOT SUCCESSFUL")
      return(0)
    }else{
      print("Writing contig to file")
      sangeranalyseR::writeFastaSC(sancon,outputDir = tab$folder[x],selection = "contig")
      return(1)
    }
  }
  )
  system(paste(c("cat",paste0(tab$folder[which(res==1)],tab$name[which(res==1)],"_contig.fa"), ">",paste0(outfld,"sequences.fasta")),collapse=" "))

}

align_seqs=function(list=NULL,outfld=NULL){
  tmp=tempfile(tmpdir = outfld)
  tmp2=tempfile(tmpdir = outfld)
  system(paste0("cat ",paste(list,collapse=" ")," > ",tmp))
  system(paste0("mafft --auto --thread 2 --adjustdirectionaccurately ",tmp," > ",tmp2))
  system(paste0("cp ",tmp2," ",outfld,"/aligned_sequences.fasta"))
  system(paste0("rm ",tmp," ",tmp2))
}

get_dists=function(aln=NULL,fileout=NULL){
  tmp=tempfile()
  system(paste0("snp-dists -b ",aln," > ",tmp))
  tab=read.table(tmp,header=T)
  d=as.dist(tab,upper=FALSE)
  tab=melt(as.matrix(d), varnames = c("From", "To"),value.name = "SNPs")
  write.table(tab,fileout,row.names = F)
  return(tab)
}
