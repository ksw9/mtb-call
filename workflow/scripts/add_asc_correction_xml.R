## Add correction for invariant sites. 
add_asc_correction_xml <- function(xml_path, fasta_path){

id <- str_replace(basename(fasta_path), pattern = '.fa', replacement = '')
id  

# Read in original XML file
xml_string <- read_file(xml_path)

# Replace text for alignment ID first.
original <- id
replacement <- paste(id,'Original', sep = '')
xml_update1 <- sub(pattern = original, replacement = replacement, x = xml_string)

# Count number of variant nucleotides in first sequence (approximates msa).
seq1 <-  seqinr::read.fasta(fasta_path)
As_var = ifelse(is.na(table(seq1[1])['a']), 0, table(seq1[1])['a'])
Cs_var = ifelse(is.na(table(seq1[1])['c']), 0, table(seq1[1])['c'])
Gs_var = ifelse(is.na(table(seq1[1])['g']), 0, table(seq1[1])['g'])
Ts_var = ifelse(is.na(table(seq1[1])['t']), 0, table(seq1[1])['t'])

# Reference nucleotide counts
As_ref = 758552
Cs_ref = 1449998
Gs_ref = 1444614
Ts_ref = 758368

# Get invariant sites.
As = As_ref - As_var
Cs = Cs_ref - Cs_var
Gs = Gs_ref - Gs_var
Ts = Ts_ref - Ts_var

# Add correction text.
correction_text <- paste("</data>\n<data id='",id,"' spec='FilteredAlignment' filter='-' data='@", id, "Original' constantSiteWeights='",As," ",Gs," ",Cs," ",Ts,"'/>", sep = '')
correction_text

# Add correction text to file.
xml_update2 <- sub(pattern = "</data>\n", replacement = correction_text, x = xml_update1)

# Write to updated XML file.
write_lines(xml_update2, file = str_replace(xml_path, pattern = '\\.xml','_snpcor.xml'))

}
