###############################################################################
# This function corrects the lengh of the transcript 
# features based on the variant selection, 
# and the padding length. For example, 
# if we select only a RS id, it takes the output from 
# annotateGeneLength and it corrects it to 1 (if padding length =0)
# INPUT:
#   the output from the function .annotateGeneLength
# OUTPUT: 
#   a modified version of the output from .annotateGeneLength
.annotateVariationLength <- function(panel 
                                     , gene_length , utr , padding_length )
{
    full_gene_specs <- c("amplification" , "deletion" , "up" , "down" , "" 
                        , "mutation_type" , "exact_fusion" , "gene_fusion")
    # every amino acid notation is composed by three bases
    aLength <- 3
    # Every time we subtract end-start in a "bed style" fashion we need to add 1
    offSet <- 1
    # List all possible target field
    allTarget <- c(full_gene_specs, "amino_acid_position", "amino_acid_variant" 
                   , "dbSNP_rs", "genomic_variant", "genomic_position")
    #Cycle through the panel and get the variation lenght based on the choice
    variation_length <- vapply(seq_len(nrow(panel)) , function(x) {
            # browse alteration
            target <- panel[x , "exact_alteration"]
            if( ! target %in% allTarget){
              stop(paste("Mutation Specification is not correct at panel line" 
                         , x))
            }
            if(target=="amino_acid_position"){
                out <- strsplit(as.character(
                  panel[x , "mutation_specification"]) , "-") %>% 
                  unlist
                # Final amino acid - First amino acid + 1. all multiply by 3
                out <- (as.numeric(out[2])-as.numeric(out[1])+offSet)*aLength
            }
            if(target=="amino_acid_variant"){
                out <- aLength + 2*padding_length
                if(out<(padding_length*2 + offSet))
                    out <- (padding_length*2 + offSet)
            }
            if(target %in% c("dbSNP_rs" , "genomic_variant")){
                out <- offSet + 2*padding_length
                if(out<(padding_length*2 + offSet))
                    out <- (padding_length*2 + offSet)
            }
            if(target=="genomic_position"){
                out <- strsplit(as.character(
                  panel[x , "mutation_specification"]) , ":|-") %>% 
                unlist
                out <- as.numeric(out[3]) - as.numeric(out[2]) + offSet
                if(out<(padding_length*2 + offSet))
                    out <- (padding_length*2 + offSet)
            }
            if(target %in% full_gene_specs){
              if(utr){
                 chosenLength <- "cds_and_utr_len"
              } else {
                chosenLength <- "cds_len"
              }
              if(panel[x , "alteration"]!="fusion"){
                out <- gene_length[ 
                    gene_length$gene_symbol==panel[x , "gene_symbol"] 
                    , chosenLength]
              } else { 
                #it is a fusion, lets remove the __ simbol
                if(grepl("__" , panel[x , "gene_symbol"])){
                  out <- gene_length[ gene_length$gene_symbol %in% 
                      unlist(strsplit(panel[x , "gene_symbol"] , "__")) 
                                      , chosenLength] %>%
                  as.integer %>% mean %>% round
                } else {
                  out <- gene_length[ 
                      gene_length$gene_symbol==panel[x , "gene_symbol"] 
                      , chosenLength] %>%
                    as.integer
                }
              }
            }
            return(out)
        } , numeric(1))
    panel$variation_len <- as.integer(variation_length)
    return(panel)
}