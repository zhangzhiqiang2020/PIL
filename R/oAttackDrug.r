#' Function to perform drug targets combinatory attack for an igraph object
#'
#' \code{oAttackDrug} is supposed to perform drug-based combinatory attack for an igraph object, where the effect of drug targets removal is defined as the fraction of network nodes disconnected from the giant component (the largest connected component) remained. There are two types of attack analysis: removing a single drug targets (accordingly we define the attackness per a drug), and removing drug targets in combinations (combinatorial attack, often used based on the concept of combinatorial optimisation). Also attackness for a drug is defined as the effect estimated when a single drug targets node removed.
#'
#' @param drug.nodes  a list containing drugs and their drugs target nodes combined for combinatorial attack.
#' @param drugs.fixed a vector of drug names that fixed with the other drugs for combinatorial attack.
#' @param combine.num  the drug combinatorial strategies: 2, two drugs are combined; 3, three drugs are combined.
#' @param ig an object of class "igraph" with node attribute 'name'.
#'
#' @return
#' a tibble with 4 columns, including 'drug', 'frac.disconnected' (fraction of network nodes disconnected from the giant component), 'frac.removed' (fraction of network nodes removed), 'nodes.removed' (nodes removed, provided separated by ',').
#' @note none
#' @export
#' @seealso \code{\link{oAttackDrug}}
#' @include oAttackDrug.r
#' @examples
#' \dontrun{
#' drug.nodes <- pbapply::pblapply(drug_names,function(x){
#' DTT %>% transmute(drug=drug,Symbol=Symbol,phase=phase) %>%
#' group_by(drug,phase) %>% 
#' summarise(Targets=str_c(Symbol,collapse = ","),.groups = 'drop') %>% 
#' filter(drug==x,phase!=1,phase!=2) %>% 
#' pull(Targets) %>% str_split(",") %>% unlist() %>% unique() %>% 
#' intersect((subg %>% oIG2TB('nodes') %>% pull(name)))
#' })
#' names(drug.nodes) <- drug_names
#' # remove the null value 
#' id  <- sapply(drug.nodes,function(x) !identical(x,character(0)))
#' drug.nodes<- drug.nodes[id]
#' # attack by drugs
#' df <- pbapply::pblapply(1:3, function(x) {
#' oAttackDrug(subg,drug.nodes = drug.nodes,combine.num = x)
#' })
#'  }

oAttackDrug <- function(ig,drug.nodes=NULL,drugs.fixed=NULL,combine.num=1) {
    
    ##check the parameter drug.nodes 
    if(is.null(drug.nodes)){
      stop("A list of drugs with their corresponding target genes is missing")
    }
    ## check drugs.fixed 
    if(!is.null(drugs.fixed)){

      ##fixed drug nodes
      fixed.nodes <-  as.character(unlist(drug.nodes[drugs.fixed]))
      
      drug.nodes <- drug.nodes[!(names(drug.nodes) %in% drugs.fixed)]
      drug.nodes <- lapply(drug.nodes, function(x){ c(x,fixed.nodes)})
      names(drug.nodes) <- sapply(names(drug.nodes),function(x){ paste0(c(x,drugs.fixed),collapse = " + ")})
    }else{
      ## ckeck combine times 
      if (combine.num==1){
        
        drug.nodes <- drug.nodes
        
      }else{
        drug.combined <- combn(names(drug.nodes),combine.num,simplify=F)
        drug.nodes.combine <- lapply(drug.combined, function (x) { 
                             unique(as.character(unlist(sapply(x, function(m){drug.nodes[[m]]}))))})
        names(drug.nodes.combine) <-  sapply(drug.combined, function(x){paste0(x,collapse = " + ")})
        drug.nodes <- drug.nodes.combine
      }
    }
    
    ###Network Attack 
    max.comp.orig <- max(igraph::components(ig)$csize)
    n <- igraph::vcount(ig)
  
    ## drug nodes combinatorial attack
    m <- length(drug.nodes)
    max.comp.removed <- rep(max.comp.orig, m)
    nodes.removed <- rep(max.comp.orig, m)
    removed.pct <- rep(max.comp.orig, m)
    pb <- dplyr::progress_estimated(m)
    for (i in seq_len(m)) {
      pb$tick()$print()
      ind <- match(V(ig)$name, drug.nodes[[i]])
      v <- V(ig)$name[!is.na(ind)]
      nodes.removed[i] <- paste(v, collapse = ",")
      g.manual <- igraph::delete_vertices(ig, v)
      max.comp.removed[i] <- max(igraph::components(g.manual)$csize)
      removed.pct[i] <- 1 - igraph::vcount(g.manual) / n
    }
    comp.pct <- max.comp.removed / max.comp.orig
  
  
  res <- tibble::tibble(drug = names(drug.nodes), frac.disconnected = 1 - comp.pct, 
                        frac.removed = removed.pct, nodes.removed = nodes.removed)
  
  return(res)
}
