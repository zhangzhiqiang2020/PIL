#' Function to visualise enrichment results using a upset-like plot
#'
#' \code{oVisDrugAtt} is supposed to visualise enrichment results using a upset-like plot. The top is the kite plot, and visualised below is the combination matrix for overlapped genes. It returns an object of class "ggplot".
#'
#' @param data a data frame. It contains all these columns including  'value' and 'member' (separated by ',' or ', ')
#' @param color the point color
#' @param shape the point shape
#' @param size the point size
#' @param label.height.unit the unit specifying the per-row height of the combination matrix. By default, it is NULL. If specified (such as 8), it will be used to decide 'combmatrix.label.height' and 'combmatrix.label.text'
#' @param member.sortBy a data frame. It contains all the target nodes and their weights for ranking. 
#' @param drug.levels  how to define the levels of drug in the combination matrix It can be 'num' (the number of non-zeros for a member) or 'customised' (see 'levels.customised' below) or 'auto' (alphabeta order; by default)
#' @param drugs.customised the customised shown of drugs in the combination matrix. 
#'
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{oVisDrugAtt}}
#' @include oVisDrugAtt.r
#' @examples
#'  \dontrun{
#' set.seed(85)
#' data <- rbind(df_best, df_single) %>% distinct() %>% 
#'     arrange(i,frac.disconnected) %>%
#'     select(drug,frac.disconnected, nodes.removed) %>%
#'     rename_with(~c("drug","value","member")) 
#'     p<- data %>% oVisDrugAtt()
#'     }

oVisDrugAtt <- function(data, member.sortBy = NULL, color = "cyan4", shape = 18, size = 2, drugs.customised = NULL) {

  i <- NULL

 levels <- NULL
 df <- NULL
  if (!is.null(drugs.customised)) {
     df <-  data %>% filter(drug %in% drugs.customised) %>% mutate(i=purrr::map_int(drug, ~ stringr::str_split(.x, "+|+ ", simplify = TRUE) %>% length())) %>% arrange(value)
  }else{
    df <- data %>%  mutate(i=purrr::map_int(drug, ~ stringr::str_split(.x, "[+]", simplify = TRUE) %>% length())) %>% arrange(value)
  }

  ########################
  gp<- df %>% mutate(drug=fct_inorder(drug)) %>% 
    ggplot(aes(x=drug,y=value)) + 
    geom_col(fill="steelblue", color="transparent", width=0.1, alpha=0.1) +
    geom_point(aes(color=as.factor(i)),shape =shape, size = size, alpha=1) + 
    guides(color='none') +
    ylab('Effect on the crosstalk measured \nby fraction of nodes disconnected\nupon node removal') +
    xlab('') + 
    theme_classic() + 
    theme(axis.text.x =  element_blank()) 
  
  ##########################
  ##########################

  if (!is.null(member.sortBy)) {
    
    df<- df %>% tidyr::separate_rows(member,sep = ",") %>%mutate(tag=1,member=fct_relevel(member,member.sortBy)) 
    
  } else{
     
    df<- df  %>% tidyr::separate_rows(member,sep = ",") %>%
      mutate(tag=1) 
  }
  ########################
gp_matrix<- df %>% 
             ggplot(aes(x=fct_inorder(drug),y=member)) + 
                  geom_point(aes(fill=tag),color="steelblue",size=2) +
                  geom_line(aes(group=drug),color = "steelblue",linewidth = 0.5) +
                  guides(fill="none") + 
                  labs(x="",y="") +
                  theme_linedraw() +
                  theme(panel.grid= element_line(color = "gray",linewidth = 0.1),
                         axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) 
                 
  gp <- gp/gp_matrix + plot_layout(ncol = 1,heights = c(4,2) )           
 
  gp
}
