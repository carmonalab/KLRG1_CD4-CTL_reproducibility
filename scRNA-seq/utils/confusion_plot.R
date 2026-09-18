plot.confusion <- function(object = NULL,
                           var.1 = NULL,
                           var.2 = NULL,
                           relative = FALSE,
                           useNA = "ifany",
                           type = "tile",
                           xlab = NULL,
                           ylab = NULL,
                           plot.title = "") {
   
   if (is.null(object)) {
      stop("Please provide a Seurat object or its metadata")
   }
   
   if (inherits(object, "Seurat")) {
      data <- object@meta.data
   } else if (is.data.frame(object)) {
      data <- object
   } else {
      stop("Please provide a Seurat object or its metadata")
   }
   
   legend.title <- "Number of cells"
   data <- as.data.frame(data)
   
   
   if (is.null(xlab)) {
      xlab <-  var.1
   }
   if (is.null(ylab)) {
      ylab <-  var.2
   }
   
   vars <- c(var.1, var.2)
   if (any(is.null(vars))) {
      stop("Please provide 2 cell type classification columns in metadata as var.1 and var.2")
   }
   
   if (any(!vars %in% names(data))) {
      stop("Not all classification variables: ",
           paste(vars, collapse = ", "),
           " were found in metadata")
   }
   
   # Get counts for every group
   pz <- table(data[[var.1]],
               data[[var.2]],
               useNA = useNA)
   
   if (relative) {
      pz <- prop.table(pz, margin = 1) %>%
         round(digits = 2)
      legend.title <- "Proportion of cells"
   }
   
   
   if (tolower(type) == "tile") {
      plot <- pz %>%
         as.data.frame() %>%
         ggplot2::ggplot(ggplot2::aes(Var1, Var2,
                                      fill = Freq)) +
         ggplot2::geom_tile() +
         ggplot2::scale_fill_gradient(name = legend.title,
                                      low = "#FFFFC8",
                                      high = "#7D0025") +
         ggplot2::labs(x = xlab,
                       y = ylab,
                       title = plot.title) +
         ggplot2::geom_label(ggplot2::aes(label = ifelse(Freq > 0,
                                                         round(Freq, 1),
                                                         NA)),
                             color = "white",
                             alpha = 0.6,
                             fill = "black")
      
   } else if (grepl("bar|col", type, ignore.case = TRUE)) {
      plot <- pz %>%
         as.data.frame() %>%
         ggplot2::ggplot(ggplot2::aes(Var1, Freq,
                                      fill = Var2)) +
         ggplot2::labs(x = xlab,
                       y = legend.title,
                       fill = ylab,
                       title = plot.title) +
         ggplot2::geom_col()
   }
   
   plot <- plot +
      theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         hjust = 1,
                                                         vjust = 1))
   return(plot)
}
