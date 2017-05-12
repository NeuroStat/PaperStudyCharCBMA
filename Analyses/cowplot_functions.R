####################
#### TITLE:     Functions of the cowplot package.
#### Contents:
####
#### Source Files: //Meta\ Analyis/R\ Code/Studie_FixRan/FixRanStudy.git/Development
#### First Modified: 03/01/2017
#### Notes:
#################

##
###############
### Notes
###############
##

# The cowplot package is not completely up and running on my machine to date (03/01/2017).
# However, I was able to pull the functions that I want to use from the Github page (https://github.com/wilkelab/cowplot).
# Here I collect them so I can source in other R files to this file.

## I CLAIM NO AUTHORSHIP NOR RESPONSIBILITY FOR THIS CODE.

# *************************************************
#                     Themes
# *************************************************

#' Create the default cowplot theme
#'
#' After loading the cowplot package, this theme will be the default
#' for all graphs made with ggplot2.
#'
#' @param font_size Overall font size. Default is 14.
#' @param font_family Default font family.
#' @param line_size Default line size.
#' @return The theme.
#' @examples
#' qplot(1:10, (1:10)^2) + theme_cowplot(font_size = 15)
#' @export
theme_cowplot <- function(font_size = 14, font_family = "", line_size = .5) {
  half_line <- font_size / 2
  small_rel <- 0.857
  small_size <- small_rel * font_size

  theme_grey(base_size = font_size, base_family = font_family) %+replace%
    theme(
      rect              = element_rect(fill = "transparent", colour = NA, color = NA, size = 0, linetype = 0),
      text              = element_text(family = font_family, face = "plain", colour = "black",
                                       size = font_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = .9,
                                       margin = ggplot2::margin(), debug = FALSE),
      axis.text         = element_text(colour = "black", size = small_size),
      #axis.title        = element_text(face = "bold"),
      axis.text.x       = element_text(margin = ggplot2::margin(t = small_size / 4), vjust = 1),
      axis.text.y       = element_text(margin = ggplot2::margin(r = small_size / 4), hjust = 1),
      axis.title.x      = element_text(
        margin = ggplot2::margin(t = small_size / 2, b = small_size / 4)
      ),
      axis.title.y      = element_text(
        angle = 90,
        margin = ggplot2::margin(r = small_size / 2, l = small_size / 4),
      ),
      axis.ticks        = element_line(colour = "black", size = line_size),
      axis.line         = element_line(colour = "black", size = line_size),
      # the following two lines are not needed for ggplot2 2.2.0 or later
      axis.line.x         = element_line(colour = "black", size = line_size),
      axis.line.y         = element_line(colour = "black", size = line_size),
      legend.key        = element_blank(),
      legend.key.size   = grid::unit(1, "lines"),
      legend.text       = element_text(size = rel(small_rel)),
      #    legend.position   = c(-0.03, 1.05),
      #    legend.justification = c("left", "top"),
      panel.background  = element_blank(),
      panel.border      = element_blank(),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(size = rel(small_rel)),
      strip.background  = element_rect(fill = "grey80", colour = "grey50", size = 0),
      plot.background   = element_blank(),
      plot.title        = element_text(face = "bold",
                                       size = font_size,
                                       margin = ggplot2::margin(b = half_line)),

      complete = TRUE
    )
}


#' Create a completely empty theme
#'
#' The theme created by this function shows nothing but the plot panel. Since ggplot2 now
#' provides \code{theme_void}, which serves the same purpose, this function is now just a
#' wrapper around \code{theme_void}.
#' @param base_size Overall font size. Default is 14.
#' @param base_family Base font family.
#' @return The theme.
#' @examples
#' qplot(1:10, (1:10)^2) + theme_nothing()
#' @export
theme_nothing <- function(base_size = 14, base_family = ""){
  theme_void(base_size = base_size, base_family = base_family)
}

#' Add/modify/remove the background grid in a ggplot2 plot
#'
#' This function provides a simple way to modify the background grid in ggplot2. It
#' doesn't do anything that can't be done just the same with \code{theme()}. However, it simplifies
#' creation of the most commonly needed variations.
#' @param major Specifies along which axes you would like to plot major grid lines. Options are "xy", "x",
#'  "y", "only_minor" (disables major grid lines but allows you to switch on minor grid lines), "none".
#' @param minor Specifies along which axes you would like to plot minor grid lines. Options are "xy", "x",
#'  "y", "none".
#' @param size.major Size of the major grid lines.
#' @param size.minor Size of the minor grid lines.
#' @param colour.major Color of the major grid lines.
#' @param colour.minor Color of the minor grid lines.
#' @export
background_grid <- function(major = c("xy", "x", "y", "only_minor", "none"),
                            minor = c("xy", "x", "y", "none"),
                            size.major = 0.2, size.minor = 0.5,
                            colour.major = "grey90", colour.minor = "grey98"){

  if (major[1] == "none") return(theme(panel.grid = element_blank()))

  t <- switch( major[1],
               x = theme(panel.grid.major   = element_line(colour = colour.major,
                                                           size = size.major),
                         panel.grid.major.y = element_blank()),
               y = theme(panel.grid.major   = element_line(colour = colour.major,
                                                           size = size.major),
                         panel.grid.major.x = element_blank()),
               xy = theme(panel.grid.major = element_line(colour = colour.major,
                                                          size = size.major)),
               yx = theme(panel.grid.major = element_line(colour = colour.major,
                                                          size = size.major)),
               theme(panel.grid.major = element_blank()))

  t + switch( minor[1],
              x = theme(panel.grid.minor   = element_line(colour = colour.minor,
                                                          size = size.minor),
                        panel.grid.minor.y = element_blank()),
              y = theme(panel.grid.minor   = element_line(colour = colour.minor,
                                                          size = size.minor),
                        panel.grid.minor.x = element_blank()),
              xy = theme(panel.grid.minor = element_line(colour = colour.minor,
                                                         size = size.minor)),
              yx = theme(panel.grid.minor = element_line(colour = colour.minor,
                                                         size = size.minor)),
              theme(panel.grid.minor = element_blank()))
}


#' Add/remove the panel border in a ggplot2 plot
#'
#' This function provides a simple way to modify the panel border in ggplot2. It
#' doesn't do anything that can't be done just the same with \code{theme()}. However, it
#' saves some typing.
#' @param colour The color of the border.
#' @param size Size.
#' @param linetype Line type.
#' @param remove If \code{TRUE}, removes the current panel border.
#' @export
panel_border <- function(colour = 'gray80', size = 0.5, linetype = 1, remove = FALSE){
  if (remove){
    return(theme(panel.border = element_blank()))
  }
  theme(panel.border = element_rect(colour = colour, fill=NA, linetype = linetype, size = size))
}

# *************************************************
#                     Drawing code
# *************************************************


# ****** Internal functions used by drawing code ******
ggplot_to_gtable <- function(plot)
{
  if (methods::is(plot, "ggplot")){
    # ggplotGrob must open a device and when a multiple page capable device (e.g. PDF) is open this will save a blank page
    # in order to avoid saving this blank page to the final target device a NULL device is opened and closed here to *absorb* the blank plot

    # commenting this out to see if it was the cause of
    grDevices::pdf(NULL)
    plot <- ggplot2::ggplotGrob(plot)
    grDevices::dev.off()
    plot
  }
  else if (methods::is(plot, "gtable")){
    plot
  }
  else{
    stop('Argument needs to be of class "ggplot" or "gtable"' )
  }
}


#' Draw a line.
#'
#' This is a convenience function. It's just a thin wrapper around \code{geom_line}.
#'
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param ... Style parameters, such as \code{colour}, \code{alpha}, \code{size}, etc.
#' @export
draw_line <- function(x, y, ...){
  geom_path(data = data.frame(x, y),
            aes(x = x, y = y),
            ...)
}

#' Draw text.
#'
#' This is a convenience function to plot multiple pieces of text at the same time. It cannot
#' handle mathematical expressions, though. For those, use \code{draw_label}.
#'
#' Note that font sizes get scaled by a factor of 2.85, so sizes given here agree with font sizes used in
#' the theme. This is different from \code{geom_text} in ggplot2.
#'
#' By default, the x and y coordinates specify the center of the text box. Set \code{hjust = 0, vjust = 0} to specify
#' the lower left corner, and other values of \code{hjust} and \code{vjust} for any other relative location you want to
#' specify.
#' @param text Character or expression vector specifying the text to be written.
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param size Font size of the text to be drawn.
#' @param ... Style parameters, such as \code{colour}, \code{alpha}, \code{angle}, \code{size}, etc.
#' @export
draw_text <- function(text, x = 0.5, y = 0.5, size = 14, ...){
  geom_text(data = data.frame(text, x, y),
            aes(label = text, x = x, y = y),
            size = size / .pt, # scale font size to match size in theme definition
            ...)
}


#' Draw a text label or mathematical expression.
#'
#' This function can draw either a character string or mathematical expression at the given
#' coordinates. It works both on top of \code{ggdraw} and directly with \code{ggplot}, depending
#' on which coordinate system is desired (see examples).
#'
#' By default, the x and y coordinates specify the center of the text box. Set \code{hjust = 0, vjust = 0} to specify
#' the lower left corner, and other values of \code{hjust} and \code{vjust} for any other relative location you want to
#' specify.
#' @param label String or plotmath expression to be drawn.
#' @param x The x location of the label.
#' @param y The y location of the label.
#' @param hjust Horizontal justification
#' @param vjust Vertical justification
#' @param fontfamily The font family
#' @param fontface The font face ("plain", "bold", etc.)
#' @param colour Text color
#' @param size Point size of text
#' @param angle Angle at which text is drawn
#' @param lineheight Line height of text
#' @param alpha The alpha value of the text
#' @examples
#' p <- ggplot(mtcars, aes(mpg, disp)) + geom_line(colour = "blue") + background_grid(minor='none')
#' c <- cor.test(mtcars$mpg, mtcars$disp, method='sp')
#' label <- substitute(paste("Spearman ", rho, " = ", estimate, ", P = ", pvalue),
#'                     list(estimate = signif(c$estimate, 2), pvalue = signif(c$p.value, 2)))
#' # adding label via ggdraw, in the ggdraw coordinates
#' ggdraw(p) + draw_label(label, .7, .9)
#' # adding label directly to plot, in the data coordinates
#' p + draw_label(label, 20, 400, hjust = 0, vjust = 0)
#' @export
draw_label <- function(label, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
                    fontfamily = "", fontface = "plain", colour = "black", size = 14,
                    angle = 0, lineheight = 0.9, alpha = 1)
{
  text_par <- grid::gpar(col = colour,
                         fontsize = size,
                         fontfamily = fontfamily,
                         fontface = fontface,
                         lineheight = lineheight,
                         alpha = alpha)

  # render the label
  text.grob <- grid::textGrob(label, x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"),
                             hjust = hjust, vjust = vjust, rot = angle, gp = text_par)
  annotation_custom(text.grob, xmin = x, xmax = x, ymin = y, ymax = y)
}


#' Add a label to a plot
#'
#' This function adds a plot label to the upper left corner of a graph (or an arbitrarily specified position). It takes all the same parameters
#' as \code{draw_text}, but has defaults that make it convenient to label graphs with letters A, B, C, etc. Just like \code{draw_text()},
#' it can handle vectors of labels with associated coordinates.
#' @param label String (or vector of strings) to be drawn as the label.
#' @param x The x position (or vector thereof) of the label(s).
#' @param y The y position (or vector thereof) of the label(s).
#' @param hjust Horizontal adjustment.
#' @param vjust Vertical adjustment.
#' @param size Font size of the label to be drawn.
#' @param fontface Font face of the label to be drawn.
#' @param ... Other arguments to be handed to \code{draw_text}.
#' @export
draw_plot_label <- function(label, x=0, y=1, hjust = -0.5, vjust = 1.5, size = 16, fontface = 'bold', ...){
  draw_text(text = label, x = x, y = y, hjust = hjust, vjust = vjust, size = size, fontface = fontface, ...)
}


#' Add a label to a figure
#'
#' This function is similar to \code{draw_plot_label}, just with slightly different arguments and defaults. The main purpose of this
#' function is to add labels specifying extra information about the figure, such as "Figure 1", which is sometimes useful.
#' @param label Label to be drawn
#' @param position Position of the label, can be one of "top.left", "top", "top.right", "bottom.left", "bottom", "bottom.right". Default is "top.left"
#' @param size (optional) Size of the label to be drawn. Default is the text size of the current theme
#' @param fontface (optional) Font face of the label to be drawn. Default is the font face of the current theme
#' @param ... other arguments passed to \code{draw_plot_label}
#'
#' @examples
#'
#' p1 <- qplot(1:10, 1:10)
#' p2 <- qplot(1:10, (1:10)^2)
#' p3 <- qplot(1:10, (1:10)^3)
#' p4 <- qplot(1:10, (1:10)^4)
#'
#' # Create a simple grid
#' p <- plot_grid(p1, p2, p3, p4, align = 'hv')
#'
#' # Default font size and position
#' p + draw_figure_label(label = "Figure 1")
#'
#' # Different position and font size
#' p + draw_figure_label(label = "Figure 1", position = "bottom.right", size = 10)
#'
#' # Using bold font face
#' p + draw_figure_label(label = "Figure 1", fontface = "bold")
#'
#' # Making the label red and slanted
#' p + draw_figure_label(label = "Figure 1", angle = -45, colour = "red")
#'
#' # Labeling an individual plot
#' ggdraw(p2) + draw_figure_label(label = "Figure 1", position = "bottom.right", size = 10)
#'
#' @author Ulrik Stervbo (ulrik.stervbo [AT] gmail.com)
#' @export
draw_figure_label <- function(label, position = c("top.left", "top", "top.right", "bottom.left", "bottom", "bottom.right"), size, fontface, ...){
  # Get the position
  position <- match.arg(position)

  # Set default font size and face from the theme
  if(missing(size)){
    size <- theme_get()$text$size
  }
  if(missing(fontface)){
    fontface <- theme_get()$text$face
  }

  # Call draw_plot_label() with appropriate label positions
  switch(position,
         top.left     = draw_plot_label(label, x = 0,   y = 1, hjust = -0.1, vjust = 1.1,  size = size, fontface = fontface, ...),
         top          = draw_plot_label(label, x = 0.5, y = 1, hjust = 0,    vjust = 1.1,  size = size, fontface = fontface, ...),
         top.right    = draw_plot_label(label, x = 1,   y = 1, hjust = 1.1,  vjust = 1.1,  size = size, fontface = fontface, ...),
         bottom.left  = draw_plot_label(label, x = 0,   y = 0, hjust = -0.1, vjust = -0.1, size = size, fontface = fontface, ...),
         bottom       = draw_plot_label(label, x = 0.5, y = 0, hjust = 0,    vjust = -0.1, size = size, fontface = fontface, ...),
         bottom.right = draw_plot_label(label, x = 1,   y = 0, hjust = 1.1,  vjust = -0.1, size = size, fontface = fontface, ...)
  )
}

#' Draw a (sub)plot.
#'
#' Places a plot somewhere onto the drawing canvas. By default, coordinates run from
#' 0 to 1, and the point (0, 0) is in the lower left corner of the canvas.
#' @param plot The plot to place. Can be either a ggplot2 plot or an arbitrary gtable.
#' @param x The x location of the lower left corner of the plot.
#' @param y The y location of the lower left corner of the plot.
#' @param width Width of the plot.
#' @param height Height of the plot.
#' @export
draw_plot <- function(plot, x = 0, y = 0, width = 1, height = 1){
  g <- ggplot_to_gtable(plot) # convert to gtable if necessary
  plot.grob <- grid::grobTree(g)
  annotation_custom(plot.grob, xmin = x, xmax = x+width, ymin = y, ymax = y+height)
}

#' Draw a grob.
#'
#' Places an arbitrary grob somewhere onto the drawing canvas. By default, coordinates run from
#' 0 to 1, and the point (0, 0) is in the lower left corner of the canvas.
#' @param grob The grob to place.
#' @param x The x location of the lower left corner of the grob.
#' @param y The y location of the lower left corner of the grob.
#' @param width Width of the grob.
#' @param height Height of the grob.
#' @export
draw_grob <- function(grob, x = 0, y = 0, width = 1, height = 1){
  annotation_custom(grid::grobTree(grob), xmin = x, xmax = x+width, ymin = y, ymax = y+height)
}



#' Set up a drawing layer on top of a ggplot
#' @param plot The plot to use as a starting point. Can be either a ggplot2 plot or an arbitrary gtable.
#' @param xlim The x-axis limits for the drawing layer (default is [0, 1]).
#' @param ylim The y-axis limits for the drawing layer (default is [0, 1]).
#' @export
ggdraw <- function(plot = NULL, xlim = c(0, 1), ylim = c(0, 1)) {

  d <- data.frame(x=0:1, y=0:1) # dummy data
  p <- ggplot(d, aes_string(x="x", y="y")) + # empty plot
    scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    theme_nothing() + # with empty theme
    labs(x=NULL, y=NULL) # and absolutely no axes

  if (!is.null(plot)){
    g <- ggplot_to_gtable(plot) # convert to gtable if necessary
    plot.grob <- grid::grobTree(g)
    p <- p + annotation_custom(plot.grob)
  }
  p # return ggplot drawing layer
}

#' Align multiple plots vertically and/or horizontally
#'
#' Align multiple plots for plotting manually. Can be used to graph two separate y axis, but still doesn't work if second y axis needs to be shown.
#' @param ... List of plots to be aligned.
#' @param plotlist (optional) List of plots to display. Alternatively, the plots can be provided
#' individually as the first n arguments of the function plot_grid (see examples).
#' @param align (optional) Specifies whether graphs in the grid should be horizontally ("h") or
#'  vertically ("v") aligned. Options are "none" (default), "hv" (align in both directions), "h", and "v".
#' @examples
#' p1 <- qplot(1:10, rpois(10, lambda=15), geom="point")
#' p2 <- qplot(1:10, (1:10)^2, geom="line") + theme_nothing()
#' # manually align and plot on top of each other
#' aligned_plots <- align_plots(p1, p2, align="hv")
#' ggdraw() + draw_grob(aligned_plots[[1]]) + draw_grob(aligned_plots[[2]])
#' @export
align_plots <- function(..., plotlist = NULL, align = c("none", "h", "v", "hv")){
  plots <- c(list(...), plotlist)
  num_plots <- length(plots)

  # convert list of plots into list of grobs / gtables
  grobs <- lapply(plots, function(x) {if (!is.null(x)) ggplot_to_gtable(x) else NULL})


  #aligning graphs.
  halign <- switch(align[1],
                   h = TRUE,
                   vh = TRUE,
                   hv = TRUE,
                   FALSE
  )
  valign <- switch(align[1],
                   v = TRUE,
                   vh = TRUE,
                   hv = TRUE,
                   FALSE
  )

  # calculate the maximum widths and heights over all graphs, and find out whether
  # they can be aligned if necessary
  if (valign)
  {
    num_widths <- unique(lapply(grobs, function(x){length(x$widths)})) # count number of unique lengths
    num_widths[num_widths==0] <- NULL # remove entry for missing graphs
    if (length(num_widths) > 1)
    {
      valign = FALSE
      warning("Graphs cannot be vertically aligned. Placing graphs unaligned.")
    }
    else
    {
      max_widths <- do.call(grid::unit.pmax, lapply(grobs, function(x){x$widths}))
    }
  }

  if (halign)
  {
    num_heights <- unique(lapply(grobs, function(x){length(x$heights)})) # count number of unique lengths
    num_heights[num_heights==0] <- NULL # remove entry for missing graphs
    if (length(num_heights) > 1)
    {
      halign = FALSE
      warning("Graphs cannot be horizontally aligned. Placing graphs unaligned.")
    }
    else
    {
      max_heights <- do.call(grid::unit.pmax, lapply(grobs, function(x){x$heights}))
    }
  }

  # now assing to all graphs
  for ( i in 1:num_plots )
  {
    if (!is.null(grobs[[i]]))
    {
      if (valign) grobs[[i]]$widths <- max_widths
      if (halign) grobs[[i]]$heights <- max_heights
    }
  }
  grobs
}




#' Arrange multiple plots into a grid
#'
#' Arrange multiple plots into a grid.
#' @param ... List of plots to be arranged into the grid. The plots can be either ggplot2 plot objects
#'              or arbitrary gtables.
#' @param plotlist (optional) List of plots to display. Alternatively, the plots can be provided
#' individually as the first n arguments of the function plot_grid (see examples).
#' @param align (optional) Specifies whether graphs in the grid should be horizontally ("h") or
#'  vertically ("v") aligned. Options are "none" (default), "hv" (align in both directions), "h", and "v".
#' @param nrow (optional) Number of rows in the plot grid.
#' @param ncol (optional) Number of columns in the plot grid.
#' @param scale (optional) Allows to set an overall scaling of each sub-plot. Can be set separately for
#'              each subplot, by giving a vector of scale values, or at once for all subplots,
#'              by giving a single value.
#' @param rel_widths (optional) Numerical vector of relative columns widths. For example, in a two-column
#'              grid, \code{rel_widths = c(2, 1)} would make the first column twice as wide as the
#'              second column.
#' @param rel_heights (optional) Numerical vector of relative columns heights. Works just as
#'              \code{rel_widths} does, but for rows rather than columns.
#' @param labels (optional) List of labels to be added to the plots. You can also set \code{labels="AUTO"} to
#'              auto-generate upper-case labels or \code{labels="auto"} to auto-generate lower-case labels.
#' @param label_size (optional) Numerical value indicating the label size. Default is 14.
#' @param hjust Adjusts the horizontal position of each label. More negative values move the label further
#'              to the right on the plot canvas. Default is -0.5.
#' @param vjust Adjusts the vertical position of each label. More positive values move the label further
#'              down on the plot canvas. Default is 1.5.
#' @param rows Deprecated. Like \code{nrow}.
#' @param cols Deprecated. Like \code{ncol}.
#' @examples
#' p1 <- qplot(1:10, 1:10)
#' p2 <- qplot(1:10, (1:10)^2)
#' p3 <- qplot(1:10, (1:10)^3)
#' p4 <- qplot(1:10, (1:10)^4)
#' # simple grid
#' plot_grid(p1, p2, p3, p4)
#' # simple grid with labels and aligned plots
#' plot_grid(p1, p2, p3, p4, labels=c('A', 'B', 'C', 'D'), align="hv")
#' # manually setting the number of rows, auto-generate upper-case labels
#' plot_grid(p1, p2, p3, nrow=3, labels="AUTO", label_size=12, align="v")
#' # missing plots in some grid locations, auto-generate lower-case labels
#' plot_grid(p1, NULL, NULL, p2, p3, NULL, ncol=2,
#'  labels="auto", label_size=12, align="v")
#' # making rows and columns of different widths/heights
#' plot_grid(p1, p2, p3, p4, align='hv', rel_heights=c(2,1), rel_widths=c(1,2))
#' @export
plot_grid <- function(..., plotlist = NULL, align = c("none", "h", "v", "hv"),
                      nrow = NULL, ncol = NULL, scale = 1, rel_widths = 1,
                      rel_heights = 1, labels = NULL, label_size = 14,
                      hjust = -0.5, vjust = 1.5,
                      cols = NULL, rows = NULL ) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  num_plots <- length(plots)

  if (!is.null(cols)){
    warning("Argument 'cols' is deprecated. Use 'ncol' instead.")
  }

  if (!is.null(rows)){
    warning("Argument 'rows' is deprecated. Use 'nrow' instead.")
  }

  # internally, this function operates with variables cols and rows instead of ncol and nrow
  if (!is.null(ncol)){
    cols <- ncol
  }
  if (!is.null(nrow)){
    rows <- nrow
  }

  # Align the plots (if specified)
  grobs <- align_plots(plotlist = plots, align=align)

  # calculate grid dimensions
  if (is.null(cols) && is.null(rows)){
    # if neither rows nor cols are given, we make a square grid
    cols <- ceiling(sqrt(num_plots))
    rows <- ceiling(num_plots/cols)
  }
  # alternatively, we know at least how many rows or how many columns we need
  if (is.null(cols)) cols <- ceiling(num_plots/rows)
  if (is.null(rows)) rows <- ceiling(num_plots/cols)

  # in general, we allow for separate scale values for each graph
  if (length(scale)==1)
    scale <- rep(scale, num_plots)

  if ("AUTO" %in% labels)
    labels <- LETTERS[1:num_plots]
  else if ("auto" %in% labels)
    labels <- letters[1:num_plots]

  # we also allow for separate hjust and vjust values for each label
  if (length(hjust)==1)
    hjust <- rep(hjust, length(labels))

  if (length(vjust)==1)
    vjust <- rep(vjust, length(labels))

  # calculate appropriate vectors of rel. heights and widths
  rel_heights <- rep(rel_heights, length.out = rows)
  rel_widths <- rep(rel_widths, length.out = cols)
  # calculate the appropriate coordinates and deltas for each row and column
  x_deltas <- rel_widths/sum(rel_widths)
  y_deltas <- rel_heights/sum(rel_heights)
  xs <- cumsum(rel_widths)/sum(rel_widths) - x_deltas
  ys <- 1 - cumsum(rel_heights)/sum(rel_heights)

  # now place all the plots
  p <- ggdraw() # start with nothing
  col_count <- 0
  row_count <- 1
  for (i in 1:(rows*cols)){
    if (i > num_plots) break

    x_delta <- x_deltas[col_count+1]
    y_delta <- y_deltas[row_count]
    # calculate width, offset, etc
    width <- x_delta * scale[i]
    height <- y_delta * scale[i]
    x_off <- (x_delta - width)/2
    y_off <- (y_delta - height)/2
    x <- xs[col_count+1] + x_off
    y <- ys[row_count] + y_off

    # place the plot
    p_next <- grobs[[i]]
    if (!is.null(p_next)){
      p <- p + draw_grob(grid::grobTree(p_next), x, y, width, height)
    }
    # place a label if we have one
    if (i <= length(labels)){
      p <- p + draw_plot_label(labels[i], x - x_off, y + height - y_off, size = label_size,
                               hjust = hjust[i], vjust = vjust[i])
    }
    # move on to next grid position
    col_count <- col_count + 1
    if (col_count >= cols){
      col_count <- 0
      row_count <- row_count + 1
    }
  }
  p
}

##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################


