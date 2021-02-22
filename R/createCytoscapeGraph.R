#' createCytoscapeGraph

# FIXME: input should be nicely formatted yaml file

createCytoscapeGraph <- function(g0, # netw_g
                                 g1, # ppi_g
                                 tissue,
                                 nodes = NULL, # subset g0
                                 n_cutoffs = 1000,
                                 netwdir = getwd(),
                                 netw_layout = "force-directed",
                                 t_sleep = 2,
                                 min_edge_transparency = 255,
                                 max_edge_transparency = 255,
                                 min_node_size = 35,
                                 max_node_size = 100,
                                 min_edge_color = col2hex("gray"),
                                 max_edge_color = sizzling_red,
                                 threshold = NULL) {

  # input for testing ---------------------------------------------------

  g0 <- netw_g
  g1 <- ppi_g
  n_cutoffs <- 1000
  netwdir <- getwd()
  netw_layout <- "force-directed"
  t_sleep <- 2
  min_edge_transparency <- 255
  max_edge_transparency <- 25
  min_node_size <- 35
  max_node_size <- 100
  min_edge_color <- col2hex("gray")
  max_edge_color <- col2hex("dark red")
  threshold <- NULL

  # misc functions ------------------------------------------------------

  sizzling_red <- "#FF595E"

  is_connected <- function(g, threshold) {
    # helper function to determine if a graph is connected or not a
    # given edge weight threshold
    # wrapper around igraph's:
    # * delete.edges - remove edges from graph
    # * is.connected - check if graph is a single component
    idx <- which(E(g)$weight <= threshold)
    filt_graph <- igraph::delete.edges(g, idx)
    return(igraph::is.connected(filt_graph))
  }

  # imports
  suppressPackageStartupMessages({
    library(RCy3)
  })

  # format graph layout string for RCy3
  cys_layout <- paste(netw_layout, "edgeAttribute=weight")

  # check that we are connected to Cytoscape
  .response <- try(RCy3::cytoscapePing(), silent = TRUE)
  if (.response != "You are connected to Cytoscape!") {
    stop("Start Cytoscape first!")
  }

  rm(list = c(".response"))

  # subset input graph: only keep input nodes
  if (is.null(nodes)) {
    # message("Not subsetting the graph.")
  } else {
    # subset graph
    idx <- match(nodes, names(V(g0)))
    subg0 <- induced_subgraph(g0, vids = V(g0)[idx])
  }

  # node size ~ hubbiness or importance in its subgraph
  adjm0 <- as.matrix(igraph::as_adjacency_matrix(subg0, attr = "weight"))
  node_importance <- apply(adjm0, 2, sum) # weighted degree centrality

  # update graph with node 'size' attribute
  subg0 <- set_vertex_attr(subg0,
    "size",
    value = node_importance[names(V(subg0))]
  )

  # remove weak edges
  # prune weak edges -- edges will be removed until the graph
  # breaks into multiple components. I have tried a matrix
  # implementation of this, but it was not faster. Thus we just
  # manually explore a range of cut offs to find the best
  # points at which the the graph decomposes.
  # Seq from min(edge.weight) to max to generate cutoffs
  # FIXME: only works for weighted graphs!
  if (is.null(threshold)) {
    min_weight <- min(E(subg0)$weight)
    max_weight <- max(E(subg0)$weight)
    cutoffs <- seq(min_weight, max_weight, length.out = n_cutoffs)
    # check if graph is connected or not at various thresholds
    checks <- sapply(cutoffs, function(threshold) {
      is_connected(subg0, threshold)
    })
    # limit is max(cutoff) at which the graph is still connected
    limit <- cutoffs[max(which(checks == TRUE))]
    if (all(checks)) {
      stop("Error thesholding graph.")
    }
    # Prune edges. NOTE: This removes all edge types.
    g <- delete.edges(subg0, which(E(subg0)$weight <= limit))
  } else {
    # Prune edges given a threshold.
    limit <- threshold
    g <- delete.edges(subg0, which(E(subg0)$weight <= limit))
  }

  # write graph to file -------------------------------------------------
  # NOTE: this is faster than sending to Cytoscape through RCy3
  myfile <- file.path(netwdir,paste0(module_name, ".gml"))
  igraph::write_graph(g, myfile, format = "gml") # this format just works!
  #warnings() # boolean converted to numeric

  # Send to Cytoscape with RCy3 ------------------------------------------
  # NOTE: any underscores in the attribute names of g are removed
  winfile <- gsub("/mnt/d/", "D:/", myfile) # must provide Windows path!
  cysnetw <- RCy3::importNetworkFromFile(winfile)
  Sys.sleep(t_sleep)
  unlink(myfile) # don't get ahead of Cytoscape!

  ## RCy3 visual style defaults ------------------------------------------
  # create a visual style
  style_name <- "mystyle"

  # DEFAULTS:
  style_defaults <- list(
    NETWORK_TITLE = "Network",
    NODE_FILL_COLOR = col2hex("gray"),
    NODE_TRANSPARENCY = 255,
    NODE_SIZE = 35,
    NODE_SHAPE = "ellipse",
    NODE_LABEL_TRANSPARENCY = 255,
    NODE_LABEL_FONT_SIZE = 12,
    NODE_LABEL_COLOR = col2hex("black"),
    NODE_BORDER_TRANSPARENCY = 255,
    NODE_BORDER_WIDTH = 5,
    NODE_BORDER_PAINT = col2hex("black"),
    NODE_TRANSPARENCY = 255,
    EDGE_STROKE_UNSELECTED_PAINT = col2hex("black"),
    EDGE_WIDTH = 2,
    NETWORK_BACKGROUND_PAINT = col2hex("white")
  )

  ## Define RCy3 mapped visual attributes -------------------------------

  # mapped properties:
  weight_range <- c(min(E(g)$weight), max(E(g)$weight))
  size_range <- c(min(V(g)$size), max(V(g)$size))
  # color_range <- c(min(V(g)$stat),max(V(g)$stat))
  # node_colors <- c("min"=col2hex("gray"),"max" = unique(V(g)$Color))

  # list of mapped params:
  # NOTE: list is formatted in this way for RCy3
  # FIXME: node color should scale from gray to module color
  # scale is proportional to how much that protein is changing...sum log(p)

  style_map <- list(
    # NODE_FILL_COLOR = mapVisualProperty("node fill color", # RCy3 attr
    # 				     "stat", # column name
    # 				     "c", # use continuous mapping
    # 				     color_range,
    # 				     node_colors),
    NODE_FILL_COLOR = mapVisualProperty(
      "node fill color", # RCy3 attr
      "Color", # column name
      "p"
    ), # use pass-through mapping
    NODE_BORDER_PAINT = mapVisualProperty(
      "node border paint", # RCy3 attr
      "NodeBorder", # column name
      "p"
      ), # use pass through mapping
    NODE_LABEL = mapVisualProperty(
      "node label",
      "protein",
      "p"
    ),
    EDGE_TRANSPARENCY = mapVisualProperty(
      "edge transparency",
      "weight",
      "c", # continuous
      weight_range,
      c(
        min_edge_transparency,
        max_edge_transparency
      )
    ),
    EDGE_STROKE_UNSELECTED_PAINT = mapVisualProperty(
      "edge stroke unselected paint",
      "weight",
      "c",
      weight_range,
      c(
        min_edge_color,
        max_edge_color
      )
    ),
    NODE_SIZE = mapVisualProperty(
      "node size",
      "size",
      "c",
      size_range,
      c(min_node_size, max_node_size)
    )
  )

  # Apply visual style to graph -------------------------------------------

  # create visual style
  RCy3::createVisualStyle(style_name,
    defaults = style_defaults,
    mappings = style_map
  )
  # apply to graph
  setVisualStyle(style_name)
  Sys.sleep(3)

  # Collect edge attributes from g1 --------------------------------------
  # add ppi edges...

  idx <- match(nodes, names(V(g1)))
  subg1 <- induced_subgraph(g1, vids = V(g1)[idx])

  # coerce graph to edge list
  # this doesnt work for weighted graphs!
  edge_list <- apply(igraph::as_edgelist(subg1, names = TRUE), 1, as.list)

  # If edge list is only of length 1, unnest it to avoid problems.
  if (length(edge_list) == 1) {
    edge_list <- unlist(edge_list, recursive = FALSE)
  }

  # Add secondary edge set to Cytoscape graph ----------------------------
  # since a Cytoscape graph can only have one  edge set, we automate the
  # process of adding a second set of edges
  if (length(edge_list) > 0) {
    add_edges <- addCyEdges(edge_list)
    # Add PPIs and set to black
    selected_edges <- selectEdges(add_edges, by.col = "SUID")
    # Set to black with edge bend
    namen <- "EDGE_STROKE_UNSELECTED_PAINT"
    setEdgePropertyBypass(
      edge.names = selected_edges$edges,
      new.values = col2hex("black"),
      visual.property = namen,
      bypass = TRUE
    )
    setEdgePropertyBypass(
      edge.names = selected_edges$edges,
      new.values = TRUE,
      visual.property = "EDGE_BEND",
      bypass = TRUE
    )
  } else {
    # something
  }

  #  clean-up
  clearSelection()
  Sys.sleep(t_sleep)

  # apply graph layout
  layoutNetwork(netw_layout)
  Sys.sleep(t_sleep)
  fitContent()
  Sys.sleep(t_sleep)

  # Mask color of non-significant nodes.
  sig <- names(V(g))[V(g)$SigProt == 1]
  ns <- names(V(g))[names(V(g)) %notin% sig]

  if (length(ns) > 0) {
    setNodeColorBypass(ns, new.colors = col2hex("gray"))
  }

  # Free up some memory.
  suppressMessages({
    cytoscapeFreeMemory()
  })

  # return the threshold used to prune edges
  if (is.null(threshold)) {
    return(limit)
  }
} # EOF
